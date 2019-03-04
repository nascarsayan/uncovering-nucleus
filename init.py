import networkx as nx
import csv
import matplotlib.pyplot as plt
from sys import argv
import json
import os
import re
from collections import Counter
import plotly.offline as py
import plotly.graph_objs as go

efile = './dummyDataset/edges.csv'
plotlykey = 'LuM1q5YsmH2ItxMWeE8I'
if len(argv) > 1:
  efile = argv[1]
corenum = None
if len(argv) > 2:
  corenum = int(argv[2])
dataset = os.path.splitext(os.path.basename(efile))[0]
outputDir = './out/%s' % (efile)


def dumpData(data, fname):
  with open(fname, 'w') as fp:
    json.dump(data, fp)


def readEdges(fname):
  rows = []
  if fname.endswith('.txt'):
    with open(fname, 'r') as txtfile:
      for line in txtfile:
        if ((line.strip())[0]).isdigit():
          break
      line = txtfile.readline()
      delim = '[ \t]+'
      if ',' in line:
        delim = ','
      rows.append((list(
          map(lambda x: int(x.strip()), re.split(delim, line.strip()))))[:2])
      for line in txtfile:
        try:
          rows.append((list(
              map(lambda x: int(x.strip()), re.split(delim,
                                                     line.strip()))))[:2])
        except Exception as _:
          pass
  else:
    with open(fname, 'r') as csvfile:
      csvr = csv.reader(csvfile)
      for row in csvr:
        rows.append(list(map(int, row)))

  return rows


def drawGraph(G, figname='graph.png'):
  nx.draw(G, with_labels=True, font_weight='bold')
  plt.savefig(figname)
  plt.clf()


def lineGraphs(y, label, title='line-graph', filename='line.png', clf=True):
  plt.ylim(0, max(list(map(lambda ye: max(ye), y))))
  xe = list(range(len(y[0])))
  for i in range(0, 12, 2):
    plt.plot(xe, y[i], '.-', label='B=%.1f' % (label[i]))
  plt.legend()
  plt.title(title)
  plt.xlabel('G_k (k-core)')
  plt.ylabel('NI (G_k, dep(i, β))')
  plt.savefig(filename)
  if clf:
    plt.clf()


def plotlyGraphs(y, label, title='line-graph', filename='line.png', clf=True):
  layout = go.Layout(
    title=title,
    xaxis=dict(
        title='G_k (k-core)s',
    ),
    yaxis=dict(
        title='NI (G_k, dep(i, β))',
    )
)
  x = list(range(len(y[0])))
  trace = list(map(lambda i: go.Scatter(x=x, y=y[i], name=label[i]), list(range(len(y)))))
  fig = go.Figure(data=trace, layout=layout)
  py.plot(fig, filename=filename)


def getDep(G, cores, outputDir):
  betas = list(map(lambda x: x / 10, list(range(11))))
  nodes = G.nodes()
  edges = G.edges()
  k_max = max(cores.values())
  print('k_max = %d' %(k_max))
  ni = []
  aggni = []
  for beta in betas:
    print('+++ β = %0.1f' % (beta))
    depBeta = [{v: 0 for v in nodes}]
    for k in range(1, k_max + 1):
      depBetaK = {}
      for v in nodes:
        depBetaK[v] = depBeta[k - 1][v]
        if k > cores[v]:
          continue
        nbrs = list(
            filter(lambda u: cores[u] == k, list(nx.all_neighbors(G, v))))
        depBetaK[v] += len(nbrs) + beta * sum(
            list(map(lambda u: depBeta[k - 1][u], nbrs)))
      depBeta.append(depBetaK)
    # print(list(map(lambda x: list(x.values()) , depBeta)))
    dumpData(list(map(lambda x: list(x.values()) , depBeta)), '%s/dep/%d.json' % (outputDir, beta * 10))
    sumDep = sum(depBeta[-1].values())
    vk_1 = nodes
    ek_1 = edges
    nib = [0]
    for k in range(1, k_max):
      vk = list(filter(lambda u: cores[u] > k, vk_1))
      ek = list(filter(lambda e: cores[e[0]] > k and cores[e[1]] > k, ek_1))
      nibk = (1 / len(vk_1)) * (len(ek) / (len(vk) * (len(vk) - 1))) * (
          sum(list(map(lambda i: depBeta[-1][i], vk))) / sumDep)
      nib.append(nibk)
      vk_1 = vk
      ek_1 = ek
    nib += [0]
    dumpData(nib, '%s/ni/%d.json' % (outputDir, beta * 10))
    mni = nib.index(max(nib))
    print('Nuclear Index = ', mni)
    ni.append(nib)
    aggni.append(mni)
    # print('# nodes in the core with max NI = %d; # nodes in the core, with max k = %d' % (len(list(filter(lambda x: cores[x] >= mni, nodes))), len(list(filter(lambda x: cores[x] == k_max, nodes)))))
  # lineGraphs(ni, betas, '%s NuclearIndex' % (dataset),
  #            '%s/ni/graph.png' % (outputDir))
  plotlyGraphs(ni, betas, '%s NuclearIndex' % (dataset),
             '%s/ni/graph.html' % (outputDir))
  cnt = Counter(aggni)
  kc = (cnt.most_common(1))[0][0]
  chosenBs = list(
               map(lambda x: betas[x],
                   filter(lambda y: aggni[y] == kc, range(10))))
  print('@@@ k_C = %d, β chosen = %r' % (kc, chosenBs[len(chosenBs) // 2]))
  nucleus = nx.k_core(G, kc)
  print('N(GC) = %d, E(GC) = %d' %(nucleus.number_of_nodes(), nucleus.number_of_edges()))
  return kc

def plotlyCentralities(cens, G, cores, kc, folder='./'):
  for cen in cens:
    print('Calculating %s' %(cen['name']))
    vals = cen['fn'](G)
    print('Plotting')
    h1 = [vals[k] for k in vals.keys() if cores[k] < kc]
    h2 = [vals[k] for k in vals.keys() if cores[k] >= kc]
    layout = go.Layout(
        title='%s Distribution' % (cen['name']),
        barmode='stack',
    )
    data = [
        go.Histogram(x=h1, name='Outside K_C core'),
        go.Histogram(x=h2, name='Inside K_C core')
    ]
    fig = go.Figure(data=data, layout=layout)
    py.plot(fig, filename='%s/%s.html' % (folder, cen['name']))

def main():
  if not (os.path.exists(outputDir)):
    os.makedirs(outputDir)
  if not (os.path.exists('%s/dep' % outputDir)):
    os.makedirs('%s/dep' % outputDir)
  if not (os.path.exists('%s/ni' % outputDir)):
    os.makedirs('%s/ni' % outputDir)
  print('Reading edges...')
  edges = readEdges(efile)
  G = nx.Graph()
  print('Done!\nAdding edges...')
  G.add_edges_from(edges)
  print('Done!\nCalculating core values after removing self-loops...')
  G.remove_edges_from(nx.selfloop_edges(G))
  print('# nodes = %d, # edges = %d' %(G.number_of_nodes(), G.number_of_edges()))
  cores = nx.core_number(G)
  if corenum is not None:
    queryg = nx.k_core(G, corenum)
    print('kth core contains %d nodes and %d edges' % (queryg.number_of_nodes(), queryg.number_of_edges()))
  coreFreq = Counter(cores.values())
  cdf = []; s = 0
  for x in range(max(coreFreq) - 1, -1, -1):
    s = s + coreFreq[x]
    cdf.append(s)
  print(list(reversed(cdf)))
  dumpData(cores, '%s/cores.json' % (outputDir))
  print('Done!\nGetting dependency values...',)
  kc = getDep(G, cores, outputDir)
  cenf = nx.algorithms.centrality
  cens = [{
      'fn': cenf.degree_centrality,
      'name': 'DegreeCentrality'
  }, {
      'fn': cenf.closeness_centrality,
      'name': 'ClosenessCentrality'
  }, {
      'fn': cenf.betweenness_centrality,
      'name': 'BetweennessCentrality'
  }, {
      'fn': cenf.eigenvector_centrality,
      'name': 'EigenvectorCentrality'
  }]
  plotlyCentralities(cens, G, cores, kc, '%s/ni/' % (outputDir))


if __name__ == '__main__':
  main()