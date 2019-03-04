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

def dumpData(data, fname):
  with open(fname, 'w') as fp:
    json.dump(data, fp)

dataset = None

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


def main():
  if len(argv) == 1:
    return
  if len(argv) > 2:
    dataset = argv[2]
  fol = argv[1]
  fils = list(filter(lambda x: x.endswith('.json'), os.listdir(fol)))
  ni = []
  betas = []
  # print(fils)
  for idx in range(1, 11):
    fil = '%d.json' % (idx)
    # print(fil)
    if fil in fils:
      with open('%s/%s' % (fol, fil), 'r') as fp:
        ni.append(json.load(fp))
      betas.append(idx)
  print(ni)
  plotlyGraphs(ni, betas, '%s NuclearIndex' % (dataset),
              '%s/graph.html' % (fol))


if __name__ == '__main__':
  main()