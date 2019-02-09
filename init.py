import networkx as nx
import csv
import matplotlib.pyplot as plt
from sys import argv
import json
import os
import re
from collections import Counter

efile = './dummyDataset/edges.csv'
if len(argv) > 1:
  efile = argv[1]
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
  plt.savefig(filename)
  if clf:
    plt.clf()


def getDep(G, cores, outputDir):
  betas = list(map(lambda x: x / 10, list(range(11))))
  nodes = G.nodes()
  edges = G.edges()
  k_max = max(cores.values())
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
    # dumpData(depBeta, '%s/dep/%d.json' % (outputDir, beta * 10))
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
  lineGraphs(ni, betas, 'NuclearIndex', '%s/ni/graph.png' % (outputDir))
  cnt = Counter(aggni)
  kc = (cnt.most_common(1))[0][0]
  print('@@@ Aggregate Nuclear Index = %d, β chosen = %r'
        % (kc,
           list(
               map(lambda x: betas[x],
                   filter(lambda y: aggni[y] == kc, range(10))))))


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
  print('Done!\nCalculating core values...')
  G.remove_edges_from(nx.selfloop_edges(G))
  cores = nx.core_number(G)
  dumpData(cores, '%s/cores.json' % (outputDir))
  print('Done!\nGetting dependency values...',)
  getDep(G, cores, outputDir)


if __name__ == '__main__':
  main()