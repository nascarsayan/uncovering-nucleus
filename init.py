import networkx as nx
import csv
import matplotlib.pyplot as plt
from sys import argv
import json
import os
import re

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
        if not line.startswith('#'):
          break
      line = txtfile.readline()
      delim = '[ \t]+'
      if ',' in line:
        delim = ','
      rows.append(list(map(lambda x: int(x.strip()), re.split(delim, line.strip()))))
      for line in txtfile:
        try:
          rows.append(list(map(lambda x: int(x.strip()), re.split(delim, line.strip()))))
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


def getDep(G, cores, outputDir):
  betas = list(map(lambda x: x / 10, list(range(11))))
  nodes = G.nodes()
  k_max = max(cores.values())
  for beta in betas:
    print('+++ Beta = %f' %(beta))
    depBeta = [{v: 0 for v in nodes}]
    for k in range(1, k_max + 1):
      print('+++ k = %f' %(k))
      depBetaK = {}
      for v in nodes:
        depBetaK[v] = depBeta[k - 1][v]
        nbrs = list(
            filter(lambda v: cores[v] == k, list(nx.all_neighbors(G, v))))
        depBetaK[v] += len(nbrs) + beta * sum(
            list(map(lambda u: depBeta[k - 1][u], nbrs)))
      depBeta.append(depBetaK)
    dumpData(depBeta, '%s/dep/%d.json' % (outputDir, beta*10))



def main():
  if not(os.path.exists(outputDir)):
    os.makedirs(outputDir)
  if not(os.path.exists('%s/dep' % outputDir)):
    os.makedirs('%s/dep' % outputDir)
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