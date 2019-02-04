import networkx as nx
import csv
import matplotlib.pyplot as plt
from sys import argv

rootdir = './twitterDataset'
if len(argv) > 1 and argv[1] == 'test':
  rootDir = './dummyDataset'

nfile = '%s/nodes.csv' % (rootdir)
efile = '%s/edges.csv' % (rootdir)


def readCsv(fname):
  rows = []
  with open(fname, 'r') as csvfile:
    csvr = csv.reader(csvfile)
    for row in csvr:
      rows.append(list(map(int, row)))
  return rows


def drawGraph(G, figname='graph.png'):
  nx.draw(G, with_labels=True, font_weight='bold')
  plt.savefig(figname)
  plt.clf()


def getDep(G, cores):
  betas = list(map(lambda x: x / 10, list(range(11))))
  dep = []
  for beta in betas:
    print('+++ Beta = %f' %(beta))
    nodes = G.nodes()
    k_max = max(cores.values())
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
    dep.append(depBeta)
  return dep


def main():
  print('Reading nodes...')
  nodes = list(map(lambda x: x[0], readCsv(nfile)))
  print('Done!\nReading edges...')
  edges = readCsv(efile)
  print('Done!\nAdding nodes...')
  G = nx.Graph()
  G.add_nodes_from(nodes)
  print('Done!\nAdding edges...')
  G.add_edges_from(edges)
  print('Done!\nCalculating core values...')
  G.remove_edges_from(nx.selfloop_edges(G))
  cores = nx.core_number(G)
  print('Done!\nGetting dependency values...',)  
  print(getDep(G, cores))


if __name__ == '__main__':
  main()