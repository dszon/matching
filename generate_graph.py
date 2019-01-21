#!/usr/local/bin/python3
import random
import sys
import argparse
import networkx as nx
import random
from collections import Counter
import numpy as np

# ----------- parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument("--nnodes",  "-n", action="store", dest="nnodes",    default=20,  type=int,     help="number of nodes")
parser.add_argument("--degree",  "-d", action="store", dest="degree",    default=5,  type=int,     help="mean degree")
parser.add_argument("--countries",  "-c", action="store", dest="countries",    default=1,  type=int,     help="number of countries")
args = parser.parse_args()

characters = '123456789abcdefghijkmnopqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ'
uppercharacters = 'ABCDEFGHJKLMNPQRSTUVWXYZ'


def randomName(l,alpha=False):
    if alpha == False:
        s = str(random.randrange(1,10))
        s = s + ''.join([str(random.randrange(0,10)) for i in range(l-1)])
    else:
        s = random.choice(uppercharacters)
        s = s + ''.join([random.choice(uppercharacters) for i in range(l-1)])
    return s


if args.countries == 1:
    G      = nx.random_regular_graph(args.degree, args.nnodes)
else:
    G      = nx.Graph()
    country_sizes = [random.random() for c in range(args.countries)]
    sc = sum(country_sizes)
    country_sizes = [s/sc for s in country_sizes]

    for i,cs in enumerate(country_sizes):
        print('country {:d} has size {:f}'.format(i,country_sizes[i]))

    edges = []
    nationality = np.random.choice(range(args.countries),size=args.nnodes,p=country_sizes)
    for n1 in range(args.nnodes):
        for n2 in range(args.nnodes):
            if n1 != n2 and nationality[n1] != nationality[n2] and random.random() > .7:
                edges.append((n1,n2))

    G.add_edges_from(edges)

if (1):
    l = int(np.ceil(np.log(args.nnodes)/np.log(26)))
    labels = { i:randomName(l,alpha=True) for i in range(args.nnodes) }
else:
    vornamenRaw = []
    with open('../../Daten_und_Analyse/Rohdaten/Vornamen/vornamen_maennlich.txt','r') as fh:
        for line in fh.readlines():
            vornamenRaw.append(line.strip().lower())
    with open('../../Daten_und_Analyse/Rohdaten/Vornamen/vornamen_weiblich.txt','r') as fh:
        for line in fh.readlines():
            vornamenRaw.append(line.strip().lower())
    vornamenRaw = list(set(vornamenRaw))

    vornamen = vornamenRaw
    factor = 2
    while len(vornamen) < args.nnodes:
        vornamen = []
        for f in range(factor):
            vornamen.extend([v + '_{:d}'.f for v in vornamenRaw])
        factor += 1

    random.shuffle(vornamen)
    labels = { i: vornamen[i] for i in range(args.nnodes) }


G = nx.relabel_nodes(G,labels)

# check_labels = [lab for lab in labels.items()]
# print(Counter(check_labels).most_common(3))

testnodes = set()
with open('graph.txt','w') as fh:
    for edge in G.edges():
        node1 = edge[0]
        node2 = edge[1]
        space1 = ' '
        space2 = ' '
        testnodes.add(node1)
        testnodes.add(node2)
        fh.write('{:s}{:s}{:s}{:s}{:f}\n'.format(node1,space1,node2,space2,random.random()))

print(G.number_of_nodes(),'nodes')
print(G.number_of_edges(),'edges')
print('mean degree',G.number_of_edges()*2/G.number_of_nodes())
# print(len(testnodes),'node names')
