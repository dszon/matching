#!/usr/local/bin/python3
import random
import sys
import argparse
import networkx as nx
import random

# ----------- parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument("--nnodes",  "-n", action="store", dest="nnodes",    default=20,  type=int,     help="number of nodes")
parser.add_argument("--degree",  "-d", action="store", dest="degree",    default=5,  type=int,     help="mean degree")
args = parser.parse_args()


testnodes = set()
with open('graph.txt','r') as fh:
    for line in fh.readlines():
        linea = line.strip().split(' ')
        testnodes.add(linea[0])
        testnodes.add(linea[1])
print(len(testnodes),'node names')
