import random
from time import time
import numpy as np
import networkx as nx
import itertools
import matplotlib.pyplot as plt


def get_chain_g(node_num):
    g = nx.Graph()
    for i in range(node_num):
        g.add_node(i, weight=0)
        if i < node_num - 1:
            g.add_edge(i, i + 1, weight=1)
    g.add_edge(node_num - 1, 0, weight=1)
    return g

def get_lattice_g(nr, nc):
    g = nx.Graph()
    for i in range(nr*nc):
        r = i // nc
        c = i % nc
        g.add_node(i, weight=0)
        if c < nc - 1:
            g.add_edge(i, i + 1, weight=1)
        else:
            g.add_edge(i, i - (nc - 1), weight=1)
        
        if i + nc < nr * nc:
            g.add_edge(i, i + nc, weight=1)
        else:
            g.add_edge(i, i - (nr - 1)*nc, weight=1)

    return g 


def get_hamiltonian(conf, g):
    Hs = 0
    # for nd in g.nodes:
    #     confu = 2*int(conf[nd]) - 1
    #     Hs -= g.nodes[nd]['weight'] * confu

    for ed in g.edges:
        confu, confv = 2*int(conf[ed[0]]) - 1, 2*int(conf[ed[1]]) - 1
        Hs += g.adj[ed[0]][ed[1]]['weight'] * confu * confv

    return Hs


st = time()

node_num = 9
betalist = [0.2, 1, 1.3, 1.9, 2.3, 5.4, 6.8, 9.1]
# betalist = [0.2, 1] # np.linspace(0, 20, 20)

# conf_set = [-1, 1]
perms = []
for ce in range(0, 2**node_num):
    key = format(ce, '0%db'%node_num)
    perms.append(key)


# test result of one conf
# it = []
# for i in range(2 ** node_num):
#     if (i & 1) == 0:
#         it.append(1)
#     else:
#         it.append(0)
# perms.append(it)
# print(perms)

# g = get_chain_g(3)
# g = get_complete_g(node_num)
g = get_lattice_g(3, 3)
# plt.figure()
# nx.draw(g, with_labels=True)
Ising_beta = []
for beta in betalist:
    Em, sp = 0, 0
    wx, mx = [], []
    for ep in perms:
        # print(ep)
        Hs = get_hamiltonian(ep, g)
        mx.append(np.sum([2*int(x) - 1 for x in ep]))
        tw = np.exp(-beta*Hs)
        wx.append(tw)
        sp += tw
    Em = np.sum(np.array(wx) / sp * np.array(mx))
    Ising_beta.append(Em)

ed = time()
print("time cost is: ", ed - st)
print("value of Ising_beta list is: ", Ising_beta)


plt.figure()
plt.plot(betalist, Ising_beta, "ob-")
plt.ylabel('expectation of magnetism')
plt.xlabel('value of beta')
plt.title("phase transition")
plt.legend()
plt.show()

# init_phi = []
# for i in range(node_num):
#     spin = random.choice([1, -1])
#     init_phi.append(spin)