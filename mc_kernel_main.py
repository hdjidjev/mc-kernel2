import math
import json
import networkx as nx
import random
from mc_kernel import findCliqueDW

def selVarsToDrop(annealResults):
    bottomPercent = 0.1 # percent of worst variables to return
    '''file_name = 'nonParallel_0_0_20_1.json'
    tp, programming_thermalization, readout_thermalization, AT, UTC, ts = file_name.split('_')
    UTC,_ = UTC.split('.')
    file = open('DWave_results/'+file_name,"r")
    annealResults = json.load(file)
    file.close()'''

    res1 = annealResults[0]
    #qpuTime = res1[0]
    all_vectors = res1[1]
    all_energies = res1[2]
    #pairs = zip(all_energies,all_vectors)
    #all_energies = (all_energies)
    m = min(all_energies)
    all_energiesMinusMin = [x-m for x in all_energies]
    sums = [sum([all_vectors[i][j]+all_energiesMinusMin[i]/10.0/len(all_vectors[0]) for i in range(len(all_vectors))]) for j in range(len(all_vectors[0]))]
    pairs = sorted([(i, sums[i]) for i in range(len(sums))], key = lambda x: x[1])
    selSize = math.ceil(bottomPercent*len(pairs))
    selectedVars = [pairs[i][0] for i in range(selSize)]
    #print(selectedVars)
    '''for i in selectedVars:
        print([all_vectors[j][i] for j in range(len(all_vectors))])'''
    return selectedVars

def kCoreSub(G, n0, relabel=True): #construct a graph induced by n0 vertices with biggest k-core #s
    #print([G.degree(v) for v in G])
    core_numbers_dict = nx.core_number(G)
    #print(core_numbers_dict)
    core_numbers_list = [k for k, v in sorted(core_numbers_dict.items(), key=lambda item: (item[1]+2*random.random(),G.degree(item[0])+2*random.random()), reverse=True)]
    #print(core_numbers_list)
    g = G.subgraph(core_numbers_list[:n0])
    if relabel:
        d = dict(zip(g,range(len(g))))
        inverse_d = {d[k]:k for k in d}
        g1 = nx.relabel_nodes(g, d) # relable to 0,1,2,...
    return g1, inverse_d

'''def extend(cl, cliqueNodes, used, iter): # construct next subgraph
    count = {v:0 for v in G.nodes()} # number of neighbors in cl for each vertex
    for v in cl:
        for w in G[v]:
            if w not in used or random.random()<threshold:
                count[w] += 1 # the # of neighbors in sub
    for v in cliqueNodes:
        count[v] = N # prioritize these nodes
        #count[v] += 1 # prioritize these nodes
    if random.random()>threshold:
        select = [k for k, v in sorted(count.items(), key=lambda item: item[1]+2*random.random(), reverse=True)]
        select = select[:n0] # top n0 vertices wrt count
    else:
        select = random.choices(list(count.keys()), weights=list(count.values()), k=n0)
        print("Random selection of sub\n")
    sub = G.subgraph(select) # use k-core numbering that prioritizes cliqueNodes
    return sub'''

def MC_iterative():
    N=100
    p = 0.8
    n0 = 40 # max graph size for solving MC
    itNumber = 20
    num_reads = 200
    #threshold = 0 #was 0.7 # probability threshold for adding used node to clgraph
    G = nx.erdos_renyi_graph(N, p)
    opt=nx.max_weight_clique(G,weight=None) # best clique
    MCsize=opt[1] # best clique size
    sub, inverse_d = kCoreSub(G, n0)
    cl = nx.max_weight_clique(sub,weight=None) # first clique
    Csize = cl[1] # first clique size
    #cl = cl[0]  # first clique
    print("N: ", N, " p: ", p, "init:",len(cl), end=' ')
    annealResults = findCliqueDW(sub, num_reads = num_reads)
    #selVarsToDrop(annealResults)
    for _ in range(itNumber):
        selectedVars = selVarsToDrop(annealResults)
        for v in selectedVars:
            G.remove_node(inverse_d[v])
        if len(G) == 0: break
        sub, inverse_d = kCoreSub(G, n0)
        # todo replace next 4 lines with DW anneal
        cl2 = nx.max_weight_clique(sub, weight=None)
        print(cl2[1], end=' ')
        if cl2[1]>Csize:
            Csize = cl2[1]
        annealResults = findCliqueDW(sub, num_reads = num_reads, solverName="Advantage_system4.1")
        #cl=cl2[0]
    print(', found:', Csize, ', opt:', MCsize)



if __name__ == '__main__':
    MC_iterative()
    