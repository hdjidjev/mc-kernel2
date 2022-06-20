import ast
import json
from dwave import embedding
import networkx as nx
import time
from mc_kernel import Start_DWave_connection

#hardware_edges_filename = "Advantage_system1.1_edges.txt"

def rename_Nodes(cliqueSize, problems):
	G = nx.Graph()
	cl = nx.complete_graph(cliqueSize)
	for problem in range(problems):
		for e in list(cl.edges()):
			G.add_edge(str(e[0])+"_"+str(problem), str(e[1])+"_"+str(problem))
	return G
def deep_search(G, hardware_connectivity):
	for attempt in range(2):
		print(attempt)
		emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=1000, max_no_improvement=1000, chainlength_patience=1000, timeout=25600, threads=20)
		if emb != {}:
			return emb
	return {}
def iterative_search(hardware_connectivity, cliqueSize):
	problems = 0
	saved_embedding = {}
	while (True):
		problems += 1
		G = rename_Nodes(cliqueSize, problems)
		print(-4, end=',')
		emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=2, max_no_improvement=2, chainlength_patience=2, timeout=100, threads=1)
		print('start ', problems)
		if emb == {}:
			print(-3)
			emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=5, max_no_improvement=5, chainlength_patience=5, timeout=200, threads=2)
			if emb == {}:
				print(-2)
				emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=10, max_no_improvement=10, chainlength_patience=10, timeout=400, threads=2)
				if emb == {}:
					print(-1)
					emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=20, max_no_improvement=20, chainlength_patience=20, timeout=800, threads=2)
					if emb == {}:
						print(0)
						emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=50, max_no_improvement=50, chainlength_patience=50, timeout=1600, threads=10)
						if emb == {}:
							'''print(1)
							emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=100, max_no_improvement=100, chainlength_patience=100, timeout=3200, threads=10)
							if emb == {}:
								print(2)
								emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=200, max_no_improvement=200, chainlength_patience=200, timeout=6400, threads=20)
								if emb == {}:
									print(3)
									emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=500, max_no_improvement=500, chainlength_patience=500, timeout=12800, threads=20)
									if emb == {}:
										print(4)
										emb = embedding.minorminer.find_embedding(G, hardware_connectivity, tries=1000, max_no_improvement=1000, chainlength_patience=1000, timeout=25600, threads=20)
										if emb == {}:
											print(5)
											emb = deep_search(G, hardware_connectivity)
											if emb == {}:'''
							embeddings = extract_separate_embeddings(saved_embedding, problems-1)
							assert len(embeddings) == problems-1
							return embeddings
		saved_embedding = emb
		embeddings = extract_separate_embeddings(saved_embedding, problems)
		assert len(embeddings)==problems

def extract_separate_embeddings(saved_embedding, problems):
	result = {}
	for problem in range(problems):
		clique_embedding = {}
		for variable in list(saved_embedding.keys()):
			s = variable.split("_")
			if problem==int(s[1]):
				clique_embedding[int(s[0])] = saved_embedding[variable]
		result[problem] = clique_embedding
	return result

def getHardwareGraph(solverName, fromFile=True):
	file = open("Hardware_graphs/"+solverName+".json", "r")
	edgeList = json.load(file)
	file.close()
	hardwareGraph = nx.from_edgelist(edgeList)
	return hardwareGraph

def writeHardwareGraph(solverName):
	hardwareGraph, _ = Start_DWave_connection(solverName)
	file = open("Hardware_graphs/"+solverName+".json", "w")
	json.dump(list(hardwareGraph.edges), file)
	file.close()

def generateEmbedding(cliqueSize, solverName):
	hardwareGraph = getHardwareGraph(solverName)
	start = time.time()
	embeddings = iterative_search(hardwareGraph.copy(), cliqueSize)
	file = open("results_disjoint_embeddings/"+solverName+"_clique_"+str(cliqueSize)+".json", "w")
	json.dump(list(embeddings.values()), file)
	file.close()
	end = time.time()
	diff = end-start
	print('Time='+str(diff))
	print('embeddings count  = ', len(embeddings))


if __name__ == '__main__':
	saveHardwareGraph = False
	solverName = 'Advantage_system4.1'
	if saveHardwareGraph:
		writeHardwareGraph(solverName)
	else:
		for cliqueSize in [50]:
			print(cliqueSize)
			generateEmbedding(cliqueSize, solverName)
