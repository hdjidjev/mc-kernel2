import dwave
from dwave import embedding
import ast
#import itertools
import time
import networkx as nx
from dwave.cloud import Client
import json
#import random
import dimod
from datetime import datetime
from create_random_QUBOs import maximum_clique_qubo


def get_qubit_list(embedding):
	out = []
	for a in embedding:
		out += embedding[a]
	return out


def Start_DWave_connection(solverName):
	if solverName == 'Advantage_system5.1':
		region = 'eu-central-1'
	else:
		region = 'na-west-1'
	client = Client.from_config(region=region)
	DWave_solver = client.get_solver(solverName)
	A = DWave_solver.undirected_edges
	connectivity_graph = nx.Graph(list(A))
	return connectivity_graph, DWave_solver


'''def findCliqueDWfromFile(solverName = "Advantage_system4.1", cliqueSize=40, embeddingFile="results_disjoint_embeddings/Advantage_system4.1_clique_40.json",
	UTC = 1, N_DWAVE_CALLS = 1, AT = 20, programming_thermalization = 0, readout_thermalization = 0,
	num_reads = 100, sequential = True, parallel = False):
	# AT=annealing time, UTC=uniform_torque_compensation

	Target, solver = Start_DWave_connection(solverName)

	# D-Wave parameters
	params = {"num_reads": num_reads, "annealing_time": AT, "programming_thermalization": programming_thermalization, "readout_thermalization": readout_thermalization}

	#get embedding from file
	#file = open("embeddings/Advantage_system1.1/clique"+str(CLIQUE)+"_embedding.json", "r")
	file = open(embeddingFile, "r")
	embeddings = json.load(file)
	file.close()

	non_parallel_results = {}


	combined_QUBOs = {}
	bqm_tracking = {}

	### Sequential QA
	if sequential:
		#for problem in range(len(embeddings)):
		for problem in range(1):
			file = open("random_QUBOs/"+str(problem)+".txt", "r")
			QUBO = ast.literal_eval(file.read())
			file.close()
			qubit_list = get_qubit_list(embeddings[problem])
			embedding_dict_json = embeddings[problem]
			embedding_dict = {}
			for i in embedding_dict_json:
				embedding_dict[int(i)] = embedding_dict_json[i]
			physical_subgraph = Target.subgraph(qubit_list)
			bqm = dimod.BinaryQuadraticModel.from_qubo(QUBO)
			bqm_tracking[problem] = bqm
			chain_strength_fixed = embedding.chain_strength.uniform_torque_compensation(bqm, prefactor=UTC)#Uniform Torque Compensation chain strength calculation
			embedded_qubo = embedding.embed_qubo(QUBO, embedding_dict, physical_subgraph, chain_strength=chain_strength_fixed)
			combined_QUBOs = {**combined_QUBOs, **embedded_qubo}
			all_vectors = []
			all_QPU = 0
			all_energies = []
			for _ in range(N_DWAVE_CALLS):
				while (True):
					try:
						sampleset = solver.sample_qubo(embedded_qubo, answer_mode='raw', **params)
						vectors = sampleset.samples
						energies = sampleset.energies
						break
					except:
						print("fail1", flush=True)#Connection problem
						time.sleep(1)
						continue
				embedded = dimod.SampleSet.from_samples(vectors, dimod.BINARY, energies)
				samples_unembedded = dwave.embedding.unembed_sampleset(embedded, embedding_dict, bqm)
				vectors = samples_unembedded.record.sample.tolist()
				QPU_time = sampleset['timing']['qpu_access_time']/float(1000000)#Timing is given in microseconds
				all_QPU += QPU_time
				all_vectors += vectors
				all_energies += energies
			non_parallel_results[problem] = [all_QPU, all_vectors, all_energies]
		file_name = str(programming_thermalization)+"_"+str(readout_thermalization)+"_"+str(AT)+"_"+str(UTC)
		file = open("DWave_results/nonParallel_"+file_name+".json", "w")
		json.dump(non_parallel_results, file)
		file.close()

	### Parallel QA
	if parallel:
		all_vectors = []
		all_QPU = 0
		all_energies = []
		for rep_solve in range(N_DWAVE_CALLS):
			while (True):
				try:
					sampleset = solver.sample_qubo(combined_QUBOs, answer_mode='raw', **params)
					vectors = sampleset.samples
					energies = sampleset.energies
					break
				except:
					print("fail", flush=True)
					time.sleep(1)
					continue
			QPU_time = sampleset['timing']['qpu_access_time']/float(1000000)
			all_QPU += QPU_time
			all_vectors += vectors
			all_energies += energies
			#parallel_results["QPU_time"] = all_QPU
		file_name = str(programming_thermalization)+"_"+str(readout_thermalization)+"_"+str(AT)+"_"+str(UTC)
		file = open("DWave_results/parallel_"+file_name+".json", "w")
		json.dump([all_QPU, all_vectors, all_energies], file)
		file.close()'''


def findCliqueDW(graph, solverName="Advantage_system4.1", cliqueSize=40,
				 UTC=1, N_DWAVE_CALLS=1, AT=20, programming_thermalization=0, readout_thermalization=0,
				 num_reads=100, sequential=True, parallel=False):
	# AT=annealing time, UTC=uniform_torque_compensation

	Target, solver = Start_DWave_connection(solverName)
	embeddingFile = "results_disjoint_embeddings/" + \
		solverName + "_clique_" + str(cliqueSize) + ".json"

	# D-Wave parameters
	params = {"num_reads": num_reads, "annealing_time": AT, "programming_thermalization":
			  programming_thermalization, "readout_thermalization": readout_thermalization}

	# get embedding from file
	#file = open("embeddings/Advantage_system1.1/clique"+str(CLIQUE)+"_embedding.json", "r")
	file = open(embeddingFile, "r")
	embeddings = json.load(file)
	file.close()

	non_parallel_results = {}

	combined_QUBOs = {}
	bqm_tracking = {}

	# Sequential QA
	if sequential:
		# for problem in range(len(embeddings)):
		for problem in range(1):
			QUBO = maximum_clique_qubo(graph)
			bqm = dimod.BinaryQuadraticModel.from_qubo(QUBO)
			bqm_tracking[problem] = bqm
			qubit_list = get_qubit_list(embeddings[problem])
			embedding_dict_json = embeddings[problem]
			embedding_dict = {}
			for i in embedding_dict_json:
				embedding_dict[int(i)] = embedding_dict_json[i]
			physical_subgraph = Target.subgraph(qubit_list)
			chain_strength_fixed = embedding.chain_strength.uniform_torque_compensation(
				bqm, prefactor=UTC)  # Uniform Torque Compensation chain strength calculation
			embedded_qubo = embedding.embed_qubo(
				QUBO, embedding_dict, physical_subgraph, chain_strength=chain_strength_fixed)
			combined_QUBOs = {**combined_QUBOs, **embedded_qubo}
			all_vectors = []
			all_QPU = 0
			all_energies = []
			for _ in range(N_DWAVE_CALLS):
				while (True):
					try:
						sampleset = solver.sample_qubo(
							embedded_qubo, answer_mode='raw', **params)
						vectors = sampleset.samples
						energies = sampleset.energies
						break
					except:
						print("fail1", flush=True)  # Connection problem
						time.sleep(1)
						continue
				embedded = dimod.SampleSet.from_samples(
					vectors, dimod.BINARY, energies)
				samples_unembedded = dwave.embedding.unembed_sampleset(
					embedded, embedding_dict, bqm)
				vectors = samples_unembedded.record.sample.tolist()
				# Timing is given in microseconds
				QPU_time = sampleset['timing']['qpu_access_time'] / \
					float(1000000)
				all_QPU += QPU_time
				all_vectors += vectors
				all_energies += energies
			non_parallel_results[problem] = [
				all_QPU, all_vectors, all_energies]
		dt = datetime.now()
		ts = datetime.timestamp(dt)  # timestamp
		file_name = str(programming_thermalization)+"_" + \
			str(readout_thermalization)+"_"+str(AT)+"_"+str(UTC)+"_"+str(ts)
		file = open("DWave_results/nonParallel_"+file_name+".json", "w")
		json.dump(non_parallel_results, file)
		file.close()
		return non_parallel_results

	# Parallel QA
	if parallel:
		all_vectors = []
		all_QPU = 0
		all_energies = []
		for _ in range(N_DWAVE_CALLS):
			while (True):
				try:
					sampleset = solver.sample_qubo(
						combined_QUBOs, answer_mode='raw', **params)
					vectors = sampleset.samples
					energies = sampleset.energies
					break
				except:
					print("fail", flush=True)
					time.sleep(1)
					continue
			QPU_time = sampleset['timing']['qpu_access_time']/float(1000000)
			all_QPU += QPU_time
			all_vectors += vectors
			all_energies += energies
			#parallel_results["QPU_time"] = all_QPU
		file_name = str(programming_thermalization)+"_" + \
			str(readout_thermalization)+"_"+str(AT)+"_"+str(UTC)
		file = open("DWave_results/parallel_"+file_name+".json", "w")
		json.dump([all_QPU, all_vectors, all_energies], file)
		file.close()


if __name__ == '__main__':
	# findCliqueDW()
	print()
