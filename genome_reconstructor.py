#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import gc

import fast_debruijn_graph_builder as fdgb
import data_io as dio

class Node:
	def __init__(self, id, name, sequence):
		self.id = id
		self.name = name
		self.sequence = sequence
		self.length = len(sequence)
		self.adjacencies = []
		self.reverse_adjacencies = []
		self.weight = -1
		self.label = False

class Graph:
	def __init__(self):
		self.number_of_nodes = 0
		self.number_of_edges = 0
		self.nodes = []
		self.node_name_dictionary = {}
		self.max_label = 0
		self.min_label = 0
		
	def add_node(self, node_name, node_seq):
		self.nodes.append(Node(self.number_of_nodes, node_name, node_seq))
		self.node_name_dictionary[node_name] = self.number_of_nodes
		self.number_of_nodes += 1
		
	def add_edge(self, node_s, node_t):
		if node_s in self.node_name_dictionary and node_t in self.node_name_dictionary:
			node_s_id = self.node_name_dictionary[node_s]
			node_t_id = self.node_name_dictionary[node_t]
			self.nodes[node_s_id].adjacencies.append(node_t_id)
			self.nodes[node_t_id].reverse_adjacencies.append(node_s_id)
			self.number_of_edges += 1
		
	def read_graph_from_asqg(self, filename):
		with open(filename, 'r') as datafile:
			for line in datafile:
				data = re.split(r'\s', line)
				if data[0] == "VT":
					self.add_node(data[1], data[2])
				if data[0] == "ED":
					self.add_edge(data[1], data[2])
						
	def construct_assembly_ordering_labels(self):
		# algorithm assumes that graph
		# 	is not empty and
		#	has only one component and
		# 	has no cycles, i.e. implies a partial order
		
		queue = [[self.nodes[0].id, 0]]
		while (len(queue) > 0):
			current_data = queue[0]
			queue.pop(0)
			current_node_id = current_data[0]
			#print current_node_id
			#print self.nodes[current_node_id].adjacencies
			for node_id in self.nodes[current_node_id].adjacencies:
				if self.nodes[node_id].label == False:
					start_label = current_data[1] + self.nodes[current_node_id].length
					if start_label > self.max_label:
						self.max_label = start_label
					self.nodes[node_id].label = start_label
					queue.append([node_id, start_label])
			for node_id in self.nodes[current_node_id].reverse_adjacencies:
				if self.nodes[node_id].label == False:
					start_label = current_data[1] - self.nodes[node_id].length
					if start_label < self.min_label:
						self.min_label = start_label
					self.nodes[node_id].label = start_label
					queue.append([node_id, start_label])
					
	def print_nodes_sorted_by_label(self):
		sorted_nodes = sorted(self.nodes, key=lambda x: x.label)
		print [str(n.id)+" ("+n.name+"): "+str(n.label) for n in sorted_nodes]

	def get_partition_of_sequences(self, number_of_parts, verbose=False):
		sorted_nodes = sorted(self.nodes, key=lambda x: x.label)
		label_div = self.max_label-self.min_label
		part_size = label_div/(number_of_parts+1)
		
		parts = []
		for i in range(number_of_parts):
			current_start = self.min_label+i*(part_size)
			if i == number_of_parts-1:
				current_end = self.max_label
			else:
				current_end = self.min_label+(i+2)*(part_size)
			parts.append([node for node in sorted_nodes if node.label >= current_start and node.label <= current_end])
			if verbose:
				print str(i)+": "+str(current_start)+" - "+str(current_end)+" : "+str(len(parts[-1]))+" sequences"
		return parts

	def write_node_sequences_to_file(self, filename, nodes):
		sequence_string = ""
		for n in nodes:
			sequence_string += n.sequence + "\n"

		outputfile = file(filename, 'w')
		outputfile.write(sequence_string)
		
g = Graph()
g.read_graph_from_asqg("Output/corona-largereads-asbk-1/corona-largereads-asbk-1-i_8000_1_[250]_1-0_30.asqg")
g.construct_assembly_ordering_labels()
#g.print_nodes_sorted_by_label()

nodesets = g.get_partition_of_sequences(100)

i = 0
for nodeset in nodesets:
	if i == 20:
		readfile_name = "reconstruct_sequence_reads_"+str(i)
		g.write_node_sequences_to_file(readfile_name, nodeset)

		reads = dio.get_reads_from_file(filename = readfile_name)
		debruijn = fdgb.GraphData(reads, 20)
		# delete reads and kmers to save ram:
		debruijn.reads = []
		debruijn.kmers = []
		# run garbage collector:
		gc.collect()
		debruijn.contract_unique_overlaps(verbose=False)
		debruijn.remove_parallel_sequences()
		
		debruijn.get_asqg_output(filename = "reconstruct_test_"+str(i)+".asqg")
		debruijn.get_csv_output(filename = "reconstruct_test_"+str(i)+".csv")
	i += 1
	
