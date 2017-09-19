#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import gc
import os

import fast_debruijn_graph_builder as fdgb
import data_io as dio
import manjasDefinitionen as md
import sampleReads

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
	# reconstructed graph from asqg-file
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
		
def test_with_read_from_asqg():
	for p in [5, 10, 15]:
		for k in [30,40]:
			for k2 in [13,15,17]:
				g = Graph()
				g.read_graph_from_asqg("Output/corona-largereads-asbk-1/corona-largereads-asbk-1-i_8000_1_[250]_"+str(p)+"-0_"+str(k)+".asqg")
				g.construct_assembly_ordering_labels()
				#g.print_nodes_sorted_by_label()
				testcasename = "test-reconstruct_p"+str(p)+"_k"+str(k)+"_k"+str(k2)+"_"
				
				nodesets = g.get_partition_of_sequences(100)
				
				i = 0
				for nodeset in nodesets:
					if i == 20:
						path="Output/test/"
						readfile_name = testcasename+"reads_"+str(i)
						g.write_node_sequences_to_file(path+readfile_name, nodeset)
				
						reads = dio.get_reads_from_file(filename = path+readfile_name)
						debruijn = fdgb.GraphData(reads, k2)
						# delete reads and kmers to save ram:
						debruijn.reads = []
						debruijn.kmers = []
						# run garbage collector:
						gc.collect()
						debruijn.contract_unique_overlaps(verbose=False)
						debruijn.remove_parallel_sequences()
						
						debruijn.get_asqg_output(filename = path+testcasename+str(i)+"_preremove"+".asqg")
						debruijn.get_csv_output(filename = path+testcasename+str(i)+"_preremove"+".csv")
				
						debruijn.remove_insignificant_sequences(verbose=False)
						
						debruijn.get_asqg_output(filename = path+testcasename+str(i)+"_postremove"+".asqg")
						debruijn.get_csv_output(filename = path+testcasename+str(i)+"_postremove"+".csv")
						
						debruijn.contract_unique_overlaps(verbose=False)
				
						debruijn.get_asqg_output(filename = path+testcasename+str(i)+"_postcontract"+".asqg")
						debruijn.get_csv_output(filename = path+testcasename+str(i)+"_postcontract"+".csv")
				
						debruijn.remove_tips()
						debruijn.contract_unique_overlaps(verbose=False)
						
						debruijn.get_asqg_output(filename = path+testcasename+str(i)+"_posttipremoval"+".asqg")
						debruijn.get_csv_output(filename = path+testcasename+str(i)+"_posttipremoval"+".csv")
					i += 1

def test_from_scratch():
	path = "Output/test/"
	output_dir = path
	set_of_viruses = [md.v1]
	readlength = 4000
	num_reads = 50
	num_different_viruses = 1
	epr = 0
	epi = 15
	
	k1 = 30
	k2 = 17
	
	number_of_parts = 30
	
	casename = "l"+str(readlength)+"_n"+str(num_reads)+"_e"+str(epi)
	
	readfilename = "test_read.txt"
	if not os.path.isfile(readfilename):
		sampleReads.samplereads(output_filename	= readfilename,	read_length = readlength, set_of_viruses = set_of_viruses[:num_different_viruses], number_of_reads = [num_reads], replace_error_percentage = epr, indel_error_percentage = epi)
	
	reads = dio.get_reads_from_file(filename = readfilename)
	
	debruijn = fdgb.GraphData(reads, k1)
	# delete reads and kmers to save ram:
	debruijn.reads = []
	debruijn.kmers = []
	# run garbage collector:
	gc.collect()
	debruijn.contract_unique_overlaps(verbose=False)
	debruijn.remove_parallel_sequences()
	debruijn.construct_assembly_ordering_labels()
	
	debruijn.get_asqg_output(filename = output_dir+"/"+casename+"_base.asqg")
	debruijn.get_csv_output(filename = output_dir+"/"+casename+"_base.csv")
	
	parts = debruijn.get_partition_of_sequences(number_of_parts, verbose=False)
	
	debruijn = []
	gc.collect()
	
	#for p in range(number_of_parts):
	p = 1
	
	part_seqs = [s.sequence for s in parts[p]]
	#part_reads = [reads[i] for i in parts[p]]
	debruijn_part = fdgb.GraphData(part_seqs, k2)
	
	postfix = "_p"+str(p)+"_preremove"
	
	debruijn_part.get_asqg_output(filename = output_dir+"/"+casename+postfix+".asqg")
	debruijn_part.get_csv_output(filename = output_dir+"/"+casename+postfix+".csv")
	
	debruijn_part.remove_insignificant_sequences(verbose=False)
	debruijn_part.remove_single_sequence_components
	
	postfix = "_p"+str(p)+"_postremove"
	debruijn_part.get_asqg_output(filename = output_dir+"/"+casename+postfix+".asqg")
	debruijn_part.get_csv_output(filename = output_dir+"/"+casename+postfix+".csv")
	
	debruijn_part.contract_unique_overlaps(verbose=False)
	
	postfix = "_p"+str(p)+"_postcontract"
	debruijn_part.get_asqg_output(filename = output_dir+"/"+casename+postfix+".asqg")
	debruijn_part.get_csv_output(filename = output_dir+"/"+casename+postfix+".csv")
	
	debruijn_part.remove_tips()
	debruijn_part.contract_unique_overlaps(verbose=False)
	
	postfix = "_p"+str(p)+"_posttipremoval"
	debruijn_part.get_asqg_output(filename = output_dir+"/"+casename+postfix+".asqg")
	debruijn_part.get_csv_output(filename = output_dir+"/"+casename+postfix+".csv")
	
	debruijn_part.write_sequences_to_file(filename = output_dir+"/"+casename+"reads_p"+str(p)+".txt")

test_from_scratch()
