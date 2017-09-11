#!/usr/bin/python
# -*- coding: utf-8 -*-

import re

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
					self.nodes[node_id].label = start_label
					queue.append([node_id, start_label])
			for node_id in self.nodes[current_node_id].reverse_adjacencies:
				if self.nodes[node_id].label == False:
					start_label = current_data[1] - self.nodes[node_id].length
					self.nodes[node_id].label = start_label
					queue.append([node_id, start_label])
					
	def print_nodes_sorted_by_label(self):
		sorted_nodes = sorted(self.nodes, key=lambda x: x.label)
		print [str(n.id)+" ("+n.name+"): "+str(n.label) for n in sorted_nodes]
		
g = Graph()
g.read_graph_from_asqg("Output/general_absk_1/bvdv_sample_250_1_[500]_0-5_16.asqg")
g.construct_assembly_ordering_labels()
g.print_nodes_sorted_by_label()