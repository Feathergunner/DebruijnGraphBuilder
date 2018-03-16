#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import matplotlib.pyplot as plt

def compute_distribution(data, verbose=False):
	if len(data) == 0:
		return [0]
	else:
		if verbose:
			print ("Max(data) is " + str(max(data)))
		dist = [0]*(max(data)+1)
		for d in data:
			dist[d] += 1
		return dist
	
class CaseRestriction:
	# the parameters have to be either -1 or a list of values
	# a specific set of values is compatible with this restriction,
	# if the restriction is either -1 or the values is included in the list
	def __init__(self, k_values=-1, readlengths=-1, error_percentages=-1, number_of_reads=-1, num_genomes=1):
		self.k_values = k_values
		self.readlengths = readlengths
		self.error_percentages = error_percentages
		self.number_of_reads = number_of_reads
		self.number_of_genomes = num_genomes
		
	def check_compatibility(self, k, readlength, error_percentage, number_of_reads):
		is_compatible = True
		if not (self.k_values == -1 or k in self.k_values):
			is_compatible = False
		if not (self.readlengths == -1 or readlength in self.readlengths):
			is_compatible = False
		if not (self.error_percentages == -1 or error_percentage in self.error_percentages):
			is_compatible = False
		if not (self.number_of_reads == -1 or number_of_reads in self.number_of_reads):
			is_compatible = False
		return is_compatible
		
	def construct_casename(self):
		casename = "ng("+str(self.number_of_genomes)+")_k("
		if self.k_values == -1:
			casename += "-1"
		else:
			for i in range(len(self.k_values)-1):
				casename += str(int(self.k_values[i]))+"-"
			casename += str(int(self.k_values[-1]))
		casename += ")_rl("
		if self.readlengths == -1:
			casename += "-1)_nr("
		else:
			for i in range(len(self.readlengths)-1):
				casename += str(self.readlengths[i])+"-"
			casename += str(self.readlengths[-1])
		casename += ")_nr("
		if self.number_of_reads == -1:
			casename += "-1"
		else:
			for i in range(len(self.number_of_reads)-1):
				casename += str(self.number_of_reads[i])+"-"
			casename += str(self.number_of_reads[-1])
		casename += ")_ep("
		if self.error_percentages == -1:
			casename += "-1"
		else:
			for i in range(len(self.error_percentages)-1):
				#ep_string = "".join(re.split(r'\.',str(self.error_percentages[i])))
				ep_string = "".join(re.split(r'-',str(self.error_percentages[i])))
				casename += ep_string+"-"
			#casename += "".join(re.split(r'\.',str(self.error_percentages[-1])))
			casename += "".join(re.split(r'-',str(self.error_percentages[-1])))
		casename += ")"
		return casename

class GraphData:
	def __init__(self,
				 dirname = "",
				 casename = "",
				 error_percentage = -1,
				 readlength = -1,
				 num_of_reads = -1,
				 k_value = -1,
				 nodes = []):
				 
		self.dirname = dirname
		self.casename = casename
		self.error_percentage = error_percentage
		self.readlength = readlength
		self.num_of_reads = num_of_reads
		self.k_value = k_value
		self.nodes = nodes
		self.num_of_nodes = len(self.nodes)
		self.num_of_edges = 0
		self.node_sequence_lengths = []
		self.node_adjacencies = {}
		self.node_reverse_adjacencies = {}
		self.node_name_dictionary = {}
		self.coverage_depths = []
		self.number_of_components = 0
		self.component_sizes = []
		
	def get_data_from_file(self, filename):
		with open(filename, 'r') as datafile:
			for line in datafile:
				data = re.split(r'\s', line)
				if data[0] == "VT":
					self.add_node(data[1], data[2])
				if data[0] == "ED":
					self.add_edge(data[1], data[2])
					
	def get_data_from_csv(self, csvfilename):
		with open(csvfilename, 'r') as datafile:
			for line in datafile:
				data = re.split(r',', line)
				if not data[0] == "Node_Label":
					self.coverage_depths.append(int(data[2]))		
		
	def add_node(self, node_name, node_seq):
		self.nodes.append(node_seq)
		self.node_sequence_lengths.append(len(node_seq))
		self.node_adjacencies[self.num_of_nodes] = []#[node_name] = []
		self.node_reverse_adjacencies[self.num_of_nodes] = []
		self.node_name_dictionary[node_name] = self.num_of_nodes
		self.num_of_nodes += 1
		# reset number of components:
		self.number_of_components = 0
	
	def add_edge(self, node_s, node_t):
		if node_s in self.node_name_dictionary and node_t in self.node_name_dictionary:
			# check is necessary, because a bug (?) left an edge from a self-inverse-sequence somewhere
			self.num_of_edges += 1
			node_s_id = self.node_name_dictionary[node_s]
			node_t_id = self.node_name_dictionary[node_t]
			
			self.node_adjacencies[node_s_id].append(node_t_id)
			self.node_reverse_adjacencies[node_t_id].append(node_s_id)
		# reset number of components:
		self.number_of_components = 0
			
	def get_degree_distribution(self):
		node_degrees = {}
		for s in self.node_adjacencies:
			for t in self.node_adjacencies[s]:
				if s in node_degrees:
					node_degrees[s] += 1
				else:
					node_degrees[s] = 1
				if t in node_degrees:
					node_degrees[t] += 1
				else:
					node_degrees[t] = 1
			
		node_degree_list = []
		for node in node_degrees:
			node_degree_list.append(node_degrees[node])
			
		return node_degree_list
		
	def get_avg_seq_length(self):
		return float(sum(self.node_sequence_lengths))/float(self.num_of_nodes)
		
	def get_avg_coverage_depth(self):
		return float(sum(self.coverage_depths))/float(self.num_of_nodes)
		
	def get_avg_component_size(self):
		num_comps = self.get_number_of_components()
		return float(sum(self.component_sizes))/float(num_comps)
		
	def get_maximum_component_size(self):
		num_comps = self.get_number_of_components()
		return max(self.component_sizes)
		
	def get_number_of_components(self):
		if self.number_of_components > 0:
			return self.number_of_components
		
		self.number_of_components = 0
		queue = []
		visited_nodes = [False]*self.num_of_nodes
		
		all_nodes_visited = False
		while not all_nodes_visited:
			node_index = 0
			while node_index < self.num_of_nodes and visited_nodes[node_index]:
				node_index += 1
			if node_index == self.num_of_nodes:
				all_nodes_visited = True
			else:
				self.number_of_components += 1
				queue.append(node_index)
				visited_nodes[node_index] = True
				current_component_size = 0
				while len(queue) > 0:
					current_node = queue.pop(0)
					current_component_size += 1
					for adj_node_id in self.node_adjacencies[current_node]:
						if not visited_nodes[adj_node_id]:
							queue.append(adj_node_id)
							visited_nodes[adj_node_id] = True
					for adj_node_id in self.node_reverse_adjacencies[current_node]:
						if not visited_nodes[adj_node_id]:
							queue.append(adj_node_id)
							visited_nodes[adj_node_id] = True
				self.component_sizes.append(current_component_size)
		return self.number_of_components

class GraphAnalyzer:
	def __init__(self, sourcedir):
		self.sourcedir = sourcedir
		self.graphdatas = []
	
	def get_data(self, verbose=False, case_restrict=CaseRestriction()):
		datapath = "Output/"+self.sourcedir
		
		for file in os.listdir(datapath):
			filenameparts = re.split(r'\.', file)
			filename = filenameparts[0]
			if len(filenameparts) > 1 and filenameparts[-1] == "asqg":
				casedata = re.split(r'_', filename)
				k_value = int(casedata[-1])
				# DIFFERENT VERSIONS OF FILENAMES:
				if len(casedata) == 6:
					readlength = int(casedata[1])
					error_percentage = float(".".join(re.split(r'-',casedata[4])))
					num_of_reads_data = re.split(r'[\[,\]]', casedata[3])
					num_of_reads = int(num_of_reads_data[1])#*(len(num_of_reads_data)-2)
				else:
					readlength = int(casedata[2])
					error_percentage = float(".".join(re.split(r'-',casedata[5])))
					num_of_reads_data = re.split(r'[\[,\]]', casedata[4])
					num_of_reads = int(num_of_reads_data[1])#*(len(num_of_reads_data)-2)
				if case_restrict.check_compatibility(k_value, readlength, error_percentage, num_of_reads):
					if verbose:
						print ("Open file "+datapath+"/"+file)
					with open(datapath+"/"+file, 'r') as datafile:
						gd = GraphData()
						gd.k_value = k_value
						gd.readlength = readlength
						gd.num_of_reads = num_of_reads
						gd.error_percentage = error_percentage
						
						gd.casename = filename			
						gd.dirname = datapath
						
						for line in datafile:
							data = re.split(r'\s', line)
							if data[0] == "VT":
								gd.add_node(data[1], data[2])
							if data[0] == "ED":
								gd.add_edge(data[1], data[2])
						
						self.graphdatas.append(gd)
					
	def get_data_from_specific_filenames(self, filenames, num_of_reads, readlengths, error_percentages, k_values, verbose=False):
		# the lists:
		#	filenames
		#	num_of_reads
		#	readlengths
		#	error_percentages
		#	k_values
		# have to contain the respective data in the order of the files defined by the list filenames
		datapath = "Output/"+self.sourcedir
		
		for i in range(len(filenames)):
			if verbose:
				print ("Open file "+datapath+"/"+filenames[i])
			with open(datapath+"/"+filenames[i]+".asqg", 'r') as datafile:
				gd = GraphData()
				gd.k_value = k_values[i]
				gd.readlength = readlengths[i]
				gd.num_of_reads = num_of_reads[i]
				gd.error_percentage = error_percentages[i]
				
				gd.casename = filenames[i]	
				gd.dirname = datapath
				
				for line in datafile:
					data = re.split(r'\s', line)
					if data[0] == "VT":
						gd.add_node(data[1], data[2])
					if data[0] == "ED":
						gd.add_edge(data[1], data[2])
				
				self.graphdatas.append(gd)
					
	def lineplot(self, data, x_axis, style='-', axis=0, legend_pos=0, verbose=False):
		if not (data in ["num_of_nodes", "num_of_edges", "num_of_components", "avg_seq_lengths"] and x_axis in ["k_value", "rel_k_value", "readlength", "error_percentage"]):
			print ("Error! Wrong Specifier!")
		else:
			x_values = []
			y_values = []
			label_data = []
			y_label = ""
			for gd in self.graphdatas:
				this_label = ""
				y_label = ""
				if data == "num_of_nodes":
					y = gd.num_of_nodes
					this_label += "nodes_"
					y_label += "nodes/edges"
				elif data == "num_of_edges":
					y = gd.num_of_edges
					this_label += "edges_"
					y_label += "nodes/edges"
				elif data == "num_of_components":
					y = gd.get_number_of_components()
					this_label += "components_"
					y_label += "components"
				elif data == "avg_seq_lengths":
					y = gd.get_avg_seq_length()
					this_label += "seqlengths_"
					y_label += "average sequence length"
				
				if x_axis == "k_value":
					x = gd.k_value
					this_label += "rl="+str(gd.readlength)+"_ep="+str(gd.error_percentage)+"_nr="+str(gd.num_of_reads)
				elif x_axis == "readlength":
					x = gd.readlength
					this_label += "k="+str(gd.k_value)+"_ep="+str(gd.error_percentage)+"_nr="+str(gd.num_of_reads)
				elif x_axis == "error_percentage":
					x = gd.error_percentage
					this_label += "k="+str(gd.k_value)+"_rl="+str(gd.readlength)+"_nr="+str(gd.num_of_reads)
				elif x_axis == "rel_k_value":
					x = float(gd.k_value)/gd.readlength
					this_label += "rl="+str(gd.readlength)+"_ep="+str(gd.error_percentage)+"_nr="+str(gd.num_of_reads)
				
				if verbose:
					print ("Add tuple ("+str(x)+","+str(y)+") to data")
				
				if not this_label in label_data:
					label_data.append(this_label)
					x_values.append([])
					y_values.append([])
					
				i = 0
				while not this_label == label_data[i]:
					i += 1
				x_values[i].append(x)
				y_values[i].append(y)
				
			
			for i in range(len(x_values)):
				x_sorted = sorted(x_values[i])
				y_sorted = [y for (x,y) in sorted(zip(x_values[i], y_values[i]))]
				if axis == 0:
					plt.plot(x_sorted, y_sorted, label=label_data[i], linestyle=style, linewidth=3)
				else:
					axis.plot(x_sorted, y_sorted, label=label_data[i], linestyle=style,  linewidth=3)
			
			if axis == 0:
				plt.xlabel(x_axis)
				plt.ylabel(y_label)
				plt.legend(loc=legend_pos)
			else:
				axis.set_xlabel(x_axis)
				axis.set_ylabel(y_label)
				axis.legend(loc=legend_pos)

	def distribution_lineplot(self, data, style='-', axis=0, legend_pos=0, verbose=False):
		if not (data in ["seq_length", "degree_dist"]):
			print ("Error! Wrong Specifier!")
		else:
			x_values = []
			y_values = []
			label_data = []
			for gd in self.graphdatas:
				this_label = ""
				y_label = "Frequency"
				x_label = ""
				if data == "seq_length":
					y = compute_distribution(gd.node_sequence_lengths)
					this_label += "sequence-lengths_"
					x_label += "Sequence-length"
				elif data == "degree_dist":
					y = compute_distribution(gd.get_degree_distribution())
					this_label += "degree_distribution"
					x_label += "Node-degree"
					
				this_label += "rl="+str(gd.readlength)+"_ep="+str(gd.error_percentage)+"_nr="+str(gd.num_of_reads)+"_k="+str(gd.k_value)
				
				if not this_label in label_data:
					label_data.append(this_label)
					x_values.append(range(len(y)))
					y_values.append(y)
			
			for i in range(len(x_values)):
				if axis == 0:
					plt.plot(x_values[i], y_values[i], label=label_data[i], linestyle=style, linewidth=3)
				else:
					axis.plot(x_values[i], y_values[i], label=label_data[i], linestyle=style,  linewidth=3)
			
			if axis == 0:
				plt.xlabel(x_label)
				plt.ylabel(y_label)
				plt.legend(loc=legend_pos)
			else:
				axis.set_xlabel(x_label)
				axis.set_ylabel(y_label)
				axis.legend(loc=legend_pos)