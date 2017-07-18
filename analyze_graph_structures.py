#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import matplotlib.pyplot as plt

def compute_distribution(data, verbose=False):
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
	def __init__(self, k_values=-1, readlengths=-1, error_percentages=-1, number_of_reads=-1):
		self.k_values = k_values
		self.readlengths = readlengths
		self.error_percentages = error_percentages
		self.number_of_reads = number_of_reads
		
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
		casename = "k("
		if self.k_values == -1:
			casename += "-1"
		else:
			for i in range(len(self.k_values)-1):
				casename += str(self.k_values[i])+"-"
			casename += str(self.k_values[-1])
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
				casename += str(self.error_percentages[i])+"-"
			casename += str(self.error_percentages[-1])
		casename += ")"
		return casename

class GraphData:
	def __init__(self):
		self.dirname = ""
		self.casename = ""
		self.error_percentage = -1
		self.readlength = -1
		self.num_of_reads = -1
		self.k_value = -1
		self.nodes = []
		self.num_of_nodes = 0
		self.num_of_edges = 0
		self.node_sequence_lengths = []
		self.node_adjacencies = {}
		
	def add_node(self, node_name, node_seq):
		self.num_of_nodes += 1
		self.nodes.append(node_name)
		self.node_sequence_lengths.append(len(node_seq))
		self.node_adjacencies[node_name] = []
	
	def add_edge(self, node_s, node_t):
		self.num_of_edges += 1
		if node_s not in self.node_adjacencies:
			self.node_adjacencies[node_s] = [node_t]
		else:
			self.node_adjacencies[node_s].append(node_t)
			
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

class GraphAnalyzer:
	def __init__(self, sourcedir):
		self.sourcedir = sourcedir
		self.graphdatas = []
	
	def get_data(self, verbose=False, case_restrict=CaseRestriction()):
		datapath = "Output/"+sourcedir
		for file in os.listdir(datapath):
			filenameparts = re.split(r'\.', file)
			filename = filenameparts[0]
			if len(filenameparts) > 1 and filenameparts[-1] == "asqg":
				casedata = re.split(r'_', filename)
				k_value = int(casedata[-1])
				readlength = int(casedata[2])
				error_percentage = float(".".join(re.split(r'-',casedata[5])))
				num_of_reads_data = re.split(r'[\[,\]]', casedata[4])
				num_of_reads = int(num_of_reads_data[1])*(len(num_of_reads_data)-2)
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
					
	def lineplot(self, data, x_axis, style='-'):
		if not (data in ["num_of_nodes", "num_of_edges"] and x_axis in ["k_value", "readlength", "error_percentage"]):
			print ("Error! Wrong Specifier!")
		else:
			x_values = []
			y_values = []
			label_data = []
			for gd in self.graphdatas:
				this_label = ""
				if data == "num_of_nodes":
					y = gd.num_of_nodes
					this_label += "nodes_"
				elif data == "num_of_edges":
					y = gd.num_of_edges
					this_label += "edges_"
				
				if x_axis == "k_value":
					x = gd.k_value
					this_label += "rl="+str(gd.readlength)+"_ep="+str(gd.error_percentage)+"_nr="+str(gd.num_of_reads)
				elif x_axis == "readlength":
					x = gd.readlength
					this_label += "k="+str(gd.k_value)+"_ep="+str(gd.error_percentage)+"_nr="+str(gd.num_of_reads)
				elif x_axis == "error_percentage":
					x = gd.error_percentage
					this_label += "k="+str(gd.k_value)+"_rl="+str(gd.readlength)+"_nr="+str(gd.num_of_reads)
				
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
				plt.plot(x_sorted, y_sorted, label=label_data[i], linestyle=style)
				
			plt.xlabel(x_axis)
			plt.ylabel("total number")
			plt.legend()