#!/usr/bin/python
# -*- coding: utf-8 -*-

import random
import networkx as nx
import matplotlib.pyplot as plt

def merge_intervals(list_of_intervals, k_value=1):
	num_of_intervals = len(list_of_intervals)
	sorted_intervals = sorted(list_of_intervals, key=lambda x: x[0])
	merged_intervals = []
	
	current_interval = 0
	current_interval_start = -1
	current_interval_end = -1
	while current_interval < num_of_intervals:
		if current_interval_start < 0:
			current_interval_start = sorted_intervals[current_interval][0]
			current_interval_end = sorted_intervals[current_interval][1]
		else:
			if sorted_intervals[current_interval][0] <= current_interval_end-k_value+1:
				current_interval_end = max(current_interval_end, sorted_intervals[current_interval][1])
			else:
				merged_intervals.append([current_interval_start, current_interval_end])
				current_interval_start = sorted_intervals[current_interval][0]
				current_interval_end = sorted_intervals[current_interval][1]

		current_interval += 1
		
	merged_intervals.append([current_interval_start, current_interval_end])
	return merged_intervals

def export_to_asqg(V, E, k_value, filename = "random_debruijn_graph.asqg"):
	print ("Writing asqg-file ...")
	headline = "HT\t\n"
	asqg_vertices = ""
	vertex_count = 0
	for v in V:
		sequence = "".join(["A"]*(V[v][1]-V[v][0]+1))
		#print V[v]
		#print sequence
		asqg_vertices += "VT\t"+v+"\t"+sequence+"\n"
		vertex_count += 1
	
	asqg_edges = ""
	for e in E:
		seq_1_length = V[e[0]][1]-V[e[0]][0]
		seq_2_length = V[e[1]][1]-V[e[1]][0]
		asqg_edges += "ED\t"+e[0]+" "+e[1]+" "+str(seq_1_length-k_value+1)+" "+str(seq_1_length-1)+" "+str(seq_1_length)+" 0 "+str(k_value-2)+" "+str(seq_2_length)+" 0 0 "+str(1)+"\n"
    
	print filename
	outputfile = file(filename, 'w')
	outputfile.write(headline)
	outputfile.write(asqg_vertices)
	outputfile.write(asqg_edges)
	
def create_random_debruijn_graph(genome_length = 10000, read_number = 1000, read_length = 1000, error_percentage = 10.0, k_value = 30):	
	V_debruijn = {}
	V_connect = {}
	
	correct_sequences = []
	
	for n_i in range(read_number):
	
		s = random.choice(range(genome_length - read_length))
		t = s + read_length
		#print (str(s)+" - " + str(t))
		
		number_of_currently_correct_bases = 0
		last_node = -1
		next_node_start = -1
		
		for r_i in range(read_length):
			if random.random() < error_percentage/100:
				#print ("Error in read " + str(n_i) + " at position " + str(s+r_i))
				if next_node_start < 0:
					next_node_start = s
				if number_of_currently_correct_bases >= k_value:
					# case: error after correct intervall
					correct_sequences.append([next_node_start, s+r_i-1])
					next_node_start = s + r_i
				if read_length - r_i < k_value:
					# case: tip at end of read
					#print ("Created new end-tip")
					V_debruijn[str(n_i)+"_"+str(next_node_start)] = [next_node_start, t]
					if next_node_start not in V_connect:
						V_connect[next_node_start] = {"in": [], "out": []}
					V_connect[next_node_start]["out"].append(str(n_i)+"_"+str(next_node_start))
					r_i = read_length
				number_of_currently_correct_bases = 0
				
			else:
				number_of_currently_correct_bases += 1
				if number_of_currently_correct_bases == k_value and next_node_start > 0:
					#print ("Created new front-tip or bubble")
					this_node_start = next_node_start
					next_node_start = s + r_i-k_value + 1
					V_debruijn[str(n_i)+"_"+str(this_node_start)] = [this_node_start, s+r_i-k_value]
					if not this_node_start == s:
						if this_node_start not in V_connect:
							V_connect[this_node_start] = {"in": [], "out": []}
						V_connect[this_node_start]["out"].append(str(n_i)+"_"+str(this_node_start))
					
					if next_node_start not in V_connect:
						V_connect[next_node_start] = {"in": [], "out": []}
					V_connect[next_node_start]["in"].append(str(n_i)+"_"+str(this_node_start))
					
		if next_node_start < 0:
			next_node_start = s
			#print ("Created new error-free read")
		
		if number_of_currently_correct_bases >= k_value:
			# add rest of read:
			#print ("Create Sequence of rest of read")
			correct_sequences.append([next_node_start, t])
			
	#print ("Sequences pre-merge:")
	#print (correct_sequences)
	#print V_connect
	
	sorted_correct_sequences = merge_intervals(correct_sequences, k_value)
	sorted_connections = sorted([vc for vc in V_connect])
	#print ("Pre debruijn-construction:")
	#print sorted_correct_sequences
	#print sorted_connections
	#print V_debruijn
	
	for correct_seq in sorted_correct_sequences:
		current_s = correct_seq[0]
		current_t = correct_seq[1]
		for vc in sorted_connections:
			if vc-1 in range(current_s, current_t+1):
				node_name = "c_"+str(current_s)
				# create debruijn-node:
				V_debruijn[node_name] = [current_s, vc-1]
				# add connections:
				if current_s in V_connect:
					V_connect[current_s]["out"].append(node_name)
				V_connect[vc]["in"].append(node_name)
				
				current_s = vc
		if current_s < current_t:
			node_name = "c_"+str(current_s)
			# create debruijn-node:
			V_debruijn[node_name] = [current_s, current_t]
			# add connections:
			if current_s in V_connect:
				V_connect[current_s]["out"].append(node_name)
			if current_t in V_connect:
				V_connect[current_t]["in"].append(node_name)
			
	
	sorted_debruijn_nodes = sorted([[name, V_debruijn[name]] for name in V_debruijn], key=lambda x: x[1][0])
	sorted_connections = sorted([vc for vc in V_connect])
	#print ("Post debruijn-construction:")
	#print sorted_debruijn_nodes
	#print sorted_correct_sequences
	#print V_connect
		
	# Construct edges:
	E_debruijn = []
	for vc_id in V_connect:
		vc = V_connect[vc_id]
		for e_s in vc["in"]:
			for e_t in vc["out"]:
				E_debruijn.append((e_s, e_t))
				
	#print ("Post edge-construction:")
	#print E_debruijn
	
	return [V_debruijn, E_debruijn]
		
[V, E] = create_random_debruijn_graph(genome_length = 30000, read_number = 1000, read_length = 10000, error_percentage = 15.0, k_value = 30)

export_to_asqg(V, E, 1)

'''
print V
print E

G = nx.DiGraph()
G.add_nodes_from(V)
G.add_edges_from(E)

nx.draw(G)
plt.show()
'''