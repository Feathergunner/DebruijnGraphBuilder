#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import gc
import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from sklearn import metrics

import data_io as dio
import fast_debruijn_graph_builder as fdgb

def construct_network_graph(filename_input):
	nodes = []
	edges = []

	with open(filename_input) as inputfile:
		lines = inputfile.readlines()
	for l in lines:
		data = re.split(r'\t', l)
		if data[0] == "VT":
			# case sequence data:
			nodes.append(int(re.split(r'_',data[1])[1]))			
		elif data[0] == "ED":
			# case overlap data:
			overlap_data = re.split(r'\s',data[1])
			source_seq_id = int(re.split(r'_',overlap_data[0])[1])
			target_seq_id = int(re.split(r'_',overlap_data[1])[1])
			if source_seq_id > 0 and target_seq_id > 0:
				edges.append((source_seq_id, target_seq_id))
	G = nx.Graph()
	G.add_nodes_from(nodes)
	G.add_edges_from(edges)
	
	nx.draw(G)
	plt.show()

	adj_mat = nx.to_numpy_matrix(G)
	
	sc = SpectralClustering(2, affinity='precomputed', n_init=100)
	sc.fit(adj_mat)
	
	print('spectral clustering')
	print(sc.labels_)

def get_sequences_by_params(filename_input, filename_output, min_weight=1, number_of_parts=1, overlap=1, partition_method="equidist", verbose=False):
	with open(filename_input) as inputfile:
		lines = inputfile.readlines()
	i = 0
	num_of_relevant_sequences = 0
	n = len(lines)
	min_label = False
	max_label = False
	
	filename_output += "_w"+str(min_weight)
	
	part_sequence_filenames = []
	
	sequences = []
	for l in lines:
		if i > 0:
			if i%10000 == 0:
				print str(i)+"/"+str(n)+" - "+str(num_of_relevant_sequences)
			data = re.split(r',',l)
			if verbose:
				print data
			if int(data[2]) > min_weight:
				label = int(data[3])
				sequences.append([data[1], label, data[2]])
				
				if not min_label or label < min_label:
					min_label = label
				if not max_label or label > max_label:
					max_label = label
				num_of_relevant_sequences += 1
		i+=1
	
	if number_of_parts > 1 and partition_method in ["equidist", "equisize"]:
		overlap = min(overlap, number_of_parts-1)
		if partition_method == "equidist":
			partition_size = (max_label-min_label)*overlap/number_of_parts
			partition_start_diff = (max_label-min_label-partition_size)/number_of_parts
			
			for i in range(number_of_parts):
				partname = filename_output+"_p"+str(i)+".txt"
				part_sequence_filenames.append(partname)
				part_seqs = []
				for s in sequences:
					if s[1] > i*partition_start_diff and s[1] < i*partition_start_diff+partition_size:
						part_seqs.append(s[0]+","+str(s[2])+"\n")
				outputfile = file(partname, 'w')
				for seq in part_seqs:
					outputfile.write(seq)
					
		elif partition_method == "equisize":
			partition_size = num_of_relevant_sequences*overlap/(number_of_parts-1)
			partition_start_diff = num_of_relevant_sequences/(number_of_parts-1)
			seq_sorted_by_label = sorted(sequences, key=lambda x :x[1])
			part_seqs = [[] for i in range(number_of_parts)]
			number_of_seq = 0
			#print ("number_of_parts: "+str(number_of_parts))
			for s in seq_sorted_by_label:
				first_part = number_of_seq / partition_start_diff
				for p_i in range(first_part, first_part+overlap):
					if p_i < number_of_parts:
						part_seqs[p_i].append(s[0]+","+str(s[2])+"\n")
				number_of_seq += 1
			for i in range(number_of_parts):
				partname = filename_output+"_p"+str(i)+".txt"
				part_sequence_filenames.append(partname)
				outputfile = file(partname, 'w')
				for seq in part_seqs[i]:
					outputfile.write(seq)
			
	else:
		if partition_method not in ["equidist", "equisize"]:
			print ("Error! No correct method for partitioning specified. No partitioning was applied.")
		part_seqs = []
		for s in sequences:
			part_seqs.append(s[0]+","+str(s[2])+"\n")
		part_sequence_filenames.append(filename_output+"_p0.txt")
		outputfile = file(filename_output+"_p0.txt", 'w')
		for seq in part_seqs:
			outputfile.write(seq)
			
	return part_sequence_filenames
	
def reconstruct_part(reads, filename_output, k, minweight=1, minlengthfactor=1, allow_recursion=False, saveall=False):
	print ("Reconstruct part "+filename_output+" ...")
	
	reduce_to_single_path = False
	
	debruijn = fdgb.GraphData(reads, k, load_weights=False, verbose=False)
	# delete reads and kmers to save ram:
	reads = []
	debruijn.reads = []
	# debruijn.kmers = []
	# run garbage collector:
	gc.collect()
	debruijn.remove_parallel_sequences()
	debruijn.contract_unique_overlaps(verbose=False)
	
	#debruijn.remove_insignificant_sequences(minimal_weight=minweight, verbose=False)
	#debruijn.remove_tips()
	#debruijn.contract_unique_overlaps(verbose = False)
	#debruijn.remove_single_sequence_components()
	#debruijn.reduce_to_single_largest_component()
	#debruijn.construct_assembly_ordering_labels()
	
	if reduce_to_single_path:
		debruijn.reduce_to_single_largest_component()
		debruijn.construct_assembly_ordering_labels(verbose = True)
		debruijn.reduce_to_single_path_max_weight(verbose = False)
		debruijn.contract_unique_overlaps(verbose = False)
		#debruijn.remove_short_sequences(length_bound_by_multiple_of_k=minlengthfactor)
		debruijn.construct_assembly_ordering_labels()
	
	if saveall:
		debruijn.get_asqg_output(filename = filename_output+".asqg")
		debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)
	
	return debruijn.get_label_span()
		
def reconstruct_parts(list_of_inputfiles, filename_output_base, k, minweight, minlengthfactor, allow_recursion=False, saveall=False, delete_parts=True):
	files_of_parts = []

	partcounter = 0
	totalparts = len(list_of_inputfiles)
	for filename in list_of_inputfiles:
	    #if filename == list_of_inputfiles[0]:
		reads = dio.get_reads_from_file(filename)
		if delete_parts:
			os.remove(filename)
		if len(reads)>0:
			filename_output = filename_output_base+"_k"+str(k)+"_p"+str(partcounter)
			partcounter += 1
			graph_span = reconstruct_part(reads, filename_output, k, minweight, minlengthfactor, allow_recursion, saveall=saveall)
			if (allow_recursion):# and graph_span > 100*k):
				part_sequence_filenames = get_sequences_by_params(filename_input=filename_output+".csv", filename_output=filename_output+"_w"+str(minweight), min_weight=minweight, number_of_parts=100, overlap=10, verbose=True)
				part_parts = reconstruct_parts_lv2(part_sequence_filenames, filename_output_base=filename_output+"_w"+str(minweight), k=15)
				
				merged_parts_seqfilename = reconstruct_merge(filename_output_base=filename_output_base+"_w"+str(minweight)+"_k"+str(17), files_to_merge=part_parts, merge_k=11, saveall=saveall)
				filename_output = merged_parts_seqfilename
			files_of_parts.append(filename_output+"_seqsonly.txt")
			
	return files_of_parts
	
def reconstruct_parts_lv2(list_of_inputfiles, filename_output_base, number_of_parts, k, saveall=False):
	files_of_parts = []
	for i in range(number_of_parts):
		filename_input = filename_input_base+"_p"+str(i)+".txt"
		filename_output = filename_output_base+"_k"+str(k)+"_p"+str(i)
		
		reads = dio.get_reads_from_file(filename_input)
		if len(reads)>0:
			
			debruijn = fdgb.GraphData(reads, k, verbose=False)
			# delete reads and kmers to save ram:
			reads = []
			debruijn.reads = []
			debruijn.kmers = []
			# run garbage collector:
			gc.collect()
			debruijn.contract_unique_overlaps(verbose=False)
			debruijn.remove_parallel_sequences()
			
			debruijn.remove_single_sequence_components()
			
			debruijn.construct_assembly_ordering_labels()
			debruijn.reduce_to_single_path_max_weight()
			debruijn.contract_unique_overlaps(verbose = False)
			
			if saveall:
				debruijn.get_asqg_output(filename = filename_output+".asqg")
				debruijn.get_csv_output(filename = filename_output+".csv")
			debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)
			files_of_parts.append(filename_output+"_seqsonly.txt")
	return files_of_parts

def reconstruct_merge(filename_output_base, files_to_merge, merge_k, number_of_parts, delete_parts=True):
	reads = []
	for file in files_to_merge:
		if os.path.isfile(file):
			reads += dio.get_reads_from_file(filename = file)
			if delete_parts:
				os.remove(file)
		else:
			print "Error! File doesent exist: " + file
	
	reduce_to_single_path = True
	
	debruijn = fdgb.GraphData(reads, merge_k)
	# delete reads and kmers to save ram:
	debruijn.reads = []
	debruijn.kmers = []
	# run garbage collector:
	gc.collect()
	
	debruijn.remove_parallel_sequences(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	
	#debruijn.remove_single_sequence_components()
	#debruijn.reduce_to_single_largest_component()
	#debruijn.remove_tips()
	debruijn.construct_assembly_ordering_labels(verbose = False)
	
	filename_output = filename_output_base+"_p"+str(number_of_parts)+"_merged_k"+str(merge_k)
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)
	
	if reduce_to_single_path:
		debruijn.reduce_to_single_path_max_weight(verbose = False)
		debruijn.contract_unique_overlaps(verbose = False)
		debruijn.construct_assembly_ordering_labels(verbose = False)
		filename_output += "_singlepath"
		debruijn.get_asqg_output(filename = filename_output+".asqg")
		debruijn.get_csv_output(filename = filename_output+".csv")
		debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)
	
	return filename_output+"_seqsonly.txt"

def singlestep_assembly():
	num_of_reads = 5000
	# reads = dio.get_reads_from_fastq_file("Data/hcov229e_only.fq", num_of_reads=10000)
	data_dir = "Output/corona_allreads"
	k=40
	for i in range(20):
		read_start = i*3500
		reads = dio.get_reads_from_fastq_file("Data/hcov229e_only.fq", num_of_reads=num_of_reads, first_read=read_start)
		
		filename_output = data_dir+"/corona_realreads_singlestep_numreads"+str(num_of_reads)+"_s"+str(read_start)+"_k"+str(k)

		debruijn = fdgb.GraphData([reads], k)
		# delete reads to save ram:
		reads = []
		debruijn.reads = []
		#debruijn.kmers = []
		# run garbage collector:
		gc.collect()
	
		debruijn.remove_parallel_sequences(verbose = False)
		debruijn.contract_unique_overlaps(verbose = False)
	
		debruijn.remove_single_sequence_components()
		debruijn.construct_assembly_ordering_labels(verbose = False)
	
		debruijn.get_asqg_output(filename = filename_output+".asqg")
		debruijn.get_csv_output(filename = filename_output+".csv")
	
		debruijn.reduce_to_single_path_max_weight(verbose = False)
		debruijn.contract_unique_overlaps(verbose = False)
		debruijn.construct_assembly_ordering_labels(verbose = False)
		filename_output += "_singlepath"
		debruijn.get_asqg_output(filename = filename_output+".asqg")
		debruijn.get_csv_output(filename = filename_output+".csv")
		debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)

def singlestep_assembly_test(nr=1000, rs=0, k=40):
	data_dir = "Output/corona_allreads"
	reads = dio.get_reads_from_fastq_file("Data/hcov229e_only.fq", num_of_reads=nr, first_read=rs)
		
	filename_output = data_dir+"/crs_n"+str(nr)+"_k"+str(k)+"_s"+str(rs)

	debruijn = fdgb.GraphData([reads], k)
	# delete reads to save ram:
	reads = []
	debruijn.reads = []
	# run garbage collector:
	gc.collect()
	
	debruijn.remove_parallel_sequences(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	
	debruijn.construct_assembly_ordering_labels(verbose = False)
	
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)
	
	debruijn.remove_single_sequence_components()
	debruijn.remove_tips()
	debruijn.remove_parallel_sequences()
	debruijn.remove_insignificant_sequences(2)
	debruijn.contract_unique_overlaps()
	debruijn.construct_assembly_ordering_labels(verbose = False)
	
	filename_output += "_reduced"
	
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)
	
	debruijn.reduce_to_single_path_max_weight(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.construct_assembly_ordering_labels(verbose = False)
	filename_output += "_singlepath"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)

	
def reconstruction_pipeline():
	data_dir = "Output/corona_allreads"
	sourcefilename = data_dir+"/corona_realreads_n-1_k40.csv"
	outputfilename = data_dir+"/corona_realreads_k40_equisizeparts"

	if not os.path.exists(data_dir):
		os.makedirs(data_dir)
	'''
	# mininum weight of sequences in original graph that are considered:
	w = 50#30#10
	# number of parts for partitioning:
	p = 1000#500#2000
	# number of overlapping parts (local redundancy: larger for better reconstruction of total graph from parts):
	o = 2
	# k for debruijn-graphs of parts:
	k1 = 15#19#23
	# k for reconstruction, should be smaller than k1 to ensure that parts can be merged:
	k2 = 13#17#21
	# minimum weight of sequences in partition-graphs:
	w2 = 200#50
	# lower bound to the length of sequences in partition-graphs as multiple of k:
	f = 1#1.5
	'''
	# mininum weight of sequences in original graph that are considered:
	w = 5#30#10
	# number of parts for partitioning:
	p = 500#500#2000
	# number of overlapping parts (local redundancy: larger for better reconstruction of total graph from parts):
	o = 2
	# k for debruijn-graphs of parts:
	k1 = 23#19#23
	# k for reconstruction, should be smaller than k1 to ensure that parts can be merged:
	k2 = 21#17#21
	# minimum weight of sequences in partition-graphs:
	w2 = 2#50
	# lower bound to the length of sequences in partition-graphs as multiple of k:
	f = 1
	
	# get partitioning of sequences with minimum weight, as defined by param w, p, o:
	part_sequence_filenames = get_sequences_by_params(filename_input=sourcefilename, filename_output=outputfilename, min_weight=w, number_of_parts=p, overlap=o, partition_method="equisize")
	
	# build debruijn graphs for each part and merge them together:
	parts = reconstruct_parts(part_sequence_filenames, filename_output_base=outputfilename+"_w"+str(w), k=k1, minweight=w2, minlengthfactor=f, allow_recursion=False, saveall=False)
	if len(parts) > 1:
		reconstruct_merge(filename_output_base=outputfilename+"_w"+str(w)+"_k"+str(k1), files_to_merge=parts, merge_k=k2, number_of_parts=p)

'''
longest read: 56510
'''
	
if __name__ == '__main__':
	#data_dir = "Output/corona_allreads"
	#construct_network_graph(data_dir+"/corona_realreads_k40_w50_k17_p500_merged_k15.asqg")
	#construct_network_graph(data_dir+"/corona_realreads_n-1_k40.asqg")
	
	#reconstruction_pipeline()
	singlestep_assembly_test(nr=1000, rs=56000, k=40)
	
	'''
	data_dir = "Output/corona_allreads"
	inputfiles = []
	for i in range(3500, 70000, 3500):
		inputfiles.append(data_dir+"/corona_realreads_singlestep_numreads5000_s"+str(i)+"_k40_singlepath_seqsonly.txt")
	reconstruct_merge(data_dir+"/singlestep_results_merge", inputfiles, merge_k=50, number_of_parts=len(inputfiles), delete_parts=False)
	'''
