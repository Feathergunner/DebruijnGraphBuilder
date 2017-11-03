#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import gc
import os

import data_io as dio
import fast_debruijn_graph_builder as fdgb

def get_sequences_by_params(filename_input, filename_output, min_weight=1, number_of_parts=1, overlap=1):
	with open(sourcefilename) as inputfile:
		lines = inputfile.readlines()
	i = 0
	j = 0
	n = len(lines)
	min_label = False
	max_label = False
	
	filename_output += "_w"+str(min_weight)
	
	sequences = []
	for l in lines:
		if i > 0:
			if i%10000 == 0:#i>0 and i<100:
				print str(i)+"/"+str(n)+" - "+str(j)
			data = re.split(r',',l)
			if int(data[2]) > min_weight:
				label = int(data[3])
				sequences.append([data[1], label, data[2]])
				
				if not min_label or label < min_label:
					min_label = label
				if not max_label or label > max_label:
					max_label = label
				j += 1
		i+=1
	
	if number_of_parts > 1:
		partition_size = (max_label-min_label)/(number_of_parts-overlap)
		partition_start_diff = (max_label-min_label)/(number_of_parts+overlap)
		
		for i in range(number_of_parts):
			partname = filename_output+"_p"+str(i)+".txt"
			part_seqs = []
			for s in sequences:
				if s[1] > i*partition_start_diff and s[1] < i*partition_start_diff+partition_size:
					part_seqs.append(s[0]+","+str(s[2])+"\n")
			outputfile = file(partname, 'w')
			for seq in part_seqs:
				outputfile.write(seq)
	else:
		part_seqs = []
		for s in sequences:
			part_seqs.append(s[0]+","+str(s[2])+"\n")
		outputfile = file(partname, 'w')
		for seq in part_seqs:
			outputfile.write(seq)
		
def reconstruct_parts(filename_input_base, filename_output_base, number_of_parts, k1, minweight, minlengthfactor):
	for i in range(number_of_parts):
		filename_input = filename_input_base+"_p"+str(i)+".txt"
		filename_output = filename_output_base+"_k"+str(k1)+"_p"+str(i)
		
		reads = dio.get_reads_from_file(filename_input)
		if len(reads)>0:
			
			debruijn = fdgb.GraphData(reads, k1, verbose=False)
			# delete reads and kmers to save ram:
			reads = []
			debruijn.reads = []
			debruijn.kmers = []
			# run garbage collector:
			gc.collect()
			debruijn.contract_unique_overlaps(verbose=False)
			debruijn.remove_parallel_sequences()
	
			debruijn.remove_insignificant_sequences(minimal_weight=minweight)
			debruijn.contract_unique_overlaps(verbose = False)
			debruijn.remove_short_sequences(length_bound_by_multiple_of_k=minlengthfactor)
			debruijn.contract_unique_overlaps(verbose = False)
			
			debruijn.remove_single_sequence_components()
			
			#debruijn.construct_assembly_ordering_labels()
			#debruijn.reduce_to_single_path_max_weight()
			#debruijn.contract_unique_overlaps(verbose = False)
			
			debruijn.get_asqg_output(filename = filename_output+".asqg")
			debruijn.get_csv_output(filename = filename_output+".csv")
			debruijn.write_sequences_to_file(filename = filename_output+"_seqsonly.txt", addweights=True)

def reconstruct_merge(filename_output_base, number_of_parts, k1, k2):
	reads = []
	for i in range(number_of_parts):
		seqfilename = filename_output_base+"_k"+str(k1)+"_p"+str(i)+"_seqsonly.txt"
		if os.path.isfile(seqfilename):
			reads += dio.get_reads_from_file(filename = seqfilename)
	
	debruijn = fdgb.GraphData(reads, k2)
	# delete reads and kmers to save ram:
	debruijn.reads = []
	debruijn.kmers = []
	# run garbage collector:
	gc.collect()
	
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	#debruijn.get_asqg_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_merged_step1.asqg")
	#debruijn.get_csv_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_merged_step1.csv")
	#debruijn.remove_tips()
	#debruijn.construct_assembly_ordering_labels()
	#debruijn.reduce_to_single_path_max_weight()
	#debruijn.contract_unique_overlaps(verbose = False)
	
	debruijn.get_asqg_output(filename = filename_output_base+"_k"+str(k1)+"_p"+str(number_of_parts)+"_merged_k"+str(k2)+".asqg")
	debruijn.get_csv_output(filename = filename_output_base+"_k"+str(k1)+"_p"+str(number_of_parts)+"_merged_k"+str(k2)+".csv")
	
	
if __name__ == '__main__':

	data_dir = "Output/corona_allreads"
	sourcefilename = data_dir+"/corona_realreads_n-1_k40.csv"
	outputfilename = data_dir+"/corona_realreads_k40"

	if not os.path.exists(data_dir):
		os.makedirs(data_dir)
	# mininum weight of sequences in original graph that are considered:
	w = 5
	# number of parts for partitioning:
	p = 1000
	# number of overlapping parts (local redundancy: larger for better reconstruction of total graph from parts):
	o = 10
	# k for debruijn-graphs of parts:
	k1 = 20
	# k for reconstruction, should be smaller than k1 to ensure that parts can be merged:
	k2 = 17
	# minimum weight of sequences in partition-graphs:
	w2 = 50
	# lower bound to the length of sequences in partition-graphs as multiple of k:
	f = 2
	
	# get partitioning of sequences with minimum weight, as defined by param w, p, o:
	get_sequences_by_params(filename_input=sourcefilename, filename_output=outputfilename, min_weight=w, number_of_parts=p, overlap=o)
	
	# build debruijn graphs for each part and merge them together:
	reconstruct_parts(filename_input_base=outputfilename+"_w"+str(w), filename_output_base=outputfilename+"_w"+str(w), number_of_parts=p, k1=k1, minweight=w2, minlengthfactor=f)
	reconstruct_merge(filename_output_base=outputfilename+"_w"+str(w), number_of_parts=p, k1=k1, k2=k2)