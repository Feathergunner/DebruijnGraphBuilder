#!/usr/bin/python
# -*- coding: utf-8 -*-

import random

import sampleReads
import data_io as dio
import fast_debruijn_graph_builder as fdgb
import debruijn_graph_builder as dgb
import manjasDefinitionen as md

sampleReads.read_genomes()

k_values = [15, 17, 20, 25]
readlength_settings = [50, 100, 250, 500]
num_different_viruses = 2
viruses = [md.v1, md.v5]
num_reads_settings = [[5000, 5000], [2500, 2500], [1000, 1000], [500, 500]]
error_percentages = [0, 0.1, 0.5, 1, 5, 10, 15]
	
for i in range(len(readlength_settings)):
	readlength = readlength_settings[i]
	num_reads = num_reads_settings[i]
	for ep in error_percentages:
		casename = "bvdv_sample_"+str(readlength)+"_"+str(num_different_viruses)+"_["
		for n_r_i in range(len(num_reads)-1):
			casename += str(num_reads[n_r_i])+","
		casename += str(num_reads[-1])+"]_"+str(ep)
		
		readfilename = "Data/reads_"+casename+".txt"
		for k in k_values:
			casename_k = casename+"_"+str(k)
			
			print ("Working on case " + casename_k)
			
			sampleReads.samplereads(output_filename			= readfilename,
									read_length				= readlength,
									set_of_viruses			= viruses,
									number_of_reads			= num_reads,
									avg_error_percentage	= ep)
		
			reads = dio.get_reads_from_file(filename = readfilename)
			debruijn = fdgb.GraphData(reads, k)
			debruijn.contract_unique_overlaps()
			debruijn.remove_parallel_sequences()
			
			debruijn.get_asqg_output(filename = "Output/"+casename_k+".asqg")
			debruijn.get_csv_output(filename = "Output/"+casename_k+".csv")
	