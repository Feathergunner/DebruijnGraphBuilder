#!/usr/bin/python
# -*- coding: utf-8 -*-

import re

import sampleReads
import data_io as dio
import fast_debruijn_graph_builder as fdgb
import manjasDefinitionen as md

sampleReads.read_genomes()

k_relative_settings = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
readlength_settings = [50, 100, 250, 500, 1000]
number_of_reads_settings = [[500,500],[250,250],[100,100],[50,50],[25,25]]
coverage_factors = [1,5,10,15,20] 
error_percentages = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]
set_of_viruses = [md.v1, md.v5]
num_different_viruses = 2

for cf in coverage_factors:
	for i in range(len(readlength_settings)):
		readlength = readlength_settings[i]
		num_reads = cf*number_of_reads_settings[i]
		
		for error_percentage in error_percentages:
			ep_string = "-".join(re.split(r'\.',str(error_percentage)))
			casename_gen = "bvdv_sample_"+str(readlength)+"_"+str(num_different_viruses)+"_["
			for n_r_i in range(len(num_reads)-1):
				casename_gen += str(num_reads[n_r_i])+","
			casename_gen += str(num_reads[-1])+"]_"+ep_string
			readfilename = "Data/reads_"+casename_gen+".txt"
		
			sampleReads.samplereads(output_filename			= readfilename,
									read_length				= readlength,
									set_of_viruses			= set_of_viruses,
									number_of_reads			= num_reads,
									avg_error_percentage	= error_percentage)
									
			for k_rel in k_relative_settings:
				k = int(readlength * k_rel)
				casename = casename_gen + "_"+str(k)
			
				reads = dio.get_reads_from_file(filename = readfilename)
				debruijn = fdgb.GraphData(reads, k)
				debruijn.contract_unique_overlaps()
				debruijn.remove_parallel_sequences()
			
				debruijn.get_asqg_output(filename = "Output/"+casename+".asqg")
				debruijn.get_csv_output(filename = "Output/"+casename+".csv")