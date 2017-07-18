#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os.path
import gc

import sampleReads
import data_io as dio
import fast_debruijn_graph_builder as fdgb
import manjasDefinitionen as md

sampleReads.read_genomes()

k_absolute_settings = [10,12,14,16,18,20,25,30]
readlength_settings = [50, 100, 250, 500, 1000]
number_of_reads_settings = [[500,500],[250,250],[100,100],[50,50],[25,25]]
coverage_factors = [1,5,10,15,20]
error_percentages = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]
set_of_viruses = [md.v1, md.v5]
num_different_viruses = 2

output_dir = "Output/general_absk"

for cf in coverage_factors:
	for i in range(len(readlength_settings)):
		readlength = readlength_settings[i]
		num_reads = [cf*nr for nr in number_of_reads_settings[i]]
		
		for error_percentage in error_percentages:
			ep_string = "-".join(re.split(r'\.',str(error_percentage)))
			casename_gen = "bvdv_sample_"+str(readlength)+"_"+str(num_different_viruses)+"_["
			for n_r_i in range(len(num_reads)-1):
				casename_gen += str(num_reads[n_r_i])+","
			casename_gen += str(num_reads[-1])+"]_"+ep_string
			readfilename = "Data/reads_"+casename_gen+".txt"
			
			if not os.path.isfile(readfilename):
				sampleReads.samplereads(output_filename			= readfilename,
										read_length				= readlength,
										set_of_viruses			= set_of_viruses,
										number_of_reads			= num_reads,
										avg_error_percentage	= error_percentage)
			else:
				print ("Reads already exist!")
									
			for k in k_absolute_settings:
				# run garbage collector:
				gc.collect()
				
				casename = casename_gen + "_"+str(k)
				
				if not os.path.isfile(output_dir+"/"+casename+".asqg"):
					print ("Working on case "+casename)
					try:
						reads = dio.get_reads_from_file(filename = readfilename)
						debruijn = fdgb.GraphData(reads, k)
						# delete reads and kmers to save ram:
						debruijn.reads = []
						debruijn.kmers = []
						# run garbage collector:
						gc.collect()
						debruijn.contract_unique_overlaps()
						debruijn.remove_parallel_sequences()
					
						debruijn.get_asqg_output(filename = output_dir+"/"+casename+".asqg")
						debruijn.get_csv_output(filename = output_dir+"/"+casename+".csv")
					except:
						pass
					
				else:
					print ("Data already exists!")
