#!/usr/bin/python
# -*- coding: utf-8 -*-

import random

import sampleReads
import data_io as dio
import fast_debruijn_graph_builder as fdgb
import manjasDefinitionen as md

sampleReads.read_genomes()

k_values = [13,14,15,16,17,20]
readlengths_settings = [100,200,1000,5000,10000]
num_different_viruses = 2
viruses = [md.v1, md.v5]
num_reads = [5000, 5000]
error_percentage = 0

for readlength in readlengths_settings:
	casename_gen = "bvdv_sample_"+str(readlength)+"_"+str(num_different_viruses)+"_["
	for n_r_i in range(len(num_reads)-1):
		casename_gen += str(num_reads[n_r_i])+","
	casename_gen += str(num_reads[-1])+"]_"+str(error_percentage)
	readfilename = "Data/reads_"+casename_gen+".txt"
	
	sampleReads.samplereads(output_filename			= readfilename,
							read_length				= readlength,
							set_of_viruses			= viruses,
							number_of_reads			= num_reads,
							avg_error_percentage	= error_percentage)

	for k in k_values:
		casename = casename_gen + "_"+str(k)
		
		reads = dio.get_reads_from_file(filename = readfilename)
		debruijn = fdgb.GraphData(reads, k)
		debruijn.contract_unique_overlaps()
		debruijn.remove_parallel_sequences()
		
		debruijn.get_asqg_output(filename = "Output/"+casename+".asqg")
		debruijn.get_csv_output(filename = "Output/"+casename+".csv")
