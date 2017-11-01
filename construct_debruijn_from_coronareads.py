#!/usr/bin/python
# -*- coding: utf-8 -*-

import gc
import os

import data_io as dio
import fast_debruijn_graph_builder as fdgb

readfilename = "Data/hcov229e_only.fq"
output_dir = "Output/corona_realparts"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

k = 40
num_reads = 2000
start = 1
	
casename = "corona_realreads_s_"+str(start)+"n_"+str(num_reads)+"_k"+str(k)

#if not os.path.isfile(output_dir+"/"+casename+".asqg"):
reads = dio.get_reads_from_fastq_file(readfilename, num_of_reads=num_reads, first_read=start)
debruijn = fdgb.GraphData([reads], k, alphabet={"A":"U", "U":"A", "G":"C", "C":"G"}, verbose=False)
# delete reads and kmers to save ram:
reads = []
debruijn.reads = []
# run garbage collector:
gc.collect()
debruijn.contract_unique_overlaps(verbose=False)
debruijn.remove_parallel_sequences()
debruijn.remove_single_sequence_components()
debruijn.get_asqg_output(filename = output_dir+"/"+casename+".asqg")
debruijn.get_csv_output(filename = output_dir+"/"+casename+".csv")
