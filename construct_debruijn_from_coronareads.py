#!/usr/bin/python
# -*- coding: utf-8 -*-

import gc

import data_io as dio
import fast_debruijn_graph_builder as fdgb

readfilename = "Data/2017-09-05_coronavirus.fq"
output_dir = "Output"

for k in [30, 40]:
	for num_reads in [2000, 5000]:
		casename = "corona_realreads_n"+str(num_reads)+"_k"+str(k)
		
		reads = dio.get_reads_from_fastq_file(readfilename, num_reads)
		
		debruijn = fdgb.GraphData([reads], k, alphabet={"A":"U", "U":"A", "G":"C", "C":"G"}, verbose=False)
		#print reads
		#print len(reads)
		# delete reads and kmers to save ram:
		reads = []
		debruijn.reads = []
		debruijn.kmers = []
		# run garbage collector:
		gc.collect()
		debruijn.contract_unique_overlaps(verbose=True)
		debruijn.remove_parallel_sequences()
		
		debruijn.get_asqg_output(filename = output_dir+"/"+casename+".asqg")
		debruijn.get_csv_output(filename = output_dir+"/"+casename+".csv")