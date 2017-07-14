#!usr/bin/python

import timeit

import debruijn_graph_builder as dgb
import fast_debruijn_graph_builder as fdgb
import data_io as dio

dna = dio.genereate_dna(length=1000)
reads, alignment = dio.genereate_reads(dna, coverage=20, avg_read_length=50, remove_pct=0, mutation_pct=0.2, mutation_alphabet=["A","C","G","T"], both_directions=True)

#reads = dio.get_reads_from_file(filename = "Data/reads_bvdv_sample_50_2_[5000,5000]_0.txt")
k = 16

def measure_runtime():
	start_fdgb = timeit.default_timer()
	
	debruijn = fdgb.GraphData(reads, k, verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	
	stop_fdgb = timeit.default_timer()
	
	start_dgb = timeit.default_timer()
	
	debruijn = dgb.GraphData(reads, k, verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	
	stop_dgb = timeit.default_timer()
	
	print ("dgb: " + str(stop_dgb - start_dgb))
	print ("fdgb: " + str(stop_fdgb - start_fdgb))

def test_tip_removal():
	debruijn = fdgb.GraphData(reads, k, verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.get_asqg_output(filename="test_1_pre_tip_removal")
	debruijn.remove_tips(verbose = True)
	debruijn.get_asqg_output(filename="test_2_post_tip_removal_pre_unify")
	debruijn.remove_parallel_sequences(verbose = False)
	debruijn.get_asqg_output(filename="test_3_post_tip_removal_post_unify")
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.get_asqg_output(filename="test_4_post_tip_removal_post_contract")
	
test_tip_removal()