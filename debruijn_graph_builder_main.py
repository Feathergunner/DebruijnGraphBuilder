#!usr/bin/python

import timeit

import debruijn_graph_builder as dgb
import fast_debruijn_graph_builder as fdgb
import data_io as dio
import sampleReads as sr
import manjasDefinitionen as md

dna = dio.genereate_dna(length=500)
dio.write_dna_to_file("Output/test/genome_dna_test.txt", dna)

#reads, alignment = dio.genereate_reads(dna, coverage=500, avg_read_length=50, remove_pct=0, mutation_pct=0.0, mutation_alphabet=["A","C","G","T"], both_directions=False, verbose=False)
#reads = dio.get_reads_from_file(filename = "Data/reads_bvdv_sample_50_2_[5000,5000]_0.txt")

sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads.txt", read_length=50, length_stddev=0, set_of_viruses=["dna_test"], number_of_reads=[500], replace_error_percentage=0.0, indel_error_percentage=0.2, inverted_reads=False)

reads = dio.get_reads_from_file("Output/test/testreads.txt")

#print reads

k = 20

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
	
	
def test_reconstruction_1():
	dna = dio.genereate_dna(length=500)
	dio.write_dna_to_file("Output/test/genome_dna_test.txt", dna)
	
	#reads, alignment = dio.genereate_reads(dna, coverage=500, avg_read_length=50, remove_pct=0, mutation_pct=0.0, mutation_alphabet=["A","C","G","T"], both_directions=False, verbose=False)
	#reads = dio.get_reads_from_file(filename = "Data/reads_bvdv_sample_50_2_[5000,5000]_0.txt")
	
	sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads.txt", read_length=50, length_stddev=0, set_of_viruses=["dna_test"], number_of_reads=[500], replace_error_percentage=0.0, indel_error_percentage=0.2, inverted_reads=False)
	
	reads = dio.get_reads_from_file("Output/test/testreads.txt")
	
	#print reads
	
	k = 20
	
	debruijn = fdgb.GraphData(reads, k, verbose = False)
	debruijn.construct_assembly_ordering_labels()
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	
	debruijn.get_asqg_output(filename="Output/test/recon2_pre_partition_k-"+str(k))
	
	parts = debruijn.get_partition_of_sequences(10, verbose=True)
	
	k2 = 20
	#print parts
	p = 0
	part_reads = [reads[i] for i in parts[p]]
	#print [len(p[0]) for p in part_reads]
	print len([r for r in reads if r not in part_reads])
	print len([r for r in part_reads])
	#print reads
	
	debruijn = fdgb.GraphData(part_reads, k2, verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	
	debruijn.get_asqg_output(filename="Output/test/recon2_part"+str(p)+"_k-"+str(k)+"-"+str(k2))
	
def test_reconstruction_2():
	readlength = 10000
	number_of_reads = 100
	sampleReads.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads.txt", read_length = readlength, set_of_viruses=[md.cv], number_of_reads = num_reads, replace_error_percentage = 0.0, indel_error_percentage = 15.0, inverted_reads=False)
	
	reads = dio.get_reads_from_file("Output/test/testreads.txt")
	
	k = 30
	
	debruijn = fdgb.GraphData(reads, k, verbose = False)
	debruijn.construct_assembly_ordering_labels()
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	
	debruijn.get_asqg_output(filename="Output/test/recon_pre_partition_k-"+str(k))
	
	parts = debruijn.get_partition_of_sequences(100, verbose=True)
	
	k2 = 20
	#print parts
	p = 0
	part_reads = [reads[i] for i in parts[p]]
	#print [len(p[0]) for p in part_reads]
	print len([r for r in reads if r not in part_reads])
	print len([r for r in part_reads])
	#print reads
	
	debruijn = fdgb.GraphData(part_reads, k2, verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	
	debruijn.get_asqg_output(filename="Output/test/recon_part"+str(p)+"_k-"+str(k)+"-"+str(k2))
	
