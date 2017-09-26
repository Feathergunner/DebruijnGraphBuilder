#!usr/bin/python

import timeit
import os
import gc
import re

import debruijn_graph_builder as dgb
import fast_debruijn_graph_builder as fdgb
import data_io as dio
import sampleReads as sr
import manjasDefinitionen as md

'''
dna = dio.genereate_dna(length=500)
dio.write_dna_to_file("Output/test/genome_dna_test.txt", dna)

reads, alignment = dio.genereate_reads(dna, coverage=500, avg_read_length=50, remove_pct=0, mutation_pct=0.0, mutation_alphabet=["A","C","G","T"], both_directions=False, verbose=False)
#reads = dio.get_reads_from_file(filename = "Data/reads_bvdv_sample_50_2_[5000,5000]_0.txt")

sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads.txt", read_length=50, length_stddev=0, set_of_viruses=["dna_test"], number_of_reads=[500], replace_error_percentage=0.0, indel_error_percentage=0.2, inverted_reads=False)

reads = dio.get_reads_from_file("Output/test/testreads.txt")

#print reads

k = 20
'''

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
	
def test_reconstruction_3():
	dna = dio.genereate_dna(length=500)
	reads, alignment = dio.genereate_reads(dna, coverage=100, avg_read_length=50, remove_pct=0, mutation_pct=0.5, mutation_alphabet=["A","C","G","T"], both_directions=False, verbose=False)
	k = 20
	#print reads
	debruijn = fdgb.GraphData(reads, k, verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	debruijn.remove_single_sequence_components()
	#debruijn.contract_unique_overlaps(verbose = False)
	debruijn.construct_assembly_ordering_labels()
	debruijn.get_asqg_output(filename="Output/test/test_partition_1.asqg")
	debruijn.get_csv_output(filename="Output/test/test_partition_1.csv")
	
	debruijn.remove_insignificant_sequences()
	debruijn.remove_single_sequence_components()
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.get_asqg_output(filename="Output/test/test_partition_2.asqg")
	debruijn.get_csv_output(filename="Output/test/test_partition_2.csv")
	
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.get_asqg_output(filename="Output/test/test_partition_3.asqg")
	debruijn.get_csv_output(filename="Output/test/test_partition_3.csv")
	
	reconstructed = debruijn.get_relevant_sequences()[0]
	min_len = min(len(reconstructed), len(dna))
	print ("the shorter squence is "+str(min_len)+" symbols long")
	for i in range(min_len):
		if not dna[i] == reconstructed[i]:
			print "error at position "+str(i)
	
def test_reconstruction_4():
	import dataset_settings as ds
	setting = ds.largereads_test_recons_2
	
	k_absolute_settings = setting["k_absolute_settings"]
	coverage_factors = setting["coverage_factors"]
	readlength_settings = setting["readlength_settings"]
	number_of_reads_settings = setting["number_of_reads_settings"]
	error_percentages = setting["error_percentages"]
	num_different_viruses = setting["num_different_viruses"]
	set_of_viruses = setting["set_of_viruses"]
	output_dir = "Output/"+setting["name"]+"-i"
	name = setting["name"]
	
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	if not os.path.exists(output_dir+"/reads"):
		os.makedirs(output_dir+"/reads")
	
	cf = coverage_factors[0]
	readlength = readlength_settings[0]
	num_reads = cf*number_of_reads_settings[0]
	error_percentage = error_percentages[0]
	k = k_absolute_settings[0]
	epr = 0
	epi = error_percentage
	
	ep_string = "-".join(re.split(r'\.',str(error_percentage)))
	casename_gen = name+"_"+str(readlength)+"_"+str(num_different_viruses)+"_["+str(num_reads)+"]_"+ep_string
	readfilename = output_dir + "/reads/" + casename_gen + ".txt"
	
	sr.samplereads(output_filename	= readfilename,	read_length = readlength, set_of_viruses = set_of_viruses[:num_different_viruses], number_of_reads = [num_reads], replace_error_percentage = epr, indel_error_percentage = epi)
	gc.collect()
					
	casename = casename_gen + "_"+str(k)
	reads = dio.get_reads_from_file(filename = readfilename)
	debruijn = fdgb.GraphData(reads, k)
	# delete reads and kmers to save ram:
	debruijn.reads = []
	debruijn.kmers = []
	# run garbage collector:
	gc.collect()
	debruijn.contract_unique_overlaps(verbose=False)
	debruijn.remove_parallel_sequences()
	
	debruijn.get_asqg_output(filename = output_dir+"/"+casename+".asqg")
	debruijn.get_csv_output(filename = output_dir+"/"+casename+".csv")
	
	number_of_parts = 30
	k2 = 20
	debruijn.construct_assembly_ordering_labels(verbose = False)
	parts = debruijn.get_partition_of_sequences(number_of_parts, verbose = False)
	
	reconstructed_sequences = []
	
	for part_id in range(number_of_parts):
		#for part_id in range(3,4):
		part_sequences = parts[part_id]
		print len(part_sequences)
		#print part_sequences
		sequences_as_reads_string = "\n".join([seq.sequence for seq in part_sequences])
		outfile = file(output_dir+"/"+casename+"_p"+str(part_id)+"_reads.txt", 'w')
		outfile.write(sequences_as_reads_string)
		
		'''	
		debruijn_part = fdgb.GraphData([[seq.sequence for seq in part_sequences]], k2)
		# delete reads and kmers to save ram:
		debruijn_part.reads = []
		debruijn_part.kmers = []
		# run garbage collector:
		gc.collect()
		debruijn_part.contract_unique_overlaps(verbose = False)
		debruijn_part.remove_parallel_sequences(verbose = False)
		
		debruijn_part.remove_single_sequence_components()
		debruijn_part.construct_assembly_ordering_labels()
		# debruijn_part.remove_insignificant_sequences()
		# debruijn_part.reduce_to_single_path_max_weight()
		debruijn_part.contract_unique_overlaps(verbose = False)
		
		reconstructed_sequences.append(debruijn_part.get_relevant_sequences()[0])
		
		debruijn_part.get_asqg_output(filename = output_dir+"/"+casename+"_p"+str(part_id)+".asqg")
		debruijn_part.get_csv_output(filename = output_dir+"/"+casename+"_p"+str(part_id)+".csv")
		'''
	
	debruijn = []
	debruijn_part = []
	gc.collect()
	
	debruijn_recons = fdgb.GraphData([reconstructed_sequences], k2)
	debruijn_recons.contract_unique_overlaps()
	debruijn_recons.get_asqg_output(filename = output_dir+"/"+casename+"_recons.asqg")
	debruijn_recons.get_csv_output(filename = output_dir+"/"+casename+"_recons.csv")
		
test_reconstruction_4()
