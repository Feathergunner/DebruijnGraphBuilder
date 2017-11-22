#!usr/bin/python
# -*- coding: utf-8 -*-

import timeit
import os
import gc
import re
import sys
import math

#import debruijn_graph_builder as dgb
import fast_debruijn_graph_builder as fdgb
import veryfast_debruijn_graph_builder as vfdgb
#import parallel_fast_debruijn_graph_builder as pfdgb
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
	dna = dio.genereate_dna(length=1000)
	dio.write_dna_to_file("Output/test/genome_dna_test.txt", dna)
	if not os.path.isfile("Output/test/testreads.txt"):
		sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads.txt", read_length=50, length_stddev=0, set_of_viruses=["dna_test"], number_of_reads=[1], replace_error_percentage=0.0, indel_error_percentage=0.0, inverted_reads=False)
	
	#reads = dio.get_reads_from_file("Output/test/testreads.txt")
	reads = dio.get_reads_from_file("Output/corona_realreads/corona_realreads_high_evidence_sequences_w10_p6.txt")
	k = 30
	
	start_fdgb = 0
	stop_fdgb = 0
	start_vfdgb = 0
	stop_vfdgb = 0
	start_pfdgb = 0
	stop_pfdgb = 0
	
	gc.set_debug(gc.DEBUG_STATS & gc.DEBUG_COLLECTABLE & gc.DEBUG_UNCOLLECTABLE)
	
	start_fdgb = timeit.default_timer()
	debruijn = fdgb.GraphData(reads, k, verbose = False)
	print ("delete uneccesary data ...")
	debruijn.reads = []
	#debruijn.kmers = []
	debruijn.kmer_dict = {}
	debruijn.print_memory_usage()
	print ("Size of tracked objects pre-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	#print ("All tracked objects: ")
	#for obj in gc.get_objects():
	#	print ("\t" + str(obj) + "\n\tsize: "+str(sys.getsizeof(obj))+"\n")
	gc.collect()
	print ("Size of tracked objects post-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.print_memory_usage()
	print ("Size of tracked objects pre-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	gc.collect()
	print ("Size of tracked objects post-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	debruijn.remove_parallel_sequences(verbose = False)
	stop_fdgb = timeit.default_timer()
	debruijn.print_memory_usage()
	print ("Size of tracked objects pre-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	gc.collect()
	print ("Size of tracked objects post-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	debruijn.get_asqg_output(filename="Output/test/fdgb_test.asqg")
	debruijn.get_csv_output(filename="Output/test/fdgb_test.csv")

	#debruijn = 0
	#gc.collect()
	
	#start_pfdgb = timeit.default_timer()
	#debruijn = pfdgb.GraphData(reads, k, verbose = False)
	#debruijn.remove_parallel_sequences(verbose = False)
	#debruijn.contract_unique_overlaps_parallel_master_process(number_of_processes=4, verbose = False)
	#stop_pfdgb = timeit.default_timer()
	#debruijn.get_asqg_output(filename="Output/test/pfdgb_test.asqg")
	
	debruijn = 0
	gc.collect()

	print ("")
	start_vfdgb = timeit.default_timer()
	debruijn = vfdgb.GraphData(reads, k, verbose = False)
	print ("delete uneccesary data ...")
	debruijn.kmers = []
	debruijn.kmer_dict = {}
	debruijn.print_memory_usage()
	print ("Size of tracked objects pre-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	gc.collect()
	print ("Size of tracked objects post-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.print_memory_usage()
	print ("Size of tracked objects pre-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	gc.collect()
	print ("Size of tracked objects post-collection: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	debruijn.remove_parallel_sequences(verbose = False)
	stop_vfdgb = timeit.default_timer()
	debruijn.print_memory_usage()
	print ("Size of tracked objects : " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	debruijn.get_asqg_output(filename="Output/test/vfdgb_test.asqg")
	debruijn.get_csv_output(filename="Output/test/vfdgb_test.csv")

	print ("fdgb: " + str("%.2f" % (stop_fdgb - start_fdgb)))
	print ("vfdgb: " + str("%.2f" % (stop_vfdgb - start_vfdgb)))
	print ("pfdgb: " + str("%.2f" % (stop_pfdgb - start_pfdgb)))

'''
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
	setting = ds.largereads_test_recons_3
	
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
	parts = debruijn.get_partition_of_sequences(number_of_parts, overlap=4, verbose = False)
	
	reconstructed_sequences = []
	
	for part_id in range(number_of_parts):
		#for part_id in range(3,4):
		part_sequences = parts[part_id]
		print len(part_sequences)
		#print part_sequences
		sequences_as_reads_string = "\n".join([seq.sequence for seq in part_sequences])
		outfile = file(output_dir+"/"+casename+"_p"+str(part_id)+"_reads.txt", 'w')
		outfile.write(sequences_as_reads_string)
		
	debruijn = []
	debruijn_part = []
	gc.collect()
	
	debruijn_recons = fdgb.GraphData([reconstructed_sequences], k2)
	debruijn_recons.contract_unique_overlaps()
	debruijn_recons.get_asqg_output(filename = output_dir+"/"+casename+"_recons.asqg")
	debruijn_recons.get_csv_output(filename = output_dir+"/"+casename+"_recons.csv")
		
def test_recons_from_sequences():
	read_dir = "Output/largereads-recons-test2-i"
	read_basename = "largereads-recons-test3_5000_1_[7500]_5-0_30_"
	k = 18
	num_of_sets = 30
	
	for i in range(num_of_sets):
		seqfilename = read_dir+"/"+read_basename+"p"+str(i)+"_reads.txt"
		reads = dio.get_reads_from_file(filename = seqfilename)
		debruijn = fdgb.GraphData(reads, k)
		# delete reads and kmers to save ram:
		debruijn.reads = []
		debruijn.kmers = []
		# run garbage collector:
		gc.collect()
		
		debruijn.contract_unique_overlaps(verbose = False)
		debruijn.remove_parallel_sequences(verbose = False)
		
		debruijn.get_asqg_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step1.asqg")
		debruijn.get_csv_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step1.csv")
		
		debruijn.remove_insignificant_sequences()
		debruijn.contract_unique_overlaps(verbose = False)
		debruijn.remove_single_sequence_components()
		debruijn.construct_assembly_ordering_labels()
		
		debruijn.get_asqg_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step2.asqg")
		debruijn.get_csv_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step2.csv")
		
		debruijn.reduce_to_single_path_max_weight()
		debruijn.contract_unique_overlaps(verbose = False)
		
		debruijn.get_asqg_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step3.asqg")
		debruijn.get_csv_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step3.csv")
		
		debruijn.remove_short_sequences()
		
		debruijn.get_asqg_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step4.asqg")
		debruijn.get_csv_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step4.csv")
		debruijn.write_sequences_to_file(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step4.txt")
		
def test_recons_merge():
	read_dir = "Output\largereads-recons-test2-i"
	read_basename = "largereads-recons-test3_5000_1_[7500]_5-0_30_"
	num_of_sets = 30
	k = 18
	
	reads = []
	for i in range(num_of_sets):
		seqfilename = read_dir+"/"+read_basename+"_k"+str(k)+"_p"+str(i)+"_step4.txt"
		reads += dio.get_reads_from_file(filename = seqfilename)
		
	debruijn = fdgb.GraphData(reads, k)
	# delete reads and kmers to save ram:
	debruijn.reads = []
	debruijn.kmers = []
	# run garbage collector:
	gc.collect()
	
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.remove_parallel_sequences(verbose = False)
	debruijn.get_asqg_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_merged_step1.asqg")
	debruijn.get_csv_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_merged_step1.csv")
	#debruijn.remove_tips()
	debruijn.construct_assembly_ordering_labels()
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps(verbose = False)
	
	debruijn.get_asqg_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_merged_step2.asqg")
	debruijn.get_csv_output(filename = read_dir+"/"+read_basename+"_k"+str(k)+"_merged_step2.csv")
	
'''

def test_assembly_ordering():
	gl = 2000
	rl = 500
	nr = 500
	ep = 5.0

	dna = dio.genereate_dna(length=gl)
	dio.write_dna_to_file("Output/test/genome_dna_test.txt", dna)
	if not os.path.isfile("Output/test/testreads_assembly.txt"):
		sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads_assembly.txt", read_length=rl, length_stddev=0, set_of_viruses=["dna_test"], number_of_reads=[nr], replace_error_percentage=0.0, indel_error_percentage=ep, inverted_reads=False)
	
	'''
	dna2 = dio.genereate_dna(length=gl)
	dio.write_dna_to_file("Output/test/genome_dna_test2.txt", dna)
	if not os.path.isfile("Output/test/testreads_assembly2.txt"):
		sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads_assembly2.txt", read_length=rl, length_stddev=0, set_of_viruses=["dna_test2"], number_of_reads=[nr], replace_error_percentage=0.0, indel_error_percentage=ep, inverted_reads=False)
	'''
		
	reads = dio.get_reads_from_file("Output/test/testreads_assembly.txt")
	#reads += dio.get_reads_from_file("Output/test/testreads_assembly2.txt")
	
	filename_output = "Output/test/assembly_ordering_test"
	k = 30
	
	debruijn = fdgb.GraphData(reads, k)
	# delete reads and kmers to save ram:
	debruijn.reads = []
	debruijn.kmers = []
	# run garbage collector:
	gc.collect()
	
	debruijn.remove_parallel_sequences(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	
	debruijn.remove_single_sequence_components()
	#debruijn.reduce_to_single_largest_component()
	#debruijn.remove_tips()
	debruijn.construct_assembly_ordering_labels(verbose = False)
	
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	
	debruijn.reduce_to_single_path_max_weight(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.construct_assembly_ordering_labels(verbose = False)
	filename_output += "_singlepath"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")

def construct_consensus_from_multiple_parts():
	size_of_parts = 100
	k = 40
	filename = "Data/hcov229e_only.fq"
	filepath_output = "Output/corona_recons_multiparts"
	filename_output_base = filepath_output+"/crm_partsize"+str(size_of_parts)
	if not os.path.exists(filepath_output):
		os.makedirs(filepath_output)

	readpartition = dio.get_read_partition_by_readlength(filename = filename, size_of_parts=size_of_parts)
	n = len(readpartition)
	p = 3*(n/4)
	# consider only second half of partition (parts with longer reads)
	for rp in readpartition[p:]:
		minreadlength = min([x[0] for x in rp])
		k = get_adaptive_k(minreadlength)
		p += 1
		filename_output = filename_output_base+"_k"+str(k)+"_p"+str(p)
		read_ids = [x[1] for x in rp]
		reads = dio.get_reads_from_fastq_file_by_length(filename = filename, read_ids = read_ids)
		
		debruijn = fdgb.GraphData([reads], k)

		# delete reads and kmers to save ram:
		debruijn.reads = []
		debruijn.kmers = []
		# run garbage collector:
		gc.collect()
	
		debruijn.remove_parallel_sequences(verbose = False)
		debruijn.contract_unique_overlaps(verbose = False)
		
		debruijn.remove_single_sequence_components()
		debruijn.construct_assembly_ordering_labels(verbose = False)

		debruijn.get_asqg_output(filename = filename_output+".asqg")
		debruijn.get_csv_output(filename = filename_output+".csv")
		debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt", addweights=True)
		
		debruijn.reduce_to_single_path_max_weight(verbose = False)
		debruijn.contract_unique_overlaps(verbose = False)
		debruijn.construct_assembly_ordering_labels(verbose = False)
		filename_output += "_singlepath"
		debruijn.get_asqg_output(filename = filename_output+".asqg")
		debruijn.get_csv_output(filename = filename_output+".csv")
		debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt")

def construct_consensus_from_part(k, read_ids, readfile, filepath_output, filename_output):
	if not os.path.exists(filepath_output):
		os.makedirs(filepath_output)

	reads = dio.get_reads_from_fastq_file_by_length(filename = readfile, read_ids = read_ids)
	
	debruijn = fdgb.GraphData([reads], k)
	# delete reads and kmers to save ram:
	debruijn.reads = []
	debruijn.kmers = []
	# run garbage collector:
	gc.collect()

	debruijn.remove_parallel_sequences(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	
	debruijn.remove_single_sequence_components()
	debruijn.construct_assembly_ordering_labels(verbose = False)

	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt", addweights=True)
	
	debruijn.reduce_to_single_path_max_weight(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.construct_assembly_ordering_labels(verbose = False)
	filename_output += "_singlepath"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt")
	
def exp_construct_consensus_from_specific_part():
	size_of_parts = 50
	k=50
	readfilename = "Data/hcov229e_only.fq"
	readpartition = dio.get_read_partition_by_readlength(filename = readfilename, size_of_parts=size_of_parts)
	# get read ids of 50 largest reads:
	read_ids = [x[1] for x in readpartition[-1]]
	filepath_output = "Output/corona_recons_multiparts"
	filename_output = filepath_output+"/crm_partsize"+str(size_of_parts)+"_k"+str(k)+"_p"+str(len(readpartition))
	
	construct_consensus_from_part(k=k, read_ids = read_ids, readfile = readfilename, filepath_output = filepath_output, filename_output = filename_output)
		
def get_adaptive_k(readlength):
	'''
	if readlength < 100:
		return 25
	elif readlength < 200:
		return 30
	elif readlength < 500:
		return 35
	elif readlength < 1000:
		return 40
	elif readlength < 3000:
		return 45
	else:
		return 50
	'''
	return int(math.log(readlength, 2.1)*4)

if __name__ == '__main__':
	#test_reconstruction_4()
	#test_recons_from_sequences()
	#test_recons_merge()
	
	#measure_runtime()
	
	#test_assembly_ordering()
	#exp_construct_consensus_from_specific_part()
	construct_consensus_from_multiple_parts()