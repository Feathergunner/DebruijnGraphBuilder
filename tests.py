#!usr/bin/python
# -*- coding: utf-8 -*-

import meta
import data_gen as dgen
import sampleReads as sr
import data_io as dio
import fast_debruijn_graph_builder as fdgb

import timeit
import sys
import gc
import os
import re
import copy

def test_basic_functionality(new_dataset=False, data_dir="Output/test", verbose=False):
	print ("Starting basic test of de bruijn graph construction")
	if not os.path.exists(data_dir):
		os.mkdir(data_dir)
	if new_dataset or not os.path.isfile(data_dir+"/test_dna.txt"):
		dna = dgen.generate_dna(length=3000)
		dgen.write_dna_to_file(data_dir+"/test_dna.txt", dna)
	else:
		with open(data_dir+"/test_dna.txt") as inputfile:
			lines = inputfile.readlines()
		dna = lines[0]
	
	if new_dataset or not os.path.isfile(data_dir+"/test_reads.txt"):
		reads = dgen.samplereads(dna, number_of_reads=2000, replace_error_percentage=0.0, indel_error_percentage=5.0, mutation_alphabet=["A","C","G","T"], read_length_mean=200, read_length_stddev=0, readlength_distribution='exponential')
		with open(data_dir+"/test_reads.txt", 'w') as outf:
			for read in reads:
				outf.write(read + '\n')
	else:
		reads = dio.get_reads_from_file(filename=data_dir+"/test_reads.txt")
	k = 20
	
	start_fdgb = 0
	stop_fdgb = 0
	
	print ("Construct de bruijn graph")
	
	start_fdgb = timeit.default_timer()
	debruijn = fdgb.GraphData([reads], k, reduce_data=True, simplify_graph=True, remove_tips=False, construct_labels=True, directed_reads=False, verbose=verbose)
	stop_fdgb = timeit.default_timer()
	debruijn.print_memory_usage()
	print ("Size of tracked objects: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	print ("Running time of fdgb: " + str("%.2f" % (stop_fdgb - start_fdgb)))
	
	debruijn.get_asqg_output(filename = data_dir+"/test_basic.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_basic.csv")
	
	print ("Number of sequences: "+str(len([seq for seq in debruijn.sequences if seq.is_relevant])))
	print ("Number of overlaps: "+str(len(debruijn.overlaps)))
	
def test_basic_functionality_small(new_dataset=False, data_dir="Output/test", verbose=False):
	print ("Starting basic test of de bruijn graph construction")
	if not os.path.exists(data_dir):
		os.mkdir(data_dir)
	if new_dataset or not os.path.isfile(data_dir+"/test_dna.txt"):
		dna = dgen.generate_dna(length=200)
		dgen.write_dna_to_file(data_dir+"/test_dna.txt", dna)
	else:
		with open(data_dir+"/test_dna.txt") as inputfile:
			lines = inputfile.readlines()
		dna = lines[0]
	
	if new_dataset or not os.path.isfile(data_dir+"/test_reads.txt"):
		reads = dgen.samplereads(dna, number_of_reads=100, replace_error_percentage=0.5, indel_error_percentage=0.0, mutation_alphabet=["A","C","G","T"], read_length_mean=50, read_length_stddev=0)
		with open(data_dir+"/test_reads.txt", 'w') as outf:
			for read in reads:
				outf.write(read + '\n')
	else:
		reads = dio.get_reads_from_file(filename=data_dir+"/test_reads.txt")
	k = 17
	
	start_fdgb = 0
	stop_fdgb = 0
	
	print ("Construct de bruijn graph")
	start_fdgb = timeit.default_timer()
	debruijn = fdgb.GraphData(reads, k, reduce_data=True, simplify_graph=True, remove_tips=False, construct_labels=True, directed_reads=False, verbose=verbose)
	stop_fdgb = timeit.default_timer()
	debruijn.print_memory_usage()
	print ("Size of tracked objects: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	print ("Running time of fdgb: " + str("%.2f" % (stop_fdgb - start_fdgb)))
	
	debruijn.get_asqg_output(filename = data_dir+"/test_basic.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_basic.csv")
	
	print ("Number of sequences: "+str(len([seq for seq in debruijn.sequences if seq.is_relevant])))
	print ("Number of overlaps: "+str(len(debruijn.overlaps)))
	
	#for seq in debruijn.sequences:
	#	if seq.is_relevant:
	#		print seq.label_p
	
def test_basic_consensus_construction(new_dataset=False, data_dir="Output/test", verbose=False):
	print ("Starting basic test of de bruijn graph construction")
	if not os.path.exists(data_dir):
		os.mkdir(data_dir)
	if new_dataset or not os.path.isfile(data_dir+"/test_dna.txt"):
		dna = dgen.generate_dna(length=1000)
		dgen.write_dna_to_file(data_dir+"/test_dna.txt", dna)
	else:
		with open(data_dir+"/test_dna.txt") as inputfile:
			lines = inputfile.readlines()
		dna = lines[0]
	
	if new_dataset or not os.path.isfile(data_dir+"/test_reads.txt"):
		reads = dgen.samplereads(dna, number_of_reads=1000, replace_error_percentage=0.5, indel_error_percentage=0.0, mutation_alphabet=["A","C","G","T"], read_length_mean=50, read_length_stddev=0)
		with open(data_dir+"/test_reads.txt", 'w') as outf:
			for read in reads:
				outf.write(read + '\n')
	else:
		reads = dio.get_reads_from_file(filename=data_dir+"/test_reads.txt")
	k = 17
	
	start_fdgb = 0
	stop_fdgb = 0
	print ("Construct de bruijn graph")
	start_fdgb = timeit.default_timer()
	debruijn = fdgb.GraphData(reads, k, reduce_data=True, simplify_graph=True, remove_tips=False, construct_labels=False, directed_reads=True, verbose=verbose)
	
	debruijn.get_asqg_output(filename = data_dir+"/test_1_base.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_1_base.csv")
	
	debruijn.remove_tips()
	
	debruijn.get_asqg_output(filename = data_dir+"/test_2_posttipremoval.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_2_posttipremoval.csv")
	
	debruijn.remove_insignificant_sequences()
	debruijn.contract_unique_overlaps()
	debruijn.remove_single_sequence_components()
	
	debruijn.get_asqg_output(filename = data_dir+"/test_3_postsequenceremoval.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_3_postsequenceremoval.csv")
	
	debruijn.remove_insignificant_overlaps()
	debruijn.contract_unique_overlaps()
	debruijn.remove_single_sequence_components()
	
	debruijn.get_asqg_output(filename = data_dir+"/test_4_postoverlapremoval.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_4_postoverlapremoval.csv")
	
	debruijn.reduce_to_single_largest_component()
	debruijn.construct_assembly_ordering_labels()
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps()
	
	debruijn.get_asqg_output(filename = data_dir+"test_5_consensus.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_5_consensus.csv")
	debruijn.write_sequences_to_file(filename = data_dir+"/test_5_consensus.fasta", asfasta = True)
	
	stop_fdgb = timeit.default_timer()
	debruijn.print_memory_usage()
	print ("Size of tracked objects: " +str(sum(sys.getsizeof(i) for i in gc.get_objects())/1000000.0) + "MB")
	print ("Running time of fdgb: " + str("%.2f" % (stop_fdgb - start_fdgb)))
	
	print ("Number of sequences: "+str(len([seq for seq in debruijn.sequences if seq.is_relevant])))
	print ("Number of overlaps: "+str(len(debruijn.overlaps)))

def test_reduce_to_largest_single_sequence():
	reads = ["ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTTAAAAACCCCCGGGGGTTTTT", "TGCATTGGCCAATTTGGGCCCAAA"]
	debruijn = fdgb.GraphData(reads, 15, reduce_data=True, simplify_graph=True, remove_tips=True, construct_labels=False, directed_reads=True, verbose=False)
	
	debruijn.print_all_sequences()
	
	debruijn.remove_single_sequence_components()
	debruijn.reduce_to_single_largest_component()
	
	debruijn.print_all_sequences()

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
	#debruijn.reduce_to_single_largest_component()
	#debruijn.remove_tips()
	#debruijn.construct_assembly_ordering_labels(verbose = False)
	
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	
	debruijn.reduce_to_single_path_max_weight(verbose = False)
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.construct_assembly_ordering_labels(verbose = False)
	filename_output += "_singlepath"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")

def test_assembly_ordering_cycle(data_dir="Output/test"):
	if not os.path.exists(data_dir):
		os.mkdir(data_dir)
	dnastring ="TAATGGGGAATGTAGCCGTCAGTCTTAGACTGTAAAGCTAAAACACACAGAAGAATCTTTTTCAAAGCCGAACTCACTGTCGGCGGTTAGGTTGCCGCGTCTTACCAAGTTGCTGTCTACACGCGCACAATACCGAGTGATCATGGCCGGATTCGCATGCCAACTTCACTAGCGCATGAGTAACATGATGGCGTAGCGTGGGCGAAGCACTCCGTGGCGGGACGGGCACCATGGCACATAAGAGAAACGGACATCGGAAGCAGGACCCCCTAATGGGGAATGTAGCCGTCAGTCTTAGAC"
	dna = [c for c in dnastring]
	reads = dgen.samplereads(dna, 100, replace_error_percentage=0.5, read_length_mean=100)
	debruijn = fdgb.GraphData([reads], k=19, directed_reads=True, load_weights=False, reduce_data=True, simplify_graph=True, construct_labels=False, remove_tips=True)
	
	debruijn.get_asqg_output(filename = data_dir+"/test_labels_cycle_1.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_labels_cycle_1.csv")
	
	debruijn.construct_assembly_ordering_labels(do_second_iteration = False, verbose=2)
	
	debruijn.get_asqg_output(filename = data_dir+"/test_labels_cycle_2.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_labels_cycle_2.csv")

def test_spectral_partitioning():
	#readpartition = dio.get_read_partition_by_readlength(filename = "Data/2017-09-05_coronavirus.fq", size_of_parts=500)
	#rp = [r[1] for r in readpartition[0]]
	reads = dio.get_reads_from_fastq_file(filename = "Data/2017-09-05_coronavirus.fq", num_of_reads = 500)
	k = 25
	
	print (len(reads))

	debruijn = fdgb.GraphData([reads], k, alphabet={"A":"U", "C":"G", "G":"C", "U":"A"})
	
	debruijn.reduce_to_single_largest_component()
	debruijn.construct_assembly_ordering_labels(verbose = False)
	
	debruijn.get_asqg_output(filename = "Output/test/mincuttest_precut.asqg")
	
	#c = debruijn.get_components()
	#print ("number of components: "+str(len(c)))
	debruijn.compute_mincut()
	
	debruijn.get_asqg_output(filename = "Output/test/mincuttest_postcut.asqg")
	
	debruijn.reduce_every_component_to_single_path_max_weight()
	debruijn.contract_unique_overlaps(verbose = False)
	debruijn.get_asqg_output(filename = "Output/test/mincuttest_postcut_postreduce.asqg")
	debruijn.get_csv_output(filename = "Output/test/mincuttest_postcut_postreduce.csv")
	debruijn.write_sequences_to_file(filename = "Output/test/mincuttest_postcut_postreduce_sequences.txt")
	
def minimal_test_spectral_partitioning():
	#'''
	dna_1 = dio.genereate_dna(length=1000)
	dna_2 = dio.genereate_dna(length=1000)
	dna_3 = dio.genereate_dna(length=100)
	dna_1+=dna_2[-25:]+dna_3
	dio.write_dna_to_file("Output/test/genome_dna_test_1.txt", dna_1)
	#dio.write_dna_to_file("Output/test/genome_dna_test_2.txt", dna_2)
	#if not os.path.isfile("Output/test/testreads.txt"):
	#sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads_1.txt", read_length=100, length_stddev=0, set_of_viruses=["dna_test_1"], number_of_reads=[300], replace_error_percentage=1.0, indel_error_percentage=0.0, inverted_reads=False)
	sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads_1.txt", read_length=300, length_stddev=0, set_of_viruses=["dna_test_1"], number_of_reads=[1000], replace_error_percentage=3.0, indel_error_percentage=0.0, inverted_reads=False)
	sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads_2.txt", read_length=300, length_stddev=0, set_of_viruses=["dna_test_2"], number_of_reads=[1000], replace_error_percentage=3.0, indel_error_percentage=0.0, inverted_reads=False)
	#'''
	
	reads_1 = dio.get_reads_from_file("Output/test/testreads_1.txt")
	reads_2 = dio.get_reads_from_file("Output/test/testreads_2.txt")
	#reads_2 = []
	k = 25
	
	debruijn = fdgb.GraphData(reads_1 + reads_2, k)
	
	#debruijn.reduce_to_single_largest_component()
	#debruijn.construct_assembly_ordering_labels(verbose = False)
	
	debruijn.get_asqg_output(filename = "Output/test/mincuttest_precut.asqg")
	debruijn.partition_graph_into_components_of_clusters(verbose=False)
	
	debruijn.get_asqg_output(filename = "Output/test/mincuttest_postcut.asqg")
	
def test_laplacian_construction():
	'''
	dna_small = dio.genereate_dna(length=100)
	dio.write_dna_to_file("Output/test/genome_dna_test_small.txt", dna_small)
	sr.samplereads(input_filedir="Output/test/", output_filename="Output/test/testreads_small.txt", read_length=30, length_stddev=0, set_of_viruses=["dna_test_small"], number_of_reads=[100], replace_error_percentage=0.2, indel_error_percentage=0.0, inverted_reads=False)
	'''
	reads = dio.get_reads_from_file("Output/test/testreads_small.txt")
	k = 15
	
	debruijn = fdgb.GraphData(reads, k)
	# delete reads to save ram:
	reads = []
	
	debruijn.get_asqg_output(filename = "Output/test/laplaciantest_precut.asqg")
	
	c = debruijn.get_components()
	print ("number of components: "+str(len(c)))
	debruijn.compute_mincut()
	
	debruijn.get_asqg_output(filename = "Output/test/laplaciantest_postcut.asqg")
	
def test_clustercut_on_quasispecies(number_of_base_dnas=1, dna_length=5000, number_of_variations=5, num_reads_per_dna=2000):
	filename_output = "Output/test/test_clustercut_on_quasispecies"

	dna_partsize_mean = dna_length/4
	dna_partsize_variation = dna_length/100
	mean_readlength = 100

	dnas = []
	reads = []
	for dna_i in range(number_of_base_dnas):
		dnas += dgen.generate_set_of_related_dnas(length=dna_length, mean_partsize=dna_partsize_mean, variation_partsize=dna_partsize_variation, number_of_dnas=number_of_variations)
		
	print ("Number of generated dnas: "+str(len(dnas)))
	
	for dna_i in range(len(dnas)):
		meta.print_progress(dna_i, len(dnas)-1, front_string="Generate reads. Progress: ")
		filename = "Output/test/test_clustercut_on_quasispecies_dna_"+str(dna_i)
		dgen.write_dna_to_file(filename, dnas[dna_i])
		reads += dgen.samplereads(dna=dnas[dna_i], number_of_reads=num_reads_per_dna, replace_error_percentage=0.0, indel_error_percentage=15.0, mutation_alphabet=["A","C","G","T"], read_length_mean=mean_readlength, read_length_stddev=0, readlength_distribution='exponential')
		
	print ("Number of generated reads: "+str(len(reads)))
	
	k = 30
	debruijn = fdgb.GraphData([reads], k, remove_tips=True)
    
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt", addweights=True)
	
	debruijn.partition_graph_into_components_of_clusters(verbose=False)
	filename_output += "_divided"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	
	debruijn.reduce_every_component_to_single_path_max_weight(verbose=False)
	filename_output += "_singlepath"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt")

def test_exponential_readlengths():
	dna = dgen.generate_dna(length=30000)
	reads = dgen.samplereads(dna, number_of_reads=10000, read_length_mean=1000, readlength_distribution='exponential')
	meta.get_readlength_distribution(reads, 200)
	
def test_hubreads():
	outdir = "./Output/test/hubreads"
	
	if not os.path.exists(outdir):
		os.mkdir(outdir)

	dna = dgen.generate_dna(length=300)
	reads = dgen.samplereads(dna, number_of_reads=150, replace_error_percentage=0.1, indel_error_percentage=0.0, mutation_alphabet=["A","C","G","T"], read_length_mean=50, read_length_stddev=0, readlength_distribution='gaussian')
	
	debruijn_master = fdgb.GraphData([reads], k=17, directed_reads=True, load_weights=False, reduce_data=False, simplify_graph=True, construct_labels=False, remove_tips=False)
	
	debruijn_master.get_asqg_output(filename = outdir+"/base.asqg")
	debruijn_master.get_csv_output(filename = outdir+"/base.csv")
	debruijn_master.write_sequences_to_file(filename = outdir+"/base_seqs.txt", asfasta = False)
	
	hubreads = debruijn_master.get_hubreads_by_adjacent_sequences(verbose=False)
	
	dio.write_reads_to_file(hubreads, outdir+"/hubreads.txt")
	
	debruijn_hubreads_basic = fdgb.GraphData([hubreads], k=17, directed_reads=True, load_weights=False, reduce_data=False, simplify_graph=True, construct_labels=False, remove_tips=False)
	#debruijn_hubreads_basic.print_all_reads()
	debruijn_hubreads_basic.get_asqg_output(filename = outdir+"/hubreads_basic.asqg")
	debruijn_hubreads_basic.get_csv_output(filename = outdir+"/hubreads_basic.csv")
	
	debruijn_hubreads_weights = fdgb.GraphData([hubreads], k=17, directed_reads=True, load_weights=True, reduce_data=False, simplify_graph=True, construct_labels=False, remove_tips=False)
	#debruijn_hubreads_weights.print_all_reads()
	debruijn_hubreads_weights.get_asqg_output(filename = outdir+"/hubreads_weights.asqg")
	debruijn_hubreads_weights.get_csv_output(filename = outdir+"/hubreads_weights.csv")
	
def test_hubpaths():
	dna = dgen.generate_dna(length=3000)
	reads = dgen.samplereads(dna, number_of_reads=1000, replace_error_percentage=5.0, indel_error_percentage=0.0, mutation_alphabet=["A","C","G","T"], read_length_mean=100, read_length_stddev=5, readlength_distribution='gaussian')
	k = 15
	filename_output_base = "Output/test/hubpath/test_hubpath"
	
	for i in range(10):
		filename_output = filename_output_base + "_" + str(k)# + "_" + str(i)
		if i > 0:
			debruijn = fdgb.GraphData([hubpaths], k, remove_tips = True	, construct_labels = False)
			debruijn.remove_tips(only_simply_connected_tips=False)
			#debruijn.get_asqg_output(filename = filename_output+"_hp-based.asqg")
		
		else:
			debruijn = fdgb.GraphData([reads], k, remove_tips = True, construct_labels = False)
		hubpaths = debruijn.get_hubpath_sequences_by_adjacent_sequences(verbose = False)
		debruijn.get_asqg_output(filename = filename_output+".asqg")#"_orig.asqg")
		k += 1
		
		#print ("Number of overlaps in "+filename_output+"_orig: "+str(len(debruijn.overlaps)))
		#print ("Number of hubpaths in "+filename_output+"_orig: "+str(len(hubpaths)))

def test_position_labels():
	dna = dgen.generate_dna(length=1000)
	reads = dgen.samplereads(dna, number_of_reads=1000, replace_error_percentage=0.0, indel_error_percentage=0.1, mutation_alphabet=["A","C","G","T"], read_length_mean=100, read_length_stddev=5, readlength_distribution='gaussian')
	k = 20
	debruijn = fdgb.GraphData([reads], k, construct_labels=False)
	print len(debruijn.sequences)
	debruijn.construct_assembly_ordering_labels(compute_position_labels=False, verbose=2);

def test_multisized_debruijn_graph(dna_length = 5000, n_reads = 1000, readlength = 1000, indel_error = 15.0, min_k = 15, max_k = 25, verbose=False):
	dna = dgen.generate_dna(length=dna_length)
	reads = dgen.samplereads(dna, number_of_reads=n_reads, replace_error_percentage=0.0, indel_error_percentage=indel_error, mutation_alphabet=["A","C","G","T"], read_length_mean=readlength, read_length_stddev=5, readlength_distribution='exponentioal')
	filename_output_base = "/hubreadbased_dna_"+str(dna_length)+"_nr"+str(n_reads)+"_nl"+str(readlength)+"_ir"+str(indel_error)+"_mink"+str(min_k)+"_maxk"+str(max_k)
	
	data_dir_base = "Output/test/multisized"
	if not os.path.exists(data_dir_base):
		os.mkdir(data_dir_base)
		
	dio.write_sequences_to_fasta([''.join(dna)], data_dir_base+"/dna.fasta")
	filename_output_base = data_dir_base + filename_output_base
	
	hubreads = []
	remaining_reads = reads
	
	for k in range(min_k, max_k+1):
		print ("Construct hubreads for k="+str(k))
		this_iterations_reads = hubreads + remaining_reads
		debruijn = fdgb.GraphData([this_iterations_reads], k, reduce_data=False, remove_tips = True, construct_labels = False)
		debruijn.get_asqg_output(filename = filename_output_base+"_k"+str(k)+".asqg")
		
		hubreads = debruijn.get_hubreads_by_adjacent_sequences(verbose = False)
		#hubreads += new_hubreads
		
		debruijn.remove_insignificant_overlaps(minimal_evidence=2)
		remaining_reads_ids = debruijn.get_relevant_reads(verbose = True)
		remaining_reads = [this_iterations_reads[i] for i in remaining_reads_ids]
		
		debruijn.reduce_every_component_to_single_path_max_weight(verbose = False)
		debruijn.remove_single_sequence_loops()
		debruijn.reduce_to_single_largest_component()
		debruijn.get_asqg_output(filename = filename_output_base+"_k"+str(k)+"_reduced.asqg")
		debruijn.get_csv_output(filename = filename_output_base+"_k"+str(k)+"_reduced.csv")
		debruijn.write_sequences_to_file(filename = filename_output_base+"_k"+str(k)+"s_reduced.fasta", asfasta = True)
		if verbose:
			print hubreads	
	
	print ("Construct debruijn-graph from hubreads:")
	debruijn = fdgb.GraphData([hubreads], max_k, remove_tips = True)
	debruijn.get_asqg_output(filename = filename_output_base+"_hubreads.asqg")
	debruijn.remove_insignificant_overlaps(minimal_evidence=2)
	debruijn.remove_single_sequence_loops()
	debruijn.construct_assembly_ordering_labels()
	debruijn.reduce_every_component_to_single_path_max_weight(verbose = False)
	debruijn.get_asqg_output(filename = filename_output_base+"_hubreads_reduced.asqg")
	debruijn.get_csv_output(filename = filename_output_base+"_hubreads_reduced.csv")
	debruijn.write_sequences_to_file(filename = filename_output_base+"_hubreads_reduced.fasta", asfasta = True)

def test_multisized_on_corona_reads():
	read_id_partition_by_size = dio.get_read_partition_by_readlength(filename="Data/hcov229e_only.fq", size_of_parts=100)
	reads = dio.get_reads_from_fastq_file(read_ids = read_id_partition_by_size[-1], filename="Data/hcov229e_only.fq")
	filename_output_base = "Output/test/multisized_corona/cm_largest100_"
	
	min_k = 40
	max_k = 50
	
	hubreads = []
	remaining_reads = reads
	
	for k in range(min_k, max_k+1):
		print ("Construct hubreads for k="+str(k))
		this_iterations_reads = hubreads + remaining_reads
		debruijn = fdgb.GraphData([this_iterations_reads], k, reduce_data=False, remove_tips = True, construct_labels = False)
		debruijn.get_asqg_output(filename = filename_output_base+"_k"+str(k)+".asqg")
		
		hubreads = debruijn.get_hubreads_by_adjacent_sequences(verbose = False)
		
		debruijn.remove_insignificant_overlaps(minimal_evidence=2)
		remaining_reads_ids = debruijn.get_relevant_reads(verbose = True)
		remaining_reads = [this_iterations_reads[i] for i in remaining_reads_ids]
		
		debruijn.reduce_every_component_to_single_path_max_weight(verbose = False)
		debruijn.remove_single_sequence_loops()
		debruijn.reduce_to_single_largest_component()
		debruijn.get_asqg_output(filename = filename_output_base+"_k"+str(k)+"_reduced.asqg")
		debruijn.get_csv_output(filename = filename_output_base+"_k"+str(k)+"_reduced.csv")
		
	print ("Construct debruijn-graph from hubreads:")
	debruijn = fdgb.GraphData([hubreads], max_k, remove_tips = True)
	debruijn.get_asqg_output(filename = filename_output_base+"_hubreads.asqg")
	debruijn.remove_insignificant_overlaps(minimal_evidence=2)
	debruijn.remove_single_sequence_loops()
	debruijn.construct_assembly_ordering_labels()
	debruijn.reduce_every_component_to_single_path_max_weight(verbose = False)
	debruijn.get_asqg_output(filename = filename_output_base+"_hubreads_reduced.asqg")
	debruijn.get_csv_output(filename = filename_output_base+"_hubreads_reduced.csv")
	debruijn.write_sequences_to_file(filename = filename_output_base+"_hubreads_reduced.fasta", asfasta = True)
	
def test_read_parttition_by_length_distribution():
	read_id_partition_by_distribution = dio.get_read_partition_by_lengthdistribution(filename="Data/hcov229e_only.fq", size_of_parts=500, verbose=True)
	reads_1 = dio.get_reads_from_fastq_file(read_ids = [r[1] for r in read_id_partition_by_distribution[-1]], filename="Data/hcov229e_only.fq")
	print ("first set of reads: total number of reads: "+str(len(reads_1)))
	print ("lengths of reads:")
	print (sorted([len(r) for r in reads_1]))
	reads_2 = dio.get_reads_from_fastq_file(read_ids = [r[1] for r in read_id_partition_by_distribution[-2]], filename="Data/hcov229e_only.fq")
	print ("second set of reads: total number of reads: "+str(len(reads_2)))
	print ("lengths of reads:")
	print (sorted([len(r) for r in reads_2]))
	
def test_removal_of_insignificant_nodes_and_edges():
	filename_output_base = "Output/test/test_deconvolute_"
	
	dna = dgen.generate_dna(2000)
	
	k = 15
	readlength = 300
	num_of_reads = 1000
	error_percentage = 15.0
	ep_string = "".join(re.split(r'\.',str("%2.2f" % error_percentage)))
	casename = "_rl"+str(readlength)+"_nr"+str(num_of_reads)+"_er"+ep_string+"_k"+str(k)
	
	reads = dgen.samplereads(dna, number_of_reads=num_of_reads, replace_error_percentage=error_percentage, read_length_mean=readlength)
	debruijn = fdgb.GraphData([reads], k, directed_reads = True, load_weights=False, construct_labels=False, remove_tips=True)
	
	debruijn.get_asqg_output(filename = filename_output_base+casename+"_base.asqg")
	debruijn.get_csv_output(filename = filename_output_base+casename+"_base.csv")
	
	debruijn2 = copy.deepcopy(debruijn)
	
	debruijn.remove_insignificant_overlaps()
	debruijn.get_asqg_output(filename = filename_output_base+casename+"_r-ov.asqg")
	debruijn.get_csv_output(filename = filename_output_base+casename+"_r-ov.csv")
	
	debruijn.remove_single_sequence_components()
	debruijn.get_asqg_output(filename = filename_output_base+casename+"_r-ov_r-ssc.asqg")
	debruijn.get_csv_output(filename = filename_output_base+casename+"_r-ov_r-ssc.csv")
	
	debruijn2.remove_insignificant_sequences()
	debruijn2.get_asqg_output(filename = filename_output_base+casename+"_r-s.asqg")
	debruijn.get_csv_output(filename = filename_output_base+casename+"_r-s.csv")
	
	debruijn2.remove_single_sequence_components()
	debruijn2.get_asqg_output(filename = filename_output_base+casename+"_r-s_r-ssc.asqg")
	debruijn.get_csv_output(filename = filename_output_base+casename+"_r-s_r-ssc.csv")
	
def compare_different_read_partitions():
	k = 50
	filename_output_base = "Output/test/compare_read_partitions_"
	read_id_partition_by_size = dio.get_read_partition_by_readlength(filename="Data/hcov229e_only.fq", size_of_parts=50)
	read_id_partition_by_distribution = dio.get_read_partition_by_lengthdistribution(filename="Data/hcov229e_only.fq", size_of_parts=500, verbose=True)
	
	read_parts_by_size_ids = [x[1] for x in read_id_partition_by_size[-1]]
	reads_by_size = dio.get_reads_from_fastq_file(filename = "Data/hcov229e_only.fq", read_ids = read_parts_by_size_ids)
	debruijn_by_size = fdgb.GraphData([reads_by_size], k)
	debruijn_by_size.get_asqg_output(filename = filename_output_base+"by_size.asqg")
	
	read_parts_by_distribution_ids = [x[1] for x in read_id_partition_by_distribution[-1]]
	reads_by_distribution = dio.get_reads_from_fastq_file(filename = "Data/hcov229e_only.fq", read_ids = read_parts_by_distribution_ids)
	debruijn_by_distribution = fdgb.GraphData([reads_by_distribution], k)
	debruijn_by_distribution.get_asqg_output(filename = filename_output_base+"by_distribution.asqg")
	
def test_overlap_weight_distribution():
	dna = dgen.generate_dna(length=5000)
	
	reads = dgen.samplereads(dna, number_of_reads=1000, replace_error_percentage=2.0, read_length_mean=100)

	debruijn = fdgb.GraphData([reads], 17, reduce_data=True, simplify_graph=True, remove_tips=False, construct_labels=True, directed_reads=True)
	print debruijn.compute_overlap_evidence_distribution()
	
	debruijn.remove_tips()
	print debruijn.compute_overlap_evidence_distribution()

def test_check_if_graph_decomposes():
	data_dir="Output/test"
	dna = dgen.generate_dna(length=5000)
	reads = dgen.samplereads(dna, number_of_reads=1000, indel_error_percentage=15.0, read_length_mean=1000)
	debruijn = fdgb.GraphData([reads], 15, reduce_data=True, simplify_graph=True, remove_tips=True, construct_labels=False, directed_reads=True)
	
	debruijn.get_asqg_output(filename = data_dir+"/test_check_decomp_1.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_check_decomp_1.csv")
	
	print debruijn.compute_overlap_evidence_distribution()
	debruijn.remove_insignificant_overlaps(2)
	debruijn.remove_tips()
	debruijn.contract_unique_overlaps()
	debruijn.remove_single_sequence_components()
	
	debruijn.get_asqg_output(filename = data_dir+"/test_check_decomp_2.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_check_decomp_2.csv")
	
	print debruijn.compute_overlap_evidence_distribution()
	debruijn.remove_low_evidence_overlaps_until_graph_decomposes(verbose=False)
	
	debruijn.get_asqg_output(filename = data_dir+"/test_check_decomp_3.asqg")
	debruijn.get_csv_output(filename = data_dir+"/test_check_decomp_3.csv")
	
if __name__ == '__main__':
	#test_basic_functionality(new_dataset=True)
	#test_basic_functionality_small(new_dataset=False)
	test_basic_consensus_construction(new_dataset=True)
	
	#test_assembly_ordering_cycle()
	
	#test_clustercut_on_quasispecies(number_of_base_dnas=3, dna_length=5000, number_of_variations=1, num_reads_per_dna=5000)
	#test_exponential_readlengths()
	#test_hubpaths()
	#test_multisized_debruijn_graph(dna_length = 2000, n_reads = 100, readlength = 400, indel_error = 5.0)
	#test_multisized_debruijn_graph(dna_length = 5000, n_reads = 100, readlength = 1000, indel_error = 10.0)
	#test_multisized_debruijn_graph(dna_length = 5000, n_reads = 100, readlength = 1000, indel_error = 15.0)
	#test_multisized_on_corona_reads()
	#test_read_parttition_by_length_distribution()
	#test_position_labels()
	#compare_different_read_partitions()
	
	#test_removal_of_insignificant_nodes_and_edges()
	#test_hubreads()
	#test_check_if_graph_decomposes()
	#test_reduce_to_largest_single_sequence()