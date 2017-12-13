#!usr/bin/python
# -*- coding: utf-8 -*-

import meta
import data_gen as dgen
import data_io as dio
import fast_debruijn_graph_builder as fdgb

import gc

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
		reads += dgen.samplereads(dna=dnas[dna_i], number_of_reads=num_reads_per_dna, replace_error_percentage=0.0, indel_error_percentage=5.0, mutation_alphabet=["A","C","G","T"], read_length_mean=mean_readlength, read_length_stddev=0, readlength_distribution='exponential')
		
	print ("Number of generated reads: "+str(len(reads)))
	
	k = 30
	debruijn = fdgb.GraphData([reads], k, remove_tips=True)
    
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt", addweights=True)
	
	debruijn.partition_graph_into_components_of_clusters(verbose=True)
	filename_output += "_divided"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	
	debruijn.reduce_every_component_to_single_path_max_weight(verbose=True)
	filename_output += "_singlepath"
	debruijn.get_asqg_output(filename = filename_output+".asqg")
	debruijn.get_csv_output(filename = filename_output+".csv")
	debruijn.write_sequences_to_file(filename = filename_output+"_sequences.txt")

def test_exponential_readlengths():
	dna = dgen.generate_dna(length=30000)
	reads = dgen.samplereads(dna, number_of_reads=10000, read_length_mean=1000, readlength_distribution='exponential')
	meta.get_readlength_distribution(reads, 200)

if __name__ == '__main__':
	test_clustercut_on_quasispecies(number_of_base_dnas=5, dna_length=10000, number_of_variations=1, num_reads_per_dna=3000)
	#test_exponential_readlengths()