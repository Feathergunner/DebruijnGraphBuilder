#!usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import subprocess
from multiprocessing import Process
import time

import data_gen as dgen
import data_io as dio
import fast_debruijn_graph_builder as fdgb
import apply_blast as abl
import construct_stat_plots as csp

def experiment_iterative_low_coverage_removal_singlecase(	outputdir,
															dna,
															num_reads				= 1000,
															readlength				= 1000,
															error_percentage		= 15.0,
															k						= 13,
															saveparts				= True,
															uniform_coveragedepth	= True,
															name_suffix				= "",
															logfile					= False):
	# initilize:
	if not os.path.exists(outputdir):
		os.mkdir(outputdir)
	if not os.path.exists(outputdir+"/reads"):
		os.mkdir(outputdir+"/reads")
	
	ep_string = "".join(re.split(r'\.',str("%2.2f" % error_percentage)))
	casename =  "itlowcovrem_rl"+str(readlength)+"_nr"+str(num_reads)+"_ei"+ep_string+"_k"+str(k)

	if not name_suffix == "":
		casename += "_"+name_suffix
		
	if logfile:
		sys.stdout = open(outputdir+"/log_"+casename, 'w')
		
	# generate reads:
	reads = dgen.samplereads(dna, num_reads, indel_error_percentage=error_percentage, read_length_mean=readlength, uniform_coveragedepth=uniform_coveragedepth)
	dio.write_reads_to_file(reads, outputdir+"/reads/reads_"+casename+".txt")
	
	debruijn = fdgb.GraphData([reads], k=k, directed_reads=True, load_weights=False, reduce_data=True, simplify_graph=True, construct_labels=False, remove_tips=False)
	# basic reduction:
	debruijn.remove_tips()
	debruijn.remove_insignificant_overlaps(2, keep_relevant_tips=True) # <- removes all overlaps with coverage 1
	debruijn.remove_tips()
	#debruijn.contract_unique_overlaps()
	debruijn.remove_single_sequence_components()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_1_base.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+casename+"_1_base.csv")
	
	# reduction step 1: remove low coverage overlaps
	debruijn.remove_low_evidence_overlaps_until_graph_decomposes(relative_component_size_bound=0.01)
	debruijn.reduce_to_single_largest_component()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_2_ovred.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+casename+"_2_ovred.csv")
	
	# reduction step 2: remove low coverage sequences:
	debruijn.remove_low_coverage_sequences_until_graph_decomposes(relative_component_size_bound=0.01)
	debruijn.reduce_to_single_largest_component()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_3_seqred.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+casename+"_3_seqred.csv")
	
	# construct a consensus sequence:
	debruijn.construct_assembly_ordering_labels()
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps()
	debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_4_singlepath.asqg")
	debruijn.get_csv_output(filename = outputdir+"/"+casename+"_4_singlepath.csv")
	debruijn.write_sequences_to_file(filename = outputdir+"/"+casename+"_4_singlepath.fasta", asfasta = True)

def experiment_iterative_low_coverage_removal(	outputdir,
												dna_length				= 5000,
												num_reads  				= 1000,
												readlength 		 		= 1000,
												error_percentage		= 15.0,
												k						= 13,
												saveparts				= True,
												create_new_dna			= False,
												uniform_coveragedepth	= True):
	if not os.path.exists(outputdir):
		os.mkdir(outputdir)

	if create_new_dna or not os.path.isfile(outputdir+"/dna.txt"):
		# generate dna:
		dna = dgen.generate_dna(dna_length)
		# write dna to fasta:
		dio.write_genome_to_file(dna, outputdir+"/dna.txt")
		dio.write_sequences_to_fasta([''.join(dna)], outputdir+"/dna.fasta")
	else:
		dna = [c for c in dio.get_genome_from_file(outputdir+"/dna.txt")]
	
	experiment_iterative_low_coverage_removal_singlecase(outputdir, dna, num_reads=num_reads, readlength=readlength, error_percentage=error_percentage, k=k, saveparts=saveparts, uniform_coveragedepth=uniform_coveragedepth)
	
def create_dataset(name, dimension_of_set=1, dna_length=5000, error_rate=15.0, num_of_reads=[250, 500, 750, 1000], k_lengths=[13,15,17,19]):
	outputdir = "Output/"+name	
	if os.path.exists(outputdir):
		print ("Cannot create a dataset with the name '"+name+"'! There already exists a dataset of the same name!")
		return
	os.mkdir(outputdir)
	os.mkdir(outputdir+"/reads")
	os.mkdir(outputdir+"/plots")
	
	print ("Generate dna ...")
	# generate dna:
	dna = dgen.generate_dna(dna_length)
	# write dna to fasta:
	dio.write_genome_to_file(dna, outputdir+"/dna.txt")
	dio.write_sequences_to_fasta([''.join(dna)], outputdir+"/dna.fasta")
	
	print ("Run experiments ...")
	threads = []
	for i in range(dimension_of_set):
		for nr in num_of_reads:
			for k in k_lengths:
				p = Process(target=experiment_iterative_low_coverage_removal_singlecase, args=(outputdir, dna, nr, 1000, error_rate, k, True, True, str(i+1), True))
				threads.append(p)
				p.start()
		
	for p in threads:
		p.join()
		
	print ("Compute BLAST ratings ...")
	abl.init_blast_db(outputdir+"/dna.fasta")
	# ensure that blast finishes:
	time.sleep(2.0)
	abl.compute_blast_results(outputdir, outputdir+"/dna.fasta")
	
	# draw plots:
	print ("Construct statistic plots ...")
	csp.construct_heatmaps_3g_lcfr(datadir=outputdir, basename="itlowcovrem", outputdir=outputdir+"/plots", number_of_reads=num_of_reads, k_values=k_lengths, error_rate=error_rate, dna_length=dna_length, readlength=1000, dimension_of_set=dimension_of_set)

if __name__ == "__main__":
	dim = 1
	num_reads = 1000
	k_lengths = []
	name = "itlowcovrem"
	ucd = True
	fullset = False
	testset = False
	for arg in sys.argv[1:]:
		arg_data = re.split(r'=', arg)
		if arg_data[0] == "nr":
			num_reads = int(arg_data[1])
		elif arg_data[0] == "k":
			k_lengths.append(int(arg_data[1]))
		elif arg_data[0] == "d":
			dim = int(arg_data[1])
		elif arg_data[0] == "name":
			name += "_"+arg_data[1]
		elif arg_data[0] == "ucd":
			ucd = True
		elif arg_data[0] == "noucd":
			ucd = False
		elif arg_data[0] == "set":
			fullset = True
		elif arg_data[0] == "testset":
			testset = True
			
	if k_lengths == []:
		k_lengths = [15]
	
	if not fullset and not testset:
		for k in k_lengths:
			experiment_iterative_low_coverage_removal("Output/"+name, num_reads=num_reads, k=k, uniform_coveragedepth=ucd)
	
	elif testset:
		create_dataset(name, dimension_of_set=dim, error_rate=5.0, num_of_reads=[50, 100, 150, 200])
	else:
		create_dataset(name, dimension_of_set=dim)