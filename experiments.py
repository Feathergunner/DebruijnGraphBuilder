#!usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import subprocess
from multiprocessing import Process
import time

import datasets_experiments as dse
import data_gen as dgen
import data_io as dio
import dbg_consensus_construction as cc

def experiment_consensus_singlecase(algorithm,
									outputdir,
									dna,
									num_reads,
									readlength,
									error_percentage,
									k,
									error_type			= "replace",
									uniform_coverage	= True,
									name_suffix			= "",
									saveparts			= True,
									logfile				= False):
	# algorithm: 		"simplecons", "locofere", "covref"
	# outputdir: 		string
	# dna: 				string
	# num_reads: 		int > 0
	# readlength: 		int > 0
	# error_percentage:	float from [0.0, 100.0]
	# k: 				int > 0
	# error_type: 		"replace", "indel"
	# uniform_coverage: boolean
	# name_suffix: 		string
	# saveparts: 		boolean
	# logfile: 			boolean
										
	# initialize:
	if not os.path.exists(outputdir):
		os.mkdir(outputdir)
	if not os.path.exists(outputdir+"/reads"):
		os.mkdir(outputdir+"/reads")
	
	ep_string = "".join(re.split(r'\.',str("%2.2f" % error_percentage)))
	if error_type == "indel":
		ep_string = "ei"+ep_string
	else:
		ep_string = "er"+ep_string
	casename =  algorithm+"_rl"+str(readlength)+"_nr"+str(num_reads)+"_"+ep_string+"_k"+str(k)

	if not name_suffix == "":
		casename += "_"+name_suffix
	if logfile:
		# redirect stdout to file:
		sys.stdout = open(outputdir+"/log_"+casename, 'w')
		
	# generate reads:
	if error_type == "indel":
		err = 0.0
		eir = error_percentage
	else:
		err = error_percentage
		eir = 0.0
	reads = dgen.samplereads(dna, num_reads, replace_error_percentage=err,  indel_error_percentage=eir, read_length_mean=readlength, uniform_coveragedepth=uniform_coverage)
	dio.write_reads_to_file(reads, outputdir+"/reads/reads_"+casename+".txt")
	
	# construct consensus:
	if algorithm == "simplecons":
		cc.simplecons(reads, k, outputdir, casename, saveparts)
	elif algorithm == "locofere":
		cc.cons_locofere(reads, k, outputdir, casename, saveparts)
	elif algorithm == "covref":
		print ("Error! Algorithm '"+algorithm+"' is not yet implemented...")
	else:
		print ("Error! Algorithm '"+algorithm+"' unknown!")
		
def create_dataset(	algorithm,
					focus,
					name,
					dimension_of_set	= 1,
					dna_length			= 5000,
					error_rates			= [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
					num_of_reads		= [2000],
					readlength			= 50,
					k_lengths			= [11,13,15,17,19,21],
					error_type			= "replace",
					uniform_coverage	= True,
					scope				= "all",
					saveparts			= True):
	# algorithm: 		"simplecons", "locofere" or "covref"
	# focus: 			"err_vs_k" or "rl_vs_k"
	# name: 			sting
	# dimension_of_set: int > 0
	# dna_length: 		int > 0
	# error_rates:		list of floats from [0.0, 100.0]
	# num_of_reads:		list of ints > 0
	# readlength: 		int > 0
	# k_lengths:		list of ints > 0
	# error_type: 		"replace" or "indel"
	# uniform_coverage: boolean
	# scope:			"all", "data" or "stat"
	# saveparts: 		boolean
	
	outputdir = "Output/"+name+"_d"+str(dimension_of_set)
	
	if scope == "all" or scope == "data":
		if os.path.exists(outputdir):
			print ("Error! Cannot create a dataset with the name '"+name+"'! There already exists a dataset of the same name.")
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
			for error_rate in error_rates:
				for k in k_lengths:
					for nr in num_of_reads:
						p = Process(target=experiment_consensus_singlecase, args=(algorithm, outputdir, dna, nr, readlength, error_rate, k, error_type, uniform_coverage, str(i+1), saveparts,  True))
						threads.append(p)
						p.start()
			
		for p in threads:
			p.join()
		
	if scope == "all" or scope == "stat":
		import apply_blast as abl
		import construct_stat_plots as csp
		
		if not os.path.exists(outputdir):
			print ("Error! Data with the name '"+name+"' does not exist!")
			return
		
		print ("Compute BLAST ratings ...")
		abl.init_blast_db(outputdir+"/dna.fasta")
		# ensure that blast finishes:
		time.sleep(2.0)
		abl.compute_blast_results(outputdir, outputdir+"/dna.fasta")
		
		# draw plots:
		print ("Construct statistic plots ...")
		if focus == "err_vs_k":
			csp.construct_heatmaps_cons_2g(datadir=outputdir, basename=algorithm, outputdir=outputdir+"/plots", number_of_reads=num_of_reads[0], k_values=k_lengths, error_rates=error_rates, dna_length=dna_length, readlength=readlength, dimension_of_set=dimension_of_set)
		elif focus == "rl_vs_k":
			csp.construct_heatmaps_cons_3g(datadir=outputdir, basename=algorithm, outputdir=outputdir+"/plots", number_of_reads=num_of_reads, k_values=k_lengths, error_rate=error_rates[0], dna_length=dna_length, readlength=readlength, dimension_of_set=dimension_of_set)
		else:
			print ("Error! Focus '"+focus+"' unknown!")
			
def create_dataset_from_setting(algorithm,
								setting,
								dimension_of_set,
								scope):
	# algorithm: 		"simplecons", "locofere" or "covref"
	# setting:			a dictionary from datasets_experiments.py
	# dimension_of_set: int > 0
	# scope:			"all", "data" or "stat"
	
	numrreads = setting["numbers_of_reads"]
	readlengths = setting["readlengths"][0]
	k_lengths = setting["kmer_lengths"]
	error_rates = setting["error_rates"]
	error_type = setting["error_type"]
	name = algorithm+"_"+setting["name"]
	
	if len(error_rates) == 1:
		focus = "rl_vs_k"
	else:
		focus = "err_vs_k"		
	
	create_dataset(algorithm, focus, name, dimension_of_set, 5000, error_rates, numrreads, readlengths, k_lengths, error_type, True, scope, True)
	
def experiments(params):
	set_id = 0
	dim = 1
	num_reads = []
	k_lengths = []
	error_rate = []
	name = "defaultexperiment"
	ucd = True
	fullset = False
	testset = 0
	only_data = False
	only_stat = False
	scope = "all"
	algo = "simplecons"
	
	for arg in params:
		arg_data = re.split(r'=', arg)
		if arg_data[0] == "nr":
			num_reads.append(int(arg_data[1]))
		elif arg_data[0] == "k":
			k_lengths.append(int(arg_data[1]))
		elif arg_data[0] == "d":
			dim = int(arg_data[1])
		elif arg_data[0] == "e":
			error_rate.append(float(arg_data[1]))
		elif arg_data[0] == "name":
			name = arg_data[1]
		elif arg_data[0] == "algo":
			if arg_data[1] == "sc":
				algo = "simplecons"
			elif arg_data[1] == "lcfr":
				algo = "locofere"
			elif arg_data[1] == "cr":
				algo = "covref"
			else:
				print ("Error! Parameter '"+arg_data[1]+"' unknown!")
			
		elif arg_data[0] == "ucd":
			ucd = True
		elif arg_data[0] == "noucd":
			ucd = False
			
		elif arg_data[0] == "set":
			set_id = int(arg_data[1])
		elif arg_data[0] == "test":
			testset = 2
		elif arg_data[0] == "testgen2":
			testset = 2
		elif arg_data[0] == "testgen3":
			testset = 3
		elif arg_data[0] == "onlydata":
			scope = "data"
		elif arg_data[0] == "onlystat":
			scope = "stat"
		else:
			print ("Error! Parameter '"+arg_data[0]+"' unkown!")
			return
			
	if set_id == 0 and testset == 0:
		if k_lengths == []:
			k_lengths = [15]
		if num_reads == []:
			num_reads = [2000]
		if error_rate == []:
			error_rate = [0.1]
			
		outputdir = "Output/"+str(name)
		if not os.path.exists(outputdir):
			os.mkdir(outputdir)
	
		if not os.path.isfile(outputdir+"/dna.txt"):
			# generate dna:
			dna = dgen.generate_dna(5000)
			# write dna to fasta:
			dio.write_genome_to_file(dna, outputdir+"/dna.txt")
			dio.write_sequences_to_fasta([''.join(dna)], outputdir+"/dna.fasta")
		else:
			dna = [c for c in dio.get_genome_from_file(outputdir+"/dna.txt")]
		for k in k_lengths:
			for nr in num_reads:
				for er in error_rate:
					experiment_consensus_singlecase(algo, outputdir, dna, nr, 50, er, k)
					
	elif testset > 0:
		if testset == 2:
			create_dataset_from_setting("simplecons", dse.testset_2g, dimension_of_set=dim, scope=scope)
		if testset == 3:
			create_dataset_from_setting("locofere", dse.testset_3g, dimension_of_set=dim, scope=scope)
	elif set_id in range(1, len(dse.allsettings)+1):
		setting = dse.allsettings[set_id-1]
		
		create_dataset_from_setting(algo, setting, dimension_of_set=dim, scope=scope)
		
if __name__ == "__main__":
	experiments(sys.argv[1:])
	