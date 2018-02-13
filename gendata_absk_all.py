#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
import gc

import sampleReads
import data_io as dio
import fast_debruijn_graph_builder as fdgb
import dataset_settings as ds

sampleReads.read_genomes()

def gendata(setting, onlyreads=False, reconstruct=False, clustercut=False):
	k_absolute_settings = setting["k_absolute_settings"]
	coverage_factors = setting["coverage_factors"]
	readlength_settings = setting["readlength_settings"]
	number_of_reads_settings = setting["number_of_reads_settings"]
	error_percentages = setting["error_percentages"]
	num_different_viruses = setting["num_different_viruses"]
	set_of_viruses = setting["set_of_viruses"]
	output_dir = "Output/"+setting["name"]
	name = setting["name"]
	
	if setting["error_type"] == "replace":
		name += "-r"
	else:
		name += "-i"
					
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	if not os.path.exists(output_dir+"/reads"):
		os.makedirs(output_dir+"/reads")
	
	for cf in coverage_factors:
		for i in range(len(readlength_settings)):
			readlength = readlength_settings[i]
			num_reads = [cf*number_of_reads_settings[i]]*num_different_viruses#[cf*nr for nr in number_of_reads_settings[i]]
			
			for error_percentage in error_percentages:
				ep_string = "-".join(re.split(r'\.',str(error_percentage)))
				casename_gen = name+"_"+str(readlength)+"_"+str(num_different_viruses)+"_["
				for n_r_i in range(len(num_reads)-1):
					casename_gen += str(num_reads[n_r_i])+","
				casename_gen += str(num_reads[-1])+"]_"+ep_string
				#readfilename = "Data/reads_"+casename_gen+".txt"
				readfilename = output_dir + "/reads/" + casename_gen + ".txt"
				
				if setting["error_type"] == "replace":
					epr = error_percentage
					epi = 0
				else:
					epr = 0
					epi = error_percentage
				
				if not os.path.isfile(readfilename):
					sampleReads.samplereads(output_filename = readfilename, read_length = readlength, set_of_viruses = set_of_viruses[:num_different_viruses], number_of_reads = num_reads,	replace_error_percentage = epr, indel_error_percentage = epi)
				else:
					print ("Reads already exist!")

				if not onlyreads:
					for k in k_absolute_settings:
						# run garbage collector:
						gc.collect()
						
						casename = casename_gen + "_"+str(k)
						
						if not os.path.isfile(output_dir+"/"+casename+".asqg"):
							print ("Working on case "+casename)
							try:
								reads = dio.get_reads_from_file(filename = readfilename)
								debruijn = fdgb.GraphData(reads, k, directed_reads=True, load_weights=False)
								
								debruijn.get_asqg_output(filename = output_dir+"/"+casename+".asqg")
								debruijn.get_csv_output(filename = output_dir+"/"+casename+".csv")
								
								if clustercut:
									debruijn.partition_graph_into_components_of_clusters(verbose=True)
									casename += "_divided"
									debruijn.get_asqg_output(filename = output_dir+"/"+casename++".asqg")
									debruijn.get_csv_output(filename = output_dir+"/"+casename++".csv")
								
								if reconstruct:
									debruijn.remove_insignificant_sequences()
									debruijn.remove_single_sequence_components()
									debruijn.construct_assembly_ordering_labels(verbose = True)
									debruijn.reduce_to_single_path_max_weight()
									debruijn.contract_unique_overlaps(verbose = False)
									casename+="_reconstructed"
									debruijn.get_asqg_output(filename = output_dir+"/"+casename+".asqg")
									debruijn.get_csv_output(filename = output_dir+"/"+casename+".csv")
							
							except:
								pass
							
						else:
							print ("Data on case "+casename+" already exists!")

#gendata(ds.bvdv_absk_1)
#gendata(ds.bvdv_absk_2)
#gendata(ds.bvdv_absk_4)

#gendata(ds.reads_for_sebastian_corona, onlyreads=True)
gendata(ds.bvdv_large_absk_2_clustercut, clustercut=True)
