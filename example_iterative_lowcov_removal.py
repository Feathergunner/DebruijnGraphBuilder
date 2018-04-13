#!usr/bin/env python
# -*- coding: utf-8 -*-

import data_gen as dgen
import data_io as dio
import fast_debruijn_graph_builder as fdgb

import os
import re
import sys

def experiment_iterative_low_coverage_removal(outputdir_base, dna_length=5000, num_reads=1000, readlength=1000, error_percentage=15.0, k=13, saveparts=True, create_new_dna=False, uniform_coveragedepth=True):		
	if not os.path.exists(outputdir_base):
		os.mkdir(outputdir_base)
	dirnum = 1
	while os.path.exists(outputdir_base+"/"+str(dirnum)):
		dirnum+=1
	outputdir = outputdir_base+"/"+str(dirnum)
	os.mkdir(outputdir)
	os.mkdir(outputdir+"/reads")

	if create_new_dna or not os.path.isfile(outputdir+"/dna.txt"):
		# generate dna:
		dna = dgen.generate_dna(dna_length)
		# write dna to fasta:
		dio.write_genome_to_file(dna, outputdir+"/dna.txt")
		dio.write_sequences_to_fasta([''.join(dna)], outputdir+"/dna.fasta")
	else:
		dna = [c for c in dio.get_genome_from_file(outputdir+"/dna.txt")]
		
	# generate reads:
	reads = dgen.samplereads(dna, num_reads, indel_error_percentage=error_percentage, read_length_mean=readlength, uniform_coveragedepth=uniform_coveragedepth)
	
	# initilize:
	ep_string = "".join(re.split(r'\.',str("%2.2f" % error_percentage)))
	casename = "rl"+str(readlength)+"_nr"+str(num_reads)+"_ei"+ep_string+"_k"+str(k)
	if uniform_coveragedepth:
		casename = "itlowcovrem_ucd_"+casename
	else:
		casename = "itlowcovrem_"+casename
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
	debruijn.remove_low_evidence_overlaps_until_graph_decomposes()
	debruijn.reduce_to_single_largest_component()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_2_ovred.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+casename+"_2_ovred.csv")
	
	# reduction step 2: remove low coverage sequences:
	debruijn.remove_low_coverage_sequences_until_graph_decomposes()
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
	
if __name__ == "__main__":
	num_reads = 1000
	k_lengths = []
	name = "itlowcovrem"
	ucd = True
	for arg in sys.argv[1:]:
		arg_data = re.split(r'=', arg)
		if arg_data[0] == "nr":
			num_reads = int(arg_data[1])
		elif arg_data[0] == "k":
			k_lengths.append(int(arg_data[1]))
		elif arg_data[0] == "name":
			name += "_"+arg_data[1]
		elif arg_data[0] == "ucd":
			ucd = True
		elif arg_data[0] == "noucd":
			ucd = False
		elif arg_data[0] == "set":
			k_lengths = [13,15,17,19]
			
	if k_lengths == []:
		k_lengths = [15]
	
	for k in k_lengths:	
		experiment_iterative_low_coverage_removal("Output/"+name, num_reads=num_reads, k=k, uniform_coveragedepth=ucd)