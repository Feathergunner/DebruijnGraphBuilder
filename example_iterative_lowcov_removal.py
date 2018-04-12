#!usr/bin/env python
# -*- coding: utf-8 -*-

import data_gen as dgen
import data_io as dio
import fast_debruijn_graph_builder as fdgb

import os
import re
import sys

def experiment_iterative_low_coverage_removal(outputdir, dna_length=5000, num_reads=1000, readlength=1000, error_percentage=15.0, k=13, saveparts=True, create_new_dna=False):
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
		
	# generate reads:
	reads = dgen.samplereads(dna, num_reads, indel_error_percentage=error_percentage, read_length_mean=readlength)
	
	# initilize:
	ep_string = "".join(re.split(r'\.',str("%2.2f" % error_percentage)))
	casename = "itlowcovrem"+"_rl"+str(readlength)+"_nr"+str(num_reads)+"_ei"+ep_string+"_k"+str(k)
	debruijn = fdgb.GraphData([reads], k=k, directed_reads=True, load_weights=False, reduce_data=True, simplify_graph=True, construct_labels=False, remove_tips=False)
	# basic reduction:
	debruijn.remove_tips()
	debruijn.remove_insignificant_overlaps(2, remove_only_unique_tips=True) # <- removes all overlaps with coverage 1
	debruijn.remove_tips()
	#debruijn.contract_unique_overlaps()
	debruijn.remove_single_sequence_components()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_base.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+casename+"_base.csv")
	
	# reduction step 1: remove low coverage overlaps
	debruijn.remove_low_evidence_overlaps_until_graph_decomposes()
	debruijn.reduce_to_single_largest_component()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_reduced.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+casename+"_reduced.csv")
	
	# reduction step 2: remove low coverage sequences:
	debruijn.remove_low_coverage_sequences_until_graph_decomposes()
	debruijn.reduce_to_single_largest_component()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_reduced_2.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+casename+"_reduced_2.csv")
	
	# construct a consensus sequence:
	debruijn.construct_assembly_ordering_labels()
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps()
	debruijn.get_asqg_output(filename = outputdir+"/"+casename+"_singlepath.asqg")
	debruijn.get_csv_output(filename = outputdir+"/"+casename+"_singlepath.csv")
	debruijn.write_sequences_to_file(filename = outputdir+"/"+casename+"_singlepath.fasta", asfasta = True)
	
if __name__ == "__main__":
	num_reads = 1000
	k = 15
	name = "itlowcovrem"
	for arg in sys.argv[1:]:
		arg_data = re.split(r'=', arg)
		if arg_data[0] == "nr":
			num_reads = int(arg_data[1])
		elif arg_data[0] == "k":
			k = int(arg_data[1])
		elif arg_data[0] == "name":
			name += "_"+arg_data[1]
	
	experiment_iterative_low_coverage_removal("Output/"+name, num_reads=num_reads, k=k)