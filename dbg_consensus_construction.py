#!usr/bin/env python
# -*- coding: utf-8 -*-

import os

import fast_debruijn_graph_builder as fdgb

def simplecons(	reads,
				k,
				name		= "",
				outputdir	= ".",
				saveparts	= True,
				saveresult	= True,
				verbose		= False):
	debruijn = fdgb.GraphData(reads, k=k, directed_reads=True, load_weights=False, reduce_data=True, simplify_graph=True, construct_labels=False)
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_1_base.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_1_base.csv")
	
	debruijn.remove_tips_simple(verbose=verbose)
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_2_posttipremoval.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_2_posttipremoval.csv")
	
	debruijn.remove_insignificant_sequences(verbose=verbose)
	debruijn.remove_single_sequence_components(verbose=verbose)
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_3_postsequenceremoval.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_3_postsequenceremoval.csv")
	
	debruijn.reduce_to_single_largest_component(verbose=verbose)
	debruijn.greedy_construct_assembly_ordering_labels(verbose=verbose)
	debruijn.reduce_to_single_path_max_weight(verbose=verbose)
	debruijn.contract_unique_overlaps()
	
	if saveresult:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_4_singlepath.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_4_singlepath.csv")
		debruijn.write_sequences_to_file(filename = outputdir+"/"+name+"_4_singlepath.fasta", asfasta = True)
	
	return debruijn.get_relevant_sequences()[0]

# low coverage feature removal:
def cons_locofere(	reads,
					k,
					name			= "",
					outputdir		= ".",
					weightedreads	= False,
					saveparts 		= True,
					saveresult		= True,
					verbose			= False):
	debruijn = fdgb.GraphData(reads, k=k, directed_reads=True, load_weights=weightedreads, reduce_data=True, simplify_graph=True, construct_labels=False, remove_tips=False)
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_0_unsimplified.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_0_unsimplified.csv")
	
	# basic reduction:
	debruijn.remove_tips()
	
	debruijn.remove_single_sequence_components()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_1_base.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_1_base.csv")
	
	# reduction step 1: remove low coverage overlaps
	debruijn.remove_low_evidence_overlaps_until_graph_decomposes(relative_component_size_bound=0.01, verbose=verbose)
	debruijn.reduce_to_single_largest_component(verbose=verbose)
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_2_ovred.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_2_ovred.csv")
	
	# reduction step 2: remove low coverage sequences:
	debruijn.remove_low_coverage_sequences_until_graph_decomposes(relative_component_size_bound=0.01, verbose=verbose)
	debruijn.reduce_to_single_largest_component(verbose=verbose)
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_3_seqred.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_3_seqred.csv")
	
	# construct a consensus sequence:
	debruijn.construct_assembly_ordering_labels(verbose=verbose)
	debruijn.reduce_to_single_path_max_weight(verbose=verbose)
	debruijn.contract_unique_overlaps()
	
	if saveresult:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_4_singlepath.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_4_singlepath.csv")
		debruijn.write_sequences_to_file(filename = outputdir+"/"+name+"_4_singlepath.fasta", asfasta = True)

	return debruijn.get_relevant_sequences()[0]
	
def cons_covref(reads,
				number_of_parts,#	=50,
				overlap,#			=10,
				k_base,#			=25,
				k_part,#			=15,
				k_merge,#			=[13,15,17],
				outputdir,
				name,
				saveparts	= True,
				verbose   	= False):
	
	if saveparts and not os.path.exists(outputdir+"/parts"):
		os.mkdir(outputdir+"/parts")
				
	debruijn_master = fdgb.GraphData(reads, k=k_base, directed_reads=True, load_weights=False, reduce_data=False, simplify_graph=True, construct_labels=False, remove_tips=True)
	debruijn_master.reduce_to_single_largest_component()
	debruijn_master.construct_assembly_ordering_labels()
			
	debruijn_master.get_asqg_output(filename = outputdir+"/"+name+"_master.asqg")
	debruijn_master.get_csv_output(filename = outputdir+"/"+name+"_master.csv")
	
	parts = debruijn_master.get_partition_of_sequences(number_of_parts, overlap=overlap)
	
	parts_consensus_sequences = []
	parts_reduced_sequences = []
	for part_id in range(number_of_parts):
		# get sequences with their weights from this subset as reads:
		part_seqreads = [seq.sequence+","+str(seq.get_total_weight()) for seq in parts[part_id]]
		
		if len(part_seqreads) > 0:
			# construct consensus of this subset:			
			parts_consensus_sequences.append(cons_locofere(part_seqreads, k_part, name = name+"_part"+str(part_id), outputdir = outputdir+"/parts", weightedreads=True, saveparts = True, saveresult = True, verbose = False))
			
	# merge once from consensus sequences:
	for km in k_merge:
		casename_merge = name+"_km"+str(km)
		cons_locofere(parts_consensus_sequences, km, name=casename_merge, outputdir=outputdir, weightedreads=True, saveparts=True, saveresult=True, verbose=False)