#!usr/bin/env python
# -*- coding: utf-8 -*-

import fast_debruijn_graph_builder as fdgb

def simplecons(	reads,
				k,
				outputdir,
				name,
				saveparts=True):
	debruijn = fdgb.GraphData(reads, k=k, directed_reads=True, load_weights=False, reduce_data=True, simplify_graph=True, construct_labels=False)
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_1_base.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_1_base.csv")
	
	debruijn.remove_tips_simple()
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_2_posttipremoval.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_2_posttipremoval.csv")
	
	debruijn.remove_insignificant_sequences()
	debruijn.remove_single_sequence_components()
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_3_postsequenceremoval.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_3_postsequenceremoval.csv")
	
	debruijn.reduce_to_single_largest_component()
	debruijn.greedy_construct_assembly_ordering_labels()
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps()
	
	debruijn.get_asqg_output(filename = outputdir+"/"+name+"_4_singlepath.asqg")
	debruijn.get_csv_output(filename = outputdir+"/"+name+"_4_singlepath.csv")
	debruijn.write_sequences_to_file(filename = outputdir+"/"+name+"_4_singlepath.fasta", asfasta = True)

# low coverage feature removal:
def cons_locofere(	reads,
					k,
					outputdir,
					name,
					saveparts=True):
	debruijn = fdgb.GraphData(reads, k=k, directed_reads=True, load_weights=False, reduce_data=True, simplify_graph=True, construct_labels=False, remove_tips=False)
	
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_0_unsimplified.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_0_unsimplified.csv")
	
	# basic reduction:
	debruijn.remove_tips(verbose=False)
	'''
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_0b_tiprm.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_0b_tiprm.csv")
	debruijn.remove_insignificant_overlaps(2, keep_relevant_tips=True) # <- removes all overlaps with coverage 1
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_0c_ov2rm.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_0c_ov2rm.csv")
	debruijn.remove_tips(verbose=True)
	'''
	#debruijn.contract_unique_overlaps()
	debruijn.remove_single_sequence_components()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_1_base.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_1_base.csv")
	
	# reduction step 1: remove low coverage overlaps
	debruijn.remove_low_evidence_overlaps_until_graph_decomposes(relative_component_size_bound=0.01)
	debruijn.reduce_to_single_largest_component()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_2_ovred.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_2_ovred.csv")
	
	# reduction step 2: remove low coverage sequences:
	debruijn.remove_low_coverage_sequences_until_graph_decomposes(relative_component_size_bound=0.01, verbose=False)
	debruijn.reduce_to_single_largest_component()
	if saveparts:
		debruijn.get_asqg_output(filename = outputdir+"/"+name+"_3_seqred.asqg")
		debruijn.get_csv_output(filename = outputdir+"/"+name+"_3_seqred.csv")
	
	# construct a consensus sequence:
	debruijn.construct_assembly_ordering_labels()
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps()
	debruijn.get_asqg_output(filename = outputdir+"/"+name+"_4_singlepath.asqg")
	debruijn.get_csv_output(filename = outputdir+"/"+name+"_4_singlepath.csv")
	debruijn.write_sequences_to_file(filename = outputdir+"/"+name+"_4_singlepath.fasta", asfasta = True)

def cons_covref(reads,
				number_of_parts,#	=50,
				overlap,#			=10,
				k_base,#			=25,
				k_part,#			=15,
				k_merge,#			=[13,15,17],
				outputdir,
				name,
				saveparts=True):
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
			# construct de Bruijn graph of this subset:
			debruijn_part = fdgb.GraphData([part_seqreads], k_part, directed_reads=True, load_weights=True, reduce_data=True, simplify_graph=True, construct_labels=False, remove_tips=True)
								
			if saveparts:
				debruijn_part.get_asqg_output(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_base.asqg")
				debruijn_part.get_csv_output(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_base.csv")
			
			# basic reduction:
			debruijn_part.remove_insignificant_overlaps(2)
			debruijn_part.remove_tips()
			debruijn_part.contract_unique_overlaps()
			debruijn_part.remove_single_sequence_components()
			
			# further reduction:
			debruijn_part.remove_low_evidence_overlaps_until_graph_decomposes()
			debruijn_part.reduce_to_single_largest_component()
			debruijn_part.remove_low_coverage_sequences_until_graph_decomposes()
			debruijn_part.reduce_to_single_largest_component()
			
			if saveparts:
				#save graph:
				debruijn_part.get_asqg_output(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_reduced.asqg")
				debruijn_part.get_csv_output(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_reduced.csv")
				# save sequences:
				debruijn_part.write_sequences_to_file(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_sequences", addweights=True)
			
			# add sequences of this partial de Bruijn graph to set of sequences to merge
			parts_reduced_sequences += [seq.sequence+","+str(seq.get_total_weight()) for seq in debruijn_part.sequences if seq.is_relevant]
				
			# construct a consensus sequence of this de Bruijn graph:
			debruijn_part.construct_assembly_ordering_labels(verbose=1)
			debruijn_part.reduce_to_single_path_max_weight()
			debruijn_part.contract_unique_overlaps()
			
			if saveparts:
				# save graph:
				debruijn_part.get_asqg_output(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_singlepath.asqg")
				debruijn_part.get_csv_output(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_singlepath.csv")
				# save consensus sequence:
				debruijn_part.write_sequences_to_file(filename = outputdir+"/parts/"+name+"_part"+str(part_id)+"_singlepath", addweights=True)
			
			parts_consensus_sequences += debruijn_part.get_relevant_sequences()
			
	# merge once from sequences and once from consensus sequences:
	readsets = [parts_reduced_sequences, parts_consensus_sequences]
	readset_names = ["sequences", "consensus"]
	for i in range(len(readsets)):
		for km in k_merge:
			casename_merge = name+"_merge_"+readset_names[i]+"_km"+str(km)
			debruijn_merge_sequences = fdgb.GraphData([readsets[i]], km, directed_reads=True, load_weights=True, reduce_data=True, simplify_graph=True, construct_labels=False, remove_tips=True)
			debruijn_merge_sequences.remove_insignificant_overlaps(2)
			debruijn_merge_sequences.remove_tips()
			debruijn_merge_sequences.contract_unique_overlaps()
			debruijn_merge_sequences.remove_single_sequence_components()
			
			debruijn_merge_sequences.get_asqg_output(filename = outputdir+"/"+casename_merge+"_base.asqg")
			debruijn_merge_sequences.get_csv_output(filename = outputdir+"/"+casename_merge+"_base.csv")
			
			debruijn_merge_sequences.remove_low_evidence_overlaps_until_graph_decomposes()
			debruijn_merge_sequences.reduce_to_single_largest_component()
			debruijn_merge_sequences.remove_low_coverage_sequences_until_graph_decomposes()
			debruijn_merge_sequences.reduce_to_single_largest_component()
			
			debruijn_merge_sequences.get_asqg_output(filename = outputdir+"/"+casename_merge+"_reduced.asqg")
			debruijn_merge_sequences.get_csv_output(filename = outputdir+"/"+casename_merge+"_reduced.csv")
			
			debruijn_merge_sequences.construct_assembly_ordering_labels()
			debruijn_merge_sequences.reduce_to_single_path_max_weight()
			debruijn_merge_sequences.contract_unique_overlaps()
			
			debruijn_merge_sequences.get_asqg_output(filename = outputdir+"/"+casename_merge+"_singlepath.asqg")
			debruijn_merge_sequences.get_csv_output(filename = outputdir+"/"+casename_merge+"_singlepath.csv")
			debruijn_merge_sequences.write_sequences_to_file(filename = outputdir+"/"+casename_merge+"_singlepath.fasta", asfasta = True)