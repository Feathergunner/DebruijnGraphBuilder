#!usr/bin/env python
# -*- coding: utf-8 -*-

import data_gen as dgen
import data_io as dio
import fast_debruijn_graph_builder as fdgb

import os
import re
import sys

def experiment_recursive_coverage_refining(outputdir, dna_length=5000, num_reads=1000, readlength=1000, error_percentage=15.0, number_of_parts=50, overlap=10, k_base=25, k_part=15, k_merge=[13,15,17], saveparts=True, create_new_dna=False):
	if not os.path.exists(outputdir):
		os.mkdir(outputdir)
	if not os.path.exists(outputdir+"/parts"):
		os.mkdir(outputdir+"/parts")

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
	casename = "reccovref"+"_rl"+str(readlength)+"_nr"+str(num_reads)+"_ei"+ep_string+"_numparts"+str(number_of_parts)+"_overlap"+str(overlap)+"_kb"+str(k_base)+"_kp"+str(k_part)
	
	debruijn_master = fdgb.GraphData([reads], k=25, directed_reads=True, load_weights=False, reduce_data=False, simplify_graph=True, construct_labels=False, remove_tips=True)
	debruijn_master.reduce_to_single_largest_component()
	debruijn_master.construct_assembly_ordering_labels()
			
	debruijn_master.get_asqg_output(filename = outputdir+"/"+casename+"_master.asqg")
	debruijn_master.get_csv_output(filename = outputdir+"/"+casename+"_master.csv")
	
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
				debruijn_part.get_asqg_output(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_base.asqg")
				debruijn_part.get_csv_output(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_base.csv")
			
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
				debruijn_part.get_asqg_output(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_reduced.asqg")
				debruijn_part.get_csv_output(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_reduced.csv")
				# save sequences:
				debruijn_part.write_sequences_to_file(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_sequences", addweights=True)
			
			# add sequences of this partial de Bruijn graph to set of sequences to merge
			parts_reduced_sequences += [seq.sequence+","+str(seq.get_total_weight()) for seq in debruijn_part.sequences if seq.is_relevant]
				
			# construct a consensus sequence of this de Bruijn graph:
			debruijn_part.construct_assembly_ordering_labels(verbose=1)
			debruijn_part.reduce_to_single_path_max_weight()
			debruijn_part.contract_unique_overlaps()
			
			if saveparts:
				# save graph:
				debruijn_part.get_asqg_output(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_singlepath.asqg")
				debruijn_part.get_csv_output(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_singlepath.csv")
				# save consensus sequence:
				debruijn_part.write_sequences_to_file(filename = outputdir+"/parts/"+casename+"_part"+str(part_id)+"_singlepath", addweights=True)
			
			parts_consensus_sequences += debruijn_part.get_relevant_sequences()
			
	# merge once from sequences and once from consensus sequences:
	readsets = [parts_reduced_sequences, parts_consensus_sequences]
	readset_names = ["sequences", "consensus"]
	for i in range(len(readsets)):
		for km in k_merge:
			casename_merge = casename+"_merge_"+readset_names[i]+"_km"+str(km)
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
	
if __name__ == "__main__":
	num_reads = 1000
	num_parts = 50
	overlap = 10
	name = "reccovref"
	for arg in sys.argv[1:]:
		arg_data = re.split(r'=', arg)
		if arg_data[0] == "nr":
			num_reads = int(arg_data[1])
		elif arg_data[0] == "np":
			num_parts = int(arg_data[1])
		elif arg_data[0] == "ov":
			overlap = int(arg_data[1])
		elif arg_data[0] == "name":
			name += "_"+arg_data[1]
	
	experiment_recursive_coverage_refining("Output/"+name, num_reads=num_reads, number_of_parts=num_parts, overlap=overlap)