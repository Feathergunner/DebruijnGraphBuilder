#!usr/bin/python
# -*- coding: utf-8 -*-

import os

import fast_debruijn_graph_builder as fdgb
import data_io as dio
import data_gen as dgen
#import meta

#class ConsensusConstructor:
def construct_consensus_from_multiple_parts(reads, filepath_output, filename_output_base, k_all, k_part, number_of_parts, overlap, save_parts=False, verbose=False):
	if not os.path.exists(filepath_output):
		os.makedirs(filepath_output)
		
	debruijn = fdgb.GraphData([reads], k_all, reduce_data=False, load_weights=False, remove_tips=True, directed_reads=True, verbose=verbose)
	debruijn.reduce_to_single_largest_component(verbose=verbose)
	debruijn.construct_assembly_ordering_labels(verbose=verbose)
	
	debruijn.get_asqg_output(filename = filepath_output+"/"+filename_output_base+"_start.asqg")
	debruijn.get_csv_output(filename = filepath_output+"/"+filename_output_base+"_start.csv")
	
	parts = debruijn.get_partition_of_sequences(number_of_parts, overlap=4, verbose=verbose)
	
	debruijn.reduce_to_single_path_max_weight()
	debruijn.contract_unique_overlaps()
	debruijn.get_asqg_output(filename = filepath_output+"/"+filename_output_base+"_start_reduced.asqg")
	debruijn.get_csv_output(filename = filepath_output+"/"+filename_output_base+"_start_reduced.csv")
	
	parts_sequences = []
	
	for part_id in range(number_of_parts):
		part_sequence = parts[part_id]
		if len(part_sequence) > 0:
			part_sequences = [debruijn.sequences[seq.id].sequence for seq in part_sequence]
			part_read_ids = debruijn.get_read_of_sequences(part_sequence)
			part_reads = [reads[i] for i in part_read_ids]
			#if verbose:
			print ("Now considering part "+str(part_id)+" which contains "+str(len(part_sequences))+" sequences and "+str(len(part_reads))+" reads.")
			debruijn_part = fdgb.GraphData([part_sequences], k_part, load_weights=False, remove_tips=True, directed_reads=True, verbose=verbose)
			if save_parts:
				debruijn_part.get_asqg_output(filename = filepath_output+"/"+filename_output_base+"_part"+str(part_id)+".asqg")
				debruijn_part.get_csv_output(filename = filepath_output+"/"+filename_output_base+"_part"+str(part_id)+".csv")
				
			debruijn_part.reduce_to_single_path_max_weight()
			debruijn_part.contract_unique_overlaps()
			parts_sequences += debruijn_part.get_relevant_sequences()
			
			if save_parts:
				debruijn_part.get_asqg_output(filename = filepath_output+"/"+filename_output_base+"_part"+str(part_id)+"_reduced.asqg")
				debruijn_part.get_csv_output(filename = filepath_output+"/"+filename_output_base+"_part"+str(part_id)+"_reduced.csv")
		
	debruijn_recons = fdgb.GraphData([parts_sequences], k_part, load_weights=False, remove_tips=True, directed_reads=True, verbose=verbose)
	
	debruijn_recons.reduce_to_single_largest_component(verbose=verbose)
	debruijn_recons.construct_assembly_ordering_labels(verbose=verbose)
	
	debruijn_recons.get_asqg_output(filename = filepath_output+"/"+filename_output_base+"_recons.asqg")
	debruijn_recons.get_csv_output(filename = filepath_output+"/"+filename_output_base+"_recons.csv")
	
	debruijn_recons.reduce_to_single_path_max_weight(verbose=2)
	debruijn_recons.contract_unique_overlaps()
	
	debruijn_recons.get_asqg_output(filename = filepath_output+"/"+filename_output_base+"_recons_reduced.asqg")
	debruijn_recons.get_csv_output(filename = filepath_output+"/"+filename_output_base+"_recons_reduced.csv")
	
if __name__ == "__main__":
	# test:
	g_n = 5000
	r_n = 500
	num_r = 5000
	er = 15.0
	
	name = "n"+str(g_n)+"_r"+str(num_r)+"_l"+str(r_n)+"_e"+str(int(er))
	
	dna = dgen.generate_dna(length = g_n)
	with open("Output/test_recursice_reconstruction/"+name+"_dna.txt", 'w') as outf:
		outf.write("".join(dna) + '\n')
	
	reads = dgen.samplereads(dna, number_of_reads=num_r, replace_error_percentage=er, read_length_mean=r_n, read_length_stddev=0, readlength_distribution='gaussian', inverted_reads=False, verbose=False)
	with open("Output/test_recursice_reconstruction/"+name+"_reads.txt", 'w') as outf:
		for r in reads:
			outf.write(r + '\n')
	
	#reads = dio.get_reads_from_file("Output/test_recursice_reconstruction/reads.txt")

	construct_consensus_from_multiple_parts(reads, filepath_output="Output/test_recursice_reconstruction", filename_output_base=name, k_all=30, k_part=13, number_of_parts=30, overlap=4, save_parts=False, verbose=False)