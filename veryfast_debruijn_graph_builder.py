#!usr/bin/python
# -*- coding: utf-8 -*-

import random
import re
import gc
import sys

import veryfast_debruijn_graph_builder_datastructures as ds
	
class GraphData:
	def __init__(self, reads=0, k=0, verbose=False, alphabet={"A":"T", "C":"G", "G":"C", "T":"A"}, init=True):
		print ("Creating empty graph ...")
		
		self.k_value = k
		self.reads = []
		self.kmers = []
		self.kmer_dict = {}
		self.sequences = {}
		self.overlaps = {}
		self.alphabet = alphabet
		
		# Flag that marks if inverse sequences have been removed
		self.is_unified = False
		# min and max label of all sequences
		self.max_label = 0
		self.min_label = 0
		
		if not reads == 0:
			self.init_read_database(reads, verbose=verbose)
			if init:
				self.init_graph_database(verbose=verbose)
				
	def print_memory_usage(self):
		size_reads = sys.getsizeof(self.reads) #sum([sys.getsizeof(r) for r in self.reads])
		size_kmers = sys.getsizeof(self.kmers) + sys.getsizeof(self.kmer_dict) #sum([sys.getsizeof(k) for k in self.kmers]) + sum([sys.getsizeof(k) + sys.getsizeof(self.kmer_dict[k]) for k in self.kmer_dict])
		size_sequences = sys.getsizeof(self.sequences) #sum([sys.getsizeof(seq) + sys.getsizeof(self.sequences[seq]) for seq in self.sequences])
		size_overlaps = sys.getsizeof(self.overlaps)
		
		print ("total memory usage: " + str((size_reads + size_kmers + size_sequences + size_overlaps)/1000000.0))
		print ("\treads: " + str(size_reads/1000000.0))
		print ("\tkmers: " + str(size_kmers/1000000.0))
		print ("\tsequences: " + str(size_sequences/1000000.0))
		print ("\toverlaps: " + str(size_overlaps/1000000.0))
	
	def init_read_database(self, reads, verbose=False):
		if verbose:
			print ("Construct read database")
		read_id = 0
		for r in reads:
			for read in r:
				self.reads.append(ds.Read(read_id, read))
				read_id += 1
	
	def init_graph_database(self, verbose=False):
		self.print_memory_usage()
		
		# construct k-mer database:
		self.get_kmerdata_from_reads(verbose)
		
		# construct sequences from kmers:
		print ("Construct Sequences from k-mers ...")
		for kmer in self.kmers:
			if verbose:
				print ("now consider kmer with id " + str(kmer.id))# + ": " + kmer.sequence)
			seq_id = kmer.id
			seq_inv_id = kmer.id_of_inverse_kmer
			weight = len(kmer.evidence_reads)
			self.sequences[seq_id] = ds.ContigSequence(seq_id, seq_inv_id, self.kmers[kmer.id].source[0], self.kmers[kmer.id].source[1], self.kmers[kmer.id].source[2], self.k_value, [kmer.id], weight)
		
		print ("Number of constructed sequences: "+str(len([seq for seq in self.sequences])))
		self.print_memory_usage()
		
		# construct overlaps between adjacent sequences with read-evidence:
		print ("Construct overlaps ...")			
		for read in self.reads:
			for kmer_index in range(len(read.kmers)-1):
				source_kmer_id = read.kmers[kmer_index]
				target_kmer_id = read.kmers[kmer_index+1]
				self.increment_overlap(source_kmer_id, target_kmer_id, read.id, verbose=False)
				
		print ("Number of constructed overlaps: "+str(len([ov for ov in self.overlaps])))
		self.print_memory_usage()

	def get_sequence_string_of_kmer(self, kmer_id):
		read_id = self.kmers[kmer_id].source[0]
		kmer_start = self.kmers[kmer_id].source[1]
		kmer_string_origin = self.reads[read_id].sequence[kmer_start:kmer_start+self.k_value]
		if not self.kmers[kmer_id].source[2]:
			return kmer_string_origin
		else:
			return get_inverse_sequence(kmer_string_origin, self.alphabet)

	def get_sequence_string_of_sequence(self, seq_id):
		seq_string = ""
		for seq_part in self.sequences[seq_id].read_sources:
			seq_part_string_origin = self.reads[seq_part[0]].sequence[seq_part[1]:(seq_part[1]+seq_part[3])]
			if not seq_part[2]:
				seq_string += seq_part_string_origin
			else:
				seq_string += get_inverse_sequence(seq_part_string_origin, self.alphabet)
		return seq_string

	def get_relevant_sequences(self):
		# returns list of all sequences with set relevance flag
		sequences = []
		for s_key in self.sequences:
			sequences.append(self.get_sequence_string_of_sequence(s_key))
		return sequences
		
	def get_kmerdata_from_reads(self, verbose = False):
		# checks all reads and constructs kmers and inverse kmers
		print ("Get kmer-data from reads ...")

		read_index = 0
		kmer_counter = 0
		for read_index in range(len(self.reads)):
			if read_index%100 == 0 and not verbose:
				ds.print_progress(read_index, len(self.reads))
				#print ("Progress: "+str("%.2f" % ((float(read_index)/(float(len(self.reads))/100)))) + "%")
			elif verbose:
				print ("Current read: "+str(read_index)+"/"+str(len(self.reads)) + " - " + self.reads[read_index].sequence)
			current_read_sequence = self.reads[read_index].sequence
			kmer_start = 0
			while kmer_start <= (len(current_read_sequence) - self.k_value):
				if verbose:
					print ("Current kmer-start: "+str(kmer_start))
				# construct new kmer from read sequence:
				new_kmer_sequence = current_read_sequence[kmer_start:kmer_start + self.k_value]
				kmer_already_existing = False
				this_kmer_id = -1
				
				# check if kmer already exists in database:
				if new_kmer_sequence in self.kmer_dict:
					kmer_already_existing = True
					this_kmer_id = self.kmer_dict[new_kmer_sequence]
					self.kmers[this_kmer_id].add_evidence(read_index)
					inv_kmer_id = self.kmers[this_kmer_id].id_of_inverse_kmer
					self.kmers[inv_kmer_id].add_evidence(read_index)
					if verbose:
						print ("Kmer already exists")
						print ("Add read ("+str(read_index)+") evidence to kmer "+str(this_kmer_id))
				if not kmer_already_existing:
					if verbose:
						print ("Kmer does not exist in database. Add new kmer ...")
					
					# add kmer:
					self.kmers.append(ds.Kmer(kmer_counter, kmer_counter+1, [read_index], read_index, kmer_start, False))
					self.kmer_dict[new_kmer_sequence] = kmer_counter
					this_kmer_id = kmer_counter
					kmer_counter += 1
					# add inverse kmer:
					self.kmers.append(ds.Kmer(kmer_counter, kmer_counter-1, [read_index], read_index, kmer_start, True))
					self.kmer_dict[get_inverse_sequence(new_kmer_sequence, self.alphabet)] = kmer_counter
					kmer_counter += 1
					
				# add kmer to read:
				self.reads[read_index].add_kmer(this_kmer_id)
				kmer_start += 1
			read_index += 1
			
	def increment_overlap(self, source_seq_id, target_seq_id, read_evidence, consider_inverse = True, verbose = False):
		# This method adds an overlap to the database. If overlap already exists, only read-evidence is added.
		if target_seq_id not in self.sequences[source_seq_id].overlaps_out:
			if verbose:
				print ("Add overlap of tuple " + str(source_seq_id) + " - " + str(target_seq_id))
				print ("seq_1: "+self.sequences[source_seq_id].sequence)
				print ("seq_2: "+self.sequences[target_seq_id].sequence)
			# add new overlap:
			ov_id = len(self.overlaps)
			self.overlaps[ov_id] = (ds.SequenceOverlap(ov_id, source_seq_id, target_seq_id))
			self.sequences[source_seq_id].overlaps_out[target_seq_id] = ov_id
			self.sequences[target_seq_id].overlaps_in[source_seq_id] = ov_id
			
			if consider_inverse:
				# add inverse overlap:
				rev_ov_id = ov_id + 1
				source_rev_seq_id = self.sequences[target_seq_id].id_of_inverse_seq
				target_rev_seq_id = self.sequences[source_seq_id].id_of_inverse_seq
				if verbose:
					print ("Add inverse overlap: "+str(source_rev_seq_id) + " - " + str(target_rev_seq_id))
					print ("seq_1: "+self.sequences[source_rev_seq_id].sequence)
					print ("seq_2: "+self.sequences[target_rev_seq_id].sequence)
				self.overlaps[rev_ov_id] = ds.SequenceOverlap(rev_ov_id, source_rev_seq_id, target_rev_seq_id)
				self.sequences[source_rev_seq_id].overlaps_out[target_rev_seq_id] = rev_ov_id
				self.sequences[target_rev_seq_id].overlaps_in[source_rev_seq_id] = rev_ov_id
			
		else:
			ov_id = self.sequences[source_seq_id].overlaps_out[target_seq_id]

	def contract_unique_overlaps(self, verbose = False):
		# This method contracts all overlaps between adjacent sequences that form an unique path
		# (i.e., there are no other outgoing or incoming overlaps between the sequences)
		print ("Contract overlaps ...")
		
		ov_index_list = [ov_id for ov_id in self.overlaps]
		num_deleted_overlaps = 0 
		for ov_index in ov_index_list:
			#if (num_deleted_overlaps > 0 and num_deleted_overlaps%100000 == 0):
			#	gc.collect()
			if (ov_index%1000 == 0):
				ds.print_progress(ov_index, len(ov_index_list))
				#print ("Progress: "+str("%.2f" % ((float(ov_index-num_deleted_overlaps)/(float(len(self.overlaps))/100)))) + "%")
				#print (str(ov_index-num_deleted_overlaps)+"/"+str(len(self.overlaps)))
			if ov_index in self.overlaps:
				source_id = self.overlaps[ov_index].contig_sequence_1
				target_id = self.overlaps[ov_index].contig_sequence_2
				if not self.is_unified:
					source_rev_id = self.sequences[target_id].id_of_inverse_seq
					target_rev_id = self.sequences[source_id].id_of_inverse_seq

				if verbose:
					print ("consider overlap: ")
					print (self.overlaps[ov_index].print_data())
					print ("Source: ")
					print (self.sequences[source_id].print_data())
					print ("")
				if (len(self.sequences[source_id].overlaps_out) == 1) and (len(self.sequences[target_id].overlaps_in) == 1):
					# if source node has exactly one outgoing edge
					# and the target node has exactly one incoming edge, 
					# then contract edge:
					ov = self.overlaps[ov_index]
					self.contract_overlap(ov_index, verbose)
					num_deleted_overlaps += 1
		            
					# contract reverse overlap if not sequence is its own inverse:
					source_sequence_string = self.get_sequence_string_of_sequence(source_id)
					if not source_sequence_string == get_inverse_sequence(source_sequence_string, self.alphabet):
						if not self.is_unified:
							rev_ov_id = self.sequences[source_rev_id].overlaps_out[target_rev_id]
							self.contract_overlap(rev_ov_id, verbose)
							num_deleted_overlaps += 1
						
						self.sequences[source_id].id_of_inverse_seq = source_rev_id
						self.sequences[source_rev_id].id_of_inverse_seq = source_id
	
	def delete_overlap(self, overlap_id, verbose=False):
		# removes an overlap from the database and from both incident sequences
		if verbose:
			print ("Delete overlap "+str(overlap_id))
			self.overlaps[overlap_id].print_data()
		if overlap_id not in self.overlaps:
			print ("Error! Overlap doesn't exist!")
		else:
			source_id = self.overlaps[overlap_id].contig_sequence_1
			target_id = self.overlaps[overlap_id].contig_sequence_2
			self.sequences[source_id].overlaps_out.pop(target_id)
			self.sequences[target_id].overlaps_in.pop(source_id)
			self.overlaps.pop(overlap_id)
		
	def delete_sequence(self, sequence_id, verbose=False):
		# deletes sequence and deletes all existing overlaps with this sequence
		if verbose:
			print ("Removing Sequence "+str(sequence_id))
		adj_seq_out = self.sequences[sequence_id].overlaps_out.keys()
		for adj_seq in adj_seq_out:
			self.delete_overlap(self.sequences[sequence_id].overlaps_out[adj_seq], verbose)
		adj_seq_in = self.sequences[sequence_id].overlaps_in.keys()
		for adj_seq in adj_seq_in:
			self.delete_overlap(self.sequences[sequence_id].overlaps_in[adj_seq], verbose)
		
		self.sequences.pop(sequence_id, None)
			
	def contract_overlap(self, overlap_id, verbose=False):
		# contracts a specific overlap
		source_id = self.overlaps[overlap_id].contig_sequence_1
		target_id = self.overlaps[overlap_id].contig_sequence_2
		if verbose:
			print ("Contract overlap: ")
			print (self.overlaps[overlap_id].print_data())
			print ("Source: ")
			print (self.sequences[source_id].print_data())
			print ("Target: ")
			print (self.sequences[target_id].print_data())
		# combine nucleotide sequences:
		#self.sequences[source_id].sequence += self.sequences[target_id].sequence[self.k_value-1:self.sequences[target_id].get_length()]
		target_sources = self.sequences[target_id].read_sources
		sseq_last_sr = self.sequences[source_id].read_sources[-1]
		tseq_first_sr = target_sources[0]
		if ((sseq_last_sr[0] == tseq_first_sr[0]) and (sseq_last_sr[2] == tseq_first_sr[2]) and (tseq_first_sr[1] <= sseq_last_sr[0] + sseq_last_sr[3])):
			# both sequence origin from same read on overlap
			# simply append last source of sourceseq
			self.sequences[source_id].read_sources[-1][3] = tseq_first_sr[1]+tseq_first_sr[3] - sseq_last_sr[1]
		#self.sequences[source_id].read_sources.append([target_sources[0][0], target_sources[0][1]+self.k_value-1, target_sources[0][2], target_sources[0][3]-self.k_value+1])
		for further_target in target_sources[1:]:
			self.sequences[source_id].read_sources.append(further_target)
		
		# update list of kmers:
		for kmer in self.sequences[target_id].kmers:
			if kmer not in self.sequences[source_id].kmers:
				self.sequences[source_id].kmers.append(kmer)
		# update maxweight:
		if self.sequences[target_id].max_weight > self.sequences[source_id].max_weight:
			self.sequences[source_id].max_weight = self.sequences[target_id].max_weight
		
		# update outgoing overlaps: 
		self.sequences[source_id].overlaps_out = self.sequences[target_id].overlaps_out
		
		# move outgoing overlaps from target_seq to source_seq:
		for ov_target_out in self.sequences[target_id].overlaps_out:
			# check if overlap_id exists:
			if self.sequences[target_id].overlaps_out[ov_target_out] in self.overlaps:
				# add overlap at source:
				self.sequences[source_id].overlaps_out[ov_target_out] = self.sequences[target_id].overlaps_out[ov_target_out]
				# update source of overlap:
				self.overlaps[self.sequences[target_id].overlaps_out[ov_target_out]].contig_sequence_1 = source_id
			else:
				self.sequences[source_id].overlaps_out.pop(ov_target_out)
		
		# update incoming overlaps for adjacent sequences:
		for adj_seq_id in self.sequences[source_id].overlaps_out:
			# add source to list of incoming overlaps:
			self.sequences[adj_seq_id].overlaps_in[source_id] = self.sequences[adj_seq_id].overlaps_in[target_id]
			# remove target from list of incoming overlaps:
			self.sequences[adj_seq_id].overlaps_in.pop(target_id)
			
		self.sequences[target_id].overlaps_in = {}
		self.sequences[target_id].overlaps_out = {}
		self.delete_sequence(target_id)
		#self.sequences[target_id].is_relevant = False
		# Don't use delete_overlap, because incident sequences have been handled manually:
		self.overlaps.pop(overlap_id)
	
	def remove_parallel_sequences(self, verbose=False):
		# For every pair of sequence and its inverse, this method removes one if both are not in the same component of the graph
		if not self.is_unified:
			print ("Remove parallel sequences ...")
			checked_sequences = {seq_id : False for seq_id in self.sequences}
			seq_id_list = [seq_id for seq_id in self.sequences]
			for seq_id in seq_id_list:
				if seq_id in self.sequences:
					seq = self.sequences[seq_id]
					#if not seq.is_relevant:
					#	checked_sequences[seq_id] = True
					if not checked_sequences[seq_id]:
						all_adjacent_sequences = set([ov for ov in seq.overlaps_out] + [ov for ov in seq.overlaps_in])
						if len(all_adjacent_sequences) == 0:
							# case no adjacent sequences:
							checked_sequences[self.sequences[seq_id].id_of_inverse_seq] = True
							self.delete_sequence(self.sequences[seq_id].id_of_inverse_seq)
							#if self.sequences[seq_id].is_relevant:
							#	self.sequences[self.sequences[seq_id].id_of_inverse_seq].is_relevant = False
						else:
							# start bfs
							bfs_queue = [seq_id]
							while (len (bfs_queue) > 0):
								if verbose:
									print bfs_queue
								current_seq_id = bfs_queue.pop(0)
								current_inv_seq_id = self.sequences[current_seq_id].id_of_inverse_seq
								# mark this sequence and its inverse as checked:
								checked_sequences[current_seq_id] = True
								checked_sequences[current_inv_seq_id] = True
								#if self.sequences[current_seq_id].is_relevant:
									# inverse sequense is not relevant:
								self.delete_sequence(current_inv_seq_id, verbose)
		
								# add all adjacent sequences to queue:
								for adj_seq_id in self.sequences[current_seq_id].overlaps_out:
									if verbose:
										print ("consider adj_seq "+str(adj_seq_id))
									if (not checked_sequences[adj_seq_id]) and (adj_seq_id not in bfs_queue):
										if verbose:
											print ("Add to bfs")
										bfs_queue.append(adj_seq_id)
								for adj_seq_id in self.sequences[current_seq_id].overlaps_in:
									if verbose:
										print ("consider adj_seq "+str(adj_seq_id))
									if (not checked_sequences[adj_seq_id]) and (adj_seq_id not in bfs_queue):
										if verbose:
											print ("Add to bfs")
										bfs_queue.append(adj_seq_id)
			self.is_unified = True
			
	def get_components(self, verbose=False):
		# returns a partition of the relevant sequences from self.sequences, as defined by the graph-components
		components = []
		
		visited_sequences = {seq.id: False for seq in self.sequences if seq.is_relevant}
		for seq in self.sequences:
			if not visited_sequences[seq.id]:
				current_comp = []
				seq_stack = [seq.id]
				visited_sequences[seq.id] = True
				while len(seq_stack) > 0:
					current_seq_id = seq_stack[0]
					seq_stack.pop(0)
					current_comp.append(current_seq_id)
					for adj_seq in self.sequences[current_seq_id].overlaps_out:
						if not visited_sequences[adj_seq]:
							seq_stack.append(adj_seq)
							visited_sequences[adj_seq] = True
					for adj_seq in self.sequences[current_seq_id].overlaps_in:
						if not visited_sequences[adj_seq]:
							seq_stack.append(adj_seq)
							visited_sequences[adj_seq] = True					
				components.append(current_comp)
		return components				
		
	def remove_tips(self, verbose=False):
		# removes tips (single-sequence-dead-ends) from the graph
		print ("Removing tips ...")
		for seq in self.sequences:
			if verbose:
				print ("Consider sequence:")
				seq.print_data()
			# check if sequence is a tip and if sequence is shorter than 2k:
			if (len(seq.overlaps_out) + len(seq.overlaps_in) == 1) and (seq.get_length() < 2*self.k_value):
				if verbose:
					print ("Sequence is a tip, remove this sequence.")
				self.delete_sequence(seq.id, verbose)

	def remove_insignificant_sequences(self, minimal_weight=2, verbose=False):
		# removes all sequences with weight less than minimal_weight
		for seq in self.sequences:
			if seq.max_weight < minimal_weight:
				self.delete_sequence(seq.id, verbose)
	
	def remove_single_sequence_components(self, verbose=False):
		# removes all components that consist only of a single sequence,
		# i.e. all sequences without adjacencies
		# only if there is at least one larger component
		if len(self.overlaps) > 0:
			for seq_id in range(len(self.sequences)):
				if len(self.sequences[seq_id].overlaps_out) == 0 and len(self.sequences[seq_id].overlaps_in) == 0:
					self.delete_sequence(seq_id, verbose)

	def get_asqg_output(self, filename="asqg_file"):
		print ("Writing asqg-file ...")
		headline = "HT\t\n"
		
		asqg_vertices = ""
		vertex_count = 0
		for seq in self.sequences:
			asqg_vertices += "VT\tk_"+str(seq)+"\t"+self.get_sequence_string_of_sequence(seq)+"\n"
			vertex_count += 1
		
		asqg_edges = ""
		for ov_id in self.overlaps:
			ov = self.overlaps[ov_id]
			seq_1_length = self.sequences[ov.contig_sequence_1].get_length()
			seq_2_length = self.sequences[ov.contig_sequence_2].get_length()
			asqg_edges += "ED\tk_"+str(ov.contig_sequence_1)+" k_"+str(ov.contig_sequence_2)+" "+str(seq_1_length-self.k_value+1)+" "+str(seq_1_length-1)+" "+str(seq_1_length)+" 0 "+str(self.k_value-2)+" "+str(seq_2_length)+" 0 0 "+str(self.k_value-1)+"\n"

		print filename
		outputfile = file(filename, 'w')
		outputfile.write(headline)
		outputfile.write(asqg_vertices)
		outputfile.write(asqg_edges)

	def get_csv_output(self, filename="csv_file.csv"):
		print ("writing csv-file ...")
		headline = "Node_Label, Sequence, maxweight, label\n"
		data = ""
		for seq in self.sequences:
			data += "k_"+str(seq)+","+self.get_sequence_string_of_sequence(seq)+","+str(self.sequences[seq].max_weight)+","+str(self.sequences[seq].label)+"\n"
		outputfile = file(filename, 'w')
		outputfile.write(headline)
		outputfile.write(data)
		
	def write_sequences_to_file(self, filename):
		# writes only the sequence-string of relevant sequences to a file
		data = ""
		for seq in self.sequences:
			data += self.get_sequence_string_of_sequence(seq.id) + "\n"
		outputfile = file(filename, 'w')
		outputfile.write(data)

	def construct_assembly_ordering_labels(self, verbose=False):
		# updated for new structure

		# constructs a fuzzy partial ordering of all relevant sequences:
		# algorithm assumes that graph
		# 	is not empty and
		#	has only one component and
		# 	has no cycles, i.e. implies a partial order
		print "Construct assembly ordering labels..."
		
		for start_seq_id in self.sequences:
			if not self.sequences[start_seq_id].label:
				queue = [[start_seq_id, 0]]
				while (len(queue) > 0):
					current_data = queue[0]
					queue.pop(0)
					current_node_id = current_data[0]
					for seq_id in self.sequences[current_node_id].overlaps_out:
						if self.sequences[seq_id].label == False:
							start_label = current_data[1] + self.sequences[current_node_id].get_length()
							if start_label > self.max_label:
								self.max_label = start_label
							self.sequences[seq_id].label = start_label
							queue.append([seq_id, start_label])
					for seq_id in self.sequences[current_node_id].overlaps_in:
						if self.sequences[seq_id].label == False:
							start_label = current_data[1] - self.sequences[seq_id].get_length()
							if start_label < self.min_label:
								self.min_label = start_label
							self.sequences[seq_id].label = start_label
							queue.append([seq_id, start_label])
		if verbose:
			for seq in self.sequences:
				print str(seq.id) + ": " + str(seq.label)
					
	def get_partition_of_sequences(self, number_of_parts, overlap=3, verbose=False):
		# returns a partition of all read-ids based on intervals of labels
		sorted_nodes = sorted([seq for seq in self.sequences if seq.label], key=lambda x: x.label)
		label_div = self.max_label-self.min_label
		part_start_difference = label_div/(number_of_parts+overlap-1)
		part_size = part_start_difference*overlap
		
		if verbose:
			print label_div
			print part_size
		
		parts_seq = []
		for i in range(number_of_parts):
			current_start = self.min_label+(i*part_start_difference)
			if i == number_of_parts-1:
				current_end = self.max_label
			else:
				current_end = current_start + part_size #self.min_label+(i+2)*(part_size)
			this_part_sequences = [seq for seq in sorted_nodes if seq.label >= current_start and seq.label <= current_end]
			parts_seq.append(this_part_sequences)
		return parts_seq
	
	def get_read_of_sequences(self, sequences, verbose=False):
		# returns all reads that contain a specific sequence
		kmers = []
		for seq in sequences:
			kmers += seq.kmers
		reads = []
		for kmer_id in kmers:
			reads += self.kmers[kmer_id].evidence_reads
		reads = list(set(reads))
		return reads
		
	def reduce_to_single_path_max_weight(self, verbose=False):
		# greedy algo that traverses through graph by choosing following nodes with max weight, deletes everythin else
		# method assumes that graph has only one component and no cycles
		# and sequences have weight-labels
		
		print ("Greedy reduction to single path with local max weight")
		
		components = self.get_components()
		component_id = 0
		while component_id < len(components):
			#for c in components:
			c = components[component_id]
			component_id += 1
			if verbose:
				print "next component"
			start_seq_id = -1
			min_label = False
			max_weight = 0
			if verbose:
				print ("Searching for start sequence...")
			for seq_id in c:
				if verbose:
					print ("Current min_label is: "+str(max_weight))
					print ("Current weight of min_sequence is: "+str(min_label))
					print ("Consider sequence "+str(seq_id)+" with label "+str(self.sequences[seq_id].label))
				new_min_found = False
				if (start_seq_id < 0) or (self.sequences[seq_id].label < min_label-self.k_value) or (self.sequences[seq_id].label < min_label+self.k_value and self.sequences[seq_id].max_weight > max_weight):
					if verbose:
						print ("Set sequence "+str(seq_id)+" as new minimal sequence")
					start_seq_id = seq_id
					min_label = self.sequences[seq_id].label
					max_weight = self.sequences[seq_id].max_weight
			
			current_seq_id = start_seq_id
			last_seq_id = -1
			if verbose:
				print "start_id: " +str(current_seq_id)
			while len(self.sequences[current_seq_id].overlaps_out) > 0:
				if verbose:
					print ("Current seq: "+str(current_seq_id))
				branching = False
				next_sequences = []
				for target_id in self.sequences[current_seq_id].overlaps_out:
					next_sequences.append(target_id)
				max_seq_id = -1
				max_seq_weight = -1
				for seq_id in next_sequences:
					if not seq_id == current_seq_id:
						if self.sequences[seq_id].max_weight > max_seq_weight:
							max_seq_weight = self.sequences[seq_id].max_weight
							max_seq_id = seq_id
				if len(next_sequences) > 1 and max_seq_weight == 1:
					# no unique continuation, keep all sequences:
					components += [[sid] for sid in next_sequences if not sid == max_seq_id]
					branching = True
				# delete all incoming sequences except path:
				if not last_seq_id == -1:
					incoming_sequences = [seq_id for seq_id in self.sequences[current_seq_id].overlaps_in]
					for seq_id in incoming_sequences:
						if not seq_id == last_seq_id:
							self.delete_sequence(seq_id, verbose=False)
				# set next sequence and delete other outgoing sequences:
				for seq_id in next_sequences:
					if seq_id == max_seq_id:
						last_seq_id = current_seq_id
						current_seq_id = seq_id
					elif not branching:
						self.delete_sequence(seq_id, verbose=False)
			# delete incoming sequences at final sequence, if there are any:
			if not last_seq_id == -1:
				incoming_sequences = [seq_id for seq_id in self.sequences[current_seq_id].overlaps_in]
				for seq_id in incoming_sequences:
					if not seq_id == last_seq_id:
						self.delete_sequence(seq_id, verbose=False)
					
	def remove_short_sequences(self, length_bound_by_multiple_of_k=5):
		# brutally removes all sequences with length less than 5 times k_value (if not specified otherwise)
		
		for seq in self.sequences:
			if seq.get_length() < 5*self.k_value:
				self.delete_sequence(seq.id)
		
		
def get_inverse_sequence(sequence, alphabet={"A":"T", "C":"G", "G":"C", "T":"A"}):
	n = len(sequence)
	inv_sequence = [""]*n
	for char_position in range(len(sequence)):
		current_char = sequence[char_position]
		if current_char in alphabet:
			inv_sequence[n-char_position-1] = alphabet[current_char]
		else:
			print (sequence)
			print (current_char)
			print ("Error! Incorrect Alphabet!")
			break
	return ''.join(inv_sequence)
