#!usr/bin/python
# -*- coding: utf-8 -*-

import random
import re
import sys

import numpy as np
import scipy.sparse
from scipy.sparse import linalg as la
import math
import gc
import Queue

import meta
import data_io as dio
	
class Read:
	def __init__(self, read_id, sequence, weight=1):
		self.id = read_id
		self.sequence = sequence
		self.length = len(sequence)
		self.kmers = []
		self.weight = weight
		#self.paired_end_partner = -1
		
	#def add_paired_end_partner(self, partner_id):
	#	self.paired_end_partner = partner_id
		
	def add_kmer(self, kmer_id):
		self.kmers.append(kmer_id)

class Kmer:
	def __init__(self, kmer_id, inverse_id , sequence, evidence_reads, baseweight=1):
		self.id = kmer_id
		self.sequence = sequence
		self.id_of_inverse_kmer = inverse_id
		self.evidence_reads = evidence_reads
		self.baseweight = baseweight
		
	def add_evidence(self, evidence_read_id):
		self.evidence_reads.append(evidence_read_id)
		
	def get_evidence_weight(self):
		return self.baseweight + len(self.evidence_reads) - 1
		
	def update_baseweight(self, weight):
		self.baseweight += weight

class ContigSequence:
	# Nodes in the Debruijn-Graph
	def __init__(self, seq_id, inv_id, sequence, kmers, weight=1, is_relevant=True):
		self.id = seq_id
		self.id_of_inverse_seq = inv_id
		self.sequence = sequence
		self.kmers = kmers
		# the maximal read eavidence this sequence has for any subsequence
		self.max_weight = weight
		# overlaps (i.e. edges) are stored in dictionaries
		# self.overlap[other_sequence_id] = overlap_id
		self.overlaps_out = {}
		self.overlaps_in = {}
		# deleted sequences will have the following flag set to False
		self.is_relevant = is_relevant
		# used for estimation of position within assembly
		self.label_p = False
		self.label_n = False
		
	def get_evidence_weight(self):
		return len(self.kmers)
		
	def get_total_weight(self):
		return self.max_weight	

	def get_length(self):
		return len(self.sequence)
		
	def check_if_overlap_exists(self, other_sequence_id):
		# checks whether this sequence has a common instance of SequenceOverlap with another sequence
		if other_sequence_id in self.overlaps_out:
			return self.overlaps_out[other_sequence_id]
		elif other_sequence_id in self.overlaps_in:
			return self.overlaps_in[other_sequence_id]
		else:
			return -1

	def print_data(self):
		print ("Sequence "+str(self.id))
		print (self.sequence)
		print (self.kmers)
		print ("Overlaps:")
		print (self.overlaps_out)
		print (self.overlaps_in)
		print ("id of inverse sequence: "+str(self.id_of_inverse_seq))
		if self.is_relevant:
			print ("is relevant")
		else:
			print ("is not relevant")

class SequenceOverlap:
	# Edges in the Debruijn-Graph
	def __init__(self, ov_id, length, seq_1, seq_2, evidence_reads, evidence_weights = []):
		self.id = ov_id
		self.length = length
		# the incident sequences of this overlap. sequence_1 is the source, sequence_2 is the target.
		self.contig_sequence_1 = seq_1
		self.contig_sequence_2 = seq_2
		self.evidence_reads = evidence_reads
		if not len(evidence_weights) == len(evidence_reads):
			self.evidence_weights = len(evidence_reads)
		else:
			self.evidence_weights = sum(evidence_weights)
		# Flag that checks whether a sequence is relevant for the graph or if it has been deleted
		self.is_relevant = True
	
	def add_evidence(self, read_id, read_weight=1):
		if read_id >= 0:
			self.evidence_reads.append(read_id)
			self.evidence_weights += read_weight
		
	def get_evidence_weight(self):
		return self.evidence_weights

	def print_data(self):
		print ("Overlap "+str(self.id))
		print ("From seq. "+str(self.contig_sequence_1)+" to seq. "+str(self.contig_sequence_2))
		if self.is_relevant:
			print ("is relevant")
		else:
			print ("is not relevant")
	
class GraphData:
	def __init__(	self,
					reads=[],
					k=0,
					alphabet={"A":"T", "C":"G", "G":"C", "T":"A"},
					directed_reads=False,
					load_weights=True,
					reduce_data=True,
					simplify_graph=True,
					remove_tips=False,
					construct_labels=True,
					verbose=False):
		print ("Creating empty graph ...")
		
		self.k_value = k
		self.reads = []
		self.kmers = []
		self.kmer_dict = {}
		self.sequences = []
		self.overlaps = {}
		self.alphabet = alphabet
		self.directed_reads = directed_reads
		
		self.hubs = []
		self.tips_to_keep = []
		self.tips_to_delete = []
		
		# Flag that marks if inverse sequences have been removed
		if not self.directed_reads:
			self.is_unified = False
		else:
			self.is_unified = True
		# Flag that marks if unique overlaps have been contracted
		self.is_contracted = False
		# min and max label of all sequences
		self.min_label_p = 0
		self.max_label_p = 0
		self.min_sequence_p = -1
		self.max_sequence_p = -1
		self.max_path_length_p = 0

		self.min_label_n = 0
		self.max_label_n = 0
		self.min_sequence_n = -1
		self.max_sequence_n = -1
		
		self.removed_reads = []
		
		if not reads == []:
			self.init_graph_database(reads, load_weights=load_weights, verbose=verbose)
			
			if reduce_data:
				# delete data to save ram:
				self.reads = []
				self.kmers = []
				# run garbage collector:
				gc.collect()
				
			if simplify_graph:
				if not self.directed_reads:
					self.remove_parallel_sequences(verbose=verbose)
				self.contract_unique_overlaps(verbose=verbose)
				
				if remove_tips:
					self.remove_tips(remove_only_unique_tips=True, verbose=verbose)
					self.remove_single_sequence_components(verbose=verbose)
			
			if construct_labels:
				if verbose:
					v=2
				else:
					v=1
				self.greedy_construct_assembly_ordering_labels(verbose=v)

	def print_memory_usage(self):
		size_reads = sys.getsizeof(self.reads)
		size_kmers = sys.getsizeof(self.kmers) + sys.getsizeof(self.kmer_dict)
		size_sequences = sys.getsizeof(self.sequences)
		size_overlaps = sys.getsizeof(self.overlaps)
		
		print ("total memory usage: " + str((size_reads + size_kmers + size_sequences + size_overlaps)/1000000.0) + "MB")
		print ("\treads: " + str(size_reads/1000000.0) + "MB")
		print ("\tkmers: " + str(size_kmers/1000000.0) + "MB")
		print ("\tsequences: " + str(size_sequences/1000000.0) + "MB")
		print ("\toverlaps: " + str(size_overlaps/1000000.0) + "MB")
		
	def init_graph_database(self, reads, load_weights, verbose=False):
		# note that reads is a list of lists of reads, since reads can have different sources.
		if verbose:
			print ("Construct read database")
			
		read_id = 0
		number_of_reads = sum((len(r) for r in reads))
		if verbose:
			print ("Number of reads: "+str(number_of_reads))
		for r in reads:
			for read in r:
				is_correct = True
				if (read_id%1000 == 0 or read_id == number_of_reads-1):
					meta.print_progress(read_id, number_of_reads-1)
				if len(read) > 0:
					# a read can be a comma-seperated tuple of (read, readweight):
					readdata = re.split(r',',read)
					readseq = readdata[0]
					if load_weights and len(readdata) > 1:
						readweight = int(readdata[1])
					else:
						readweight = 1
					# check if read has correct alphabet:
					for c in readseq:
						if c not in self.alphabet:
							is_correct = False
							print ("Error! Character "+str(c)+" not in alphabet "+str(self.alphabet.keys()))
				else:
					is_correct = False
				if is_correct:
					self.reads.append(Read(read_id, readseq, readweight))
				read_id += 1
		if verbose:
			print ("Number of reads in database: "+str(len(self.reads)))
		
		# construct k-mer database:
		self.get_kmerdata_from_reads(verbose)
		
		# construct sequences from kmers:
		print ("Construct Sequences from k-mers ...")
		for kmer in self.kmers:
			if (kmer.id % 10000 == 0 or kmer.id == len(self.kmers)-1):
				meta.print_progress(kmer.id, len(self.kmers)-1)
			seq_id = kmer.id
			seq_inv_id = kmer.id_of_inverse_kmer
			weight = kmer.get_evidence_weight()
			if verbose:
				print ("now consider kmer with id " + str(kmer.id) + ": " + kmer.sequence)
				print ("\tseq_id = "+str(seq_id))
				print ("\tseq_inv_id = "+str(seq_inv_id))
				print ("\tweight = "+str(weight))
			self.sequences.append(ContigSequence(seq_id, seq_inv_id, kmer.sequence, [kmer.id], weight))
			
		# construct overlaps between adjacent sequences with read-evidence:
		print ("Construct overlaps ...")
		read_number = 0
		for read in self.reads:
			if (read.id%100 == 0 or read_number == len(self.reads)-1):
				meta.print_progress(read_number, len(self.reads)-1)
			for kmer_index in range(len(read.kmers)-1):
				source_kmer_id = read.kmers[kmer_index]
				target_kmer_id = read.kmers[kmer_index+1]
				self.increment_overlap(source_kmer_id, target_kmer_id, read.id, read.weight, verbose=False)
			read_number += 1
		if verbose:
			self.print_memory_usage()

	def print_all_reads(self):
		print ("Reads:")
		for read in self.reads:
			print (str(read.id) + ": "+read.sequence+ ", "+str(read.weight))
			kmer_string = "kmers: "
			for kmer_id in read.kmers:
				kmer_string += str(kmer_id) + " "
			print kmer_string

	def print_all_kmers(self):
		print ("Kmers:")
		for kmer in self.kmers:
			kmer_info_string = "("+str(kmer.id) + ") " + kmer.sequence + ": "
			for read in kmer.evidence_reads:
				kmer_info_string += str(read) + " "
			print kmer_info_string
	
	def print_all_sequences(self):
		print ("Sequences:")
		for s in self.sequences:
			if s.is_relevant:
				print ("("+str(s.id)+") "+s.sequence) + ": "
				s.print_data()
				print("")
	
	def get_relevant_sequences(self):
		# returns list of all sequences with set relevance flag
		sequences = []
		for s in self.sequences:
			if s.is_relevant:
				sequences.append(s.sequence)
		return sequences
		
	def get_kmerdata_from_reads(self, verbose = False):
		# checks all reads and constructs kmers and inverse kmers (if not self.directed_reads)
		print ("Get kmer-data from reads ...")
		read_index = 0
		kmer_counter = 0
		
		progress_step = max(1, int(math.ceil(len(self.reads)/1000)))
		for read_index in range(len(self.reads)):
			if (read_index%progress_step == 0 or read_index == range(len(self.reads)-1)) and not verbose:
				meta.print_progress(read_index, len(self.reads)-1)
			elif verbose:
				print ("Current read: "+str(read_index)+"/"+str(len(self.reads)) + " - " + self.reads[read_index].sequence)
			current_read_sequence = self.reads[read_index].sequence
			current_read_weight = self.reads[read_index].weight
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
					self.kmers[this_kmer_id].update_baseweight(current_read_weight-1)
					if not self.directed_reads:
						inv_kmer_id = self.kmers[this_kmer_id].id_of_inverse_kmer
						self.kmers[inv_kmer_id].add_evidence(read_index)
						self.kmers[inv_kmer_id].update_baseweight(current_read_weight-1)
					if verbose:
						print ("Kmer already exists")
						print ("Add read ("+str(read_index)+") evidence to kmer "+str(this_kmer_id))
				if not kmer_already_existing:
					if verbose:
						print ("Kmer does not exist in database. Add new kmer ...")
					
					if not self.directed_reads:
						inv_kmer_seq = meta.get_inverse_sequence(new_kmer_sequence, self.alphabet)
						if not inv_kmer_seq == new_kmer_sequence:
							id_of_inv_kmer = kmer_counter+1
						else:
							id_of_inv_kmer = kmer_counter
					else:
						id_of_inv_kmer = -1
						
					# add kmer:
					self.kmers.append(Kmer(kmer_counter, id_of_inv_kmer, new_kmer_sequence, [read_index], baseweight=current_read_weight))
					self.kmer_dict[new_kmer_sequence] = kmer_counter
					this_kmer_id = kmer_counter
					kmer_counter += 1
					if not self.directed_reads:
						# add inverse kmer:
						self.kmers.append(Kmer(kmer_counter, kmer_counter-1, inv_kmer_seq, [read_index], baseweight=current_read_weight))
						self.kmer_dict[inv_kmer_seq] = kmer_counter
						kmer_counter += 1
					
				if verbose:
					print ("Add kmer "+self.kmers[this_kmer_id].sequence+"("+str(this_kmer_id)+") to read "+str(read_index))
				# add kmer to read:
				self.reads[read_index].add_kmer(this_kmer_id)
				kmer_start += 1
			read_index += 1
			
	def increment_overlap(self, source_seq_id, target_seq_id, read_evidence, read_weight=1, consider_inverse = True, verbose = False):
		# This method adds an overlap to the database. If overlap already exists, only read-evidence is added.
		if target_seq_id not in self.sequences[source_seq_id].overlaps_out:
			if verbose:
				print ("Add overlap of tuple " + str(source_seq_id) + " - " + str(target_seq_id))
				print ("seq_1: "+self.sequences[source_seq_id].sequence)
				print ("seq_2: "+self.sequences[target_seq_id].sequence)
			# add new overlap:
			ov_id = len(self.overlaps)
			self.overlaps[ov_id] = (SequenceOverlap(ov_id, self.k_value-1, source_seq_id, target_seq_id, []))
			self.sequences[source_seq_id].overlaps_out[target_seq_id] = ov_id
			self.sequences[target_seq_id].overlaps_in[source_seq_id] = ov_id
			
			if not self.directed_reads and consider_inverse:
				# add inverse overlap:
				rev_ov_id = ov_id + 1
				source_rev_seq_id = self.sequences[target_seq_id].id_of_inverse_seq
				target_rev_seq_id = self.sequences[source_seq_id].id_of_inverse_seq
				if verbose:
					print ("Add inverse overlap: "+str(source_rev_seq_id) + " - " + str(target_rev_seq_id))
					print ("seq_1: "+self.sequences[source_rev_seq_id].sequence)
					print ("seq_2: "+self.sequences[target_rev_seq_id].sequence)
				self.overlaps[rev_ov_id] = SequenceOverlap(rev_ov_id, self.k_value-1, source_rev_seq_id, target_rev_seq_id, [])
				self.sequences[source_rev_seq_id].overlaps_out[target_rev_seq_id] = rev_ov_id
				self.sequences[target_rev_seq_id].overlaps_in[source_rev_seq_id] = rev_ov_id
			
		else:
			ov_id = self.sequences[source_seq_id].overlaps_out[target_seq_id]
			
		self.overlaps[ov_id].add_evidence(read_evidence, read_weight)
	
	'''
	def add_overlaps_for_sequences_with_small_insert_distance(self, max_insert_distance=1, verbose = False):
		print ("Adding overlaps to sequences with small insert distance ...")
		for sequence_id_1 in range(len(self.sequences)):
			meta.print_progress(sequence_id_1, len(self.sequences))
			if self.sequences[sequence_id_1].is_relevant:
				for sequence_id_2 in range(sequence_id_1+1, len(self.sequences)):
					if self.sequences[sequence_id_2].is_relevant:
						s1 = self.sequences[sequence_id_1].sequence
						s2 = self.sequences[sequence_id_2].sequence
						l1 = len(s1)
						l2 = len(s2)
						# note that a sequence from the same read has always position-distance = insert-distance, but we do not want to add an overlap to these sequences:
						# check if both sequence do share common evidence-reads:
						reads_seq_1 = []
						for kmer_id in self.sequences[sequence_id_1].kmers:
							reads_seq_1 += self.kmers[kmer_id].evidence_reads
						reads_seq_2 = []
						for kmer_id in self.sequences[sequence_id_2].kmers:
							reads_seq_2 += self.kmers[kmer_id].evidence_reads
						common_reads = [r for r in reads_seq_1 if r in reads_seq_2]
						if len(common_reads) == 0:
							#if not sequence_id_2 in self.sequences[sequence_id_1].overlaps_in:
							if meta.compute_insert_distance(s1[l1-self.k_value:], s2[:self.k_value], maxdist = max_insert_distance) <= max_insert_distance:
								self.increment_overlap(sequence_id_1, sequence_id_2, -1)
							#if not sequence_id_1 in self.sequences[sequence_id_2].overlaps_in:
							if meta.compute_insert_distance(s2[l2-self.k_value:], s1[:self.k_value], maxdist = max_insert_distance) <= max_insert_distance:
								self.increment_overlap(sequence_id_2, sequence_id_1, -1)
	'''

	def contract_unique_overlaps(self, verbose = 1):
		# This method contracts all overlaps between adjacent sequences that form an unique path
		# (i.e., there are no other outgoing or incoming overlaps between the sequences)
		if not self.is_unified:
			self.remove_parallel_sequences(verbose)
		
		if verbose > 0:
			print ("Contract overlaps ...")
		
		ov_index_list = [ov_id for ov_id in self.overlaps]
		num_deleted_overlaps = 0
		ov_counter = 0
		for ov_index in ov_index_list:
			if verbose > 0:
				if (ov_counter%1000 == 0 or ov_counter == len(ov_index_list)-1):
					meta.print_progress(ov_counter, len(ov_index_list)-1)
			ov_counter += 1
			if ov_index in self.overlaps:
				source_id = self.overlaps[ov_index].contig_sequence_1
				target_id = self.overlaps[ov_index].contig_sequence_2

				if verbose == 2:
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
					self.contract_overlap(ov_index, verbose == 2)
					num_deleted_overlaps += 1
		self.is_contracted = True
	
	def delete_overlap(self, overlap_id, verbose=False):
		# removes an overlap from the database and from both incident sequences
		if verbose:
			print ("Delete overlap "+str(overlap_id))
			self.overlaps[overlap_id].print_data()
		if overlap_id not in self.overlaps:
			print ("Error! Overlap does not exist!")
		else:
			source_id = self.overlaps[overlap_id].contig_sequence_1
			target_id = self.overlaps[overlap_id].contig_sequence_2
			self.sequences[source_id].overlaps_out.pop(target_id)
			self.sequences[target_id].overlaps_in.pop(source_id)
			self.overlaps.pop(overlap_id)
		
	def delete_sequence(self, sequence_id, verbose=False):
		# flags sequence as not relevant and deletes all existing overlaps with this sequence
		if verbose:
			print ("Removing Sequence "+str(sequence_id))
		adj_seq_out = self.sequences[sequence_id].overlaps_out.keys()
		for adj_seq in adj_seq_out:
			if verbose:
				print ("Delete overlap to sequence " +str(adj_seq))
			self.delete_overlap(self.sequences[sequence_id].overlaps_out[adj_seq], verbose)
		adj_seq_in = self.sequences[sequence_id].overlaps_in.keys()
		for adj_seq in adj_seq_in:
			if verbose:
				print ("Delete overlap to sequence " +str(adj_seq))
			self.delete_overlap(self.sequences[sequence_id].overlaps_in[adj_seq], verbose)
		self.sequences[sequence_id].is_relevant = False
			
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
		
		if not self.sequences[source_id].sequence == meta.get_inverse_sequence(self.sequences[target_id].sequence, self.alphabet):
			# combine nucleotide sequences:
			self.sequences[source_id].sequence += self.sequences[target_id].sequence[self.k_value-1:self.sequences[target_id].get_length()]
			if verbose:
				print ("combined sequence: " + self.sequences[source_id].sequence)
			# update outgoing overlaps:
			self.sequences[source_id].overlaps_out = self.sequences[target_id].overlaps_out
			# update list of kmers:
			for kmer in self.sequences[target_id].kmers:
				if kmer not in self.sequences[source_id].kmers:
					self.sequences[source_id].kmers.append(kmer)
			# update maxweight:
			if self.sequences[target_id].get_total_weight() > self.sequences[source_id].get_total_weight():
				self.sequences[source_id].max_weight = self.sequences[target_id].get_total_weight()
			
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
				
			# Don't use delete_overlap or delete_sequence, because incident sequences have been handled manually:
			self.sequences[target_id].overlaps_in = {}
			self.sequences[target_id].overlaps_out = {}
			self.sequences[target_id].is_relevant = False
			self.overlaps.pop(overlap_id)
			
		else:
			if verbose:
				print ("Not contracting sequence with its own inverse:")
				print ("Source: ")
				print (self.sequences[source_id].print_data())
				print ("Target: ")
				print (self.sequences[target_id].print_data())
	
	def remove_parallel_sequences(self, verbose=False):
		# For every pair of sequence and its inverse, this method removes one s.t. of each pair of isomorphic components, one component remains intact.
		# note: only works correct if there is no sequence that is its own inverse!
		if not self.is_unified:
			print ("Remove parallel sequences ...")
			checked_sequences = [False for seq in self.sequences]
			for seq_id in range(len(self.sequences)):
				seq = self.sequences[seq_id]
				if not seq.is_relevant:
					checked_sequences[seq_id] = True
				elif not checked_sequences[seq_id]:
					all_adjacent_sequences = set([ov for ov in seq.overlaps_out] + [ov for ov in seq.overlaps_in])
					if len(all_adjacent_sequences) == 0:
						# case no adjacent sequences:
						checked_sequences[self.sequences[seq_id].id_of_inverse_seq] = True
						if self.sequences[seq_id].is_relevant:
							self.sequences[self.sequences[seq_id].id_of_inverse_seq].is_relevant = False
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
							if self.sequences[current_seq_id].is_relevant:
								# inverse sequense is not relevant:
								self.delete_sequence(current_inv_seq_id, verbose)
	
								# add all adjacent sequences to queue:
								for adj_seq_id in self.sequences[current_seq_id].overlaps_out:
									if verbose:
										print ("consider adj_seq "+str(adj_seq_id))
									if (not checked_sequences[adj_seq_id]) and (adj_seq_id not in bfs_queue) and (self.sequences[adj_seq_id].is_relevant):
										if verbose:
											print ("Add to bfs")
										bfs_queue.append(adj_seq_id)
								for adj_seq_id in self.sequences[current_seq_id].overlaps_in:
									if verbose:
										print ("consider adj_seq "+str(adj_seq_id))
									if (not checked_sequences[adj_seq_id]) and (adj_seq_id not in bfs_queue) and (self.sequences[adj_seq_id].is_relevant):
										if verbose:
											print ("Add to bfs")
										bfs_queue.append(adj_seq_id)
			self.is_unified = True
			
	def get_components(self, verbose=False):
		# returns a partition of the relevant sequences from self.sequences, as defined by the graph-components
		components = []
		
		visited_sequences = {seq.id: False for seq in self.sequences if seq.is_relevant}
		for seq in self.sequences:
			if seq.is_relevant:
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
		
	def remove_tips(self, only_simply_connected_tips=True, maximum_tip_weight=-1, remove_only_unique_tips=True, verbose=False):
		# removes tips (single-sequence-dead-ends) from the graph
		# if only_simple_connected_tips==True, only tips will be removed that are connected to exactly one other node
		# if maximum_tip_weight > 0, only tips with weight < maximum_tip_weight will be removed.
		# if remove_only_unique_tips==True, a tip will only be removed if the node it is connected to as other connections in the same orientation, which have to be non-tips.
		
		print ("Removing tips ...")
		num_of_removed_tips = -1
		iteration = 0
		while (num_of_removed_tips != 0):
			iteration += 1
			print ("Current iteration of tip removal: "+str(iteration))
			if remove_only_unique_tips:
				# initialize list of tips to keep:
				self.get_tip_classification()
			num_of_removed_tips = 0
			for seq in self.sequences:
				if seq.is_relevant:
					#if verbose:
					#	print ("Consider sequence:")
					#	seq.print_data()
					# check if sequence is a tip and if sequence is shorter than 2k:
					if (len(seq.sequence) < 2*self.k_value):
						if seq.get_total_weight() < maximum_tip_weight or maximum_tip_weight <= 0:
							is_tip = False
							if (only_simply_connected_tips and (len(seq.overlaps_out) + len(seq.overlaps_in) == 1)):
								# (sequence is a simply connected tip)
								is_tip = True
							elif (len(seq.overlaps_out) == 0 or len(seq.overlaps_in) == 0) and (len(seq.overlaps_in) + len(seq.overlaps_out) > 0):
								# (sequence is a tip)
								is_tip =True							
							if is_tip:
								if verbose:
									print ("Consider tip sequence "+str(seq.id))
								remove_tip = False
								if (not remove_only_unique_tips) or (seq.id not in self.tips_to_keep):
									remove_tip = True											
								if remove_tip:
									if verbose:
										print ("Sequence is a tip, remove this sequence.")
									# remove sequence:
									self.delete_sequence(seq.id, verbose)
									num_of_removed_tips += 1
									# mark reads of tip:
									tip_reads = self.get_read_of_sequences([seq])
									self.removed_reads += tip_reads
							
			self.contract_unique_overlaps(verbose = 0)

	def unrestricted_tip_removal(self, only_simply_connected_tips=True, verbose=False):
		# removes tips of any lengths, but only if the node to which it is connected has another adjacent node in the same direction, which also has to be longer than the tip.
		
		num_of_removed_tips = -1
		while (num_of_removed_tips != 0):
			num_of_removed_tips = 0
			for seq in self.sequences:
				if seq.is_relevant:
					if verbose:
						print ("Consider sequence:")
						seq.print_data()
					# check if sequence is a tip:
					is_tip = False
					if (only_simply_connected_tips and (len(seq.overlaps_out) + len(seq.overlaps_in) == 1)):
						# (sequence is a simply connected tip)
						is_tip = True
					elif (len(seq.overlaps_out) == 0 or len(seq.overlaps_in) == 0):
						# (sequence is a tip)
						is_tip =True
					if is_tip:
						length_of_tip = len(seq.sequence)
						# check if tip is longest adjacent node:						
						number_of_adj_nodes_with_longer_adj_nodes = 0
						if len(seq.overlaps_in) > 0:
							# case 1: tip is outgoing
							for adj_seq in seq.overlaps_in:
								for adj_adj_seq in self.sequences[adj_seq].overlaps_out:
									if len(self.sequences[adj_adj_seq].sequence) > length_of_tip:
										number_of_adj_nodes_with_longer_adj_nodes += 1
						elif len(seq.overlaps_out) > 0:
							# case 2: tip is incoming
							for adj_seq in seq.overlaps_out:
								for adj_adj_seq in self.sequences[adj_seq].overlaps_in:
									if len(self.sequences[adj_adj_seq].sequence) > length_of_tip:
										number_of_adj_nodes_with_longer_adj_nodes += 1
										
						if number_of_adj_nodes_with_longer_adj_nodes == len(seq.overlaps_out) + len(seq.overlaps_in):
							# every adjacent node has another longer adjacent sequence
							self.delete_sequence(seq.id, verbose)
							num_of_removed_tips += 1
			self.contract_unique_overlaps()
		
	def remove_insignificant_sequences(self, minimal_weight=2, verbose=False):
		# removes all sequences with weight less than minimal_weight
		print ("Removing sequences with evidence_weight < "+str(minimal_weight)+" ...")
		for seq in self.sequences:
			if seq.get_total_weight() < minimal_weight:
				self.delete_sequence(seq.id, verbose)
				self.removed_reads += self.get_read_of_sequences([seq])
	
	def remove_insignificant_overlaps(self, minimal_evidence=2, do_contraction=True, keep_relevant_tips=False, verbose=False):
		# remove all overlaps with evidence_weight less than minimal_evidence
		if keep_relevant_tips:
			self.get_tip_classification()
		
		print ("Removing overlaps with evidence_weight < "+str(minimal_evidence)+" ...")
		ov_ids_to_delete = []
		for ov_id in self.overlaps:
			if self.overlaps[ov_id].get_evidence_weight() < minimal_evidence:
				if verbose:
					print ("Removing overap "+str(ov_id))
				if not keep_relevant_tips:
					ov_ids_to_delete.append(ov_id)
				else:
					if (self.overlaps[ov_id].contig_sequence_1 not in self.tips_to_keep) and (self.overlaps[ov_id].contig_sequence_2 not in self.tips_to_keep):
						ov_ids_to_delete.append(ov_id)
		for ov_id in ov_ids_to_delete:
			self.delete_overlap(ov_id)
			
		if do_contraction:
			self.contract_unique_overlaps()
	
	def remove_single_sequence_components(self, verbose=False):
		# removes all components that consist only of a single sequence,
		# i.e. all sequences without adjacencies
		# only if there is at least one larger component
		if len(self.overlaps) > 0:
			for seq_id in range(len(self.sequences)):
				if self.sequences[seq_id].is_relevant and len(self.sequences[seq_id].overlaps_out) == 0 and len(self.sequences[seq_id].overlaps_in) == 0:
					self.sequences[seq_id].is_relevant = False
					self.removed_reads += self.get_read_of_sequences([self.sequences[seq_id]])
			return True
		else:
			return False
					
	def remove_single_sequence_loops(self, do_contraction=True, verbose=False):
		ov_ids_to_delete = []
		for ov_id in self.overlaps:
			ov = self.overlaps[ov_id]
			if ov.contig_sequence_1 == ov.contig_sequence_2:
				ov_ids_to_delete.append(ov_id)
		for ov_id in ov_ids_to_delete:
			self.delete_overlap(ov_id)
		if do_contraction:
			self.contract_unique_overlaps()

	def get_asqg_output(self, filename="asqg_file"):
		print ("Writing asqg-file ...")
		headline = "HT\t\n"
		
		asqg_vertices = ""
		vertex_count = 0
		for seq in self.sequences:
			if seq.is_relevant:
				asqg_vertices += "VT\tk_"+str(seq.id)+"\t"+seq.sequence+"\n"
				vertex_count += 1
		
		asqg_edges = ""
		for ov_id in self.overlaps:
			ov = self.overlaps[ov_id]
			if ov.is_relevant:
				seq_1_length = self.sequences[ov.contig_sequence_1].get_length()
				seq_2_length = self.sequences[ov.contig_sequence_2].get_length()
				asqg_edges += "ED\tk_"+str(ov.contig_sequence_1)+" k_"+str(ov.contig_sequence_2)+" "+str(seq_1_length-self.k_value+1)+" "+str(seq_1_length-1)+" "+str(seq_1_length)+" 0 "+str(self.k_value-2)+" "+str(seq_2_length)+" 0 0 "+str(len(ov.evidence_reads))+"\n"

		print filename
		outputfile = file(filename, 'w')
		outputfile.write(headline)
		outputfile.write(asqg_vertices)
		outputfile.write(asqg_edges)

	def get_csv_output(self, filename="csv_file.csv"):
		print ("writing csv-file ...")
		headline = "Node_Label, Sequence, maxweight, label, read_ids\n"
		data = ""
		for seq in self.sequences:
			if seq.is_relevant:
				source_reads = []
				if not self.kmers == []:
					for k_id in seq.kmers:
						source_reads += [r for r in self.kmers[k_id].evidence_reads]
					source_reads = set(source_reads)
				data += "k_"+str(seq.id)+","+seq.sequence+","+str(seq.get_total_weight())+","+str(seq.label_p)+","+str(source_reads).replace(",",";")+"\n"
				#for r in source_reads:
				#	data += str(r)+" "
				#data += "\n"
		outputfile = file(filename, 'w')
		outputfile.write(headline)
		outputfile.write(data)
		
	def write_sequences_to_file(self, filename, addweights=False, asfasta=False):
		# writes only the sequence-string of relevant sequences to a file
		# if asfasta=True, writes sequences to file in fasta format, this ignores(!) the option addweights
		allseqs = []
		for seq in self.sequences:
			if seq.is_relevant:
				allseqs.append(seq.sequence)
		if asfasta:
			dio.write_sequences_to_fasta(allseqs, filename)
		else:
			data = ""
			for seq in allseqs:
				data += seq
				if addweights:
					data += ","+str(seq.get_total_weight())
				data += "\n"
				
			outputfile = file(filename, 'w')
			outputfile.write(data)
		
	def load_from_asqg(self, filename="asqg_file", verbose=False):
		# reconstructs a graph from a asqg-file
		print ("Loading from asqg-file ...")
		self.is_unified = True
		self.directed_reads = True
		with open(filename) as asqg_source:
			for line in asqg_source:
				linedata = re.split(r'\t', line)
					
				if linedata[0] == "VT":
					# case sequence data:
					sequence_id = int(re.split(r'_',linedata[1])[1])
					while len(self.sequences) < sequence_id:
						# add dummy sequences:
						self.sequences.append(ContigSequence(len(self.sequences), -1, "", [-1], is_relevant=False))
					# add actual sequence (id of inverse and evidence_kmers are unknown)
					self.sequences.append(ContigSequence(sequence_id, -1, linedata[2].strip(), [-1]))
					
				elif linedata[0] == "ED":
					# case overlap data:
					overlap_data = re.split(r'\s',linedata[1])
					source_seq_id = int(re.split(r'_',overlap_data[0])[1])
					target_seq_id = int(re.split(r'_',overlap_data[1])[1])
					if source_seq_id > 0 and target_seq_id > 0:
						# increase overlap (read_evidence is unknown)
						self.increment_overlap(source_seq_id, target_seq_id, [-1], verbose=verbose,consider_inverse = False)
					# if k is still unknown at this point, recover k from overlap_data:
					if self.k_value < 1:
						self.k_value = int(overlap_data[6])					
						
	def construct_assembly_ordering_labels(self, start_sequence = 0, do_second_iteration=True, delete_cycles=True, verbose=1):
		# constructs a fuzzy partial ordering of all relevant sequences based on directed graph structure and sequence lengths:
		# algorithm assumes that graph
		# 	is not empty and
		#	has only one component
		# edges that induce a cycle are ignored when they are encountered.
		# if graph contains many cycles, this algorithm has exponential running time!
		
		if verbose > 0:
			print "Construct longest assembly ordering labels..."
	
		# initialize:
		if not (start_sequence == 0 or self.sequences[start_sequence].is_relevant):
			if verbose > 0:
				print ("Error! Start sequence does not exist! Start with first sequence instead.")
			start_sequence = 0
		while start_sequence < len(self.sequences) and not self.sequences[start_sequence].is_relevant:
			start_sequence += 1
		
		if not start_sequence < len(self.sequences):
			print("Error! No sequence found!")
			return
		
		# reset labels
		for seq in self.sequences:
			if seq.is_relevant:
				seq.label_p = 0
		# init for start sequence
		self.min_label_p = 0
		self.max_label_p = 0
		self.min_sequence_p = start_sequence
		self.max_sequence_p = start_sequence
		self.max_path_length_p = self.sequences[start_sequence].get_length()
		self.sequences[start_sequence].label_p = 0

		if verbose >= 2:
			print ("Start with sequence " + str(start_sequence))
		
		predecessor = [-1 for seq in self.sequences]
		successor = [-1 for seq in self.sequences]
		was_visited = [False for seq in self.sequences]
		priority_queue = Queue.PriorityQueue()
		priority_queue.put((0,self.sequences[start_sequence].id))
		while not priority_queue.empty():
			current_data = priority_queue.get()
			current_node_id = current_data[1]
			current_position = -current_data[0]
			was_visited[current_node_id] = True
			
			successor_label = current_position + self.sequences[current_node_id].get_length() - (self.k_value-1)
			
			if verbose == 2:
				print ("Next sequence: " + str(current_node_id) + ", current label: " + str(current_position))
				
			next_nodes = [seq_id for seq_id in self.sequences[current_node_id].overlaps_out]
			for seq_id in next_nodes: # additional list necessary if cycles are deleted!
				if seq_id in self.sequences[current_node_id].overlaps_out:
					if not was_visited[seq_id] or self.sequences[seq_id].label_p < current_position:
						# check for cycle:
						is_cycle = False
						current_pred = current_node_id
						while self.sequences[current_pred].label_p > self.sequences[seq_id].label_p and not predecessor[current_pred] < 0 and not is_cycle:
							current_pred = predecessor[current_pred]
							if current_pred == seq_id:
								is_cycle = True
								if verbose >= 2:
									print ("Cycle Found! Ignore path to sequence "+str(seq_id))
									
								if delete_cycles:
									# fix the cycle: find and delete the overlap with minimum read evidence of the cycle:
									# find overlap with minimun evidence weight within cycle:
									min_ov_weight = -1
									min_ov_id = -1
									current_source = current_node_id
									current_target = seq_id
									while not current_source == seq_id:
										current_ov_id = self.sequences[current_source].overlaps_out[current_target]
										current_ov_weight = self.overlaps[current_ov_id].get_evidence_weight()
										if min_ov_weight < 0 or current_ov_weight < min_ov_weight:
											min_ov_weight = current_ov_weight
											min_ov_id = current_ov_id
										current_target = current_source
										current_source = predecessor[current_source]
									# delete this overlap:
									self.delete_overlap(min_ov_id, verbose = verbose>2)
									
									# for all succeeding nodes of deleted overlap: reset the position label:
									current_source = current_node_id
									current_target = seq_id
									while current_source in self.sequences[current_target].overlaps_in:
										current_target = current_source
										current_source = predecessor[current_source]
										predecessor[current_target] = -1
										was_visited[current_target] = False
										self.sequences[current_target].label_p = 0
								
						# if no cycle found, update labels:
						if not is_cycle:
							if successor_label > self.max_label_p:
								self.max_label_p = successor_label
								self.max_sequence_p = seq_id
							self.sequences[seq_id].label_p = successor_label
							if verbose >= 2:
								print ("updated label of node "+str(seq_id)+": "+str(successor_label))
							predecessor[seq_id] = current_node_id
							end_label = successor_label + self.sequences[seq_id].get_length()
							if end_label > self.max_path_length_p:
								self.max_path_length_p = end_label
							priority_queue.put((-successor_label, seq_id))
			
			# going backwards: only set labels for previously unvisited nodes
			for seq_id in self.sequences[current_node_id].overlaps_in:
				predecessor_label = current_position - self.sequences[seq_id].get_length() + (self.k_value-1)
				if not was_visited[seq_id] and not seq_id == predecessor[current_node_id]:
					if predecessor_label < self.min_label_p:
						self.min_label_p = predecessor_label
						self.min_sequence_p = seq_id
					self.sequences[seq_id].label_p = predecessor_label
					was_visited[seq_id] = True
					if verbose >= 2:
						print ("updated label of node "+str(seq_id)+": "+str(predecessor_label))
					priority_queue.put((-predecessor_label, seq_id))
		
		if verbose == 2:
			for seq in self.sequences:
				if seq.is_relevant:
					print str(seq.id) + ": " + str(seq.label_p)
		
		if do_second_iteration:
			next_start = 0
			if abs(self.max_label_p) > abs(self.min_label_p):
				next_start = self.max_sequence_p
			else:
				next_start = self.min_sequence_p
			if verbose == 2:
				print ("max sequence of first iteration: "+str(self.max_sequence_p) + " with label "+str(self.max_label_p))
				print ("min sequence of first iteration: "+str(self.min_sequence_p) + " with label "+str(self.min_label_p))
				print ("Start a second iteration of labelling, starting from sequence "+str(next_start))
			
			self.construct_assembly_ordering_labels(start_sequence=next_start, do_second_iteration=False, verbose=verbose)
		
	def greedy_construct_assembly_ordering_labels(self, start_sequence = 0, do_second_iteration=True, compute_position_labels=True, verbose=1):
		# constructs a fuzzy partial ordering of all relevant sequences based on directed graph structure and sequence lengths:
		# algorithm assumes that graph
		# 	is not empty and
		#	has only one component and
		# 	has no cycles, i.e. implies a partial order
		if verbose > 0:
			print "Construct assembly ordering labels..."
		
		if not (start_sequence == 0 or self.sequences[start_sequence].is_relevant):
			if verbose > 0:
				print ("Error! Start sequence does not exist! Start with first sequence instead.")
			start_sequence = 0
		while start_sequence < len(self.sequences) and not self.sequences[start_sequence].is_relevant:
			start_sequence += 1
		
		if not start_sequence < len(self.sequences):
			print("Error! No sequence found!")
			return
		
		# reset labels
		for seq in self.sequences:
			if seq.is_relevant:
				seq.label_p = False
		# init for start sequence
		if compute_position_labels:
			self.min_label_p = 0
			self.max_label_p = 0
			self.min_sequence_p = start_sequence
			self.max_sequence_p = start_sequence
			self.max_path_length_p = self.sequences[start_sequence].get_length()
			self.sequences[start_sequence].label_p = 0
		else:
			self.min_label_n = 0
			self.max_label_n = 0
			self.min_sequence_n = start_sequence
			self.max_sequence_n = start_sequence
			self.sequences[start_sequence].label_n = 0

		if verbose == 2:
			print ("Start with sequence " + str(start_sequence))
				
		queue = [[self.sequences[start_sequence].id, 0]]
		while (len(queue) > 0):
			current_data = queue[0]
			queue.pop(0)
			current_node_id = current_data[0]
			current_position = current_data[1]

			for seq_id in self.sequences[current_node_id].overlaps_out:
				if compute_position_labels:
					if self.sequences[seq_id].label_p == False:
						start_label = current_position + self.sequences[current_node_id].get_length()
						if start_label > self.max_label_p:
							self.max_label_p = start_label
							self.max_sequence_p = seq_id
						end_label = start_label + self.sequences[seq_id].get_length()
						if end_label > self.max_path_length_p:
							self.max_path_length_p = end_label
						self.sequences[seq_id].label_p = start_label
						queue.append([seq_id, start_label])
				else:
					if self.sequences[seq_id].label_n == False:
						start_label = current_position + 1
						if start_label > self.max_label_n:
							self.max_label_n = start_label
							self.max_sequence_n = seq_id
						self.sequences[seq_id].label_n = start_label
						queue.append([seq_id, start_label])
			for seq_id in self.sequences[current_node_id].overlaps_in:
				if compute_position_labels:
					if self.sequences[seq_id].label_p == False:
						start_label = current_position - self.sequences[seq_id].get_length()
						if start_label < self.min_label_p:
							self.min_label_p = start_label
							self.min_sequence_p = seq_id
						self.sequences[seq_id].label_p = start_label
						queue.append([seq_id, start_label])
				else:
					if self.sequences[seq_id].label_n == False:
						start_label = current_position - 1
						if start_label < self.min_label_n:
							self.min_label_n = start_label
							self.min_sequence_n = seq_id
						self.sequences[seq_id].label_n = start_label
						queue.append([seq_id, start_label])
		
		if verbose == 2:
			for seq in self.sequences:
				if seq.is_relevant:
					if compute_position_labels:
						print str(seq.id) + ": " + str(seq.label_p)
					else:
						print str(seq.id) + ": " + str(seq.label_n)
		
		if do_second_iteration:
			next_start = 0
			if compute_position_labels:
				if abs(self.max_label_p) > abs(self.min_label_p):
					next_start = self.max_sequence_p
				else:
					next_start = self.min_sequence_p
				if verbose == 2:
					print ("max sequence of first iteration: "+str(self.max_sequence_p) + " with label "+str(self.max_label_p))
					print ("min sequence of first iteration: "+str(self.min_sequence_p) + " with label "+str(self.min_label_p))
					print ("Start a second iteration of labelling, starting from sequence "+str(next_start))
			else:
				if abs(self.max_label_n) > abs(self.min_label_n):
					next_start = self.max_sequence_n
				else:
					next_start = self.min_sequence_n
				if verbose == 2:
					print ("max sequence of first iteration: "+str(self.max_sequence_n) + " with label "+str(self.max_label_n))
					print ("min sequence of first iteration: "+str(self.min_sequence_n) + " with label "+str(self.min_label_n))
					print ("Start a second iteration of labelling, starting from sequence "+str(next_start))

			self.greedy_construct_assembly_ordering_labels(start_sequence=next_start, do_second_iteration=False, compute_position_labels=compute_position_labels, verbose=verbose)
					
	def get_partition_of_sequences(self, number_of_parts, overlap=3, verbose=False):
		# returns a set of sets of parallel aligned sequencs that cover the whole graph
		# each interior position of the graph is covered by mutliple different sets of sequences, defined by the parameter overlap
		sorted_nodes = sorted([seq for seq in self.sequences if seq.label_p], key=lambda x: x.label_p)
		label_div = self.max_label_p-self.min_label_p
		part_start_difference = label_div/(number_of_parts+overlap-1)
		part_size = part_start_difference*overlap
		
		if verbose:
			print label_div
			print part_size
		
		parts_seq = []
		for i in range(number_of_parts):
			current_start = self.min_label_p+(i*part_start_difference)
			if i == number_of_parts-1:
				current_end = self.max_label_p
			else:
				current_end = current_start + part_size
			this_part_sequences = [seq for seq in sorted_nodes if seq.label_p >= current_start and seq.label_p <= current_end]
			parts_seq.append(this_part_sequences)
		return parts_seq
	
	def get_read_of_sequences(self, sequences, verbose=False):
		# returns all reads that are contained in specific set of sequences
		reads = []
		if len(self.kmers) > 0:
			kmers = []
			for seq in sequences:
				kmers += seq.kmers
			for kmer_id in kmers:
				reads += self.kmers[kmer_id].evidence_reads
			reads = list(set(reads))
		return reads
		
	def get_relevant_reads(self, consider_hubread_length=3, verbose=False):
		# get all reads of sequences that are not marked as removed
		#  (because they induced a removed tip or a removed low-weight-sequence))
		# and reads have also to be evidence for a specified minimum number of sequences (at least four be default)
		#  (all other reads (that are evidence for at most three sequences) are subsequences of hubreads)
		
		print ("Get relevant reads:")
		print ("* Checking sequences ...")
		read_appearances = {}
		i = 0
		for s in self.sequences:
			if i%100 == 0 or i == len(self.sequences)-1:
				meta.print_progress(i, len(self.sequences)-1)
			evidence_reads = self.get_read_of_sequences([s])
			for r in evidence_reads:
				if r not in read_appearances:
					read_appearances[r] = 0
				read_appearances[r] += 1
			i += 1
				
		print ("* Grabbing relevant reads ...")
		list_of_relevant_reads = []
		i = 0
		for r in read_appearances:
			if i%100 == 0 or i == len(read_appearances)-1:
				meta.print_progress(i, len(read_appearances)-1)
			if read_appearances[r] > consider_hubread_length and r not in self.removed_reads:
				list_of_relevant_reads.append(r)
			i += 1
			
		if verbose:
			print ("The original set of "+str(len(self.reads))+" reads was reduced to "+str(len(list_of_relevant_reads))+" relevant reads")
		return list_of_relevant_reads
		
	def reduce_to_single_path_max_weight(self, start_sequence=False, restrict_to_component=False, verbose=False):
		# greedy algo that traverses through graph by choosing following nodes with max weight, deletes everythin else
		# with backtracking: if a dead end is reached, the algorithm returns to the previous node until a second-best path is found.
		# i.e. depth-first search, starting from node with smallest label to node with highest label where in each step the node with largest weight is chosen.
		# method assumes that graph has only one component
		# and sequences have position- and weight-labels
		
		if verbose > 0:
			print ("DFS for path with local maximum weight")
			
		# node-labels of the dfs:
		#	0: not yet visited
		#	1: node visited, is on current path
		#	2: node visited and returned because it leads to a dead end
		dfs_labels = [0 for s in self.sequences]
		path = []
		
		if not start_sequence:
			current_seq_id = self.min_sequence_p
		else:
			current_seq_id = start_sequence
		
		if verbose == 2:
			print ("Initial sequence is "+str(current_seq_id))
			
		current_label = self.sequences[current_seq_id].label_p
		if verbose == 2:
			print ("maximum possible path length is: "+str(self.max_path_length_p-current_label))
		
		path.append(current_seq_id)
		path_finding_failed = False
		while (not path_finding_failed) and current_label+self.sequences[current_seq_id].get_length() < self.max_path_length_p - 2*self.k_value:
			# go as far as possible, and if no successor, backtrack if not at least 95% of longest path is covered
			dfs_labels[current_seq_id] = 1
			next_sequences = [target_id for target_id in self.sequences[current_seq_id].overlaps_out if dfs_labels[target_id] == 0]
			if verbose == 2:
				print ("Possible next sequences are: "+str(next_sequences))
			#self.sequences[current_seq_id].print_data()
			if len(next_sequences) > 0:
				# go to next node:
				max_seq_id = -1
				max_seq_weight = -1
				for seq_id in next_sequences:
					if not seq_id == current_seq_id:
						if self.sequences[seq_id].is_relevant and self.sequences[seq_id].get_total_weight() > max_seq_weight:
							max_seq_weight = self.sequences[seq_id].get_total_weight()
							max_seq_id = seq_id
				
					current_seq_id = max_seq_id
					current_label  = self.sequences[current_seq_id].label_p
					path.append(current_seq_id)
				if verbose == 2:
					print ("Next sequence is sequence "+str(max_seq_id)+" with label " + str(current_label) + " and weight "+str(max_seq_weight))
			else:
				# backtrack:
				if len(path) > 1:
					dfs_labels[current_seq_id] = 2
					if verbose == 2:
						print ("Dead end! Remove sequence "+str(path[-1])+" from path")
					path.pop(-1)
					current_seq_id = path[-1]
					current_label  = self.sequences[current_seq_id].label_p
					if verbose == 2:
						print ("Next sequence is sequence "+str(current_seq_id)+" with label " + str(current_label) + " and weight "+str(self.sequences[current_seq_id].get_total_weight()))
				else:
					if verbose > 0:
						print ("Error! No path found")
					path_finding_failed = True
			#print ("current path: " +str(path))
	
		dfs_labels[current_seq_id] = 1
		
		# delete unused sequences:
		for seq in self.sequences:
			if seq.is_relevant and not dfs_labels[seq.id] == 1 and (not restrict_to_component or seq.id in restrict_to_component):
				self.delete_sequence(seq.id)
		if verbose > 0:
			print ("Reduction finished")
		
	def reduce_every_component_to_single_path_max_weight(self, do_contraction=True, verbose=False):
		print ("Reducing every component to a single path with maximum local weight ...")
		# applies the reduce_to_single_path_max_weight algorithm on each component of the graph
		# i.e. computes a consensus sequence for each component.
		
		components = self.get_components()
		n = len(components)
		if verbose:
			print ("Number of components: "+str(n))
		c_i = 0
		for c in components:
			meta.print_progress(c_i, n-1)
			if verbose:
				print ("Size of this component: "+str(len(c)))
			self.construct_assembly_ordering_labels(start_sequence=c[0], verbose = 0)
			start_sequence = -1
			min_label_p = False
			for seq_id in c:
				if not min_label_p or self.sequences[seq_id].label_p < min_label_p:
					start_sequence = seq_id
					min_label_p = self.sequences[seq_id].label_p
			if verbose:
				print ("start_sequence for reduction: "+str(start_sequence)+" with label: "+str(min_label_p))
			self.reduce_to_single_path_max_weight(start_sequence = start_sequence, restrict_to_component=c, verbose=0)
			c_i += 1
		
		if do_contraction:
			self.contract_unique_overlaps(verbose = verbose)
		self.construct_assembly_ordering_labels(verbose = 0)
		
	def greedy_reduce_to_single_path_max_weight(self, verbose=False):
		# greedy algo that traverses through graph by choosing following nodes with max weight, deletes everythin else
		# method assumes that graph has only one component // and no cycles
		# and sequences have position- and weight-labels
		
		print ("Greedy reduction to single path with local max weight")
		
		components = self.get_components()
		print ("Number of components: "+str(len(components)))
		component_id = 0
		while component_id < len(components):
			c = components[component_id]
			component_id += 1
			if verbose:
				print ("next component")
			start_seq_id = self.min_sequence_p
			
			current_seq_id = start_seq_id
			last_seq_id = -1
			if verbose:
				print "start_id: " +str(current_seq_id)
			while len(self.sequences[current_seq_id].overlaps_out) > 0:
				if verbose:
					print ("Current seq: "+str(current_seq_id))
				next_sequences = []
				for target_id in self.sequences[current_seq_id].overlaps_out:
					next_sequences.append(target_id)
				max_seq_id = -1
				max_seq_weight = -1
				for seq_id in next_sequences:
					if not seq_id == current_seq_id:
						if self.sequences[seq_id].is_relevant and self.sequences[seq_id].get_total_weight() > max_seq_weight:
							max_seq_weight = self.sequences[seq_id].get_total_weight()
							max_seq_id = seq_id
				
				if verbose:
					print ("Next sequence is sequence "+str(max_seq_id)+" with weight "+str(max_seq_weight))
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
					else:#if not branching:
						self.delete_sequence(seq_id, verbose=False)
			# delete incoming sequences at final sequence, if there are any:
			if not last_seq_id == -1:
				incoming_sequences = [seq_id for seq_id in self.sequences[current_seq_id].overlaps_in]
				for seq_id in incoming_sequences:
					if not seq_id == last_seq_id:
						self.delete_sequence(seq_id, verbose=False)
		
	def reduce_to_single_largest_component(self, verbose=False):
		# deletes all sequences that are not part of the largest (by number of sequences) component
		# thus only the largest component will remain.
		
		print ("Reducing graph to largest component...")
		components = self.get_components()
		if len(components) > 1:
			max_size = 0#len(components[0])
			max_comp_id = 0
			for i in range(len(components)):
				c = components[i]
				if len(c) > max_size:
					max_size = len(c)
					max_comp_id = i
				if verbose:
					print ("Component "+str(i)+" consists of "+str(len(c))+" sequences.")
					
			if verbose:
				print ("Largest component is "+str(max_comp_id)+" with "+str(max_size)+" sequences.")
			for i in range(len(components)):
				if not i == max_comp_id:
					for seq_id in components[i]:
						self.delete_sequence(seq_id)
					
	def remove_short_sequences(self, length_bound_by_multiple_of_k=5):
		# brutally removes all sequences with length less than 5 times k_value (if not specified otherwise)
		for seq in self.sequences:
			if seq.is_relevant:
				if seq.get_length() < length_bound_by_multiple_of_k*self.k_value:
					self.delete_sequence(seq.id)
	
	def get_label_span(self):
		# get an approximation of the length of the genome
		# this depends on the quality of the labels. if graph containes cycles, this may be significantly incorrect.
		return self.max_label_p - self.min_label_p
		
	def construct_laplacian(self, list_of_nodes, verbose=False):
		# this method constructs the laplacian for a specific subset of nodes
		# maps original sequence_ids to smaller id for only relevant sequences (-> matrix position)
		seq_id_to_index = {}
		index_to_seq_id = {}
		
		# construct symmetric (i.e. undirected) adjacency matrix and diagonal matrix of vertex-degrees:
		n = 0
		x = []
		y = []
		degrees = [0.0]*len([s_id for s_id in list_of_nodes if self.sequences[s_id].is_relevant])
		for seq_s_id in list_of_nodes:
			if self.sequences[seq_s_id].is_relevant:
				if seq_s_id not in seq_id_to_index:
					seq_id_to_index[seq_s_id] = n
					index_to_seq_id[n] = seq_s_id
					n += 1
				s = seq_id_to_index[seq_s_id]
				for seq_t_id in self.sequences[seq_s_id].overlaps_out:
					if seq_t_id in list_of_nodes and self.sequences[seq_t_id].is_relevant:
						if seq_t_id not in seq_id_to_index:
							seq_id_to_index[seq_t_id] = n
							index_to_seq_id[n] = seq_t_id
							n += 1
						t = seq_id_to_index[seq_t_id]
						degrees[s] += 1
						degrees[t] += 1
						x.append(s)
						y.append(t)
						x.append(t)
						y.append(s)
		data = [1.0]*len(x)
		
		adj_mat = scipy.sparse.coo_matrix((data, (x,y)), shape = (n,n)).tocsc()
		deg_mat = scipy.sparse.diags(degrees).tocsc()
		
		# laplacian L = D - A
		laplacian = -(adj_mat-deg_mat)
		
		return laplacian, seq_id_to_index, index_to_seq_id
		
	def construct_spectral_clusters(self, laplacian, index_to_seq_id, verbose=False):
		# apply the spectral cluster algorithm on the de bruijn graph with a pre computed laplacian matrix of the graph
		part_a = []
		part_b = []
		part_c = []
		
		try:
			[w,v] = la.eigs(laplacian, which='SM', k=2)
		
			min_eigenvalue = -1
			min_eigenvalue_id = -1
			secmin_eigenvalue = -1
			secmin_eigenvalue_id = -1
			
			for i in range(len(w)):
				if min_eigenvalue_id < 0 or w[i] < min_eigenvalue:
					secmin_eigenvalue = min_eigenvalue
					secmin_eigenvalue_id = min_eigenvalue_id
					min_eigenvalue = w[i]
					min_eigenvalue_id = i
				elif secmin_eigenvalue_id < 0 or w[i] < secmin_eigenvalue:
					secmin_eigenvalue = w[i]
					secmin_eigenvalue_id = i
				
			secmin_eigenvector = v[:,i]
			
			for i in range(len(secmin_eigenvector)):
				if self.sequences[index_to_seq_id[i]].is_relevant:
					if secmin_eigenvector[i] < -10e-8:
						part_a.append(i)
					elif secmin_eigenvector[i] > 10e-8:
						part_b.append(i)
					else:
						part_c.append(i)
			
			if verbose:
				print ("Number of nodes in first part: "+str(len(part_a)))
				print ("Number of nodes in second part: "+str(len(part_b)))
				print ("Number of nodes in neither part: "+str(len(part_c)))
						
		
		except (la.ArpackNoConvergence):
			print ("No convergence during computation of eigenvalues!")
			print ("No partition computed.")
			
			secmin_eigenvalue = -1
			for i in range(laplacian.shape[0]):
				part_a.append(i)
					
		return secmin_eigenvalue, part_a, part_b, part_c
		
	def compute_mincut(self, component=-1, divide_clusters=True, verbose=False):
		# compute the spectral cut of the de bruijn graph.
		if verbose:
			print ("Do spectral clustering of the nodes into two clusters ...")
		absoulute_minimum_component_size = 1000
		
		self.remove_irrelevant_overlaps()
		
		if component < 0:
			component = [seq.id for seq in self.sequences]
		if len(component) > absoulute_minimum_component_size:
			if verbose:
				print ("Size of component: "+str(len(component)))
				print ("Component is large enough to consider a spectral cut")
			laplacian, seq_id_to_index, index_to_seq_id = self.construct_laplacian(component)
			
			secmin_eigenvalue, part_a, part_b, part_c = self.construct_spectral_clusters(laplacian, index_to_seq_id, verbose)
			
			min_part_size = max(absoulute_minimum_component_size, math.sqrt(len(component)))
			
			if len(part_a) > min_part_size and len(part_b) > min_part_size:
				# only consider a cut if parts of decomposition have significant size.
				# compute eigenvalues of parts to check if decomposition increases clustering
				laplacian_a, seq_id_to_index_a, index_to_seq_id_a = self.construct_laplacian([index_to_seq_id[i] for i in part_a])
				secmin_eigenvalue_a, part_a_a, part_b_a, part_c_a = self.construct_spectral_clusters(laplacian_a, index_to_seq_id_a, verbose)
				
				laplacian_b, seq_id_to_index_b, index_to_seq_id_b = self.construct_laplacian([index_to_seq_id[i] for i in part_b])
				secmin_eigenvalue_b, part_a_b, part_b_b, part_c_b = self.construct_spectral_clusters(laplacian_b, index_to_seq_id_b, verbose)
				
				if verbose:
					print ("lambda_2 of this component: "+str(secmin_eigenvalue))
					print ("lambda_2 of decomposition: "+str(secmin_eigenvalue_a)+","+str(secmin_eigenvalue_b))
			else:
				divide_clusters = False
			
			# only cut if decomposition significantly increases clustering-coefficient:
			# one of the second smallest eigenvalues has to decrese by a factor
			decreasing_factor_bound = 5
			if divide_clusters and (secmin_eigenvalue*decreasing_factor_bound < secmin_eigenvalue_a or secmin_eigenvalue*decreasing_factor_bound < secmin_eigenvalue_b):
				if verbose:
					print ("Cut graph according to partitions.")
				self.cut_graph_into_partitions([part_a, part_b] ,seq_id_to_index)
				return True, part_a, part_b
			else:
				if verbose:
					print ("Do not cut graph.")
				return False, part_a, part_b
		else:
			return False, component, []
			
	def partition_graph_into_components_of_clusters(self, verbose=False):
		print ("Decompose graph into components of clusters ...")
		components_to_potentially_cut = self.get_components()
		#number_of_clustered_components = 0
		while (len(components_to_potentially_cut) > 0):
			print "current number of components left to consider: " + str("%6d" % len(components_to_potentially_cut)) + "\r",
			sys.stdout.flush()
			res, part_a, part_b = self.compute_mincut(components_to_potentially_cut[0], verbose=verbose)
			components_to_potentially_cut.pop(0)
			if res:
				components_to_potentially_cut.append(part_a)
				components_to_potentially_cut.append(part_b)
		
	def cut_graph_into_partitions(self, partitions, seq_id_to_index, verbose=False):
		if verbose:
			print ("Delete all overlaps between different parts ...")
		id_to_part={}
		for i in range(len(partitions)):
			for id in partitions[i]:
				id_to_part[id]=i
		
		overlaps_to_delete = []
		for ov in self.overlaps:
			s1_id = self.overlaps[ov].contig_sequence_1
			s2_id = self.overlaps[ov].contig_sequence_2
			s1_index = seq_id_to_index[s1_id]
			s2_index = seq_id_to_index[s2_id]
			if not id_to_part[s1_index] == id_to_part[s2_index]:
				overlaps_to_delete.append(ov)
		
		for ov_id in set(overlaps_to_delete):
			self.delete_overlap(ov_id)
	
	def remove_irrelevant_overlaps(self):
		# remove overlaps from self.overlaps that are marked is_relevant==False
		ov_to_del = []
		for ov in self.overlaps:
			if not self.sequences[self.overlaps[ov].contig_sequence_1].is_relevant or not self.sequences[self.overlaps[ov].contig_sequence_2].is_relevant:
				ov_to_del.append(ov)
		for ov in ov_to_del:
			self.overlaps.pop(ov)
			
	def get_hubreads_by_adjacent_sequences(self, id_restriction=[], verbose=False):
		# constructs extended hubreads
		# these are all directed paths of length three in the contracted (simplified) debruijn graph.
		# if id_restriction is not empty, only hubreads are considered for which the id of the first sequence is in id_restriction
		# hubread gets weight = minimum of the weigths of all contained sequences
		if not self.is_contracted:
			self.contract_unique_overlaps()
			
		if verbose:
			print ("get hubreads from graph")
			
		hubreads = []
		for seq_id in range(len(self.sequences)):
			if seq_id % 1000 == 0 or seq_id == len(self.sequences)-1:
				meta.print_progress(seq_id, len(self.sequences)-1)
			seq = self.sequences[seq_id]
			if seq.is_relevant:
				if len(id_restriction) == 0 or seq.id in id_restriction:
					hp_start = seq.sequence
					hp_start_id = seq.id
					hp_start_weight = seq.max_weight
					for ov_out_seq_id in seq.overlaps_out.keys():
						hp_mid = self.sequences[ov_out_seq_id]
						hp_mid_id = hp_mid.id
						hp_mid_weight = hp_mid.max_weight
						hp_temp = meta.merge_sequences(hp_start, hp_mid.sequence, self.k_value-1)
						if len(hp_mid.overlaps_out) > 0:
							for ov_out_seq_id in hp_mid.overlaps_out.keys():
								hp_end = self.sequences[ov_out_seq_id]
								hp_end_weight = hp_end.max_weight
								hp_end_id = hp_end.id
								hubread = meta.merge_sequences(hp_temp, hp_end.sequence, self.k_value-1)
								hubread_weight = min(hp_start_weight, hp_mid_weight, hp_end_weight)
								if verbose:
									print ("hp_start: "+hp_start+"("+str(hp_start_id)+")")
									print ("hp_mid: "+hp_mid.sequence+"("+str(hp_mid_id)+")")
									print ("hp_end: "+hp_end.sequence+"("+str(hp_end_id)+")")
									print ("New hubread: "+hubread)
								hubreads.append(hubread+","+str(hubread_weight))
		return hubreads
	
	def get_hubreads_by_overlaps(self, verbose=False):
		# constructs hubreads as described in spades-paper:
		# these are all directed paths of length two in the contracted (simplified) debruijn graph.
		# hubread gets weight = minimum of the weigths of all contained sequences
		if not self.is_contracted:
			self.contract_unique_overlaps()
			
		if verbose:
			print ("get hubreads from graph")
			
		hubreads = []
		for ov_id in self.overlaps:
			seq_start = self.sequences[self.overlaps[ov_id].contig_sequence_1].sequence
			seq_end = self.sequences[self.overlaps[ov_id].contig_sequence_2].sequence
			hubread = meta.merge_sequences(seq_start, seq_end, self.k_value-1)
			hubread_weight = min(self.sequences[self.overlaps[ov_id].contig_sequence_1].max_weight, self.sequences[self.overlaps[ov_id].contig_sequence_2].max_weight)
			hubreads.append(hubread+","+hubread_weight)
						
		return hubreads

	def compute_common_genome_evidence(self, verbose=False):
		common_genome_evidence = {}

		sequences_by_node_labels = {}
		for seq in self.sequences:
			if seq.label_n not in sequences_by_node_labels:
				sequences_by_node_labels[seq.label_n] = []
			sequences_by_node_labels[seq.label_n].append(seq.id)

		for seq in self.sequences:
			seq_read_ids = self.get_read_of_sequences([seq])
			for other_seq_id in sequences_by_node_labels[seq.label_n]:
				other_seq_read_ids = self.get_read_of_sequences([self.sequences[other_seq_id]])
				for r_id in seq_read_ids:
					if r_id not in common_genome_evidence:
						common_genome_evidence[r_id] = {}
					for or_id in other_seq_read_ids:
						common_genome_evidence[r_id][or_id] = -1
		return common_genome_evidence

	def split_sequences_with_contradicting_genome_evidence(self, common_genome_evidence, verbose=False):
		# # use common_genome_evidence as computed by compute_common_genome_evidence()
		# for each sequence:
			# get all read_ids from all kmers of sequence
			# # 1) partition read_ids into minimal number of sets such that there are no pairs of read_ids with negative common_genome_evidence within any set:
			# # 2) split sequence into multiple sequences, one for each set. split incident edges accordingly.
		# # remove all edges without read_evidence (i.e. if source and target edge have no read_ids in common)
		return

	def compute_overlap_evidence_distribution(self):
		print ("Compute distribution of overlap evidences")
		overlap_weights = [self.overlaps[ov].get_evidence_weight() for ov in self.overlaps]
		overlap_weight_distribution = [0 for i in range(max(overlap_weights)+1)]
		for ovw in overlap_weights:
			overlap_weight_distribution[ovw] += 1
		
		return overlap_weight_distribution
		
	'''
	def get_bound_on_erroneous_overlap_weights(self):
		# computes the most-likely divider between the evidence weights of erroneous overlaps and correct overlaps:
		
		ovw = self.compute_overlap_evidence_distribution()
		div = 2
		while ovw[div] > ovw[div+1]:
			div += 1
			
		if div < len(ovw)/2:
			return div
		else:
			return 2
	'''
			
	def check_if_graph_decomposes_edgeremoval(self, overlaps_to_remove, relative_component_size_bound=0.05, verbose=False):
		# method checks if graph decomposes into multiple major components, if a set of overlaps is removed
		print ("Checking if graph decomposes ...")
		
		# same code as in get_components, but ignores overlaps that are defined by overlaps_to_remove:
		components = []
		visited_sequences = {seq.id: False for seq in self.sequences if seq.is_relevant}
		for seq in self.sequences:
			if seq.is_relevant:
				if not visited_sequences[seq.id]:
					current_comp = []
					seq_stack = [seq.id]
					visited_sequences[seq.id] = True
					while len(seq_stack) > 0:
						current_seq_id = seq_stack[0]
						seq_stack.pop(0)
						current_comp.append(current_seq_id)
						for adj_seq in self.sequences[current_seq_id].overlaps_out:
							if not visited_sequences[adj_seq] and not self.sequences[current_seq_id].overlaps_out[adj_seq] in overlaps_to_remove:
								seq_stack.append(adj_seq)
								visited_sequences[adj_seq] = True
						for adj_seq in self.sequences[current_seq_id].overlaps_in:
							if not visited_sequences[adj_seq] and not self.sequences[current_seq_id].overlaps_in[adj_seq] in overlaps_to_remove:
								seq_stack.append(adj_seq)
								visited_sequences[adj_seq] = True					
					components.append(current_comp)
		
		# check if largest component has siginificant size:
		if len(components) == 0:
			max_comp_size = 0
		else:
			max_comp_size = 0
			second_max_comp_size = 0
			for c in components:
				s = len(c)
				if s > max_comp_size:
					second_mac_comp_size = max_comp_size
					max_comp_size = s
				elif s > second_max_comp_size:
					second_max_comp_size = s
			#comp_sizes = [len(c) for c in components]
			#max_comp_size = max([len(c) for c in components])
		
		if verbose:
			print ("Number of Components after overlaps have been removed: "+str(len(components)))
			print ("Max component size: "+str(max_comp_size))
			print ("Second max component size: "+str(second_max_comp_size))
			#print ("Bound on comp size: "+str(relative_component_size_bound)+" * "+str(len(self.get_relevant_sequences()))+" = "+str(relative_component_size_bound*len(self.get_relevant_sequences())))
		
		if second_max_comp_size > relative_component_size_bound*max_comp_size:
			return True
		else:
			return False
			
	def check_if_graph_decomposes_seqremove(self, sequences_to_remove, relative_component_size_bound=0.05, verbose=False):
		# method checks if graph decomposes into multiple major components, if a set of overlaps is removed
		print ("Checking if graph decomposes ...")
		
		# same code as in get_components, but ignores sequences that are defined by sequences_to_remove:
		components = []
		visited_sequences = {seq.id: False for seq in self.sequences if seq.is_relevant}
		for seq in self.sequences:
			if seq.is_relevant:
				if not visited_sequences[seq.id]:
					current_comp = []
					seq_stack = [seq.id]
					visited_sequences[seq.id] = True
					while len(seq_stack) > 0:
						current_seq_id = seq_stack[0]
						seq_stack.pop(0)
						current_comp.append(current_seq_id)
						for adj_seq in self.sequences[current_seq_id].overlaps_out:
							if not visited_sequences[adj_seq] and not adj_seq in sequences_to_remove:
								seq_stack.append(adj_seq)
								visited_sequences[adj_seq] = True
						for adj_seq in self.sequences[current_seq_id].overlaps_in:
							if not visited_sequences[adj_seq] and not adj_seq in sequences_to_remove:
								seq_stack.append(adj_seq)
								visited_sequences[adj_seq] = True					
					components.append(current_comp)
		
		# check if largest component has siginificant size:
		if len(components) == 0:
			max_comp_size = 0
		else:
			max_comp_size = 0
			second_max_comp_size = 0
			for c in components:
				s = len(c)
				if s > max_comp_size:
					second_mac_comp_size = max_comp_size
					max_comp_size = s
				elif s > second_max_comp_size:
					second_max_comp_size = s
		
		if verbose:
			print ("Number of Components after overlaps have been removed: "+str(len(components)))
			print ("Max component size: "+str(max_comp_size))
			print ("Second max component size: "+str(second_max_comp_size))
		
		if second_max_comp_size > relative_component_size_bound*max_comp_size:
			return True
		else:
			return False
		
	def remove_low_evidence_overlaps_until_graph_decomposes(self, relative_component_size_bound=0.05, verbose=False):
		# remove overlaps with low evidence until graph decomposes
		print ("Remove low-evidence overlaps until graph decomposes:")
			
		if len(self.overlaps) == 0:
			# if graph is already trivial, return
			return
			
		#initilization:
		min_cov_evidence = 2
		# leave at least 10% of the initial overlaps:
		min_overlap_number = len(self.overlaps)/10
		
		# compute set of edges to remove in next iteration:
		overlaps_with_small_evidence = [ov_id for ov_id in self.overlaps if self.overlaps[ov_id].get_evidence_weight() < min_cov_evidence]
		# compute if graph decomposes:
		graph_decomposees = self.check_if_graph_decomposes_edgeremoval(overlaps_with_small_evidence, relative_component_size_bound, verbose=verbose)
		while not graph_decomposees and (len(self.overlaps)-len(overlaps_with_small_evidence)) > min_overlap_number:
			#print ("Now deleting overlaps with evidence < "+str(min_cov_ecidence))
			if verbose:
				print ("Next set of overlaps to remove:")
				for ov in overlaps_with_small_evidence:
					self.overlaps[ov].print_data()
				
			# remove the overlaps:
			self.remove_insignificant_overlaps(minimal_evidence=min_cov_evidence, keep_relevant_tips=True)
			# remove decomposed parts:
			self.reduce_to_single_largest_component()
			# remove new tips:
			self.remove_tips(maximum_tip_weight = min_overlap_number)
			# simplify graph:
			self.contract_unique_overlaps()
			#self.remove_single_sequence_components()
	
			if len(overlaps_with_small_evidence) < 10:
				# if last iteration had very little effect, increase the cutoff bound by 5
				min_cov_evidence += 5
			else:
				# normal case: increase cutoff bound by 1
				min_cov_evidence += 1
			
			# compute set of edges to remove in next iteration:
			overlaps_with_small_evidence = [ov_id for ov_id in self.overlaps if self.overlaps[ov_id].get_evidence_weight() < min_cov_evidence]
			# compute if graph decomposes:
			graph_decomposees = self.check_if_graph_decomposes_edgeremoval(overlaps_with_small_evidence, verbose=verbose)
			
		print ("Removal of low evidence overlaps stopped at min_cov_evidence = " +str(min_cov_evidence))
			
	def remove_low_coverage_sequences_until_graph_decomposes(self, relative_component_size_bound=0.05, max_coverage_depth_to_remove = 20, verbose=False):
		# remove sequences with low coverage depth until graph decomposes
		print ("Remove low-coverage sequences until graph decomposes:")
		
		if len(self.overlaps) == 0:
			# if graph is already trivial, return
			return
			
		# initilization:
		min_cov_depth = 2
		# compute partition of nodes into tips and non-tips (hubs)
		self.get_tip_classification()
		# compute set of nodes to remove in next iteration, do not consider tips:
		small_cov_depth_sequences = [seq.id for seq in self.sequences if seq.is_relevant and seq.get_total_weight() < min_cov_depth and seq.id in self.hubs]
		
			
		# increase cutoff bound untill there is at least one sequence to remove:
		while len(small_cov_depth_sequences) == 0:
			min_cov_depth += 1
			small_cov_depth_sequences = [seq.id for seq in self.sequences if seq.is_relevant and seq.get_total_weight() < min_cov_depth]
		
		# compute if graph decomposes:
		graph_decomposees = self.check_if_graph_decomposes_seqremove(small_cov_depth_sequences, relative_component_size_bound, verbose=verbose)	
		while not graph_decomposees and min_cov_depth < max_coverage_depth_to_remove:# and (len(self.overlaps)-len(overlaps_with_small_evidence)) > min_overlap_number:
			#print ("Now deleting overlaps with evidence < "+str(min_cov_ecidence))
			if verbose:
				print ("Graph does not decompose.")
				print ("Number of next sequences to remove: "+ str(len(small_cov_depth_sequences)))
				
			# remove the sequences:
			self.remove_insignificant_sequences(minimal_weight=min_cov_depth)
			# remove decomposed components:
			self.reduce_to_single_largest_component()
			# remove new tips:
			self.remove_tips()
			# simplify graph:
			self.contract_unique_overlaps()
	
			if len(small_cov_depth_sequences) < 10:
				# if last iteration had very little effect, increase the cutoff bound by 5
				min_cov_depth += 5
			else:
				# normal case: increase cutoff bound by 1
				min_cov_depth += 1
				
			# compute set of nodes to remove in next iteration, do not consider tips:
			small_cov_depth_sequences = [seq.id for seq in self.sequences if seq.is_relevant and seq.get_total_weight() < min_cov_depth and seq.id in self.hubs]
			
			# increase cutoff bound untill there is at least one sequence to remove:
			while len(small_cov_depth_sequences) == 0:
				min_cov_depth += 1
				small_cov_depth_sequences = [seq.id for seq in self.sequences if seq.is_relevant and seq.get_total_weight() < min_cov_depth]
			
			# compute if graph decomposes:
			graph_decomposees = self.check_if_graph_decomposes_seqremove(small_cov_depth_sequences, relative_component_size_bound, verbose=verbose)
			
		print ("Removal of low coverage sequences stopped at min_cov_depth = "+str(min_cov_depth))
				
	def get_tip_classification(self, verbose=False):
		#get all ids of sequences that are tips, and classify them as 'deletable', if their parent node has other adjacent nodes in the same orientation that are not tips or have significant higher coverage depth:
		
		hubs = []
		tips_deletable = []
		tips_keep = []
		
		visited_sequences = [False for seq in self.sequences]
		for i in range(len(self.sequences)):
			if not self.sequences[i].is_relevant:
				visited_sequences[i] = True
			elif not visited_sequences[i]:
				seq = self.sequences[i]
				if len(seq.overlaps_in) > 0 and len(seq.overlaps_out) > 0:
					# sequence is not a tip:
					visited_sequences[i] = True
					hubs.append(i)
				else:
					local_tips = []
					other_non_tip_exists = False
					if len(seq.overlaps_in) > 0:
						# sequence is an outgoing tip:
						for adj_seq in seq.overlaps_in:
							visited_sequences[adj_seq] = True
							# check all sequences that have edge into this tip:
							# for all these sequences, check if any outgoing sequence is not a tip or has a significant higher evidence weight:
							for adj_adj_seq in self.sequences[adj_seq].overlaps_out:
								visited_sequences[adj_adj_seq] = True
								if len(self.sequences[adj_adj_seq].overlaps_out) > 0:
									# there is another sequence that is not a tip
									other_non_tip_exists = True
									hubs.append(adj_adj_seq)
								else:
									local_tips.append(adj_adj_seq)
							
					else:
						# sequence is an incoming tip:
						for adj_seq in seq.overlaps_out:
							visited_sequences[adj_seq] = True
							# check all sequences that have edge into this tip:
							# for all these sequences, check if any outgoing sequence is not a tip or has a significant higher evidence weight:
							for adj_adj_seq in self.sequences[adj_seq].overlaps_in:
								visited_sequences[adj_adj_seq] = True
								if len(self.sequences[adj_adj_seq].overlaps_in) > 0:
									# there is another sequence that is not a tip
									other_non_tip_exists = True
									hubs.append(adj_adj_seq)
								else:
									local_tips.append(adj_adj_seq)
					
					if other_non_tip_exists:
						tips_deletable += local_tips
					elif len(local_tips) > 0:
						# mark all tips as 'deletable' that have less than 1/2 of the maximal coverage depth
						maxweight = max([self.sequences[tip_id].get_total_weight() for tip_id in local_tips])
						local_tips_to_delete = [tip_id for tip_id in local_tips if self.sequences[tip_id].get_total_weight() < maxweight/1.5]
						tips_deletable += local_tips_to_delete
						tips_keep += [tip_id for tip_id in local_tips if not tip_id in local_tips_to_delete]
						
		self.hubs = hubs
		self.tips_to_keep = tips_keep
		self.tips_to_delete = tips_deletable
		if verbose:
			print ("number of tips to keep: "+str(len(tips_keep)))
			print ("number of tips to delete: "+str(len(tips_deletable)))