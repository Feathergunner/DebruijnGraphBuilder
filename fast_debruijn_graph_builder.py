#!usr/bin/python

import random
import re
import sys

import numpy as np
import scipy.sparse
from scipy.sparse import linalg as la
import math

#import gc

'''
def print_progress(part, total):
	print ("Progress: "+str("%.2f" % ((float(part)/(float(total)/100)))) + "%")
'''
	
def print_progress(part, total, front_string="Progress:", end_string=""):
	if not total == 0:
		print front_string+" "+str("%6.2f" % ((float(part)/(float(total)/100)))) + "% "+end_string+"\r",
		if part >= total:
			print 
		sys.stdout.flush()
	
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
		if self.baseweight < weight:
			self.baseweight = weight

class ContigSequence:
	# Nodes in the Debruijn-Graph
	def __init__(self, seq_id, inv_id, sequence, kmers, weight=1, is_relevant=True):
		self.id = seq_id
		self.id_of_inverse_seq = inv_id
		self.sequence = sequence
		self.kmers = kmers
		# the maximal read eavidence this sequence has for any subsequence
		# self.max_evidence_weight = len(kmers)
		self.max_weight = weight
		# overlaps (i.e. edges) are stored in dictionaries
		# self.overlap[other_sequence_id] = overlap_id
		self.overlaps_out = {}
		self.overlaps_in = {}
		# deleted sequences will have the following flag set to False
		self.is_relevant = is_relevant
		# used for estimation of position within assembly
		self.label = False
		
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
	def __init__(self, ov_id, length, seq_1, seq_2, evidence_reads):
		self.id = ov_id
		self.length = length
		# the incident sequences of this overlap. sequence_1 is the source, sequence_2 is the target.
		self.contig_sequence_1 = seq_1
		self.contig_sequence_2 = seq_2
		self.evidence_reads = evidence_reads
		# Flag that checks whether a sequence is relevant for the graph or if it has been deleted
		self.is_relevant = True
		self.bla = 0
	
	def add_evidence(self, read_id):
		if read_id >= 0:
			self.evidence_reads.append(read_id)
		
	def get_evidence_weight(self):
		return len(self.evidence_reads)

	def print_data(self):
		print ("Overlap "+str(self.id))
		print ("From seq. "+str(self.contig_sequence_1)+" to seq. "+str(self.contig_sequence_2))
		if self.is_relevant:
			print ("is relevant")
		else:
			print ("is not relevant")
	
class GraphData:
	def __init__(self, reads=0, k=0, alphabet={"A":"T", "C":"G", "G":"C", "T":"A"}, load_weights=True, verbose=False):
		print ("Creating empty graph ...")
		
		self.k_value = k
		self.reads = []
		self.kmers = []
		self.kmer_dict = {}
		self.sequences = []
		self.overlaps = {}
		self.alphabet = alphabet
		
		# Flag that marks if inverse sequences have been removed
		self.is_unified = False
		# min and max label of all sequences
		self.min_label = 0
		self.max_label = 0
		self.min_sequence = -1
		self.max_sequence = -1
		
		if not reads == 0:
			self.init_graph_database(reads, load_weights=load_weights, verbose=verbose)

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
		if verbose:
			print ("Construct read database")
			
		read_id = 0
		number_of_reads = sum((len(r) for r in reads))
		if verbose:
			print ("Number of reads: "+str(number_of_reads))
		for r in reads:
			for read in r:
				if (read_id%1000 == 0 or read_id == number_of_reads-1):
					print_progress(read_id, number_of_reads-1)
				readdata = re.split(r',',read)
				readseq = readdata[0]
				if load_weights and len(readdata) > 1:
					readweight = int(readdata[1])
				else:
					readweight = 1
				# check if read has correct alphabet:
				is_correct = True
				for c in readseq:
					if c not in self.alphabet:
						is_correct = False
						print ("Error! Character "+str(c)+" not in alphabet "+str(self.alphabet.keys()))
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
				print_progress(kmer.id, len(self.kmers)-1)
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
			if (read.id%1000 == 0 or read_number == len(self.reads)-1):
				print_progress(read_number, len(self.reads)-1)
			for kmer_index in range(len(read.kmers)-1):
				source_kmer_id = read.kmers[kmer_index]
				target_kmer_id = read.kmers[kmer_index+1]
				self.increment_overlap(source_kmer_id, target_kmer_id, read.id, verbose=False)
			read_number += 1
		if verbose:
			self.print_memory_usage()

	def print_all_reads(self):
		print ("Reads:")
		for read in self.reads:
			print (str(read.id) + ": "+read.sequence)
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
		# checks all reads and constructs kmers and inverse kmers
		print ("Get kmer-data from reads ...")
		read_index = 0
		kmer_counter = 0
		
		progress_step = max(1, int(math.ceil(len(self.reads)/1000)))
		for read_index in range(len(self.reads)):
			if (read_index%progress_step == 0 or read_index == range(len(self.reads)-1)) and not verbose:
				print_progress(read_index, len(self.reads)-1)
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
					self.kmers[this_kmer_id].update_baseweight(current_read_weight)
					inv_kmer_id = self.kmers[this_kmer_id].id_of_inverse_kmer
					self.kmers[inv_kmer_id].add_evidence(read_index)
					self.kmers[inv_kmer_id].update_baseweight(current_read_weight)
					if verbose:
						print ("Kmer already exists")
						print ("Add read ("+str(read_index)+") evidence to kmer "+str(this_kmer_id))
				if not kmer_already_existing:
					if verbose:
						print ("Kmer does not exist in database. Add new kmer ...")
					
					inv_kmer_seq = get_inverse_sequence(new_kmer_sequence, self.alphabet)
					if not inv_kmer_seq == new_kmer_sequence:
						# add kmer:
						self.kmers.append(Kmer(kmer_counter, kmer_counter+1, new_kmer_sequence, [read_index], baseweight=current_read_weight))
						self.kmer_dict[new_kmer_sequence] = kmer_counter
						this_kmer_id = kmer_counter
						kmer_counter += 1
						# add inverse kmer:
						self.kmers.append(Kmer(kmer_counter, kmer_counter-1, inv_kmer_seq, [read_index], baseweight=current_read_weight))
						self.kmer_dict[inv_kmer_seq] = kmer_counter
						kmer_counter += 1
					else:
						# add kmer:
						self.kmers.append(Kmer(kmer_counter, kmer_counter, new_kmer_sequence, [read_index], baseweight=current_read_weight))
						self.kmer_dict[new_kmer_sequence] = kmer_counter
						this_kmer_id = kmer_counter
						kmer_counter += 1
					
				if verbose:
					print ("Add kmer "+self.kmers[this_kmer_id].sequence+"("+str(this_kmer_id)+") to read "+str(read_index))
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
			self.overlaps[ov_id] = (SequenceOverlap(ov_id, self.k_value-1, source_seq_id, target_seq_id, []))
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
				self.overlaps[rev_ov_id] = SequenceOverlap(rev_ov_id, self.k_value-1, source_rev_seq_id, target_rev_seq_id, [])
				self.sequences[source_rev_seq_id].overlaps_out[target_rev_seq_id] = rev_ov_id
				self.sequences[target_rev_seq_id].overlaps_in[source_rev_seq_id] = rev_ov_id
			
		else:
			ov_id = self.sequences[source_seq_id].overlaps_out[target_seq_id]
		self.overlaps[ov_id].add_evidence(read_evidence)
	
	def add_overlaps_for_sequences_with_small_insert_distance(self, max_insert_distance=1, verbose = False):
		print ("Adding overlaps to sequences with small insert distance ...")
		for sequence_id_1 in range(len(self.sequences)):
			print_progress(sequence_id_1, len(self.sequences))
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
							if compute_insert_distance(s1[l1-self.k_value:], s2[:self.k_value], maxdist = max_insert_distance) <= max_insert_distance:
								self.increment_overlap(sequence_id_1, sequence_id_2, -1)
							#if not sequence_id_1 in self.sequences[sequence_id_2].overlaps_in:
							if compute_insert_distance(s2[l2-self.k_value:], s1[:self.k_value], maxdist = max_insert_distance) <= max_insert_distance:
								self.increment_overlap(sequence_id_2, sequence_id_1, -1)

	def contract_unique_overlaps(self, verbose = False):
		# This method contracts all overlaps between adjacent sequences that form an unique path
		# (i.e., there are no other outgoing or incoming overlaps between the sequences)
		if not self.is_unified:
			self.remove_parallel_sequences(verbose)
		print ("Contract overlaps ...")
		
		ov_index_list = [ov_id for ov_id in self.overlaps]
		num_deleted_overlaps = 0
		ov_counter = 0
		for ov_index in ov_index_list:
			if (ov_counter%1000 == 0 or ov_counter == len(ov_index_list)-1):
				print_progress(ov_counter, len(ov_index_list)-1)
			ov_counter += 1
			if ov_index in self.overlaps:
				source_id = self.overlaps[ov_index].contig_sequence_1
				target_id = self.overlaps[ov_index].contig_sequence_2

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
		
		
		if not self.sequences[source_id].sequence == get_inverse_sequence(self.sequences[target_id].sequence, self.alphabet):
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
		# For every pair of sequence and its inverse, this method removes one if both are not in the same component of the graph
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
		
	def remove_tips(self, verbose=False):
		# removes tips (single-sequence-dead-ends) from the graph
		print ("Removing tips ...")
		num_of_removed_tips = -1
		while (num_of_removed_tips != 0):
			num_of_removed_tips = 0
			for seq in self.sequences:
				if seq.is_relevant:
					if verbose:
						print ("Consider sequence:")
						seq.print_data()
					# check if sequence is a tip and if sequence is shorter than 2k:
					if (len(seq.overlaps_out) + len(seq.overlaps_in) == 1) and (len(seq.sequence) < 2*self.k_value):
						if verbose:
							print ("Sequence is a tip, remove this sequence.")
						self.delete_sequence(seq.id, verbose)
						num_of_removed_tips += 1
			self.contract_unique_overlaps()

	def remove_insignificant_sequences(self, minimal_weight=2, verbose=False):
		# removes all sequences with weight less than minimal_weight
		for seq in self.sequences:
			if seq.get_total_weight() < minimal_weight:
				self.delete_sequence(seq.id, verbose)
	
	def remove_single_sequence_components(self, verbose=False):
		# removes all components that consist only of a single sequence,
		# i.e. all sequences without adjacencies
		# only if there is at least one larger component
		if len(self.overlaps) > 0:
			for seq_id in range(len(self.sequences)):
				if self.sequences[seq_id].is_relevant and len(self.sequences[seq_id].overlaps_out) == 0 and len(self.sequences[seq_id].overlaps_in) == 0:
					self.sequences[seq_id].is_relevant = False

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
				data += "k_"+str(seq.id)+","+seq.sequence+","+str(seq.get_total_weight())+","+str(seq.label)+","+str(source_reads)+"\n"
				#for r in source_reads:
				#	data += str(r)+" "
				#data += "\n"
		outputfile = file(filename, 'w')
		outputfile.write(headline)
		outputfile.write(data)
		
	def write_sequences_to_file(self, filename, addweights=False):
		# writes only the sequence-string of relevant sequences to a file
		data = ""
		for seq in self.sequences:
			if seq.is_relevant:
				data += seq.sequence
				if addweights:
					data += ","+str(seq.get_total_weight())
				data += "\n"
		outputfile = file(filename, 'w')
		outputfile.write(data)
		
	def load_from_asqg(self, filename="asqg_file", verbose=False):
		# reconstructs a graph from a asqg-file
		print ("Loading from asqg-file ...")
		self.is_unified = True
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
						
	def construct_assembly_ordering_labels(self, start_sequence = 0, do_second_iteration=True, verbose=False):
		# constructs a fuzzy partial ordering of all relevant sequences:
		# algorithm assumes that graph
		# 	is not empty and
		#	has only one component and
		# 	has no cycles, i.e. implies a partial order
		print "Construct assembly ordering labels..."
		
		if not (start_sequence == 0 or self.sequences[start_sequence].is_relevant):
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
				seq.label = False
		# init for start sequence
		self.min_label = 0
		self.max_label = 0
		self.min_sequence = start_sequence
		self.max_sequence = start_sequence
		self.sequences[start_sequence].label = 0
		
		if verbose:
			print ("Start with sequence " + str(start_sequence))
				
		queue = [[self.sequences[start_sequence].id, 0]]
		while (len(queue) > 0):
			current_data = queue[0]
			queue.pop(0)
			current_node_id = current_data[0]
			for seq_id in self.sequences[current_node_id].overlaps_out:
				if self.sequences[seq_id].label == False:
					start_label = current_data[1] + self.sequences[current_node_id].get_length()
					if start_label > self.max_label:
						self.max_label = start_label
						self.max_sequence = seq_id
					self.sequences[seq_id].label = start_label
					queue.append([seq_id, start_label])
			for seq_id in self.sequences[current_node_id].overlaps_in:
				if self.sequences[seq_id].label == False:
					start_label = current_data[1] - self.sequences[seq_id].get_length()
					if start_label < self.min_label:
						self.min_label = start_label
						self.min_sequence = seq_id
					self.sequences[seq_id].label = start_label
					queue.append([seq_id, start_label])
		
		if verbose:
			for seq in self.sequences:
				if seq.is_relevant:
					print str(seq.id) + ": " + str(seq.label)
		
		if do_second_iteration:
			next_start = 0
			if abs(self.max_label) > abs(self.min_label):
				next_start = self.max_sequence
			else:
				next_start = self.min_sequence
			if verbose:
				print ("max sequence of first iteration: "+str(self.max_sequence) + " with label "+str(self.max_label))
				print ("min sequence of first iteration: "+str(self.min_sequence) + " with label "+str(self.min_label))
				print ("Start a second iteration of labelling, starting from sequence "+str(next_start))
			self.construct_assembly_ordering_labels(start_sequence=next_start, do_second_iteration=False, verbose=verbose)
					
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
		
	def reduce_to_single_path_max_weight(self, start_sequence=False, restrict_to_component=False, verbose=False):
		# greedy algo that traverses through graph by choosing following nodes with max weight, deletes everythin else
		# with backtracking: if a dead end is reached, the algorithm returns to the previous node until a second-best path is found.
		# i.e. depth-first search, starting from node with smallest label to node with highest label where in each step the node with largest weight is chosen.
		# method assumes that graph has only one component
		# and sequences have position- and weight-labels
		
		print ("DFS for path with local maximum weight")
		
		'''
		components = self.get_components()
		if verbose:
			print ("Number of components: "+str(len(components)))
		if len(components) > 1:
			print ("Too many components.")
			self.reduce_to_single_largest_component()
		'''
			
		# node-labels of the dfs:
		#	0: not yet visited
		#	1: node visited, is on current path
		#	2: node visited and returned because it leads to a dead end
		dfs_labels = [0 for s in self.sequences]
		path = []
		
		if not start_sequence:
			current_seq_id = self.min_sequence
		else:
			current_seq_id = start_sequence
		current_label = self.sequences[current_seq_id].label
		path_finding_failed = False
		while (not path_finding_failed) and current_label < self.max_label-2*self.k_value:
			# go as far as possible, and if no successor, backtrack if not at least 95% of longest path is covered
			dfs_labels[current_seq_id] = 1
			next_sequences = [target_id for target_id in self.sequences[current_seq_id].overlaps_out if dfs_labels[target_id] == 0]
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
					current_label  = self.sequences[current_seq_id].label
					path.append(current_seq_id)
				if verbose:
					print ("Next sequence is sequence "+str(max_seq_id)+" with label " + str(current_label) + " and weight "+str(max_seq_weight))
			else:
				# backtrack:
				if len(path) > 1:
					dfs_labels[current_seq_id] = 2
					if verbose:
						print ("Dead end! Remove sequence "+str(path[-1])+" from path")
					path.pop(-1)
					current_seq_id = path[-1]
					current_label  = self.sequences[current_seq_id].label
					if verbose:
						print ("Next sequence is sequence "+str(current_seq_id)+" with label " + str(current_label) + " and weight "+str(self.sequences[current_seq_id].get_total_weight()))
				else:
					print ("Error! No path found")
					path_finding_failed = True
			#print ("current path: " +str(path))
		
		# delete unused sequences:
		for seq in self.sequences:
			if seq.is_relevant and not dfs_labels[seq.id] == 1 and (not restrict_to_component or seq.id in restrict_to_component):
				self.delete_sequence(seq.id)
		print ("Reduction finished")
		
	def reduce_every_component_to_single_path_max_weight(self, verbose=False):
		print ("Reducing every component to a single path with maximum local weight ...")
		components = self.get_components()
		if verbose:
			print ("Number of components: "+str(len(components)))
		for c in components:
			if verbose:
				print ("Size of this component: "+str(len(c)))
			self.construct_assembly_ordering_labels(start_sequence=c[0])
			start_sequence = -1
			min_label = False
			for seq_id in c:
				if not min_label or self.sequences[seq_id].label < min_label:
					start_sequence = seq_id
					min_label = self.sequences[seq_id].label
			if verbose:
				print ("start_sequence for reduction: "+str(start_sequence)+" with label: "+str(min_label))
			self.reduce_to_single_path_max_weight(start_sequence = start_sequence, restrict_to_component=c)
		
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
			start_seq_id = self.min_sequence
			'''
			min_label = False
			max_weight = 0
			if verbose:
				print ("Searching for start sequence...")
			for seq_id in c:
				if verbose:
					print ("Current min_label is: "+str(min_label))
					print ("Current weight of min_sequence is: "+str(max_weight))
					print ("Consider sequence "+str(seq_id)+" with label "+str(self.sequences[seq_id].label)+ " and weight "+str(self.sequences[seq_id].max_weight))
				#if (start_seq_id < 0) or (self.sequences[seq_id].label < min_label-self.k_value) or (self.sequences[seq_id].label < min_label+self.k_value and self.sequences[seq_id].max_weight > max_weight):
				#if (start_seq_id < 0) or (self.sequences[seq_id].label < self.min_label+self.k_value and self.sequences[seq_id].max_weight > max_weight):
				if (self.sequences[seq_id].label < self.min_label+self.k_value and self.sequences[seq_id].max_weight > max_weight):
					if verbose:
						print ("Set sequence "+str(seq_id)+" as new minimal sequence")
					start_seq_id = seq_id
					max_weight = self.sequences[seq_id].max_weight
			'''
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
		
	def reduce_to_single_largest_component(self):
		# deletes all sequences that are not part of the largest (by number of sequences) component
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
		return self.max_label - self.min_label
		
	def construct_laplacian(self, list_of_nodes):
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
		
		part_a = []
		part_b = []
		part_c = []
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
					
		return secmin_eigenvalue, part_a, part_b, part_c
					
	def compute_mincut(self, component=-1, divide_clusters=True, verbose=False):
		print ("Do spectral clustering of the nodes into two clusters ...")
		self.remove_irrelevant_overlaps()
		
		if component < 0:
			component = [seq.id for seq in self.sequences]
		laplacian, seq_id_to_index, index_to_seq_id = self.construct_laplacian(component)
		
		secmin_eigenvalue, part_a, part_b, part_c = self.construct_spectral_clusters(laplacian, index_to_seq_id, verbose)
		
		if len(part_a) > sqrt(len(componentlen)) and len(part_b) > sqrt(len(componentlen)):
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
		if divide_clusters and (secmin_eigenvalue*10 < secmin_eigenvalue_a or secmin_eigenvalue*10 < secmin_eigenvalue_b):
			if verbose:
				print ("Cut graph according to partitions.")
			self.cut_graph_into_partitions([part_a, part_b] ,seq_id_to_index)
			return True, part_a, part_b
		else:
			if verbose:
				print ("Do not cut graph.")
			return False, part_a, part_b
			
	def partition_graph_into_components_of_clusters(self, verbose=False):
		components_to_potentially_cut = self.get_components()
		while (len(components_to_potentially_cut) > 0):
			if verbose:
				print ("current number of components to consider: "+str(len(components_to_potentially_cut)))
			res, part_a, part_b = self.compute_mincut(components_to_potentially_cut[0], verbose=verbose)
			components_to_potentially_cut.pop(0)
			if res:
				components_to_potentially_cut.append(part_a)
				components_to_potentially_cut.append(part_b)
		
	def cut_graph_into_partitions(self, partitions, seq_id_to_index):
		print ("Delete all overlaps between different parts ...")
		id_to_part={}
		for i in range(len(partitions)):
			for id in partitions[i]:
				id_to_part[id]=i
		
		overlaps_to_delete = []
		for ov in self.overlaps:
			if not id_to_part[seq_id_to_index[self.overlaps[ov].contig_sequence_1]] == id_to_part[seq_id_to_index[self.overlaps[ov].contig_sequence_2]]:
				overlaps_to_delete.append(ov)
		
		for ov_id in set(overlaps_to_delete):
			self.delete_overlap(ov_id)
	
	def remove_irrelevant_overlaps(self):
		ov_to_del = []
		for ov in self.overlaps:
			if not self.sequences[self.overlaps[ov].contig_sequence_1].is_relevant or not self.sequences[self.overlaps[ov].contig_sequence_2].is_relevant:
				ov_to_del.append(ov)
		for ov in ov_to_del:
			self.overlaps.pop(ov)
		
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

def compute_insert_distance(sequence_1, sequence_2, maxdist = -1):
	# algorithm may not work properly for arbitrary large insert-distances,
	# but is correct if local insert-distace is <= 1
	# returns insert-distance if insert_distance <= maxdist,
	# otherwise returns maxdist + x for a x >= 1
	if not len(sequence_1) == len(sequence_2):
		return -1
	index_1 = 0
	index_2 = 0
	
	insert_distance = 0
	n = len(sequence_1)
	if maxdist < 0:
		maxdist = n
	while (insert_distance < (maxdist+1) and index_1 < n and index_2 < n):
		t1 = 0
		t2 = 0
		
		while (index_1 + t1 < n and (not sequence_1[index_1+t1] == sequence_2[index_2])):
			t1 += 1
		while (index_2 + t2 < n and (not sequence_1[index_1] == sequence_2[index_2+t2])):
			t2 += 1
			
		if t1 <= t2:
			index_1 += t1
			insert_distance += t1
		elif t2 < t1:
			index_2 += t2
			insert_distance += t2
			
		index_1 += 1
		index_2 += 1
	return insert_distance
