#!usr/bin/python

import random
import re

import gc

class Read:
	def __init__(self, r_id, sequence):
		self.id = r_id
		self.sequence = sequence
		self.length = len(sequence)
		self.kmers = []
		self.paired_end_partner = -1
		
	def add_paired_end_partner(self, partner_id):
		self.paired_end_partner = partner_id
		
	def add_kmer(self, kmer_id):
		self.kmers.append(kmer_id)

class Kmer:
	def __init__(self, k_id, inv_id , sequence, evidence_reads):
		self.id = k_id
		self.sequence = sequence
		self.id_of_inverse_kmer = inv_id
		self.evidence_reads = evidence_reads
		
	def add_evidence(self, evidence_read_id):
		self.evidence_reads.append(evidence_read_id)
		
	def get_evidence_weight(self):
		return len(self.evidence_reads)

class ContigSequence:
	def __init__(self, seq_id, inv_id, sequence, kmers, weight = 1, is_relevant = True):
		self.id = seq_id
		self.id_of_inverse_seq = inv_id
		self.sequence = sequence
		self.kmers = kmers
		# the maximal read eavidence this sequence has for any subsequence
		self.max_weight = weight
		# overlaps (aka edges) are stored in dictionaries
		# self.overlap[other_sequence_id] = overlap_id
		self.overlaps_out = {}
		self.overlaps_in = {}
		self.is_relevant = is_relevant
		self.label = False

	def get_length(self):
		return len(self.sequence)
		
	def check_if_overlap_exists(self, other_sequence_id):
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
	def __init__(self, ov_id, length, seq_1, seq_2, evidence_reads):
		self.id = ov_id
		self.length = length
		self.contig_sequence_1 = seq_1
		self.contig_sequence_2 = seq_2
		self.evidence_reads = evidence_reads
		self.is_relevant = True
	
	def add_evidence(self, read_id):
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
	def __init__(self, reads=0, k=0, verbose=False, alphabet={"A":"T", "C":"G", "G":"C", "T":"A"}):
		print ("Creating empty graph ...")
		
		self.k_value = k
		self.reads = []
		self.kmers = []
		self.kmer_dict = {}
		self.sequences = []
		self.overlaps = {}
		self.is_unified = False
		
		self.alphabet = alphabet
		
		# min and max label of sequences
		self.max_label = 0
		self.min_label = 0
		
		if not reads == 0:
			self.init_graph_database(reads, verbose=verbose)

	def init_graph_database(self, reads, verbose=False):
		if verbose:
			print ("Construct read database")
		read_id = 0
		for r in reads:
			for read in r:
				self.reads.append(Read(read_id, read))
				read_id += 1
			
		# construct k-mer database:
		self.get_kmerdata_from_reads(verbose)
		# construct sequences from kmers:
		print ("Construct Sequences from k-mers ...")
		for kmer in self.kmers:
			if verbose:
				print ("now consider kmer with id " + str(kmer.id) + ": " + kmer.sequence)
			seq_id = kmer.id
			seq_inv_id = kmer.id_of_inverse_kmer
			weight = len(kmer.evidence_reads)
			
			self.sequences.append(ContigSequence(seq_id, seq_inv_id, kmer.sequence, [kmer.id], weight))

		print ("Construct overlaps ...")			
		# construct overlaps between adjacent sequences with read-evidence:
		for read in self.reads:
			for kmer_index in range(len(read.kmers)-1):
				source_kmer_id = read.kmers[kmer_index]
				target_kmer_id = read.kmers[kmer_index+1]
				self.increment_overlap(source_kmer_id, target_kmer_id, read.id, verbose=False)

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
				for kmer_id in s.kmers:
					kmer_string = self.kmers[kmer_id].sequence + "\t"
					for read_id in self.kmers[kmer_id].evidence_reads:
						kmer_string += str(read_id)+" "
					print kmer_string
				print("")

	def get_kmerdata_from_reads(self, verbose = False):
		print ("Get kmer-data from reads ...")

		read_index = 0
		kmer_counter = 0
		for read_index in range(len(self.reads)):
			if read_index%100 == 0 and not verbose:
				print ("Current read: "+str(read_index)+"/"+str(len(self.reads)))
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
					if verbose:
						print ("Kmer already exists")
						print ("Add read ("+str(read_index)+") evidence to kmer "+str(this_kmer_id))
				if not kmer_already_existing:
					if verbose:
						print ("Kmer does not exist in database. Add new kmer ...")
					
					# add kmer:
					self.kmers.append(Kmer(kmer_counter, kmer_counter+1, new_kmer_sequence, [read_index]))
					self.kmer_dict[new_kmer_sequence] = kmer_counter
					this_kmer_id = kmer_counter
					kmer_counter += 1
					# add inverse kmer:
					self.kmers.append(Kmer(kmer_counter, kmer_counter-1, get_inverse_sequence(new_kmer_sequence, self.alphabet), [read_index]))
					self.kmer_dict[get_inverse_sequence(new_kmer_sequence, self.alphabet)] = kmer_counter
					kmer_counter += 1
					
				if verbose:
					print ("Add kmer "+self.kmers[this_kmer_id].sequence+"("+str(this_kmer_id)+") to read "+str(read_index))
				# add kmer to read:
				self.reads[read_index].add_kmer(this_kmer_id)
				kmer_start += 1
			read_index += 1
			
	def increment_overlap(self, source_seq_id, target_seq_id, read_evidence, consider_inverse = True, verbose = False):
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

	def contract_unique_overlaps(self, verbose = False):
		print ("Contract overlaps ...")
		
		ov_index_list = [ov_id for ov_id in self.overlaps]
		for ov_index in ov_index_list:
			if (ov_index%1000 == 0):
				print (str(ov_index)+"/"+str(len(self.overlaps)))
			if (ov_index%100000 == 0):
				# run garbage collector:
				gc.collect()
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
		            
					# contract reverse overlap if not sequence is its own inverse:
					if not self.sequences[source_id].sequence == get_inverse_sequence(self.sequences[source_id].sequence, self.alphabet):
						source_rev_id = self.sequences[target_id].id_of_inverse_seq
						target_rev_id = self.sequences[source_id].id_of_inverse_seq
						
						if not self.is_unified:
							rev_ov_id = self.sequences[source_rev_id].overlaps_out[target_rev_id]
							self.contract_overlap(rev_ov_id, verbose)
						
						self.sequences[source_id].id_of_inverse_seq = source_rev_id
						self.sequences[source_rev_id].id_of_inverse_seq = source_id
	
	def delete_overlap(self, overlap_id, verbose=False):
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
		if verbose:
			print ("Removing Sequence "+str(sequence_id))
		adj_seq_out = self.sequences[sequence_id].overlaps_out.keys()
		for adj_seq in adj_seq_out:
			self.delete_overlap(self.sequences[sequence_id].overlaps_out[adj_seq], verbose)
		adj_seq_in = self.sequences[sequence_id].overlaps_in.keys()
		for adj_seq in adj_seq_in:
			self.delete_overlap(self.sequences[sequence_id].overlaps_in[adj_seq], verbose)
		self.sequences[sequence_id].is_relevant = False
			
	def contract_overlap(self, overlap_id, verbose=False):
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
		if self.sequences[target_id].max_weight > self.sequences[source_id].max_weight:
			self.sequences[source_id].max_weight = self.sequences[target_id].max_weight
		
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
		self.sequences[target_id].is_relevant = False
		# Don't use delete_overlap, because incident sequences have been handled manually:
		self.overlaps.pop(overlap_id)
		
	def remove_tips(self, verbose=False):
		print ("Removing tips ...")
		for seq in self.sequences:
			if seq.is_relevant:
				if verbose:
					print ("Consider sequence:")
					seq.print_data()
				# check if sequence is a tip and if sequence is shorter than 2k:
				if (len(seq.overlaps_out) + len(seq.overlaps_in) == 1) and (len(seq.sequence) < 2*self.k_value):
					if verbose:
						print ("Sequence is a tip, needs to be removed")
					self.delete_sequence(seq.id, verbose)

	def remove_insignificant_sequences(self, minimal_weight=2, verbose=False):
		# removes all sequences with weight less than minimal_weight
		for seq in self.sequences:
			if seq.max_weight < minimal_weight:
				self.delete_sequence(seq.id, verbose)
	
	def remove_parallel_sequences(self, verbose=False):
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
		headline = "Node_Label, Sequence, maxweight\n"
		data = ""
		for seq in self.sequences:
			if seq.is_relevant:
				data += "k_"+str(seq.id)+","+seq.sequence+","+str(seq.max_weight)+"\n"
		outputfile = file(filename, 'w')
		outputfile.write(headline)
		outputfile.write(data)
		
	def write_sequences_to_file(self, filename):
		data = ""
		for seq in self.sequences:
			if seq.is_relevant:
				data += seq.sequence + "\n"
		outputfile = file(filename, 'w')
		outputfile.write(data)
		
	def load_from_asqg(self, filename="asqg_file", verbose=False):
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
						
	def construct_assembly_ordering_labels(self):
		# algorithm assumes that graph
		# 	is not empty and
		#	has only one component and
		# 	has no cycles, i.e. implies a partial order
		print "Construct assembly ordering labels..."
		
		for start_seq_id in range(len(self.sequences)):
			if not self.sequences[start_seq_id].label:
				queue = [[self.sequences[0].id, 0]]
				while (len(queue) > 0):
					current_data = queue[0]
					queue.pop(0)
					current_node_id = current_data[0]
					#print self.sequences[current_node_id].label
					#print current_node_id
					#print self.sequences[current_node_id].overlaps_out
					#print self.sequences[current_node_id].overlaps_in
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
					
	def get_partition_of_sequences(self, number_of_parts, verbose=False):
		# returns a partition of all read-ids based on intervals of labels
		sorted_nodes = sorted(self.sequences, key=lambda x: x.label)
		label_div = self.max_label-self.min_label
		part_size = label_div/(number_of_parts+1)
		
		if verbose:
			print label_div
			print part_size
		
		#parts = []
		parts_seq = []
		for i in range(number_of_parts):
			current_start = self.min_label+i*(part_size)
			if i == number_of_parts-1:
				current_end = self.max_label
			else:
				current_end = self.min_label+(i+2)*(part_size)
			this_part_sequences = [seq for seq in sorted_nodes if seq.label >= current_start and seq.label <= current_end]
			parts_seq.append(this_part_sequences)
		return parts_seq
		'''
			this_part_kmers = []
			for seq in this_part_sequences:
				this_part_kmers += seq.kmers
			this_part_reads = []
			for kmer_id in this_part_kmers:
				this_part_reads += self.kmers[kmer_id].evidence_reads
			this_part_reads = list(set(this_part_reads))#[self.kmers[kmer_id].evidence_reads for kmer_id in this_part_kmers]))
			parts.append(this_part_reads)
			if verbose:
				print str(i)+": "+str(current_start)+" - "+str(current_end)+" : "+str(len(parts[-1]))+" sequences"
		return parts
		'''
	
	def get_read_of_sequences(self, sequences, verbose=False):
		kmers = []
		for seq in sequences:
			kmers += seq.kmers
		reads = []
		for kmer_id in kmers:
			reads += self.kmers[kmer_id].evidence_reads
		reads = list(set(reads))
		return reads
		
	'''
	def merge_with_other_graph(self, other_graph, verbose=False):
		# store sequences and overlaps that were added in dictionaries that translate other_id -> this_id
		merged_sequencs = {}
		added_sequences = {}
		added_overlaps = {}
		for seq_other in other_graph.sequences:
			sequence_already_exists = False
			for seq_own in self.sequences:
				if seq_own.sequence == seq_other.sequence:
					sequence_already_exists = True
					if len(seq_other.overlaps_out) > 0 or len(seq_other.overlaps_in) > 0:
						merged_sequences[seq_other.id] = seq_own.id
						for ov_other_id in seq_other.overlaps_out:
							if other_graph.overlaps[ov_other_id].contig_sequence_2 in merged_sequences:
								new_ov_target_id = merged_sequencs[other_graph.overlaps[ov_other_id].contig_sequence_2]
							elif other_graph.overlaps[ov_other_id].contig_sequence_2 in added_sequences:
								new_ov_target_id = added_sequences[other_graph.overlaps[ov_other_id].contig_sequence_2]
							else:
								new_ov_target_id = len(self.sequences)
								new_ov_target_id_inv = new_ov_target_id+1
								self.sequences.append(ContigSequence(new_ov_target_id, new_ov_target_id_inv, other_graph.sequences[other_graph.overlaps[ov_other_id].contig_sequence_2], [], 1))
								self.sequences.append(ContigSequence(new_ov_target_id_inv, new_ov_target_id, get_inverse_sequence(other_graph.sequences[other_graph.overlaps[ov_other_id].contig_sequence_2], self.alphabet), [], 1))
								added_sequences[other_graph.overlaps[ov_other_id].contig_sequence_2] = new_ov_target_id
								added_sequences[other_graph.sequences[other_graph.overlaps[ov_other_id].contig_sequence_2].id_of_inverse_seq] = new_ov_target_id_inv
							self.increase_overlap(seq_own, new_ov_target_id, -1)
							
						for ov_other_id in seq_other.overlaps_in:
							if other_graph.overlaps[ov_other_id].contig_sequence_1 in merged_sequences:
								new_ov_target_id = merged_sequencs[other_graph.overlaps[ov_other_id].contig_sequence_1]
							elif other_graph.overlaps[ov_other_id].contig_sequence_1 in added_sequences:
								new_ov_target_id = added_sequences[other_graph.overlaps[ov_other_id].contig_sequence_1]
							else:
								new_ov_target_id = len(self.sequences)
								new_ov_target_id_inv = new_ov_target_id+1
								self.sequences.append(ContigSequence(new_ov_target_id, new_ov_target_id_inv, other_graph.sequences[other_graph.overlaps[ov_other_id].contig_sequence_1], [], 1))
								self.sequences.append(ContigSequence(new_ov_target_id_inv, new_ov_target_id, get_inverse_sequence(other_graph.sequences[other_graph.overlaps[ov_other_id].contig_sequence_1], self.alphabet), [], 1))
								added_sequences[other_graph.overlaps[ov_other_id].contig_sequence_1] = new_ov_target_id
								added_sequences[other_graph.sequences[other_graph.overlaps[ov_other_id].contig_sequence_1].id_of_inverse_seq] = new_ov_target_id_inv	
							self.increase_overlap(new_ov_target_id, seq_own, -1)
			if not sequence_already_exists:
				# TO DO: add sequence and overlaps
	'''
		
	def reduce_to_single_path_max_weight(self, verbose=False):
		# method assumes that graph has only one component and no cycles
		# and sequences have weight-labels
		start_seq_id = -1
		min_label = False
		for seq in self.sequences:
			if seq.is_relevant and (min_label == False or seq.label < min_label):
				start_seq = seq.id
				min_label = seq.label
			
		current_seq_id = start_seq_id
		while len(self.sequences[current_seq_id].overlaps_out) > 0:
			next_sequences = []
			for ov_id in self.sequences[current_seq_id].overlaps_out:
				next_sequences.append(self.overlaps[ov_id].contig_sequence_2)
			max_seq_id = -1
			max_seq_weight = -1
			for seq_id in next_sequences:
				if self.sequences[seq_id].is_relevant and self.sequences[seq_id].max_weight > max_seq_weight:
					max_seq_weight = self.sequences[seq_id].max_weight
					max_seq_id = seq_id
			for seq_id in next_sequences:
				if seq_id == max_seq_id:
					current_seq_id = seq_id
				else:
					self.delete_sequence(seq_id)
		
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
