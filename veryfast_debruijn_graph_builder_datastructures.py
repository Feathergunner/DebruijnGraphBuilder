#!usr/bin/python
# -*- coding: utf-8 -*-

def print_progress(part, total):
	print ("Progress: "+str("%.2f" % ((float(part)/(float(total)/100)))) + "%")

class Read:
	def __init__(self, read_id, sequence):
		self.id = read_id
		self.sequence = sequence
		self.length = len(sequence)
		self.kmers = []
		
	def add_kmer(self, kmer_id):
		self.kmers.append(kmer_id)

class Kmer:
	def __init__(self, kmer_id, inverse_id, evidence_reads, reference_read, read_start_pos, seq_is_inversed):
		self.id = kmer_id
		self.id_of_inverse_kmer = inverse_id
		self.evidence_reads = evidence_reads
		self.source = [reference_read, read_start_pos, seq_is_inversed]
		
	def add_evidence(self, evidence_read_id):
		self.evidence_reads.append(evidence_read_id)
		
	def get_evidence_weight(self):
		return len(self.evidence_reads)

class ContigSequence:
	# Nodes in the Debruijn-Graph
	def __init__(self, seq_id, inv_id, reference_read, read_start, seq_is_inversed, length, kmers, weight = 1):
		self.id = seq_id
		self.id_of_inverse_seq = inv_id

		# a read_id the start position and length of this sequence for sequence reconstruction
		self.read_sources=[[reference_read, read_start, seq_is_inversed, length]]

		self.kmers = kmers
		# the maximal read eavidence this sequence has for any subsequence
		self.max_weight = weight
		# overlaps (i.e. edges) are stored in dictionaries
		# self.overlap[other_sequence_id] = overlap_id
		self.overlaps_out = {}
		self.overlaps_in = {}
		# used for estimation of position within assembly
		self.label = False

	def get_length(self):
		return sum([r[2] for r in self.read_sources])
		
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
		print (self.kmers)
		print ("Overlaps:")
		print (self.overlaps_out)
		print (self.overlaps_in)
		print ("id of inverse sequence: "+str(self.id_of_inverse_seq))
		print ("Label: "+str(self.label))
		print ("weiht "+str(self.max_weight))

class SequenceOverlap:
	# Edges in the Debruijn-Graph
	def __init__(self, ov_id, seq_1, seq_2):
		self.id = ov_id
		# the incident sequences of this overlap. sequence_1 is the source, sequence_2 is the target.
		self.contig_sequence_1 = seq_1
		self.contig_sequence_2 = seq_2

	def print_data(self):
		print ("Overlap "+str(self.id))
		print ("From seq. "+str(self.contig_sequence_1)+" to seq. "+str(self.contig_sequence_2))