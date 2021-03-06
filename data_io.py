#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import math

def write_genome_to_file(dna, filename):
	outputfile = file(filename, 'w')
	outputfile.write("".join(dna))
	
def write_reads_to_file(reads, filename):
	outputfile = file(filename, 'w')
	for r in reads:
		outputfile.write(r+"\n")

def get_genome_from_file(filename):
	# assumes that every line starting with an uppercase alphabetic character contains a section of a single genome, with all lines in order of the genome.
	with open(filename) as fh:
		lines = fh.readlines()
	dna = []
	for line in lines:
		if re.match(r'[A-Z]', line[0]):
			dna+=line.strip()
	return ''.join(dna)

def get_reads_from_file(filename="samplereads.txt"):
	with open(filename) as fh:
		reads = fh.readlines()
	for read_index in range(len(reads)):
		reads[read_index] = re.sub(r"\s", '', reads[read_index])
	return reads
	
def get_reads_from_fastq_file_by_partition(filename, read_ids=-1):
	# if read_ids < -1, will return all reads from file,
	# if read_ids is a nonempty list, only the reads specified by indices of read_ids will be returned
	status = 0
	reads = []
	n = 0
	with open(filename) as fh:
		readdata = fh.readlines()
		for line in readdata:
			if status == 1 and not line[0] == "@":
				if read_ids<0 or n in read_ids:
					reads.append(line.strip())
				status = 0
				n += 1
			if line[0] == "@":
				status = 1
	return reads

def get_reads_from_fastq_file(filename="fastqreads.fq", num_of_reads=-1, first_read=1):
	# if num_of_reads > 0, only the specified number of reads will be read from the file
	# always all reads before first_read are ignored
	status = 0
	reads = []
	n = 1
	with open(filename) as fh:
		readdata = fh.readlines()
		for line in readdata:
			if status == 1 and not line[0] == "@":
				if n >= first_read:
					reads.append(line.strip())
				status = 0
				n += 1
				if num_of_reads > 0 and n >= first_read+num_of_reads:
					break
			if line[0] == "@":
				status = 1
	return reads

def write_asqg_file(kmers, contig_seqs, edges, k, filename="asqg_file"):
	print ("Writing asqg-file ...")
	headline = "HT\t\n"
	asqg_vertices = ""
	vertex_count = 0
	for kmer in kmers:
		asqg_vertices += "VT\tk_"+str(kmer)+"\t"+contig_seqs[kmer]+"\n"
		vertex_count += 1
	asqg_edges = ""
	for source in edges:
		for target in edges[source]:
			asqg_edges += "ED\tk_"+str(source)+" k_"+str(target)+" 1 "+str(k-1)+" "+str(k)+" 0 "+str(k-2)+" "+str(k)+" 0 0\n"
	outputfile = file(filename, 'w')
	outputfile.write(headline)
	outputfile.write(asqg_vertices)
	outputfile.write(asqg_edges)
	
def write_sequences_to_fasta(sequences, filename="fasta.fasta"):
	# this method assumes that the parameter "sequences" is a list
	# if it is a single string, every character will be interpreted as a unique sequence!
	#print ("Writing fasta-file ...")
	fasta_string = ""
	seq_counter = 0
	for seq in sequences:
		seq_counter += 1
		fasta_string += ">SEQUENCE_"+str(seq_counter)+"\n"
		fasta_string += seq
		fasta_string += "\n\n"
		
	outputfile = file(filename, 'w')
	outputfile.write(fasta_string)

def print_alignment(dna, alignments):
	print ("Alignment of reads:")
	print ('  '+' |  '.join(dna))
	read_index = 0
	alignment_position = 0
	read_align_string = ""
	while read_index < len(alignments):
		while alignment_position < alignments[read_index][0]:
			read_align_string += "   . "
			alignment_position += 1
		while alignment_position <= alignments[read_index][1]:
			read_align_string += "{:4d} ".format(read_index)
			alignment_position += 1
		read_index += 1
		if alignment_position >= len(dna):
			print read_align_string
			read_align_string = ""
			alignment_position = 0

def get_readlengths(filename, verbose=False):
	if verbose:
		print ("Get readlengths from file "+filename)
	reads = get_reads_from_fastq_file(filename)
	n = len(reads)
	if verbose:
		print ("Number of reads: "+str(n))
	readlengths = [0]*n
	for i in range(n):
		readlengths[i] = [len(reads[i]),i]
	return readlengths

def get_read_subset_by_readlength(filename, minlength, maxlength):
	readlengths = get_readlengths(filename)
	read_ids = []
	for r in readlengths:
		if r[0] >= minlength and r[0] <= maxlength:
			read_ids.append(r[1])

	return read_ids

def get_read_partition_by_readlength(filename, number_of_parts=-1, size_of_parts=-1, max_readlength=-1, verbose=False):
	readlengths = get_readlengths(filename)
	readlengths_sorted = sorted(readlengths, key=lambda x: x[0])
	readlengths_sorted_reduced = [x for x in readlengths_sorted if (x[0] <= max_readlength or max_readlength == -1)]
	
	if number_of_parts*size_of_parts>0:
		print ("Error! Exactly one of the parameters number_of_parts and size_of_parts has to be specified!")
		return [i[1] for i in readlengths]

	parts = []
	if number_of_parts > 0:
		number_of_parts = min(number_of_parts, len(readlengths_sorted_reduced))
		size_of_parts = len(readlengths_sorted_reduced)/number_of_parts
	elif size_of_parts > 0:
		size_of_parts = min(size_of_parts, len(readlengths_sorted_reduced))
		number_of_parts = len(readlengths_sorted_reduced)/size_of_parts

	if verbose:
		print ("size_of_parts = "+str(size_of_parts))
		print ("number_of_parts = "+str(number_of_parts))

	for i in range(number_of_parts):
		if i == number_of_parts-1:
			partend = len(readlengths_sorted_reduced)-1
		else:
			partend = size_of_parts*(i+1)
		parts.append(readlengths_sorted_reduced[size_of_parts*i:partend])

	return parts

def get_read_partition_by_lengthdistribution(filename, number_of_parts=-1, size_of_parts=-1, verbose=False):
	readlengths = get_readlengths(filename)
	readlengths_sorted = sorted(readlengths, key=lambda x: x[0])
	n = len(readlengths_sorted)
	
	if number_of_parts*size_of_parts>0:
		print ("Error! Exactly one of the parameters number_of_parts and size_of_parts has to be specified!")
		return [i[1] for i in readlengths]
		
	tmp_parts = []
	if number_of_parts > 0:
		size_of_parts = int(math.ceil(n/float(number_of_parts)))
	elif size_of_parts > 0:
		number_of_parts = int(math.ceil(n/float(size_of_parts)))

	if verbose:
		print ("size_of_parts = "+str(size_of_parts))
		print ("number_of_parts = "+str(number_of_parts))
	
	parts = [[] for i in range(number_of_parts)]
	for i in range(size_of_parts):
		for j in range(number_of_parts):
			k = j + (i*number_of_parts)
			if k < n:
				parts[number_of_parts-j-1].append(readlengths_sorted[n-k-1])
	return parts
		