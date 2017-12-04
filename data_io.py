#!usr/bin/python

import random
import re

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
	
def get_random_mutation(nucleotide, alphabet=['A', 'C', 'G', 'T']):
	possible_mutations = [c for c in alphabet if not c == nucleotide]
	return random.choice(possible_mutations)

def genereate_dna(length=100, alphabet=["A","C","G","T"]):
	print ("Generating dna sequence ...")
	dna = []
	for i in range(length):
		dna.append(random.choice(alphabet))
	return dna
	
def write_dna_to_file(filename, dna):
	with open(filename, 'w') as outf:
		outf.write(''.join(dna) + '\n')

def genereate_reads(dna, coverage=50, avg_read_length=50, mutation_pct=0.1, mutation_alphabet=["A","C","G","T"], remove_pct=5, variation=0, paired_end=False, both_directions=False, verbose=False):
	print ("Generating reads ...")
	#print ("mutation_pct = "+str(mutation_pct) + ", remove_pct = "+str(remove_pct))
	dna_length = len(dna)
	reads = []
	alignment = []
	for i in range(coverage):
		if verbose:
			print ("iteration: " + str(i))
		cur_index = 0
		while cur_index < dna_length:
			#next_index = min(cur_index+avg_read_length+random.choice(range(-10,10)), dna_length)
			if variation < 0:
				variation = 1+avg_read_length/5
			read_lengt_variation = int(random.gauss(mu=0, sigma=variation))
			if verbose:
				print "read_length_variation: "+str(read_lengt_variation)
			next_index = min(cur_index+avg_read_length+read_lengt_variation, dna_length)
			new_read = dna[cur_index:next_index]
			for s_i in range(len(new_read)):
				if random.choice(range(10000)) < mutation_pct*100:
					if verbose:
						print ("mutation at read "+str(s_i))
					new_read[s_i] = random.choice(mutation_alphabet)
			if random.choice(range(100)) > remove_pct:
				if both_directions and random.choice(range(100)) > 50:
					reads.append([get_inverse_sequence(new_read)])
				else:
					reads.append([''.join(new_read)])
				alignment.append([cur_index, next_index-1])
			cur_index = next_index

	### TO DO: generate paired-end information

	return reads, alignment

def get_reads_from_file(filename="samplereads.txt"):
	with open(filename) as fh:
		reads = fh.readlines()
	for read_index in range(len(reads)):
		reads[read_index] = [re.sub(r"\s", '', reads[read_index])]
	return reads
	
def get_reads_from_fastq_file_by_length(read_ids, filename="fastqreads.fq"):
	# if num_of_reads > 0, only the specified number of reads will be read from the file.
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
	# if num_of_reads > 0, only the specified number of reads will be read from the file.
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

def get_readlengths(filename):
	print ("Get readlengths from file "+filename)
	reads = get_reads_from_fastq_file(filename)
	n = len(reads)
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

def get_read_partition_by_readlength(filename, number_of_parts=-1, size_of_parts=-1, verbose=False):
	readlengths = get_readlengths(filename)
	readlengths_sorted = sorted(readlengths, key=lambda x: x[0])
	
	if number_of_parts*size_of_parts>0:
		print ("Error! Exactly one of the parameters number_of_parts and size_of_parts has to be specified!")
		return [i[1] for i in readlengths]

	parts = []
	if number_of_parts > 0:
		size_of_parts = len(readlengths)/number_of_parts
	elif size_of_parts > 0:
		number_of_parts = len(readlengths)/size_of_parts

	if verbose:
		print ("size_of_parts = "+str(size_of_parts))
		print ("number_of_parts = "+str(number_of_parts))

	for i in range(number_of_parts):
		if i == number_of_parts-1:
			partend = len(readlengths)-1
		else:
			partend = size_of_parts*(i+1)
		parts.append(readlengths_sorted[size_of_parts*i:partend])

	return parts
