#!/usr/bin/env python
# -*- coding: utf-8 -*-

from random import randint, gauss

import data_io

def read_genomes(input_filename="Data/bvdv_selected_filtered.linsi.aln", output_filename="Data/bvdv_genomes.txt"):
	with open(input_filename) as fh:
		lines = fh.readlines()
	
	genomes = []
	sequences = {}
	for line in lines:
		if (line != ''):
			if (line[0] == '>'):
				gname = line.strip()[1:]
				genomes.append(gname)
				sequences[gname] = ''
			else:
				sequences[gname] += line.strip()
	
	with open(output_filename, 'w') as outf:
		for genome_name in genomes:
			outf.write(sequences[genome_name]+"\n")

def samplereads(input_filename="Data/bvdv_genomes.txt",
				output_filename="samplereads.txt",
				read_length=50,
				length_stddev=0,
				set_of_viruses=[],
				number_of_reads=[5000,5000],
				avg_error_percentage=0,
				inverted_reads=False):
		
	alphabet = ['a','c','g','t']

	sequences = []
	with open(input_filename) as inputfile:
		lines = inputfile.readlines()
	
	for line in lines:
		sequences.append(line.strip())
		
	closegaps = True
	numreads = sum(number_of_reads)
	if len(set_of_viruses) == 0:
		set_of_viruses = range(len(number_of_reads))
	reads = []
	
	# construct reads:
	for curgenomeind_id in range(len(set_of_viruses)):
		
		curgenomeind = set_of_viruses[curgenomeind_id]
		curnumreads = number_of_reads[curgenomeind_id]
		
		for n in range(curnumreads):
			
			readlen = read_length + int(gauss(mu=0, sigma=length_stddev))
			ind = randint(0, len(sequences[curgenomeind])-readlen)
			sampleread = sequences[curgenomeind][ind:ind+readlen]
			while ('-' in sampleread):
				
				if closegaps:
					nreplace = sampleread.count('-')
					if (ind + nreplace > len(sequences[curgenomeind])):
						ind = randint(0, numreads-readlen)
						sampleread = sequences[curgenomeind][ind:ind+readlen]
					else:
						added = sequences[curgenomeind][ind+readlen:ind+readlen+nreplace]
						if ('-' in added):
							ind = randint(0, numreads-readlen)
							sampleread = sequences[curgenomeind][ind:ind+readlen]
						else:
							sampleread = sampleread.replace('-','') + added
					
				else:
					ind = randint(0, numreads-readlen)
					sampleread = sequences[curgenomeind][ind:ind+readlen]
			
			for i in range(len(sampleread)):
				if randint(0,10000) < int(100*avg_error_percentage):
					sampleread = sampleread[:i]+data_io.get_random_mutation(sampleread[i], alphabet)+sampleread[(i+1):]
					
			if (inverted_reads and randint(0, 100)>50):
				reads.append(data_io.get_inverse_sequence(sampleread.upper()))
			else:
				reads.append(sampleread.upper())
	
	with open(output_filename, 'w') as outf:
		for read in reads:
			outf.write(read + '\n')
		
#samplereads()