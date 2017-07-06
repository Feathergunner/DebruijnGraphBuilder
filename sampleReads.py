#!/usr/bin/env python
# -*- coding: utf-8 -*-

from random import randint, gauss

from data_io import get_inverse_sequence

def samplereads(input_filename="Data/bvdv_selected_filtered.linsi.aln",
				output_filename="samplereads.txt",
				read_length=50,
				length_stddev=0,
				inverted_reads=False):

	# read input genomes:
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
	
	closegaps = True
	numreads = 10000
	reads = []
	
	fractions = [0.5, 0.5, 0, 0, 0]
	
	# construct reads:s
	for curgenomeind in range(len(fractions)):
		
		curgenome = genomes[curgenomeind]
		curnumreads = int(numreads * fractions[curgenomeind])
		
		for n in range(curnumreads):
			
			readlen = read_length + int(gauss(mu=0, sigma=length_stddev))
			ind = randint(0, len(sequences[curgenome])-readlen)
			sampleread = sequences[curgenome][ind:ind+readlen]
			while ('-' in sampleread):
				
				if closegaps:
					nreplace = sampleread.count('-')
					if (ind + nreplace > len(sequences[curgenome])):
						ind = randint(0, numreads-readlen)
						sampleread = sequences[curgenome][ind:ind+readlen]
					else:
						added = sequences[curgenome][ind+readlen:ind+readlen+nreplace]
						if ('-' in added):
							ind = randint(0, numreads-readlen)
							sampleread = sequences[curgenome][ind:ind+readlen]
						else:
							sampleread = sampleread.replace('-','') + added
					
				else:
					ind = randint(0, numreads-readlen)
					sampleread = sequences[curgenome][ind:ind+readlen]
			if (inverted_reads and randint(0, 100)>50):
				reads.append(get_inverse_sequence(sampleread.upper()))
			else:
				reads.append(sampleread.upper())
	
	with open(output_filename, 'w') as outf:
		for read in reads:
			outf.write(read + '\n')
	
samplereads()