#!/usr/bin/env python
# -*- coding: utf-8 -*-

from random import randint

with open('bvdv_selected_filtered.linsi.aln') as fh:
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
readlen = 50
reads = []

fractions = [0.5, 0.5, 0, 0, 0]


for curgenomeind in range(len(fractions)):
	
	curgenome = genomes[curgenomeind]
	curnumreads = int(numreads * fractions[curgenomeind])
	
	for n in range(curnumreads):
		
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
		reads.append(sampleread.upper())

with open('samplereads.txt','w') as outf:
	for read in reads:
		outf.write(read + '\n')
