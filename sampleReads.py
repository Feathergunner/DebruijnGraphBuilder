#!/usr/bin/env python
# -*- coding: utf-8 -*-

from random import randint, gauss
import re
import os.path

import data_io
import manjasDefinitionen as md

def read_genomes(input_filename="Data/bvdv_selected_filtered.linsi.aln", output_filedir="Data/genomes/"):
	with open(input_filename) as fh:
		lines = fh.readlines()
	
	alphabet = ['a','c','g','t']
	
	genomes = []
	sequences = {}
	for line in lines:
		if (line != ''):
			if (line[0] == '>'):
				gname = line.strip()[1:]
				gname_short = re.split(r'\s', gname)[0]
				genomes.append(gname_short)
				sequences[gname_short] = ''
			else:
				i = 0
				while i < len(line):
					if line[i] not in alphabet:
						line = line[:i]+line[(i+1):]
					else:
						i += 1
				sequences[gname_short] += line.strip()
	
	for genome_name in genomes:
		output_filename = output_filedir+"genome_"+genome_name+".txt"
		with open(output_filename, 'w') as outf:
			outf.write(sequences[genome_name]+"\n")

def samplereads(input_filedir="Data/genomes/",
				output_filename="samplereads.txt",
				read_length=50,
				length_stddev=0,
				set_of_viruses=[md.v1, md.v5],
				number_of_reads=[5000,5000],
				avg_error_percentage=0,
				inverted_reads=False):
				
	print ("Sample reads from genomes ...")
		
	alphabet = ['a','c','g','t']
		
	numreads = sum(number_of_reads)
	if len(set_of_viruses) == 0:
		set_of_viruses = range(len(number_of_reads))
	reads = []
	
	# construct reads:
	for curgenomeind_id in range(len(set_of_viruses)):
		
		curgenomename = set_of_viruses[curgenomeind_id]
		curnumreads = number_of_reads[curgenomeind_id]
		
		input_filename = input_filedir+"genome_"+curgenomename+".txt"
		if not os.path.isfile(input_filename):
			print ("Error! Wrong name of genome!")
			
		with open(input_filename) as inputfile:
			lines = inputfile.readlines()
		genome = lines[0]
		
		print ("Constructing reads from genome " + str(curgenomename))
		
		for n in range(curnumreads):
			
			readlen = read_length + int(gauss(mu=0, sigma=length_stddev))
			ind = randint(0, len(genome)-readlen)
			sampleread = genome[ind:ind+readlen]
			'''
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
			'''
			
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
		
#read_genomes()
#samplereads()