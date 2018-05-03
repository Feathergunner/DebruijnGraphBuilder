#!usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import subprocess

path_to_blast_binary = "../blast/bin"

def init_blast_db(dna_fasta_filename):
	cmd = subprocess.Popen([path_to_blast_binary+"\makeblastdb", "-dbtype", "nucl", "-in", dna_fasta_filename])
	
def compute_blast_results(data_dir_sub, dna_fasta_filename, new_data=False):
	for file in os.listdir(data_dir_sub):
		file_name_parts = re.split(r'\.',file)
		if file_name_parts[-1] == "fasta":
			filename = ".".join(file_name_parts[:-1])
			if not filename == "dna":
				if not os.path.isfile(data_dir_sub+"/"+filename+"_blastresults.out") or new_data:
					print ("Consider file "+filename)
					querystring = data_dir_sub+"/"+file
					outputfile = data_dir_sub+"\\"+filename+"_blastresults.out"
					fullcommand = path_to_blast_binary+"\\blastn -query " +querystring+ " -db " + dna_fasta_filename + " -out " + outputfile
					fullcommand = fullcommand.replace("/","\\")
					commandlist = fullcommand.split(' ')
					#print fullcommand
					subprocess.call(commandlist, cwd=".\\")

if __name__ == "__main__":
	data_dir_base = "Output"
	data_dir_sub = data_dir_base+"/example_simplecons"
	
	dna_fasta_filename = data_dir_sub+"/dna.fasta"
	
	init_blast_db(dna_fasta_filename)
	compute_blast_results(data_dir_sub, dna_fasta_filename)