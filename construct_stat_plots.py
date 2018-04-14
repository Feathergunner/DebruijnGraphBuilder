#!usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import seaborn as sns

import graph_analyzer as ga

def construct_heatmaps_3g_lcfr(	datadir,
								basename,
								outputdir,
								number_of_reads		= [500, 1000, 2000, 5000],
								k_values			= [13,15,17,19],
								error_rate			= 15.0,
								dna_length			= 5000,
								readlength			= 1000,
								dimension_of_set	= 1):
	#datadir = "Output/examples_for_thesis/chapter_thirdgen/itlowcovrem_exp_thesis"
	#datadir = "..\\tex\\tex_thesis\\img\\examples_thirdgen_reconstruction_lowcovremove"
	#datadir = "Output/itlowcovrem_exp_0413"
	#basename = "itlowcovrem"
	#outputdir = "Output/itlowcovrem_exp_0413"
	
	#datadir = "Output/examples_for_thesis/example_thirdgen_reconstruct_edgeremove_test0412"
	#outputdir = "Output/examples_for_thesis/example_thirdgen_reconstruct_edgeremove_test0412"
	#basename = "thirdgen_reconstruct"
	
	# initialize:
	ep_string = "".join(re.split(r'\.',str("%2.2f" % error_rate)))			
	casename_gen_base = basename+"_rl"+str(readlength)
	
	seqlengths = np.zeros((len(number_of_reads), len(k_values)))
	covdepths = np.zeros((len(number_of_reads), len(k_values)))
	blast_identity_ratings = np.zeros((len(number_of_reads), len(k_values)))
	blast_num_of_gaps = np.zeros((len(number_of_reads), len(k_values)))
	blast_correct_fraction = np.zeros((len(number_of_reads), len(k_values)))
	
	for i in range(dimension_of_set):
		for k_i in range(len(k_values)):
			for n_i in range(len(number_of_reads)):
				k = k_values[k_i]
				nr = number_of_reads[n_i]
				#filename_base = datadir+"/"+casename_gen_base+"_nr"+str(nr)+"_ei"+ep_string+"_k"+str(k)+"_rcsb0.1_singlepath"
				filename_base = datadir+"/"+casename_gen_base+"_nr"+str(nr)+"_ei"+ep_string+"_k"+str(k)+"_"+str(i+1)+"_4_singlepath"
				#filename_base = datadir+"/"+casename_gen_base+"_nr"+str(nr)+"_ei"+ep_string+"_newtipremoval_test_k"+str(k)+"_rcsb0.5_singlepath"
				
				gd = []
				gd = ga.GraphData(error_percentage=error_rate, readlength=readlength, num_of_reads=nr, k_value=k, nodes=[])		
				gd.get_data_from_file(filename_base+".asqg")
				gd.get_data_from_csv(filename_base+".csv")
				
				if gd.num_of_nodes > 0:
					seqlengths[n_i][k_i] += gd.get_avg_seq_length()
					covdepths[n_i][k_i] += gd.get_avg_coverage_depth()
				else:
					print ("no nodes in graph of case: er:"+str(error_rate)+", k:"+str(k))
					seqlengths[n_i][k_i] += 0
					covdepths[n_i][k_i] += 0
					
				blastfile = filename_base+"_blastresults.out"
				#print outputfile
				if not os.path.isfile(blastfile):
					print ("Error! No Blast Results found! File '"+blastfile+"' is missing!")
					# works only if only one sequence:
					dnadbstring = datadir+"/dna.fasta"
					querystring = filename_base+".fasta"
					fullcommand = path_to_blast_binary+"\\blastn -query " +querystring+ " -db " + dnadbstring + " -out " + blastfile
					fullcommand = fullcommand.replace("/","\\")
					commandlist = fullcommand.split(' ')
					subprocess.call(commandlist, cwd=".\\")
				
				num_correct_identities = 0
				num_gaps = 0
				blast_identity_rating = 0
				with open(blastfile, 'r') as blastresultfile:
					linenumber = 0
					no_hits_found = False
					for line in blastresultfile:
						linenumber += 1
						if "No hits found" in line:
							no_hits_found = True
							num_correct_identities = 0.0
							num_gaps = 0.0
							break;
						if not no_hits_found and "Identities" in line:
							# format of this line (mind the single leading whitespace!)
							#  Identities = 4776/4779 (99%), Gaps = 0/4779 (0%)
							data = re.split(r'\s', line)
							blast_identity_rating = (float)(re.split(r'/', data[3])[0])/(float)(re.split(r'/', data[3])[1])
							num_correct_identities = (float)(re.split(r'/', data[3])[0])/(float)(dna_length)
							num_gaps = (int)(re.split(r'/', data[7])[0])
							
				blast_identity_ratings[n_i][k_i] += 100*blast_identity_rating
				blast_num_of_gaps[n_i][k_i] += num_gaps
				blast_correct_fraction[n_i][k_i] += 100*num_correct_identities
				
	seqlengths /= dimension_of_set
	covdepths /= dimension_of_set
	blast_identity_ratings /= dimension_of_set
	blast_num_of_gaps /= dimension_of_set
	blast_correct_fraction /= dimension_of_set
	
	outputfilename = casename_gen_base
	cb = True
	makesquare = True
	additional_arg={}#{"orientation": "horizontal"} # changes colorbar orientation
	#xlabel_desc = "Kmer Length"
	xlabel_desc = "Kmer length"
	ylabel_desc = "Number of reads"
	
	plt.clf()
	ax = sns.heatmap(seqlengths, cmap="gnuplot", linewidths=1, square=makesquare, yticklabels=number_of_reads, xticklabels=k_values, annot=True, fmt='.0f', cbar_kws=additional_arg, robust=True, cbar=cb, vmin=4000, vmax=5000)
	ax.invert_yaxis()
	ax.set(ylabel=ylabel_desc, xlabel=xlabel_desc)
	plt.title("Length of reconstructed Sequences")#in Nucleotide Bases")
	#plt.show()
	plt.savefig(outputdir+"\\"+outputfilename+"_seqlength.png", dpi=500)
	
	plt.clf()
	ax = sns.heatmap(covdepths, cmap="gnuplot", linewidths=1, square=makesquare, yticklabels=number_of_reads, xticklabels=k_values, annot=True, fmt='.0f', cbar_kws=additional_arg, robust=True, cbar=cb, vmin=0, vmax=30)
	ax.invert_yaxis()
	ax.set(ylabel=ylabel_desc, xlabel=xlabel_desc)
	plt.title("Coverage Depths of reconstructed Sequences")
	#plt.show()
	plt.savefig(outputdir+"\\"+outputfilename+"_covdepth.png", dpi=500)
	
	plt.clf()
	ax = sns.heatmap(blast_identity_ratings, cmap="gnuplot", linewidths=1, square=makesquare, yticklabels=number_of_reads, xticklabels=k_values, annot=True, annot_kws={"size": 8}, fmt='.2f', cbar_kws=additional_arg, robust=True, cbar=cb, vmin=99, vmax=100)
	ax.invert_yaxis()
	ax.set(ylabel=ylabel_desc, xlabel=xlabel_desc)
	plt.title("BLAST identity rating in %")
	#plt.show()
	plt.savefig(outputdir+"\\"+outputfilename+"_identityrating.png", dpi=500)
	
	plt.clf()
	ax = sns.heatmap(blast_correct_fraction, cmap="gnuplot", linewidths=1, square=makesquare, yticklabels=number_of_reads, xticklabels=k_values, annot=True, annot_kws={"size": 8}, fmt='.2f', cbar_kws=additional_arg, robust=True, cbar=cb, vmin=80, vmax=100)
	ax.invert_yaxis()
	ax.set(ylabel=ylabel_desc, xlabel=xlabel_desc)
	plt.title("Total fraction of correct nucleotide bases in %")
	#plt.show()
	plt.savefig(outputdir+"\\"+outputfilename+"_correctbasess.png", dpi=500)
	
	plt.clf()
	ax = sns.heatmap(blast_num_of_gaps, cmap="gnuplot_r", linewidths=1, square=makesquare, yticklabels=number_of_reads, xticklabels=k_values, annot=True, fmt='.1f', cbar_kws=additional_arg, robust=True, cbar=cb, vmin=0, vmax=20)
	ax.invert_yaxis()
	ax.set(ylabel=ylabel_desc, xlabel=xlabel_desc)
	plt.title("Number of gaps in reconstructed sequence")
	#plt.show()
	plt.savefig(outputdir+"\\"+outputfilename+"_numofgaps.png", dpi=500)
	
	print ("wrote images to "+outputdir+"\\"+outputfilename)