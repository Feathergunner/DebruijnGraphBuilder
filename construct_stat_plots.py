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

def construct_consensus_heatmaps(	outputdir,
									outputfilename,
									data,
									datanames,
									xticks,
									yticks,
									xlabels,
									ylabels,
									titles,
									cellformats,
									cbbounds,
									cmaps):
									
	matplotlib.rc('xtick', labelsize=10)
	if (len(ylabels) >= 8):
		matplotlib.rc('ytick', labelsize=8)
	else:
		matplotlib.rc('ytick', labelsize=10)
			
	xlabel_desc = "Kmer Length"# in Nucleotide Bases"
	ylabel_desc = "Error rate in %"
	makesquare = True
	cb = True
	additional_arg={"orientation": "horizontal"} # changes colorbar orientation
	dpi_setting=500
	
	cbar_position = "bottom"
	cbar_size = "5%"
	cbar_pad = 0.5
		
	for i in range(len(data)):
		plt.clf()
		plt.title(titles[i])
		ax = plt.gca()
		divider = make_axes_locatable(ax)
		cax = divider.append_axes(cbar_position, size=cbar_size, pad=cbar_pad)
		if not cbbounds[i][0] == cbbounds[i][1]:
			sns.heatmap(data[i],
						cmap = cmaps[i],
						linewidths = 1,
						square = makesquare,
						yticklabels = yticks,
						xticklabels = xticks,
						annot = True,
						annot_kws = {"size": 8},
						fmt = cellformats[i],
						cbar_kws = additional_arg,
						robust = True,
						cbar = cb,
						ax = ax,
						cbar_ax = cax,
						vmin=cbbounds[i][0],
						vmax=cbbounds[i][1])
		else:
			sns.heatmap(data[i],
						cmap = cmaps[i],
						linewidths = 1,
						square = makesquare,
						yticklabels = yticks,
						xticklabels = xticks,
						annot = True,
						annot_kws = {"size": 8},
						fmt = cellformats[i],
						cbar_kws = additional_arg,
						robust = True,
						cbar = cb,
						ax = ax,
						cbar_ax = cax)
					
		ax.invert_yaxis()
		ax.set(ylabel=ylabels, xlabel=xlabels)
		plt.savefig(outputdir+"\\"+outputfilename+"_"+datanames[i]+".png", dpi=dpi_setting)

def construct_heatmaps_cons_3g(	datadir,
								basename,
								outputdir,
								number_of_reads		= [500, 750, 1000, 1500],
								k_values			= [11,13,15,17,19,21],
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
					print ("no nodes in graph of case: nr:"+str(nr)+", k:"+str(k))
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
	
	construct_consensus_heatmaps(	outputdir = outputdir,
									outputfilename = casename_gen_base,
									data = [seqlengths, covdepths, blast_identity_ratings, blast_num_of_gaps, blast_correct_fraction],
									datanames = ["seqlength", "covdepth", "identityrating", "numofgaps", "correctbasess"],
									xticks = k_values,
									yticks = number_of_reads,
									xlabels = "Kmer length",
									ylabels = "Number of reads",
									titles = ["Average length of reconstructed sequences", "Average coverage depth of sequences", "Average BLAST identity rating in %", "Average number of gaps", "Total fraction of correct bases in %"],
									cellformats = ['.1f', '.1f', '.2f', '.2f', '.1f'],
									cbbounds = [[4500, 5000], [1, 30], [99, 100], [0, 20], [90, 100]],
									cmaps = ["gnuplot", "gnuplot", "gnuplot", "gnuplot_r", "gnuplot"])
	
def construct_heatmaps_cons_2g(	datadir,
								basename,
								outputdir,
								number_of_reads		= 2000,
								k_values			= [13,15,17,19],
								error_rates			= [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
								dna_length			= 5000,
								readlength			= 50,
								error_type			= "replace",
								dimension_of_set	= 1):
	# initialize:
	casename_gen_base = basename+"_rl"+str(readlength)+"_nr"+str(number_of_reads)
	
	avg_seqlengths = np.zeros((len(error_rates), len(k_values)))
	avg_covdepths = np.zeros((len(error_rates), len(k_values)))
	blast_identity_ratings = np.zeros((len(error_rates), len(k_values)))
	blast_correct_fraction = np.zeros((len(error_rates), len(k_values)))
	blast_num_of_gaps = np.zeros((len(error_rates), len(k_values)))
	
	for i in range(dimension_of_set):
		for i_er in range(len(error_rates)):
			for i_k in range(len(k_values)):
				er = error_rates[i_er]
				k = k_values[i_k]
				ep_string = "".join(re.split(r'\.',str("%2.2f" % er)))
				if error_type == "indel":
					ep_string = "ei"+ep_string
				else:
					ep_string = "er"+ep_string
				filename_base = datadir+"/"+casename_gen_base+"_"+ep_string+"_k"+str(k)+"_"+str(i+1)+"_4_singlepath"
				
				gd = []
				gd = ga.GraphData(error_percentage=er, readlength=readlength, num_of_reads=number_of_reads, k_value=k, nodes=[])		
				gd.get_data_from_file(filename_base+".asqg")
				gd.get_data_from_csv(filename_base+".csv")
				
				if gd.num_of_nodes > 0:
					avg_seqlengths[i_er][i_k] += gd.get_avg_seq_length()
					avg_covdepths[i_er][i_k] += gd.get_avg_coverage_depth()
					
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
								
					blast_identity_ratings[i_er][i_k] += 100*blast_identity_rating
					blast_num_of_gaps[i_er][i_k] += num_gaps
					blast_correct_fraction[i_er][i_k] += 100*num_correct_identities
				else:
					print ("no nodes in graph of case: er:"+str(er)+", k:"+str(k))
					avg_seqlengths[i_er][i_k] += 0
					avg_covdepths[i_er][i_k] += 0
					blast_identity_ratings[i_er][i_k] += 0
					blast_correct_fraction[i_er][i_k] += 0
					blast_num_of_gaps[i_er][i_k] += 0
	
	avg_seqlengths /= dimension_of_set
	avg_covdepths /= dimension_of_set
	blast_identity_ratings /= dimension_of_set
	blast_num_of_gaps /= dimension_of_set
	blast_correct_fraction /= dimension_of_set
				
	outputfilename = casename_gen_base
	if not error_type == "replace":
		outputfilename += "_indel"
		
	construct_consensus_heatmaps(	outputdir = outputdir,
									outputfilename = outputfilename,
									data = [avg_seqlengths, avg_covdepths, blast_identity_ratings, blast_num_of_gaps, blast_correct_fraction],
									datanames = ["seqlength", "covdepth", "identityrating", "numofgaps", "correctbasess"],
									xticks = k_values,
									yticks = error_rates,
									xlabels = "Kmer length",
									ylabels = "Error rate in %",
									titles = ["Average length of reconstructed sequences", "Average coverage depth of sequences", "Average BLAST identity rating in %", "Average number of gaps", "Total fraction of correct bases in %"],
									cellformats = ['.1f', '.1f', '.2f', '.2f', '.1f'],
									cbbounds = [[4500, 5000], [1, 30], [99, 100], [0, 20], [90, 100]],
									cmaps = ["gnuplot", "gnuplot", "gnuplot", "gnuplot_r", "gnuplot"])
									
def construct_heatmaps_dbg(	datadir,
							basename,
							outputdir,
							number_of_reads		= 2000,
							k_values			= [13,15,17,19],
							error_rates			= [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
							dna_length			= 5000,
							readlength			= 50,
							error_type			= "replace",
							dimension_of_set	= 1):
	
	# initialize:
	casename_gen_base = basename+"_rl"+str(readlength)+"_nr"+str(number_of_reads)
	
	num_nodes = np.zeros((len(error_rates), len(k_values)))
	num_edges = np.zeros((len(error_rates), len(k_values)))
	frac_edgespernode = np.zeros((len(error_rates), len(k_values)))
	num_comps = np.zeros((len(error_rates), len(k_values)))
	avg_seqlengths = np.zeros((len(error_rates), len(k_values)))
	avg_covdepths = np.zeros((len(error_rates), len(k_values)))
	avg_compsizes = np.zeros((len(error_rates), len(k_values)))
	max_compsizes = np.zeros((len(error_rates), len(k_values)))
	
	for i in range(dimension_of_set):
		for i_er in range(len(error_rates)):
			for i_k in range(len(k_values)):
				er = error_rates[i_er]
				k = k_values[i_k]
				ep_string = "".join(re.split(r'\.',str("%2.2f" % er)))
				if error_type == "indel":
					ep_string = "ei"+ep_string
				else:
					ep_string = "er"+ep_string
				filename_base = datadir+"/"+casename_gen_base+"_"+ep_string+"_k"+str(k)+"_"+str(i+1)
				
				gd = []
				gd = ga.GraphData(error_percentage=er, readlength=readlength, num_of_reads=number_of_reads, k_value=k, nodes=[])		
				gd.get_data_from_file(filename_base+".asqg")
				gd.get_data_from_csv(filename_base+".csv")
				
				if gd.num_of_nodes > 0:
					avg_seqlengths[i_er][i_k] += gd.get_avg_seq_length()
					avg_covdepths[i_er][i_k] += gd.get_avg_coverage_depth()
					num_nodes[i_er][i_k] += gd.num_of_nodes
					num_edges[i_er][i_k] += gd.num_of_edges
					num_comps[i_er][i_k] += gd.get_number_of_components()
					frac_edgespernode[i_er][i_k] += float(gd.num_of_edges)/float(gd.num_of_nodes)
					avg_compsizes[i_er][i_k] += gd.get_avg_component_size()
					max_compsizes[i_er][i_k] += float(gd.get_maximum_component_size())/float(gd.num_of_nodes)
		
				else:
					print ("no nodes in graph of case: er:"+str(er)+", k:"+str(k))
					avg_seqlengths[i_er][i_k] += 0
					avg_covdepths[i_er][i_k] += 0
					num_nodes[i_er][i_k] += 0
					num_edges[i_er][i_k] += 0
					num_comps[i_er][i_k] += 0
					frac_edgespernode[i_er][i_k] += 0
					avg_compsizes[i_er][i_k] += 0
					max_compsizes[i_er][i_k] += 0
	
	avg_seqlengths /= dimension_of_set
	avg_covdepths /= dimension_of_set
	num_nodes /= dimension_of_set
	num_edges /= dimension_of_set
	num_comps /= dimension_of_set
	frac_edgespernode /= dimension_of_set
	avg_compsizes /= dimension_of_set
	max_compsizes /= dimension_of_set
				
	outputfilename = casename_gen_base
	if not error_type == "replace":
		outputfilename += "_indel"
		
	construct_consensus_heatmaps(	outputdir = outputdir,
									outputfilename = outputfilename,
									data = [avg_seqlengths, avg_covdepths, num_nodes, num_edges, num_comps, frac_edgespernode, avg_compsizes, max_compsizes],
									datanames = ["seqlength", "covdepth", "numnodes", "numedges", "numcomps", "edgepernode", "avgcompsize", "maxcompsize"],
									xticks = k_values,
									yticks = error_rates,
									xlabels = "Kmer length",
									ylabels = "Error rate in %",
									titles = ["Average length of sequences", "Average coverage depth of sequences", "Average number of nodes", "Average number of edges", "Average number of components", "fraction of edges per node", "Average component size", "Fraction of nodes in largest component"],
									cellformats = ['.1f', '.1f', '.1f', '.1f', '.1f', '.2f', '.1f', '.2f'],
									cbbounds = [[0, 0], [1, 30], [0, 0], [0, 0], [0, 0], [0.9, 1.4], [0,0], [0.8, 1.0]],
									cmaps = ["gnuplot", "gnuplot", "gnuplot_r", "gnuplot_r", "gnuplot_r", "gnuplot", "gnuplot", "gnuplot"])
	
