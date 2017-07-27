#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import graph_analyzer as ga
import stochastic as st

sourcedirs_absk = ["general_absk"]#, "general_absk_wot"]
sourcedirs_relk = ["general_relk"]#, "general_relk_wot"]

def multiple_lineplot_nodes_edges_components(graph_analyzer, data, filename):
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	fig.set_size_inches(20,20)
	graph_analyzer.lineplot("num_of_nodes", data, axis=ax1, legend_pos=1)
	graph_analyzer.lineplot("num_of_edges", data, style='--', axis=ax1, legend_pos=1)
	graph_analyzer.lineplot("num_of_components", data, style=':', axis=ax2, legend_pos=2)
	fig.savefig(filename, dpi=80)
	plt.close()

def plots_absk(sourcedirs):
	readlength_settings = [50, 100, 250, 500, 1000]
	number_of_reads_settings = [1000, 500, 200, 100, 50]
	coverage_factors = [1,5,10,15,20]
	for sourcedir in sourcedirs:
		for rl in readlength_settings:
			for i in range(len(coverage_factors)):
				nr = number_of_reads_settings[i]*coverage_factors[i]
				gaga = ga.GraphAnalyzer(sourcedir)
				gacr = ga.CaseRestriction(readlengths=[rl], number_of_reads=[nr], error_percentages=[0, 0.1, 0.25, 0.5, 1, 2, 5], k_values=[12,14,16,18,20,25,30])
				gaga.get_data(case_restrict=gacr, verbose=True)
				
				if len(gaga.graphdatas) > 0:
					multiple_lineplot_nodes_edges_components(gaga, "error_percentage", "Output/"+sourcedir+"/plots/"+gacr.construct_casename()+"_error_percentage.png")
					multiple_lineplot_nodes_edges_components(gaga, "k_value", "Output/"+sourcedir+"/plots/"+gacr.construct_casename()+"_k_value.png")
				
					fig = plt.figure()
					ax1 = fig.add_subplot(111)
					ax2 = ax1.twinx()
					fig.set_size_inches(20,20)
					data = "error_percentage"
					gaga.lineplot("avg_seq_lengths", data)
					fig.savefig("Output/"+sourcedir+"/plots/"+gacr.construct_casename()+"_"+data+"_seqlengths.png")
					plt.close()
					
	
		for ep in [0, 0.1, 0.25, 0.5, 1, 2, 5]:
			for i in range(len(coverage_factors)):
				nrs = [nr*coverage_factors[i] for nr in number_of_reads_settings]
				gaga = ga.GraphAnalyzer(sourcedir)
				for j in range(len(readlength_settings)):
					gacr = ga.CaseRestriction(error_percentages=[ep], number_of_reads=[nrs[j]], readlengths=[readlength_settings[j]])
					gaga.get_data(case_restrict=gacr, verbose=True)
					
				if len(gaga.graphdatas) > 0:
					multiple_lineplot_nodes_edges_components(gaga, "k_value", "Output/"+sourcedir+"/plots/k(-1)_ep("+str(ep)+")_coveragefactor("+str(coverage_factors[i])+")_"+data+".png")

def plots_relk(sourcedirs):
	k_rel_values = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
	for sourcedir in sourcedirs:
		for nr in [1000, 5000, 10000, 15000, 20000]:
			gaga = ga.GraphAnalyzer(sourcedir)
			rl = 50
			k_values = [rl*kr for kr in k_rel_values]
			gacr = ga.CaseRestriction(readlengths=[rl], number_of_reads=[nr], error_percentages=[0, 0.1, 0.25, 0.5, 1, 2, 5], k_values=k_values)
			gaga.get_data(case_restrict=gacr, verbose=True)
			
			if len(gaga.graphdatas) > 0:
				multiple_lineplot_nodes_edges_components(gaga, "error_percentage", "Output/"+sourcedir+"/plots/"+gacr.construct_casename()+"_error_percentage.png")
				multiple_lineplot_nodes_edges_components(gaga, "rel_k_value", "Output/"+sourcedir+"/plots/"+gacr.construct_casename()+"_rel_k_value.png")

		readlength_settings = [50, 100, 150, 200]
		number_of_reads_settings = [1000,500,333,250]
		coverage_factor = [1,5,10]
		for ep in [0, 0.1, 0.25, 0.5, 1, 2, 5]:
			for i in range(len(coverage_factor)):
				print ("Case ep = "+str(ep)+", coverage_factor = "+str(coverage_factor[i]))
				nrs = [nr*coverage_factor[i] for nr in number_of_reads_settings]
				gaga = ga.GraphAnalyzer(sourcedir)
				for j in range(len(readlength_settings)):
					k_values = [readlength_settings[j]*kr for kr in k_rel_values]
					gacr = ga.CaseRestriction(error_percentages=[ep], number_of_reads=[nrs[j]], readlengths=[readlength_settings[j]], k_values=k_values)
					gaga.get_data(case_restrict=gacr, verbose=True)
				
				if len(gaga.graphdatas) > 0:
					multiple_lineplot_nodes_edges_components(gaga, "rel_k_value", "Output/"+sourcedir+"/plots/k(-1)_ep("+str(ep)+")_coveragefactor("+str(coverage_factor[i])+")_rel_k_value.png")
				
plots_relk(sourcedirs_relk)
plots_absk(sourcedirs_absk)