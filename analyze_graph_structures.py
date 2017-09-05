#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import graph_analyzer as ga
import stochastic as st

file_extension = "pdf"
figure_size = 30 # 20

sourcedirs_absk = ["general_absk"]#, "general_absk_wot"]
sourcedirs_relk = ["general_relk"]#, "general_relk_wot"]

def multiple_lineplot_nodes_edges_components(graph_analyzer, data, filename, y1=["num_of_nodes", "num_of_edges"], y2="num_of_components"):
	fig = plt.figure(figsize=(figure_size,figure_size))
	if not y2 == "":
		ax1 = fig.add_subplot(111)
		ax2 = ax1.twinx()
	#fig.set_size_inches(20,20)
	graph_analyzer.lineplot(y1[0], data, axis=ax1, legend_pos=1)
	if len(y1) > 1:
		graph_analyzer.lineplot(y1[1], data, style='--', axis=ax1, legend_pos=1)
	if not y2 == "":
		graph_analyzer.lineplot(y2, data, style=':', axis=ax2, legend_pos=2)
	fig.savefig(filename+"."+file_extension, dpi=80, format=file_extension)
	plt.close()

def plots_absk(sourcedirs, readlength_settings=[50, 100, 250, 500, 1000], number_of_reads_settings=[500, 250, 100, 50, 25], coverage_factors=[1,5,10,15,20], eps=[0, 0.1, 0.25, 0.5, 1, 2, 5], ks=[12,14,16,18,20,25,30], num_genomes=1):
	for sourcedir in sourcedirs:
		for i in range(len(readlength_settings)):
			rl = readlength_settings[i]
			for cov in coverage_factors:
				nr = number_of_reads_settings[i]*cov
				gaga = ga.GraphAnalyzer(sourcedir)
				gacr = ga.CaseRestriction(readlengths=[rl], number_of_reads=[nr], error_percentages=eps, k_values=ks, num_genomes=num_genomes)
				print ("Case "+gacr.construct_casename())
				gaga.get_data(case_restrict=gacr, verbose=True)
				
				if len(gaga.graphdatas) > 0:
					multiple_lineplot_nodes_edges_components(gaga, "error_percentage", filename="Output/"+sourcedir+"/plots/num_nec_"+gacr.construct_casename()+"_error_percentage")
					multiple_lineplot_nodes_edges_components(gaga, "k_value", filename="Output/"+sourcedir+"/plots/num_nec_"+gacr.construct_casename()+"_k_value")
				
					multiple_lineplot_nodes_edges_components(gaga, "error_percentage", filename="Output/"+sourcedir+"/plots/avgseqlengths_"+gacr.construct_casename()+"_error_percentage", y1=["avg_seq_lengths"], y2 = "num_of_components")
	
		for ep in eps:
			for i in range(len(coverage_factors)):
				nrs = [nr*coverage_factors[i] for nr in number_of_reads_settings]
				gaga = ga.GraphAnalyzer(sourcedir)
				for j in range(len(readlength_settings)):
					gacr = ga.CaseRestriction(error_percentages=[ep], number_of_reads=[nrs[j]], readlengths=[readlength_settings[j]])
					gaga.get_data(case_restrict=gacr, verbose=True)
					
				if len(gaga.graphdatas) > 0:
					multiple_lineplot_nodes_edges_components(gaga, "k_value", "Output/"+sourcedir+"/plots/num_nec_coverage_"+gacr.construct_casename()+"_k_value")

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
				#multiple_lineplot_nodes_edges_components(gaga, "error_percentage", "Output/"+sourcedir+"/plots/"+gacr.construct_casename()+"_error_percentage")
				multiple_lineplot_nodes_edges_components(gaga, "rel_k_value", "Output/"+sourcedir+"/plots/num_nec_"+gacr.construct_casename()+"_rel_k_value")

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
					multiple_lineplot_nodes_edges_components(gaga, "rel_k_value", "Output/"+sourcedir+"/plots/num_nec_coverage_"+gacr.construct_casename()+"_rel_k_value")
					
def plot_distribution(sourcedirs, k_values=[12,14,16,18,20,25,30], rl=50, nr=5000, eps=[0, 0.1, 0.25, 0.5, 1, 2, 5], num_genomes=1):
	for sourcedir in sourcedirs:
		for k_value in k_values:
			#print k_value
			gaga = ga.GraphAnalyzer(sourcedir)
			gacr = ga.CaseRestriction(readlengths=[rl], number_of_reads=[nr], error_percentages=eps[1:], k_values=[k_value], num_genomes=num_genomes)
			gaga.get_data(case_restrict=gacr, verbose=True)
			
			if len(gaga.graphdatas) > 0:
				fig = plt.figure()
				fig.set_size_inches(figure_size,figure_size)
				gaga.distribution_lineplot(data="degree_dist", legend_pos=1)
				
				#fig.savefig("Output/"+sourcedir+"/plots/distributions_deg_k("+str(k_value)+")_ep("+str(eps)+")_nr("+str(nr)+")_rl("+str(rl)+")."+file_extension, dpi=80)
				fig.savefig("Output/"+sourcedir+"/plots/distributions_deg_"+gacr.construct_casename()+"_k_value."+file_extension, dpi=80)
				plt.close()
				
				fig = plt.figure()
				fig.set_size_inches(figure_size,figure_size)
				gaga.distribution_lineplot(data="seq_length", legend_pos=2)
				#fig.savefig("Output/"+sourcedir+"/plots/distributions_seqlength_k("+str(k_value)+")_ep("+str(eps)+")_nr("+str(nr)+")_rl("+str(rl)+")."+file_extension, dpi=80)
				fig.savefig("Output/"+sourcedir+"/plots/distributions_seqlength_"+gacr.construct_casename()+"_k_value."+file_extension, dpi=80)
				plt.close()
				
		for ep in eps:
			#print ep
			gaga = ga.GraphAnalyzer(sourcedir)
			gacr = ga.CaseRestriction(readlengths=[rl], number_of_reads=[nr], error_percentages=[ep], k_values=k_values, num_genomes=num_genomes)
			gaga.get_data(case_restrict=gacr, verbose=True)
			
			if len(gaga.graphdatas) > 0:
				fig = plt.figure()
				fig.set_size_inches(figure_size,figure_size)
				gaga.distribution_lineplot(data="degree_dist", legend_pos=1)
				#fig.savefig("Output/"+sourcedir+"/plots/distributions_deg_k("+str(k_values)+")_ep("+str(ep)+")_nr("+str(nr)+")_rl("+str(rl)+")."+file_extension, dpi=80)
				fig.savefig("Output/"+sourcedir+"/plots/distributions_deg_"+gacr.construct_casename()+"_error_percentage."+file_extension, dpi=80)
				plt.close()
				
				fig = plt.figure()
				fig.set_size_inches(figure_size,figure_size)
				gaga.distribution_lineplot(data="seq_length", legend_pos=2)
				#fig.savefig("Output/"+sourcedir+"/plots/distributions_seqlength_k("+str(k_values)+")_ep("+str(ep)+")_nr("+str(nr)+")_rl("+str(rl)+")."+file_extension, dpi=80)
				fig.savefig("Output/"+sourcedir+"/plots/distributions_seqlength_"+gacr.construct_casename()+"_error_percentage."+file_extension, dpi=80)
				plt.close()
				
def plot_avg_seqlength(sourcedirs, num_genomes=1):
	eps = [0.1, 0.25, 0.5, 1, 2, 5]
	ks = [12,14,16,18,20,25,30]
	for [rl,nr] in [[50,500], [50,2500], [50,5000], [100, 2500]]:
		gaga = ga.GraphAnalyzer(sourcedirs[0])
		gacr = ga.CaseRestriction(readlengths=[rl], number_of_reads=[nr], error_percentages=eps, k_values=ks, num_genomes=num_genomes)
		gaga.get_data(case_restrict=gacr, verbose=True)
		if len(gaga.graphdatas) > 0:
			multiple_lineplot_nodes_edges_components(gaga, "error_percentage", filename="Output/"+sourcedirs[0]+"/plots/avgseqlengths_"+gacr.construct_casename()+"_error_percentage", y1=["avg_seq_lengths"], y2="num_of_components")
				
def plot_distribution_from_setting(setting):
	for cf in setting["coverage_factors"]:
		for i in range(len(setting["readlength_settings"])):
			plot_distribution(
				sourcedirs = [setting["name"]],
				k_values = setting["k_absolute_settings"],
				rl = setting["readlength_settings"][i],
				nr = cf * setting["number_of_reads_settings"][i],
				eps = setting["error_percentages"],
				num_genomes = setting["num_different_viruses"])

def plots_absk_from_setting(setting):
	plots_absk(
		sourcedirs = [setting["name"]],
		readlength_settings = setting["readlength_settings"],
		number_of_reads_settings = setting["number_of_reads_settings"],
		coverage_factors = setting["coverage_factors"],
		eps = setting["error_percentages"],
		ks = setting["k_absolute_settings"],
		num_genomes = setting["num_different_viruses"])
				
#plots_relk(["general_relk_1"])
#plots_absk(["general_absk_2"], num_genomes=2)
#plot_avg_seqlength(["general_absk_2"], num_genomes=2)
#plot_distribution(["general_absk_2"], rl=50, nr=5000, num_genomes=2)
#plot_distribution(["general_absk_2"], rl=100, nr=2500, num_genomes=2)

#plot_distribution(sourcedirs_absk,[12,14,16,18,20,25,30])
#plot_distribution(sourcedirs_relk,[10,15,20,25,30,35,40])

import dataset_settings as ds
plot_distribution_from_setting(ds.corona_large_absk_1)
#plots_absk_from_setting(ds.corona_large_absk_1)