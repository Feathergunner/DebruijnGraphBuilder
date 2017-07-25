#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import graph_analyzer as ga
import stochastic as st

sourcedirs_absk = ["general_absk"]#, "general_absk_wot"]
aourcedirs_relk = ["general_relk"]#, "general_relk_wot"]


def plots_absk(sourcedirs):
	for sourcedir in sourcedirs:
		for nr in [1000, 5000, 10000, 15000, 20000]:
			gaga = ga.GraphAnalyzer(sourcedir)
			gacr = ga.CaseRestriction(readlengths=[50], number_of_reads=[nr], error_percentages=[0, 0.1, 0.25, 0.5, 1, 2, 5], k_values=[12,14,16,18,20,25,30])
			gaga.get_data(case_restrict=gacr, verbose=True)
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twinx()
			fig.set_size_inches(20,20)
			
			data = "error_percentage"
			gaga.lineplot("num_of_nodes", data, axis=ax1, legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_edges", data, style='--', axis=ax1, legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_components", data, style=':', axis=ax2, legend_pos=2)
			fig.savefig("Output/"+sourcedir+"/"+gacr.construct_casename()+"_"+data+".png")
			plt.close()
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twinx()
			fig.set_size_inches(20,20)
		
			data = "k_value"
			gaga.lineplot("num_of_nodes", data, axis=ax1, legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_edges", data, axis=ax1, style='--', legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_components", data, style=':', axis=ax2, legend_pos=2)
			fig.savefig("Output/"+sourcedir+"/"+gacr.construct_casename()+"_"+data+".png")
			plt.close()
	
		readlength_settings = [50, 100, 250, 500, 1000]
		number_of_reads_settings = [1000, 500, 200, 100, 50]
		coverage_factor = [1,5,10]
		for ep in [0, 0.1, 0.25, 0.5, 1, 2, 5]:
			for i in range(len(coverage_factor)):
				nrs = [nr*coverage_factor[i] for nr in number_of_reads_settings]
				gaga = ga.GraphAnalyzer(sourcedir)
				for j in range(len(readlength_settings)):
					gacr = ga.CaseRestriction(error_percentages=[ep], number_of_reads=[nrs[j]], readlengths=[readlength_settings[j]])
					gaga.get_data(case_restrict=gacr, verbose=True)
				
				data = "k_value"
				#fig, ax1 = plt.subplots()
				fig = plt.figure()
				ax1 = fig.add_subplot(111)
				ax2 = ax1.twinx()
				fig.set_size_inches(20,20)
				#plt.figure(figsize=(40,40), dpi=80)
				
				gaga.lineplot("num_of_nodes", data, axis=ax1, legend_pos=1)
				#plt.gca().set_prop_cycle(None)
				gaga.lineplot("num_of_edges", data, style='--', axis=ax1, legend_pos=1)
				#plt.gca().set_prop_cycle(None)
				gaga.lineplot("num_of_components", data, style=':', axis=ax2, legend_pos=2)
				
				fig.savefig("Output/"+sourcedir+"/k(-1)_ep("+str(ep)+")_coveragefactor("+str(coverage_factor[i])+")_"+data+".png", dpi=80)
				plt.close()

def plots_relk(sourcedirs):
	k_rel_values = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
	for sourcedir in sourcedirs:
		for nr in [1000, 5000, 10000, 15000, 20000]:
			gaga = ga.GraphAnalyzer(sourcedir)
			k_values = [50*kr for kr in k_rel_values]
			gacr = ga.CaseRestriction(readlengths=[50], number_of_reads=[nr], error_percentages=[0, 0.1, 0.25, 0.5, 1, 2, 5], k_values=k_values)
			gaga.get_data(case_restrict=gacr, verbose=True)
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twinx()
			fig.set_size_inches(20,20)
			
			data = "error_percentage"
			gaga.lineplot("num_of_nodes", data, axis=ax1, legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_edges", data, style='--', axis=ax1, legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_components", data, style=':', axis=ax2, legend_pos=2)
			fig.savefig("Output/"+sourcedir+"/"+gacr.construct_casename()+"_"+data+".png")
			plt.close()
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twinx()
			fig.set_size_inches(20,20)
		
			data = "rel_k_value"
			gaga.lineplot("num_of_nodes", data, axis=ax1, legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_edges", data, axis=ax1, style='--', legend_pos=1)
			#plt.gca().set_prop_cycle(None)
			gaga.lineplot("num_of_components", data, style=':', axis=ax2, legend_pos=2)
			fig.savefig("Output/"+sourcedir+"/"+gacr.construct_casename()+"_"+data+".png")
			plt.close()
	
		readlength_settings = [50, 100, 150, 200]
		number_of_reads_settings = [500,250,166,125]
		coverage_factor = [1,5,10]
		for ep in [0, 0.1, 0.25, 0.5, 1, 2, 5]:
			for i in range(len(coverage_factor)):
				nrs = [nr*coverage_factor[i] for nr in number_of_reads_settings]
				gaga = ga.GraphAnalyzer(sourcedir)
				for j in range(len(readlength_settings)):
					k_values = [readlength_settings[j]*kr for kr in k_rel_values]
					gacr = ga.CaseRestriction(error_percentages=[ep], number_of_reads=[nrs[j]], readlengths=[readlength_settings[j]])
					gaga.get_data(case_restrict=gacr, verbose=True)
				
				data = "rel_k_value"
				#fig, ax1 = plt.subplots()
				fig = plt.figure()
				ax1 = fig.add_subplot(111)
				ax2 = ax1.twinx()
				fig.set_size_inches(20,20)
				#plt.figure(figsize=(40,40), dpi=80)
				
				gaga.lineplot("num_of_nodes", data, axis=ax1, legend_pos=1)
				#plt.gca().set_prop_cycle(None)
				gaga.lineplot("num_of_edges", data, style='--', axis=ax1, legend_pos=1)
				#plt.gca().set_prop_cycle(None)
				gaga.lineplot("num_of_components", data, style=':', axis=ax2, legend_pos=2)
				
				fig.savefig("Output/"+sourcedir+"/k(-1)_ep("+str(ep)+")_coveragefactor("+str(coverage_factor[i])+")_"+data+".png", dpi=80)
				plt.close()
				
plots_relk(aourcedirs_relk)