#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os.path
import gc

import sampleReads
import data_io as dio
import fast_debruijn_graph_builder as fdgb
import manjasDefinitionen as md

sampleReads.read_genomes()

check_mdv5 = [md.v5]
set_of_viruses_corona = [md.cv, md.v1]
set_of_viruses_2 = [md.v1, md.v5]
set_of_viruses_4 = [md.v1, md.v7, md.v3, md.v5]

setting_check_mdv5 = {
	"k_absolute_settings" : [16,20,25,30],
	"readlength_settings" : [50],
	"number_of_reads_settings" : [5000],
	"coverage_factors" : [1],
	"error_percentages" : [0.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : check_mdv5,
	"output_dir" : "Output/checkmdv5"}

setting_absk_1 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
	"readlength_settings" : [50, 100, 250, 500, 1000],
	"number_of_reads_settings" : [500, 250, 100, 50, 25],
	"coverage_factors" : [1, 5, 10, 15, 20],
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_of_viruses_2,
	"output_dir" : "Output/general_absk_1"}
	
setting_absk_2 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
	"readlength_settings" : [50, 100, 250, 500, 1000],
	"number_of_reads_settings" : [500, 250, 100, 50, 25],
	"coverage_factors" : [1, 5, 10, 15, 20],
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],#, 10.0, 15.0. 20.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_of_viruses_2,
	"output_dir" : "Output/general_absk_2"}
	
setting_absk_4 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
	"readlength_settings" : [50, 100, 250, 500, 1000],
	"number_of_reads_settings" : [500, 250, 100, 50, 25],
	"coverage_factors" : [1, 5, 10, 15, 20],
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],#, 10.0, 15.0. 20.0],
	"num_different_viruses" : 4,
	"set_of_viruses" : set_of_viruses_4,
	"output_dir" : "Output/general_absk_4"}
	
setting_corona_absk_1 = {
	"k_absolute_settings" : [14,16,18,20,25,30],
	"readlength_settings" : [50, 100],
	"number_of_reads_settings" : [500, 250],
	"coverage_factors" : [1, 5, 10, 15, 20],
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_of_viruses_corona,
	"output_dir" : "Output/general_corona_absk_1"}

setting_corona_vs_bvdv = {
	"k_absolute_settings" : [14,16,18,20,25,30],
	"readlength_settings" : [50, 100],
	"number_of_reads_settings" : [500, 250],
	"coverage_factors" : [1, 5, 10, 15, 20],
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_of_viruses_corona,
	"output_dir" : "Output/general_corona_vs_bvdv"}

def gendata(setting):
	k_absolute_settings = setting["k_absolute_settings"]
	coverage_factors = setting["coverage_factors"]
	readlength_settings = setting["readlength_settings"]
	number_of_reads_settings = setting["number_of_reads_settings"]
	error_percentages = setting["error_percentages"]
	num_different_viruses = setting["num_different_viruses"]
	set_of_viruses = setting["set_of_viruses"]
	output_dir = setting["output_dir"]
	for cf in coverage_factors:
		for i in range(len(readlength_settings)):
			readlength = readlength_settings[i]
			num_reads = [cf*number_of_reads_settings[i]]*num_different_viruses#[cf*nr for nr in number_of_reads_settings[i]]
			
			for error_percentage in error_percentages:
				ep_string = "-".join(re.split(r'\.',str(error_percentage)))
				casename_gen = "bvdv_sample_"+str(readlength)+"_"+str(num_different_viruses)+"_["
				for n_r_i in range(len(num_reads)-1):
					casename_gen += str(num_reads[n_r_i])+","
				casename_gen += str(num_reads[-1])+"]_"+ep_string
				readfilename = "Data/reads_"+casename_gen+".txt"
				
				if not os.path.isfile(readfilename):
					sampleReads.samplereads(output_filename			= readfilename,
											read_length				= readlength,
											set_of_viruses			= set_of_viruses[:num_different_viruses],
											number_of_reads			= num_reads,
											avg_error_percentage	= error_percentage)
				else:
					print ("Reads already exist!")
										
				for k in k_absolute_settings:
					# run garbage collector:
					gc.collect()
					
					casename = casename_gen + "_"+str(k)
					
					if not os.path.isfile(output_dir+"/"+casename+".asqg"):
						print ("Working on case "+casename)
						try:
							reads = dio.get_reads_from_file(filename = readfilename)
							debruijn = fdgb.GraphData(reads, k)
							# delete reads and kmers to save ram:
							debruijn.reads = []
							debruijn.kmers = []
							# run garbage collector:
							gc.collect()
							debruijn.contract_unique_overlaps(verbose=False)
							debruijn.remove_parallel_sequences()
						
							debruijn.get_asqg_output(filename = output_dir+"/"+casename+".asqg")
							debruijn.get_csv_output(filename = output_dir+"/"+casename+".csv")
						except:
							pass
						
					else:
						print ("Data already exists!")

gendata(setting_corona_absk_1)
gendata(setting_corona_vs_bvdv)