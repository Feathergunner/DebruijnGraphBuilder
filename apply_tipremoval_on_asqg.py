#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os.path

import fast_debruijn_graph_builder as fdgb

source_directory_name = "Output/general_absk"
target_directory_name = "Output/general_absk_wot"

for filename in os.listdir(source_directory_name):
	filenameparts = re.split(r'\.', filename)
	if len(filenameparts) > 1 and filenameparts[-1] == "asqg":
		if not filename in os.listdir(target_directory_name):
			print ("Working on case "+filenameparts[0])
			try:
				debruijn = fdgb.GraphData()
				debruijn.load_from_asqg(source_directory_name+"/"+filename, verbose=True)
				debruijn.remove_tips()
				debruijn.contract_unique_overlaps()
				debruijn.get_asqg_output(target_directory_name+"/"+filename)
				debruijn.get_csv_output(target_directory_name+"/"+filenameparts[0]+".csv")
			except:
				pass