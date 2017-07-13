#!usr/bin/python

import os
import re
import subprocess

path = "Output/errorimpact"

for datafile in os.listdir(path):
	filenameparts = re.split(r'\.',datafile)
	if filenameparts[-1] == "asqg":
		print datafile
		if len(filenameparts) == 2:
			casename = filenameparts[0]
		else:
			casename = ".".join(filenameparts[0:-1])
		print casename
		
		command_bandage_get_image = "../../../../code/Bandage/Bandage image ../../uni/Master/code/debruijn_graph_builder/"+path+"/"+datafile+" ../../uni/Master/code/debruijn_graph_builder/"+path+"/"+casename+".png"
		print command_bandage_get_image
		subprocess.call(command_bandage_get_image.split())
		
