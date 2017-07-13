#!usr/bin/python

import os
import re
import subprocess

datapath = "Output/errorimpact"
bandagepath = "../Bandage"

for datafile in os.listdir(datapath):
	filenameparts = re.split(r'\.',datafile)
	if filenameparts[-1] == "asqg":
		print datafile
		if len(filenameparts) == 2:
			casename = filenameparts[0]
		else:
			casename = ".".join(filenameparts[0:-1])
		print casename
		
		command_bandage_get_image = bandagepath+"/Bandage image "+datapath+"/"+datafile+" "+datapath+"/"+casename+".png"
		print command_bandage_get_image
		subprocess.call(command_bandage_get_image.split())
		
