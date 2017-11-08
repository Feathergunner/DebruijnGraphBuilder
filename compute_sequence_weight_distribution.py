#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import math
import matplotlib.pyplot as plt


def compute_weight_distribution(filename_input, num_buckets = 100, logscale=False):
	with open(filename_input) as inputfile:
		lines = inputfile.readlines()
	weights = []
	minweight = -1
	maxweight = 0
	for i in range(len(lines)):
		if i > 0:
			data = re.split(r',', lines[i])
			w = int(data[2])
			weights.append(w)
			if minweight < 0 or w < minweight:
				minweight = w
			if w > maxweight:
				maxweight = w

	bucket_size = (maxweight/num_buckets)+1
	weight_distirbution_buckets = [0]*num_buckets
	for w in weights:
		weight_distirbution_buckets[w/bucket_size] += 1

	buckets = [bucket_size*i for i in range(num_buckets)]
	if logscale:
		return [buckets, [math.log10(k+1) for k in weight_distirbution_buckets]]
	else:
		return [buckets, weight_distirbution_buckets]

if __name__ == '__main__':
	data_dir = "Output/corona_allreads"
	#inputfile = data_dir+"/corona_realreads_n-1_k40.csv"
	inputfile = data_dir+"/corona_realreads_k40_w50_k15_p1000_merged_k13.csv"
	#inputfile = data_dir+"/corona_realreads_k40_w5_k25_p1000_merged_k23.csv"

	number_of_buckets = 300
	logscale = False

	wd = compute_weight_distribution(inputfile, number_of_buckets, logscale)
	#print wd
	plt.plot(wd[0], wd[1])
	if logscale:
		plt.ylabel("log(#sequences)")
	else:
		plt.ylabel("#sequences")
	plt.xlabel("weight of sequence")
	plt.show()
