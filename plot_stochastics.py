#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import stochastic as st

ep = np.arange(0,0.7, 0.1)#np.arange(0.05, 0.4, 0.05)

ks = range(10,51,1)
g = 10000
n = 5000
r = 50

def plot_expected_correct_kmers():
	x = ks
	ys = []
	for e in ep:
		ys.append([max(0, st.expected_number_of_correct_kmers_at_position(g, r, n, k, e)) for k in ks])
			
	for i in range(len(ep)):
		plt.plot(x,ys[i], label = "correct kmers, error_percentage = "+str(ep[i]))
			
	plt.legend(loc=1)
	plt.xlabel("k")
	plt.ylabel("Expected number of k-mers at each position")
	plt.show()

def plot_average_number_of_kmers_with_errors():
	x = range(12,26,1)
	ys = []
	for e in ep:
		ys.append([max(0, st.average_number_of_kmers_with_exactly_identical_errors(g, r, n, k, e)) for k in x])

	for i in range(len(ep)):
		plt.plot(x,ys[i], label = "error_percentage = "+str(ep[i]))

	plt.legend(loc=1)
	plt.xlabel("k")
	plt.ylabel("Expected avg. number of identical kmers with errors")
	plt.show()

def plot_expected_number_of_kmers():
	x = ks
	ys = []
	for e in ep:
		ys.append([max(0, st.expected_number_of_kmers(g, r, n, k, e)) for k in ks])
			
	for i in range(len(ep)):
			plt.plot(x,ys[i], label = "number of kmers, error_percentage = "+str(ep[i]))
			
	plt.legend(loc=1)
	plt.xlabel("k")
	plt.ylabel("Expected total number of different k-mers")
	plt.show()

plot_average_number_of_kmers_with_errors()