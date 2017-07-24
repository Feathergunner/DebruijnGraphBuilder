#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import stochastic as st

ep = np.arange(0.05, 0.4, 0.05)

ks = range(10,51,1)
g = 10000
n = 5000
r = 50

x = ks
ys = []
zs = []
for e in ep:
	ys.append([max(0, st.expected_number_of_correct_kmers_at_position(g, r, n, k, e)) for k in ks])
	
for i in range(len(ep)):
	plt.plot(x,ys[i], label = "correct kmers, error_percentage = "+str(ep[i]))
	
plt.legend(loc=1)
plt.xlabel("k")
plt.ylabel("Expected number of k-mers at each position")
plt.show()