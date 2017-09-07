#!/usr/bin/python
# -*- coding: utf-8 -*-
import scipy.stats as scs
import scipy.misc as scm

def compute_n(r, k, i):
	n = 0
	max_number_reached = False
	a = [0]*r
	#print a
	while not max_number_reached:
		if sum(a) == i:
			number_valid = True
			num_ajd_zeros = 0
			for j in a:
				if j == 1:
					num_ajd_zeros = 0
				else:
					num_ajd_zeros += 1
					if num_ajd_zeros >= k:
						number_valid = False
			if number_valid:
				#print a
				n += 1
		# increase a
		j = 0
		while j < r and a[j] == 1:
			a[j] = 0
			j+=1
		if j < r:
			a[j] = 1
		else:
			max_number_reached = True
	return n

r = 10
k = 2
for i in range (r/k, r):
	print i
	print compute_n(r, k, i)

def compute_T(n, k, i):
	# n: size of array
	# k: total number of 1's
	# i: (i-1) is maximal allowed number of consecutive 1's
	if k > n:
		# not possible
		return 0
	if i > k:
		# every combination is valid
		return scm.comb(k,n)
		
	if k*(i-1) > n:
		# no legal combination possible:
		return 0
	
	d = n-k
	data_T = [[-1]*k for nn in range(n)]
	# data_T[n_j][k_j] = T(n_j, k_j)
	
	for j in range(d)
		data_T[j] = 
		
	
	
	return compute_T_recursion(n, n-k, i, data_T, 0)
	
def compute_T_recursion(n, k, i, data_T, j):
	#print ("n="+str(n)+" k="+str(k)+" i="+str(i))
	#print data_T
	
	if data_T[j] > 0:
		return data_T
		
	'''
	if n < k/i:
		data_T[j] = 1
		return data_T
	
	if k == 0:
		if i == n:
			data_T[j] = 1
		else:
			data_T[j] = 0
		return data_T
	
	if k <= 0:
		data_T[j] = 0
		return data_T
		
	if n <= 0:
		data_T[j] = 0
		return data_T
	
	if k < n/i:
		data_T[j] = 0
		return data_T
	'''
	
	data_T_tmp = 0		
	for jj in range(0, i-1):
		#print ("j+jj="+str(j+jj))
		data_T = compute_T_recursion(n-j-jj-1, k-j-jj, i, data_T, j+jj)
		data_T_tmp += data_T[j+jj]
		
	data_T[j] = data_T_tmp
	return data_T
	
r = 10
k = 2
for i in range (r/k, r):
	print i
	print compute_T(r, k, i)