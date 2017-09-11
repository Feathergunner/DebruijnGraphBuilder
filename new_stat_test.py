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

def compute_T(n, k, i):
	# n: size of array
	# k: total number of 1's
	# i: (i-1) is maximal allowed number of consecutive 1's
	
	data_T = [[-1]*(k+1) for nn in range(n+1)]
	return compute_T_recursion(n, k, i, data_T)[n][k]
	
def compute_T_recursion(n, k, i, data_T):
	print "n="+str(n)
	print "k="+str(k)
	print "i="+str(i)
	if data_T[n][k] > 0:
		return data_T
	if k > n:
		#print "k>n"
		# not possible
		data_T[n][k] = 0
		return data_T
	if i > k:
		#print "i>k"
		# every combination is valid
		data_T[n][k] = int(scm.comb(n,k))
		print "tnk="+str(data_T[n][k])
		return data_T 
	
	data_T_tmp = 0		
	for j in range(i):
		data_T = compute_T_recursion(n-j-1, k-j, i, data_T)
		data_T_tmp += data_T[n-j-1][k-j]
		
	data_T[n][k] = data_T_tmp
	return data_T

r = 4
k = 2
for i in range (r/k, r+1):
	print ""
	print i
	print compute_n(r, k, i)
	print compute_T(r, r-k, i)
	#compute_T(r, k, i)
