#!/usr/bin/python
# -*- coding: utf-8 -*-

import scipy.stats as scs
import scipy.misc as scm

import math

def compute_n(r, k, i):
	# brute-force-check
	n = 0
	max_number_reached = False
	a = [0]*r
	#print a
	while not max_number_reached:
		if sum(a) == k:
			number_valid = True
			num_ajd_ones = 0
			for j in a:
				if j == 0:
					num_ajd_ones = 0
				else:
					num_ajd_ones += 1
					if num_ajd_ones >= i:
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
	#print "n="+str(n)
	#print "k="+str(k)
	#print "i="+str(i)
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
		return data_T 
	
	data_T_tmp = 0		
	for j in range(i):
		data_T = compute_T_recursion(n-j-1, k-j, i, data_T)
		data_T_tmp += data_T[n-j-1][k-j]
		
	data_T[n][k] = data_T_tmp
	return data_T
	
def compute_T_nonrecurcive(nn, kk, i):
	max_n = int(max(nn))
	max_k = int(max(kk))
	
	data_T = [[-1]*(max_k+1) for nnn in range(max_n+1)]
	for n in range(max_n+1):
		for k in range(max_k+1):
			#print str(n)+","+str(k)
			if k > n:
				data_T[n][k] = 0
			elif i > k:
				data_T[n][k] = int(scm.comb(n,k))
			else:
				data_T[n][k] = 0
				for j in range(i):
					data_T[n][k] += data_T[n-j-1][k-j]
	
	'''
	for n in nn:
		for k in kk:
			print ("T("+str(n)+","+str(k)+") = "+str(data_T[n][k]))	
	'''
	return data_T

if __name__ == '__main__':
	e = 0.15
	'''
	r = 10
	k = 5
	for i in range (r/k, r+1):
		print ""
		print i
		print compute_n(r, k, i)
		print compute_T(r, r-k, i)
		#compute_T(r, k, i)
	'''
	readlengths = [500,1000,2000,5000,10000,15000,20000]
	r = [int(n*(1-e)) for n in readlengths]
	
	for k in [20,30,40,50]:
		print ("k="+str(k))
		data_T = compute_T_nonrecurcive(readlengths, r, k)
		for j in range(len(readlengths)):
			n = readlengths[j]
			i = r[j]
			log_m_e = math.log(data_T[n][i],2)
			log_m_t = math.log(scm.comb(n,i),2)
			log_p = log_m_e-log_m_t
			p = 2**(log_p)
			print ("log_m_e = "+str(log_m_e))
			print ("log_m_t = "+str(log_m_t))
			print ("P(n="+str(n)+",k="+str(k)+") = 2^"+str(log_p)+" ~= "+str(p))