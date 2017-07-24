#!/usr/bin/python
# -*- coding: utf-8 -*-

import scipy.misc as scm

def prob_of_complete_sequence_reconstruction(genome_length, read_length, number_of_reads, k_value, error_percentage):
	g = float(genome_length)
	l = float(read_length)
	n = float(number_of_reads)
	k = float(k_value)
	e = float(error_percentage)
	prob_not_covered = ((g-l+k)/g)**n
	prob_covered = 1-(((g-l+k)/g)**n)
	prob_error = (1-(1-e)**(2*(k-1)))
	
	return 1 - (prob_not_covered + (prob_covered * prob_error) * (g-l+k))
	
def prob_kmer_is_correct(k_value, error_percentage):
	k = float(k_value)
	e = float(error_percentage)
	
	return (1-e)**k
	
def expected_kmer_coverage_idd(genome_length, read_length, number_of_reads, k_value):
	g = float(genome_length)
	l = float(read_length)
	n = float(number_of_reads)
	k = float(k_value)
	
	return ((l-k)*n*k)/g
	
def expected_number_of_correct_kmers_at_position(genome_length, read_length, number_of_reads, k_value, error_percentage):
	g = float(genome_length)
	l = float(read_length)
	n = float(number_of_reads)
	k = float(k_value)
	e = float(error_percentage)
	
	prob_correct_kmer = prob_kmer_is_correct(k_value, error_percentage)
	return (l*n/g) * prob_correct_kmer

def prob_kmer_has_specific_number_of_errors(k_value, error_percentage, number_of_errors):
	k = float(k_value)
	e = float(error_percentage)
	n = float(number_of_errors)
	
	p = (1-e)**(k-n) * e * scm.comb(k,n)
	return p
	
def prob_multiple_identical_wrong_kmers(genome_length, read_length, number_of_reads, k_value, error_percentage):
	g = float(genome_length)
	l = float(read_length)
	n = float(number_of_reads)
	k = float(k_value)
	e = float(error_percentage)
	
	p_single_error = prob_kmer_has_specific_number_of_errors(k_value, error_percentage, 1)
	p_not_identical = 1-(((1-e)**(k-1))*e*(1.0/3.0))
	
	return p_single_error * (1 - p_not_identical**((l*n/g)-1))

def expected_number_of_kmers_with_exact_same_errors(genome_length, read_length, number_of_reads, k_value, error_percentage):
	g = float(genome_length)
	l = float(read_length)
	n = float(number_of_reads)
	k = float(k_value)
	e = float(error_percentage)
	
	prob = prob_multiple_identical_wrong_kmers(genome_length, read_length, number_of_reads, k_value, error_percentage)
	
	return prob*(l*n/g)