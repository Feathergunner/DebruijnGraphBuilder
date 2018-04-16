#!usr/bin/env python
# -*- coding: utf-8 -*-		

testset_2g = {
	"kmer_lengths" : [13,15,17,19],
	"readlengths" : [50],
    "numbers_of_reads" : [2000],
	"error_type" : "replace",
	"error_rates" : [0.0, 0.1, 0.25, 0.5],
	"name" : "testset_2g"}
	
testset_3g = {
	"kmer_lengths" : [13,15,17,19],
	"readlengths" : [1000],
    "numbers_of_reads" : [100,200,300,500],
	"error_type" : "indel",
	"error_rates" : [5.0],
	"name" : "testset_3g"}

set_cons2g_vlowcov = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [50],
    "numbers_of_reads" : [1500],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "cons_2g_vlowcov"}
	
set_cons2g_vvlowcov = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [50],
    "numbers_of_reads" : [1000],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "cons_2g_vvlowcov"}
	
set_cons2g_lowcov = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [50],
    "numbers_of_reads" : [2000],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "cons_2g_lowcov"}
	
set_cons2g_highcov = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [50],
    "numbers_of_reads" : [4000],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "cons_2g_highcov"}
	
set_cons3g_r = {
	"kmer_lengths" : [13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [500,750,1000,1500],
	"error_type" : "replace",
	"error_rates" : [15.0],
	"name" : "cons_3g_replace"}
	
set_cons3g_i = {
	"kmer_lengths" : [13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [500,750,1000,1500],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"name" : "cons_3g_indel"}
		
set_cons3g_r_lowcov = {
	"kmer_lengths" : [13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [50,100,200,500],
	"error_type" : "replace",
	"error_rates" : [15.0],
	"name" : "cons_3g_replace"}
	
set_cons3g_i_lowcov = {
	"kmer_lengths" : [13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [50,100,200,500],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"name" : "cons_3g_indel"}

allsettings = [ set_cons2g_lowcov,
				set_cons2g_highcov,
				set_cons3g_r_lowcov,
				set_cons3g_i_lowcov,
				set_cons2g_vlowcov,
				set_cons2g_vvlowcov]