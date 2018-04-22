#!usr/bin/env python
# -*- coding: utf-8 -*-	

testset_2g = {
	"kmer_lengths" : [13,15,17,19],
	"readlengths" : [50],
    "numbers_of_reads" : [2000],
	"error_type" : "replace",
	"error_rates" : [0.0, 0.1, 0.25, 0.5],
	"name" : "testset_2g"}
	
testset_3g_locofere = {
	"kmer_lengths" : [13,15,17,19],
	"readlengths" : [1000],
    "numbers_of_reads" : [50,100,200],
	"error_type" : "indel",
	"error_rates" : [5.0],
	"name" : "testset_3g_lcfr"}
	
testset_3g_covref = {
	"kmer_lengths" : [17,19],
	"readlengths" : [1000],
    "numbers_of_reads" : [100,200],
	"error_type" : "indel",
	"error_rates" : [5.0],
	"k_base" : [20, 25],
	"k_part" : [13,15],
	"number_of_parts" : [10,15],#[10,25,50],
	"overlaps" : [3,5],#[5,10,20],
	"name" : "testset_3g_cr"}

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
	"readlengths" : [100],
    "numbers_of_reads" : [2000],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "cons_2g_highcov"}

set_cons3g_r_hc = {
	"kmer_lengths" : [13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [500,750,1000],#,1500],
	"error_type" : "replace",
	"error_rates" : [15.0],
	"k_base" : [21,25,39],
	"number_of_parts" : [10,15,20,50],
	"overlaps" : [3,5],#20],
	"name" : "cons_3g_replace_hc"}
	
set_cons3g_i_hc = {
	"kmer_lengths" : [13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [500,750,1000,1500],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"k_base" : 25,
	"number_of_parts" : [10,25,50],
	"overlaps" : [5,10,20],
	"name" : "cons_3g_indel_hc"}
	
set_cons3g_i_hc_rc_t = {
	"kmer_lengths" :[13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [500,750,1000],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"k_base" : [19,25,31],
	"k_part" : [13,15],
	"number_of_parts" : [10,15,20],
	"overlaps" : [3,5],
	"name" : "cons_3g_rc_indel_rangetest"}
	
set_cons3g_i_hc_rc = {
	"kmer_lengths" : [17],#,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [300,400,500,600,1000],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"k_base" : [19],
	"k_part" : [13],
	"number_of_parts" : [8,10],
	"overlaps" : [3,5],#,7],
	"name" : "cons_3g_rc_indel"}
	
set_cons3g_i_hc_rc_small = {
	"kmer_lengths" : [17],#[13,15,17,19,21],
	"readlengths" : [1000],
    "numbers_of_reads" : [200,300,400],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"k_base" : [19],
	"k_part" : [11],
	"number_of_parts" : [30,40,50],
	"overlaps" : [3,5,10],
	"name" : "cons_3g_rc_indel_small"}
	
set_cons3g_r_hc_large = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [1000],
    "numbers_of_reads" : [500,750,1000,1500],
	"error_type" : "replace",
	"error_rates" : [15.0],
	"k_base" : 25,
	"number_of_parts" : [10,25,50],
	"overlaps" : [5,10,20],
	"name" : "cons_3g_replace_hc"}
	
set_cons3g_i_hc_large = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [1000],
    "numbers_of_reads" : [500,750,1000,1500],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"k_base" : 25,
	"number_of_parts" : [10,25,50],
	"overlaps" : [5,10,20],
	"name" : "cons_3g_indel_hc"}
	
set_cons3g_r_hc_detail = {
	"kmer_lengths" : [13, 14,15,16,17,18,19],
	"readlengths" : [1000],
    "numbers_of_reads" : [500, 600, 700, 800, 900, 1000, 1100, 1200, 1300],
	"error_type" : "replace",
	"error_rates" : [15.0],
	"k_base" : 25,
	"number_of_parts" : [10,25,50],
	"overlaps" : [5,10,20],
	"name" : "cons_3g_replace_hc"}
	
set_cons3g_i_hc_detail = {
	"kmer_lengths" : [13, 14,15,16,17,18,19],
	"readlengths" : [1000],
    "numbers_of_reads" : [400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300],
	"error_type" : "indel",
	"error_rates" : [15.0],
	"name" : "cons_3g_indel_hc"}

set_dbgbasic_sr_lc = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [50],
    "numbers_of_reads" : [1000],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "set_dbgbasic_sr_lc"}
	
set_dbgbasic_sr_hc = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [100],
    "numbers_of_reads" : [2000],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "set_dbgbasic_sr_hc"}
	
set_dbgbasic_sr_vhc = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [50],
    "numbers_of_reads" : [10000],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "set_dbgbasic_sr_vhc"}
	
set_dbgbasic_lr_lc = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [500],
    "numbers_of_reads" : [100],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0],
	"name" : "set_dbgbasic_lr_lc"}
	
set_dbgbasic_lr_hc = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [1000],
    "numbers_of_reads" : [200],
	"error_type" : "replace",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0],
	"name" : "set_dbgbasic_lr_hc"}
	
set_dbgbasic_lr_hc_i = {
	"kmer_lengths" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlengths" : [1000],
    "numbers_of_reads" : [200],
	"error_type" : "indel",
	"error_rates" : [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0],
	"name" : "set_dbgbasic_lr_hc_i"}
	
allsettings = [ set_cons2g_lowcov,		# 1
				set_cons2g_highcov,		# 2
				set_cons3g_r_hc,		# 3
				set_cons3g_i_hc,		# 4
				set_cons2g_vlowcov,		# 5
				set_dbgbasic_sr_lc,		# 6
				set_dbgbasic_sr_hc,		# 7
				set_dbgbasic_lr_lc,		# 8
				set_dbgbasic_lr_hc,		# 9
				set_dbgbasic_sr_vhc,	# 10
				set_dbgbasic_lr_hc_i,	# 11
				set_cons3g_r_hc_large,	# 12
				set_cons3g_i_hc_large,	# 13
				set_cons3g_r_hc_detail, # 14
				set_cons3g_i_hc_detail,	# 15
				set_cons3g_i_hc_rc_t,	# 16
				set_cons3g_i_hc_rc,		# 17
				set_cons3g_i_hc_rc_small # 18
				]