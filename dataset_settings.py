
import manjasDefinitionen as md

set_corona_1 = [md.cv, md.v1]
set_bvdv_2 = [md.v1, md.v5]
set_bvdv_4 = [md.v1, md.v7, md.v3, md.v5]

setting_check_mdv5 = {
	"k_absolute_settings" : [16,20,25,30],
	"readlength_settings" : [50],
	"number_of_reads_settings" : [5000],
	"coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : [md.v5],
	"name" : "checkmdv5",
	"output_dir" : "Output/checkmdv5"}

setting_absk_1 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
    "readlength_settings" : [50, 100, 250, 500, 1000],
    "number_of_reads_settings" : [500, 250, 100, 50, 25],
    "coverage_factors" : [1, 5, 10, 15, 20],
	"error_type" : "replace",
    "error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
    "num_different_viruses" : 1,
    "set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-absk-1",
    "output_dir" : "Output/general_absk_1"}

setting_absk_large_1 = {
    "k_absolute_settings" : [14,16,18,20,25,30],
    "readlength_settings" : [4000, 8000, 10000],
    "number_of_reads_settings" : [50, 25, 20],
    "coverage_factors" : [1, 2, 5],
	"error_type" : "replace",
    "error_percentages" : [0.0, 1.0, 5.0, 10.0, 15.0],
    "num_different_viruses" : 1,
    "set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-absk-largereads-1",
    "output_dir" : "Output/general_absk_large_1"}
	
setting_absk_2 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
	"readlength_settings" : [50, 100, 250, 500, 1000],
	"number_of_reads_settings" : [500, 250, 100, 50, 25],
	"coverage_factors" : [1, 5, 10, 15, 20],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],#, 10.0, 15.0. 20.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-absk-2",
	"output_dir" : "Output/general_absk_2"}
	
setting_absk_4 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
	"readlength_settings" : [50, 100, 250, 500, 1000],
	"number_of_reads_settings" : [500, 250, 100, 50, 25],
	"coverage_factors" : [1, 5, 10, 15, 20],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],#, 10.0, 15.0. 20.0],
	"num_different_viruses" : 4,
	"set_of_viruses" : set_bvdv_4,
	"name" : "bvdv-absk-4",
	"output_dir" : "Output/general_absk_4"}
	
setting_corona_absk_1 = {
	"k_absolute_settings" : [14,16,18,20,25,30],
	"readlength_settings" : [50, 100],
	"number_of_reads_settings" : [5000, 2500],
	"coverage_factors" : [1, 2, 5],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_corona_1,
	"name" : "corona-absk-1",
	"output_dir" : "Output/general_corona_absk_1"}

setting_corona_vs_bvdv = {
	"k_absolute_settings" : [14,16,18,20,25,30],
	"readlength_settings" : [50, 100],
	"number_of_reads_settings" : [5000, 2500],
	"coverage_factors" : [1, 2, 5],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_corona_1,
	"name" : "corona-vs-bvdv-absk",
	"output_dir" : "Output/general_corona_vs_bvdv"}

corona_large_absk_1 = {
	"k_absolute_settings" : [14,16,18,20,25,30,40],
	"readlength_settings" : [8000, 10000],
	"number_of_reads_settings" : [250,200],
	"coverage_factors" : [1, 2, 5],
	"error_type" : "indel",
	"error_percentages" : [1.0, 5.0, 10.0, 15.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_corona_1,
	"name" : "corona-largereads-asbk-1",
	"output_dir" : "Output/general_corona_absk_large_1"}

corona_large_manyreads = {
	"k_absolute_settings" : [130,40],
	"readlength_settings" : [8000, 10000],
	"number_of_reads_settings" : [250,200],
	"coverage_factors" : [4, 5, 10],
	"error_type" : "indel",
	"error_percentages" : [5.0, 10.0, 15.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_corona_1,
	"name" : "corona-largereads-manyreads"}
	
bvdv_absk_1 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
    "readlength_settings" : [50, 100],
    "number_of_reads_settings" : [2500, 1250],
    "coverage_factors" : [1, 2, 4],
	"error_type" : "replace",
    "error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
    "num_different_viruses" : 1,
    "set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-absk-1",
    "output_dir" : "Output/bvdv_absk_1"}
	
bvdv_absk_2 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
	"readlength_settings" : [50, 100],
    "number_of_reads_settings" : [2500, 1250],
    "coverage_factors" : [1, 2, 4],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-absk-2",
	"output_dir" : "Output/bvdv_absk_2"}
	
bvdv_absk_4 = {
	"k_absolute_settings" : [10,12,14,16,18,20,25,30],
	"readlength_settings" : [50, 100],
    "number_of_reads_settings" : [2500, 1250],
    "coverage_factors" : [1, 2, 4],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"num_different_viruses" : 4,
	"set_of_viruses" : set_bvdv_4,
	"name" : "bvdv-absk-4",
	"output_dir" : "Output/bvdv_absk_4"}
	
bvdv_large_absk_2 = {
	"k_absolute_settings" : [20,30],
	"readlength_settings" : [6000, 8000],
    "number_of_reads_settings" : [200, 150],
    "coverage_factors" : [1, 2],
	"error_type" : "indel",
	"error_percentages" : [10.0, 15.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-largereads-absk-2"}

bvdv_large_absk_4 = {
	"k_absolute_settings" : [20,30],
	"readlength_settings" : [6000, 8000],
	"number_of_reads_settings" : [200, 150],
	"coverage_factors" : [1, 2],
	"error_type" : "indel",
	"error_percentages" : [10.0, 15.0],
	"num_different_viruses" : 4,
	"set_of_viruses" : set_bvdv_4,
	"name" : "bvdv-largereads-absk-4"}

bvdv_absk_1_test_recons = {
	"k_absolute_settings" : [20,25],
    "readlength_settings" : [50, 100],
    "number_of_reads_settings" : [2500, 1250],
    "coverage_factors" : [2,4],
	"error_type" : "replace",
    "error_percentages" : [0.1, 0.5, 1.0, 5.0],
    "num_different_viruses" : 1,
    "set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-absk-1-test-recons"}
	
largereads_test_recons = {
	"k_absolute_settings" : [30],
	"readlength_settings" : [5000],
	"number_of_reads_settings" : [250],
	"coverage_factors" : [1],
	"error_type" : "indel",
	"error_percentages" : [15.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_bvdv_2,
	"name" : "largereads-recons-test"}

largereads_test_recons_2 = {
	"k_absolute_settings" : [30],
	"readlength_settings" : [5000],
	"number_of_reads_settings" : [500],
	"coverage_factors" : [1],
	"error_type" : "indel",
	"error_percentages" : [5.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_bvdv_2,
	"name" : "largereads-recons-test2"}

reads_for_sebastian_corona = {
	"k_absolute_settings" : [30],
	"readlength_settings" : [2000, 8000],
	"number_of_reads_settings" : [20000, 20000],
	"coverage_factors" : [1, 10],
	"error_type" : "indel",
	"error_percentages" : [5.0, 15.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_corona_1,
	"name" : "largereads-largenamount-corona"}

reads_for_sebastian_bvdv = {
	"k_absolute_settings" : [30],
	"readlength_settings" : [2000, 8000],
	"number_of_reads_settings" : [20000, 20000],
	"coverage_factors" : [1, 10],
	"error_type" : "indel",
	"error_percentages" : [5.0, 15.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_bvdv_2,
	"name" : "largereads-largeamount-bvdv"}

bvdv_large_absk_2_clustercut = {
	"k_absolute_settings" : [30],
	"readlength_settings" : [1000,2000,4000,6000],
    "number_of_reads_settings" : [100,100],
    "coverage_factors" : [1,2],
	"error_type" : "indel",
	"error_percentages" : [5.0, 10.0, 15.0],
	"num_different_viruses" : 2,
	"set_of_viruses" : set_bvdv_2,
	"name" : "bvdv-largereads-absk-2-cuttest"}

examples_random_reads_1_smallreads_r = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [50,100,250],
    "number_of_reads_settings" : [50,100,250,500,1000,2000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_smallreads_r_smallset = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [100],
    "number_of_reads_settings" : [2000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_smallreads_r_smallset_lowcov = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [50],
    "number_of_reads_settings" : [1000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_r = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [500,1000,2000],
    "number_of_reads_settings" : [5,10,25,50,100,200],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_r_smallset = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [200],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [5.0, 10.0, 15.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_i_smallset = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [200],
    "coverage_factors" : [1],
	"error_type" : "indel",
	"error_percentages" : [5.0, 10.0, 15.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_i_highcov = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [500,1000,2000,5000],
    "coverage_factors" : [1],
	"error_type" : "indel",
	"error_percentages" : [15.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_r_highcov = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [500,1000,2000,5000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [15.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_i_lowcov = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [50,100,200,500],
    "coverage_factors" : [1],
	"error_type" : "indel",
	"error_percentages" : [15.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_r_lowcov = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [50,100,200,500],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [15.0],
	"name" : "examples_random_reads_1"}
	
examples_random_reads_1_largereads_i = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [500,1000],
    "number_of_reads_settings" : [100,200],
    "coverage_factors" : [1],
	"error_type" : "indel",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0],
	"name" : "examples_random_reads_1i"}
	
examples_random_reads_1_smallreads_hc_r = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [50],
    "number_of_reads_settings" : [5000,10000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "examples_random_reads_1"}
	
set_conscons_shortreads_lowcov = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [50],
    "number_of_reads_settings" : [2000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "conscons_sr_lc"}
	
set_conscons_shortreads_highcov = {
	"k_absolute_settings" : [13,15,17,19,21,23,25,29,33,37,41],
	"readlength_settings" : [50],
    "number_of_reads_settings" : [4000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0],
	"name" : "conscons_sr_hc"}
	
set_conscons_longreads_highcov_r = {
	"k_absolute_settings" : [13,15,17,19,21,23,25],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [200,300,500,1000],
    "coverage_factors" : [1],
	"error_type" : "replace",
	"error_percentages" : [15.0],
	"name" : "conscons_lr_hc_r"}
	
set_conscons_longreads_highcov_i = {
	"k_absolute_settings" : [13,15,17,19,21,23,25],
	"readlength_settings" : [1000],
    "number_of_reads_settings" : [200,300,500,1000],
    "coverage_factors" : [1],
	"error_type" : "indel",
	"error_percentages" : [15.0],
	"name" : "conscons_lr_hc_i"}