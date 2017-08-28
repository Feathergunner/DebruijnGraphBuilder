
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
	"name" : "bvdv_absk_1",
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
	"name" : "bvdv_absk_largereads_1",
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
	"name" : "bvdv_absk_2",
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
	"name" : "bvdv_absk_4",
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
	"name" : "corona_absk_1",
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
	"name" : "coronavsbvdv_absk",
	"output_dir" : "Output/general_corona_vs_bvdv"}

setting_absk_large_corona_1 = {
	"k_absolute_settings" : [14,16,18,20,25,30],
	"readlength_settings" : [4000, 8000, 10000],
	"number_of_reads_settings" : [50, 25,20],
	"coverage_factors" : [1, 2, 5],
	"error_type" : "replace",
	"error_percentages" : [1.0, 5.0, 10.0, 15.0],
	"num_different_viruses" : 1,
	"set_of_viruses" : set_corona_1,
	"name" : "corona_asbk_largereads_1",
	"output_dir" : "Output/general_corona_absk_large_1"}