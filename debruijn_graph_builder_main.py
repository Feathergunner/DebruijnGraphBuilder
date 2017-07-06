#!usr/bin/python

import debruijn_graph_builder as dgb
import data_io as dio

#dna = dio.genereate_dna(length=500)
#reads, alignment = dio.genereate_reads(dna, coverage=20, avg_read_length=40, remove_pct=0, mutation_pct=0.2, mutation_alphabet=["A","C","G","T"], both_directions=True)

reads = dio.get_reads_from_file(filename = "Data/samplereads.txt")
k = 14

#dio.print_alignment(dna, alignment)

debruijn = GraphData(reads, k, verbose = False)
#debruijn.get_asqg_output(filename = "graph_pre_contract")
debruijn.contract_unique_overlaps(verbose = False)
#debruijn.get_asqg_output(filename = "graph_post_contract")
debruijn.remove_parallel_sequences(verbose = False)
debruijn.get_asqg_output(filename = "graph_post_unify")
debruijn.get_csv_output()
