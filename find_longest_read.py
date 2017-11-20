#!usr/bin/python
# -*- coding: utf-8 -*-

import data_handler as dio

def find_longest_read(filename):
	reads = dio.get_reads_from_fastq_file(filename)
	i = 0
	longreads = []
	for r in reads:
		if len(r) > 10000:
			longreads.append([i, len(r)])
		i += 1
	
	return longreads
	
if __name__ == '__main__':
	filename = "Data/hcov229e_only.fq"
	r = find_longest_read(filename)
	print r