experiments.py: Run experiments on de Bruijn graphs and consensus reconstruction

python experiments.py {option=parameter}

GENERAL OPTIONS:
These options set the algorithm and other options.

Possible options are:
	algo:	sets the algorithm to use.
		algo=sc		use simple consensus construction
		algo=lcfr	use low coverage feature removal
		algo=cr		use coverage refining
		algo=dbg	only construct de Bruijn graph, no consensus construction
		
	onlydata:	only the data will be generated, no evalutating plots will be generated
	onlystat:	only plots will be generated. Needs pyplot and BLAST
		(default: both data and statistics will be generated)
	
	nothr:		compute experiments single-threaded
	
	verbose:	enable additional output
	
	overwrite:	set whether existing results will be overwritten
		(deault: nothing will be overwritten)
		overwrite=all	everything will be overwritten
		overwrite=res	only results will be overwritten, existing dna and reads will be re-used
		overwrite=add	only missing results will be added, everything else will be kept
		
	arc:		set level of archiving
		arc=none	no archive will be generated
		arc=results	an archive containing all results will be generated (default)
		arc=all		an archive with all generated data will be generated
		
	name:		give an individual name for this experiment (not possible for predefined settings)
	
	d:			set number of repetitions for each individual experiment
		(default: d=1)
	
OPTIONS CONCERNING THE EXPERIMENTAL SETTING:
You can either specify parameters by the command line or use a predefined set of parameters.

1) Specify setting by command line parameters:
For example:
Runs a single experiment with the simple consensus algorithm on 1000 reads of length 100 with 0.2% error rate and k=13:
	python experiments.py nr=1000 rl=100 k=13 e=0.2 algo=sc
	
The same option can be used multiple times, then all specified parameters will be used. For example:
Construct two de Bruijn graphs from 500 reads of length 500 with 2.5% error rate and k=13 and k=17, respectively:
	python experiments.py nr=500 rl=500 e=2.5 k=13 k=17 algo=dbg
	
Possible options are:
	nr:		set number of reads
		(default: nr=2000)
	rl:		set length of reads
		(default: rl=100)
	k:		set kmer length
		(default: k=15)
	e:		set error rate
		(default: e=0.1)
		
2) Use a predefined set (see datasets_experiments.py)
For example:
Repeat the first experiment from chapter 3 of the thesis (de Bruijn graphs from short reads, low coverage):
	python experiments.py algo=dbg d=10 set=1
		
Possible options are:
	test:	uses a special test set, depending on the algorithm
	
	set:	specify a set from datasets_experiments.py
		set=i	chooses the i-th setting from the list allsettings defined in datasets_experiments.py