Code for experiments for my master thesis
# Untangle the Wool: Constructing Consensus Sequences with De Bruijn Graphs

## The important files are the following:
- `data_gen.py`: contains methods to generate articial genome strings and to sample reads (as described in Section 3.2.1 of the thesis).
- `fast_debruijn_graph_builder.py`: contains the de Bruijn graph construction algorithm and all methods that modify de Bruijn graph (e. g. simplication, error removal, position label computation).
- `dbg_consensus_construction.py`: contains the implementations of Algorithm 9 from my thesis (simple consensus construction), Algorithm 10 from my thesis (low coverage feature removal) and Algorithm 12 (coverage rening).
- `experiments.py`: a program to automatically run experiments with the algorithms.

## How to run the experiments:
`experiments.py` can be used to easily reproduce all experiments from my thesis. For each experiment from the thesis, the respective setting (number of reads, readlength, kmer-lengths etc.) is dened in the le datasets_experiments.py. These sets can be referenced to by a command line arguments for the program experiments.py. The table below shows how to reproduce all experiments from this thesis, sorted by their section.
The parameter `algo` sets the algorithm:
- `dbg` sets Algorithm 4
- `sc` sets Algorithm 9
- `lcfr` sets Algorithm 10
- `cr` sets Algorithm 12
The parameter `set` sets the setting; in particular `set=i` chooses the i-th setting from the list `allsettings` that is dened in the file `datasets_experiments.py`. Note that one can also easily perform experiments with custom settings by defining a new setting and adding it to the list `allsettings`.
The parameter `d` sets the dimension, i. e. the number of repetitions for each experiment.

For a full list of all parameters, see `readme_experiments.txt` or run `python experiments.py help`.

## The output:
Note that by default, `experiments.py` also computes statistics and generates the respective plots.
It is necessary that BLAST is installed and the path to the BLAST executable is correctly set in the global variable `path_to_blast_binary` in the file `apply_blast.py`. To generate the plots, matplotlib and the library seaborn (http://seaborn.pydata.org/) are needed. This can be deactivated by the additional parameter `onlydata`. In contrast, the parameter `onlystat` deactivates the experiments, in this case only the graphs are generated (if the data already exists).

All generated data will be saved in a specic subfolder of ./Output. The name of this subfolder is constructed from the name of the algorithm, the name of the setting and the dimension. Also by default `experiments.py` will not overwrite exisiting data, i. e. it will not start if a directory of the same name already exists. With the parameter `overwrite` one can define different levels of overwriting.
- With `overwrite=add`, only missing results will be added.
- With `overwrite=res`, all experiments will be computed anew, but the original genome sequence and already generated reads will be reused.
- With `overwrite=all`, everything will be overwritten by new experiments.

## Recreate the experiments from my thesis:

Section (and Figure) in Thesis | Program Call
-------------------------------|---------------------------
Section 3.3.2                  | `python experiments.py algo=dbg set=1 d=10`
Section 3.3.2                  | `python experiments.py algo=dbg set=2 d=10`
Section 3.3.2                  | `python experiments.py algo=dbg set=3 d=10`
Section 3.3.2                  | `python experiments.py algo=dbg set=4 d=10`
Section 3.3.3, Figure 3.15     | `python experiments.py algo=dbg set=5 d=10`
Section 3.3.3, Figure 3.16     | `python experiments.py algo=dbg set=6 d=10`
Section 5.3.1                  | `python experiments.py algo=sc set=8 d=10`
Section 5.3.2                  | `python experiments.py algo=sc set=9 d=10`
Section 5.3.3, Figure 5.12     | `python experiments.py algo=sc set=10 d=10`
Section 5.3.3, Figure 5.14     | `python experiments.py algo=sc set=11 d=10`
Section 6.5.1                  | `python experiments.py algo=lcfr set=11 d=10`
Section 6.5.2                  | `python experiments.py algo=cr set=12 d=1`
Section 6.5.2, Figure 6.10     | `python experiments.py algo=cr set=13 d=10`
Section 6.5.2, Figure 6.11     | `python experiments.py algo=cr set=14 d=10`

