reassemble_from_23andme
=======================

Software for reassembling whole chromosomes from 23andme data.

A key limitation is that you are only considering point mutations, and only those captured by the 23andme chip. This has been tested on v3 files for now.

Note that it doesn't handle X chromosomes entirely well for men, due to inclusion of the pseudoautosomal region from the Y chromosome being listed as part of the X chromosome in the 23andme raw data file.

This software is written in Python 3 and requires BioPython to be installed for Python 3.

You need to download the chromosomes separately, from [here](ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/) and you are after only the ones named hs_ref_GRCh37.p13_chr*.fa.gz
