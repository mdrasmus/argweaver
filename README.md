ARGweaver
=========

*Sampling and manipulating genome-wide ancestral recombination graphs (ARG).*  

The ARGweaver software package contains programs and libraries for
sampling and manipulating ancestral recombination graphs (ARG). An ARG
is a rich data structure for representing the ancestry of DNA
sequences undergoing recombination.

*ARGweaver citation:*
Matthew D. Rasmussen, Adam Siepel. Genome-wide inference of ancestral
recombination graphs. 2013. arXiv:1306.5110 [q-bio.PE]


## Requirements

The following dependencies must be installed to compile and run the
ARGweaver software:

- C++ compiler (e.g. [g++](http://gcc.gnu.org))
- [Python](http://python.org)


## Install

To compile the ARGweaver commands and library use the Makefile.

```
Make
```

Once compiled, to install the ARGweaver programs (default install in
`/usr`) use:

```
make install
```

To specify your own installation path use:

```
make install prefix=$HOME/local
```

ARGweaver can also run directly from the source directory.  Simply add the
`bin/` directory to your `PATH` environment variable or create symlinks to the
scripts within `bin/` to any directory on your `PATH`.


## Quick Start

Here is a brief example of an ARG simulation and analysis.
To generate simulated data containing a set of DNA sequences and an
ARG describing their ancestry, the following command can be used:

```
arg-sim \
    -k 8 -L 100000 \
    -N 10000 -r 1.6e-8 -m 1.8e-8 \
    -o example/example
```

This will create an ARG with 8 sequences each 100kb in length evolving in
a population of effective size 10,000 (diploid), with recombination rate
1.6e-8 recombinations/site/generation and mutation rate 1.8e-8 
mutations/generation/site. The output will be stored in the following files:

```
example/example.arg   -- an ARG stored in *.arg format
example/example.sites -- sequences stored in *.sites format
```

To infer an ARG from the simulated sequences, the following command 
can be used:

```
arg-sample \
    -s example/example.sites \
    -N 10000 -r 1.6e-8 -m 1.8e-8 \
    --ntimes 20 --maxtime 200e3 -c 10 -n 100 \
    -o example/example.sample/out
```

This will use the sequences in `example/example.sites`, assume the
same population parameters as the simulation (i.e. `-N 10000 -r 1.6e-8
-m 1.8e-8`), and several sampling specific options (i.e. 20 discretized
time steps, a maximum time of 200,000 generations, a compression of 10bp
for the sequences, and 100 sampling iterations. After sampling the following
files will be generated:

```
example/example.sample/out.log
example/example.sample/out.stats
example/example.sample/out.0.smc.gz
example/example.sample/out.10.smc.gz
example/example.sample/out.20.smc.gz
...
example/example.sample/out.100.smc.gz
```

`out.log` contains a log of the sampling procedure, `out.stats` contains 
various ARG statistics (e.g. number of recombinations, ARG posterior 
probability, etc), and `out.0.smc.gz` through `out.100.smc.gz` contain 
11 samples of an ARG in *.smc file format.

To estimate the time to most recent common ancestor (TMRCA) across
these samples, the following command can be used:

```
arg-extract-tmrca example/example.sample/out.%d.smc.gz \
    > example/example.tmrca.txt
```

This will create a tab-delimited text file containing six columns:
chromosome, start, end, posterior mean TMRCA (generations),
lower 2.5 percentile TMRCA, and upper 97.5 percentile TMRCA. The first
four columns define a track of TMRCA across the genomic region in
BED file format.


## Development

The following Python libraries are needed for developing ARGweaver:

```
nose
pyflakes
pep8
```

These can be installed using

```
pip install -r requirements-dev.txt
```
