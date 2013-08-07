ArgWeaver
---------

Sampling and manipulating genome-wide ancestral recombination graphs (ARG).


## Install

To compile the ArgWeaver commands and library use the Makefile.

```
make
```

Once compiled, to install the ArgWeaver programs (default install in /usr) use:

```
make install
```

To specify your own installation path use:
    
```
make install prefix=$HOME/local
```

ArgWeaver can also run directly from the source directory.  Simply add the
bin/ directory to your PATH environment variable or create symlinks to the 
scripts within bin/ to any directory on your PATH.


## Requirements

The following dependencies must be installed to compile and run the 
ArgWeaver software:

* C++ compiler (e.g. [g++](http://gcc.gnu.org))
* [Python](http://python.org)


### Development Requirements

The following Python libraries are needed for developing ArgWeaver:

```
nose
pyflakes
pep8
scipy
```

These can be installed using 
```
pip install -r requirements-dev.txt
```
