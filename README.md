# Sigmond Analysis Suite #

This repository contains an analysis software suite for the analysis of Monte Carlo data in lattice QCD.
Currently it only supports two-point correlation functions.

## Requires 

- cmake? -> unsure if I needed to install these because my Ubuntu VM was screwy, please try without first
- c++? -> unsure if I needed to install these because my Ubuntu VM was screwy, please try without first

## Instructions for pip install ## 

```
git clone git@github.com:andrewhanlon/sigmond.git
git checkout pip
pip install .
```
With pip, pybindings will be built. 


## Old README Below ##

The code is XML driven.
A single XML input file is the only argument passed to the main program.
There is also a binary `sigmond_query` used for reading the sigmond binary files.
Finally, python bindings to the sigmond classes and functions can be created using [pybind11](https://pybind11.readthedocs.io/en/stable/).


### Directory structure ###

- build - contains the makefiles for building the main programs, sigmond_query, and the python bindings.
- doc - contains documentation for how to use sigmond (not up to date).
- source - contains the source files.
- testing - contains various testing input files.
- ensembles.xml - contains definitions of various ensembles

### How do I build sigmond? ###

- `cd build/batch`
- edit the Makefile and adjust `INSTALL_DIR` as you see fit.
- `make`
- `make install` if you want the programs `sigmond` and `sigmond_query` copied to `INSTALL_DIR` (not required of course).

### How do I build python bindings? ###

- You will need `pybind11` installed
- `cd build/pysigmond`
- edit the Makefile and adjust `INSTALL_DIR`, `PYTHON`, and `PYTHON_CONFIG`
- `make`
- `make install`

You should now have a sigmond.<extension_suffix>.so file located in your local python site-packages directory.
It's not vital that it be there, but that it be somewhere that python can find it (e.g. in your local directory).

### Ensembles file ###

Sigmond needs information about the ensemble the data was run with.
This information is provided in the file `ensembles.xml`.
Unless an esembles file is specified in the input XML passed to sigmond,
the `ensembles.xml` file in the root directory will be used.

### XML input files ###

Once you have created an input XML file for sigmond (see the documenation for more details), you can run it using
  sigmond input.xml
