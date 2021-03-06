# BLUES - Bond Locator Utilizing Electronic Structure
[![Build Status](https://travis-ci.org/quantum2classical/blues.svg?branch=master)](https://travis-ci.org/quantum2classical/blues)
[![Coverage Status](https://coveralls.io/repos/github/quantum2classical/blues/badge.svg?branch=master)](https://coveralls.io/github/quantum2classical/blues?branch=master)
[![Documentation Status](https://readthedocs.org/projects/blues/badge/?version=latest)](https://blues.readthedocs.io/en/latest/?badge=latest)

A pipeline written in Python that converts NWChem output files into InChI and SMILES strings by estimating bond orders with the help of [JANPA](http://janpa.sourceforge.net/), an open source natural population analysis program.

BLUES is the result of a collaborative capstone project between students at the [University of Washington](https://www.washington.edu) and staff scientists at the [Pacific Northwest National Lab](https://www.pnnl.gov/) as a part of the [DIRECT (Data Intensive Research Enabling Clean Technologies) Training Program](http://depts.washington.edu/uwdirect/) at the [University of Washington](https://www.washington.edu).

Documentation can be found [here](https://blues.readthedocs.io/en/latest/index.html)

## Dependencies
- rdkit
- numpy

## Installation
```
$ git clone https://github.com/quantum2classical/blues.git
$ cd blues
$ conda env create
$ source activate blues
$ python setup.py install
```

## Usage
```python
import blues

# Pass any proper Molden file to the converter object when instantiated
converter = blues.MoldenConverter("molecule.molden")

## Convert to SMILES or InChI identifier strings
# First call runs JANPA to get bond-order matrix and build molecule
smiles = converter.tosmiles()

# Second call is much faster and only needs to convert molecule to string
inchi = converter.toinchi()
```
