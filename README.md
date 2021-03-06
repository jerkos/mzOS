mzOS
====

[![Build Status](https://travis-ci.org/jerkos/mzOS.svg?branch=master)](https://travis-ci.org/jerkos/mzOS)
[![Coverage Status](https://coveralls.io/repos/jerkos/mzOS/badge.svg?branch=master&service=github)](https://coveralls.io/github/jerkos/mzOS?branch=master)
[![Code Issues](https://www.quantifiedcode.com/api/v1/project/f458ff06e0bd4eb8b2b98cee0dab6ecb/badge.svg)](https://www.quantifiedcode.com/app/project/f458ff06e0bd4eb8b2b98cee0dab6ecb)

## Small library for features annotations in  LC-MS metabolomics experiments:

* [x] isotopes detection
* [x] adducts/fragments annotation using simple heuristics (when conflicts)
* [x] database search (HMDB, LMSD) from local sqlite file (fast)
* [x] database assignment scoring based on predicted isotopic distribution and bayesian algorithm. 

Some documentation is available [here](http://jerkos.github.io/mzOS) (to be completed).

## Dependancies

* numpy
* scipy
* scikit-learn
* pandas

## Roadmap

* support more databases
* export results as a graph