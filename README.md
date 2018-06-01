# hlmm
hlmm is a python library for fitting heteroskedastic linear mixed models to genetic data. 

A heteroskedastic linear model (HLM) can model the effect
of a set of variables on the mean of a response (such as a continuous phenotype) and the 
effect of a (potentially different) set of variables of the variability of the response

A HLM models all effects on both the mean and variance of a response as fixed effects.
A heteroskedastic linear mixed model (HLMM) adds modelling of random effects of a set of variables on the mean of the response. 

Modelling random effects can make fitting a HLMM very computationally demanding,
with the number of operations scaling with the cube of sample size. 

The package hlmm provides the ability to fit a HLMM
much more quickly when the number of variables with random effects is small compared 
to the sample size. We use an algorithm whose operations scale in proportion to the sample
size multiplied by the number of random effects squared. 

# Main features:

HLM class: given data, find the parameters that maximise
the likelihood and their sampling distribution. 

HLMM class: given data, find the parameters that maximise
the likelihood and their sampling distribution. 

hlmm_chr.py: command line script that fits heteroskedastic linear models or 
heteroskedastic linear mixed models to a contiguous segment of the genome.
The script takes bed formatted genotypes as input. and can incorporate
covariates for the fixed effects on the mean and/or variance. 

# Easy installation

We recommend installing using pip (https://pip.pypa.io/en/stable/). 
At the command line, type

    'sudo pip install hlmm'
    
# Documentation

We recommend reading the documentation for this package (http://hlmm.readthedocs.io/en/master/) and 
working through the tutorial (http://hlmm.readthedocs.io/en/master/tutorial.html).

# Detailed Package Install Instructions

hlmm has the following dependencies:

python 2.7

Packages: 

- numpy
- scipy
- pysnptools

We highly recommend using a python distribution such as Anaconda (https://store.continuum.io/cshop/anaconda/). 
This will come with both numpy and scipy installed and can include an MKL-compiled distribution
for optimal speed. 

To install from source, clone the git repository, and in the directory
containing the HLMM source code, at the shell type

    'sudo python setupy.py install'

or, on the windows command prompt, type

    'python setup.py install' 
    
# Running tests

The tests directory contains scripts for testing the computation of 
the likelihoods, gradients, and maximum likelihood solutions for
both heteroskedastic linear models (test_hetlm.py) and
for heteroskedastic linear mixed models (test_hetlmm.py).
To run these tests, a further dependency is required: numdifftools. 

To run the tests, first install hlmm. Change to the tests/ directory and at the shell type

    'python test_hetlm.py'
    'python test_hetlmm.py'

Both tests should run without any failures. 