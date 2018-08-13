# RDR
Here we provide scripts for estimating heritability using some of the methods in 'Relatedness Disequilibrium Regression estimates heritability without environmental bias'
Nature Genetics (2018).

# RELT

This script performs relatedness thresholded (RELT) heritability estimation. It regresses elements of the sample
phenotypic covariance matrix onto corresponding elements of a relatedness matrix for those pairs with relatedness
below a user-set threshold (default 0.05).

The script calculates standard errors of genetic variance and heritability estimates by a procedure that takes into
account dependence between pairs. This can be computationally demanding for large sample sizes.

It takes a binary relatedness matrix as input. The matrix is formatted the same way as produced by
GCTA with the --make-grm-bin option set. The matrix is in lower-triangular order with 32 bit floating point numbers.

Associated with the relatedness matrix is an plain text id file. The first column of the id file gives the ids of the
individuals in the order that they appear in the relatedness matrix. If IDs are specified uniquely by the first column
of the GCTA grm.id file, then the GCTA grm.id file can be used.

The trait file is a plain text file with columns: FID, IID, trait1, trait2, etc.

If covariates are supplied, the trait will first be adjusted for covariates before heritability is estimated.

The script outputs a file outprefix.herit that records variance component estimates and standard errors.

Given a trait file y.txt and the output of GCTA --make-grm-bin as R.grm.bin and R.grm.id, an example usage would be

    'python RELT.py R.grm.id R.grm.bin y.txt y'

# Sib-Regression

This script performs Sib-Regression heritability estimation. It performs simple univariate regression of
the squared difference of siblings' phenotype observations onto their genetic relatedness.