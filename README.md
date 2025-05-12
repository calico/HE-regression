# HE_regression_private
These functions are used in the upcoming manuscript on the application of Haseman-Elston regression to estimate heritability (h<sup>2</sup>) and genetic correlations (r<sub>g</sub>) among traits measured in Diversity Outbred (DO) mouse populations. See [**Test_HE_Reg.pdf**](https://github.com/calico/HE_regression_private/blob/main/Test_HE_Reg.pdf) for a knit of the .Rmd file under /code.

## Introduction
Haseman-Elston regression is a method for estimating the heritability (h<sup>2</sup>) of traits or genetic correlations (r<sub>g</sub>) among pairs of traits. At a high level, this involves regressing phenotypic similarity (defined here as the product of z-normalized trait values) on kinship among pairs of individuals. The regression coefficient (β) can be used to estimate the covariance of traits as a function of kinship (⍴):

⍴ = β / (2*σ<sub>1</sub>*σ<sub>2</sub>)

where σ<sub>1</sub> is the standard deviation of the first trait and σ<sub>2</sub> is the standard deviation of the second. For a single trait, ⍴ = h<sup>2</sup>. For a pair of traits ⍴ can be used to estimate r<sub>g</sub> via a structural equation model if the heritabilities of the pair of traits are known. 

Here, we have implemented Haseman-Elston regression for use with data derived from studies in Diversity Outbred (DO) mice. However, these functions can, in theory, be applied to any set of z-score normalized phenotype data with an accompanying kinship matrix of the sort generated via R/qtl2.

## Core functionality
Functions include:
1) z-score normalization of vectors of trait values.
2) generation of combinatorial pairwise products from pairs of z-normalized trait vectors.
3) estimation of h<sup>2</sup> or r <sub>g</sub>.
4) adjustments to the standard error of the estimates as a function of trait h<sup>2</sup> and sample size.

## Code
The main [**HER_functions.R**](https://github.com/calico/HE_regression_private/blob/main/code/HER_functions.R) file is a collection of functions that can be imported in R to estimate h2 and rg for sets of z score normalized phenotypes. This repo also contains a markdown file [**Test_HE_Reg.Rmd**](https://github.com/calico/HE_regression_private/blob/main/code/Test_HE_Reg.Rmd) that serves as documentation and demonstrates the utility of each of the functions. 

## Data
Contains a set of phenotype data [**DO_Phenotypes_non_molecular.RData**](https://github.com/calico/HE_regression_private/blob/main/data/DO_Phenotypes_non_molecular.RData) (corresponding to the 240507 version of the non-molecular data) and corresponding kinship information [**kinship_loco.RData**](https://github.com/calico/HE_regression_private/blob/main/data/kinship_loco.RData) that can be used to test the functions. Functions are loaded via:
`load('/path/to/HER_functions.R')`

## Dependencies
- dplyr
- parallel

For the .Rmd specifically:
- scales
