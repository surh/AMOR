# Abundance Matrix Operations in R

[![DOI](https://zenodo.org/badge/10928/surh/AMOR.svg)](https://zenodo.org/badge/latestdoi/10928/surh/AMOR) [![Build Status](https://travis-ci.org/surh/AMOR.svg?branch=master)](https://travis-ci.org/surh/AMOR)

This package is a collection of functions that are useful for performing common tasks with abundance tables.

The original code was written mostly to deal with abundance of bacterial marker genes in metagenomic data, but
the functions can be useful for many other types of data where a table of samples by features is accompanied by
metadata on both the samples and features.

A few statistical analysis utilities are present, like GLMs and ordination methods, but still are a work in progress.

# Installation

The easiest way to install is with the devtools package. Make sure you have it installed and use

```r
devtools::install_github("surh/AMOR")
```

If you want to get the development version, you can use

```r
devtools::install_github("surh/AMOR",ref = "dev")
```

Status of the develpment version can be seen below:

[![Build Status](https://travis-ci.org/surh/AMOR.svg?branch=dev)](https://travis-ci.org/surh/AMOR)
