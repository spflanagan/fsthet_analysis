# fsthet_analysis

This repository contains the scripts and programs used in Flanagan & Jones (2017), "Constraints on the Fst-heterozygosity outlier approach".

This document outlines the contents of each directory in the fst_outliers repository.

## fsthet
fhsthet is an R package that calculates smoothed quantiles from empirical Fst-heterozygosity distributions to identify loci with extreme Fst values. The directory is the entire package and contains the R code, documentation files, and a vignette. The fsthet-manual.pdf is its manual and it is compiled in *.zip and *.tar.gz files in the fsthet/ directory.

You can install it using devtools::install_github("spflanagan/fsthet_analysis/fsthet")

## programs/fdist2
This is the code obtained from http://www.maths.bris.ac.uk/~mamab/software/. It has been modified slightly so that it could be run en masse at the command line using the run_fdist2.sh script (see below).

## programs/numerical_analysis
This directory contains the C++ program used to conduct a numerical analysis of migration and drift under Wright's island model. It contains a Visual studios project (numerical_analyiss.vcxproj) as well as the source code (numerical_analysis.cpp and random_numbers.h) and the compiled program (for Linux: numerical_analysis, for Windows: Debug/numerical_analysis.exe).
Within this directory are zipped folders containing data output from the various parameter combinations.

## scripts
The scripts/ directory contains several scripts that were used to generate the data for the paper. 
- bootstrap_on_nullnum.R runs and analyzes the fsthet analysis on the numerical analysis results without selection.
- bootstrap_on_selnum.R runs and analyzes the fsthet analysis on the numerical analysis results WITH selection.
- plot_fdist.R plots the results of running the fdist2 simulator.
- run_fdist2.sh shows how we ran the various parameter sets for fdist2.
- run_numerical_analysis.sh, run_numerical_analysis_exp.sh, and run_numerical_analysis_forfigs.sh all are various parameter set runs of the numerical analysis program (in programs/numerical_analysis)/
- testing_expectations.R is the script used to show that our numerical analysis results adhere to the expectations of the island model with drift and migration.
- testing_fst_calcs.R is a script comparing the different Fst calculation methods and is where the code for Fig. S2 is located.

