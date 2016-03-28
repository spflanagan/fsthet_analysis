#!/bin/bash
OUTDIR="../../results/numerical_analysis_testingexp"
cd ../programs/numerical_analysis

./numerical_analysis -o ${OUTDIR}/Nm0.1.d100.s100. -n 0.1 -d 100 -s 20 -r 1
./numerical_analysis -o ${OUTDIR}/Nm1.d100.s100. -n 1 -d 100 -s 20 -r 1
./numerical_analysis -o ${OUTDIR}/Nm10.d100.s100. -n 10 -d 100 -s 20 -r 1

