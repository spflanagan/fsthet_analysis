#!/bin/bash
OUTDIR="../../results/numerical_analysis"
cd ../programs/numerical_analysis
./numerical_analysis -o ${OUTDIR}/Nm0.1.d2.s2.fig. -n 0.1 -d 2 -s 2 -r 0
./numerical_analysis -o ${OUTDIR}/Nm1.d2.s2.fig. -n 1 -d 2 -s 2 -r 0
./numerical_analysis -o ${OUTDIR}/Nm10.d2.s2.fig. -n 10 -d 2 -s 2 -r 0
./numerical_analysis -o ${OUTDIR}/Nm0.1.d5.s2.fig. -n 0.1 -d 5 -s 5 -r 0
./numerical_analysis -o ${OUTDIR}/Nm1.d5.s2.fig. -n 1 -d 5 -s 5 -r 0
./numerical_analysis -o ${OUTDIR}/Nm10.d5.s2.fig. -n 10 -d 5 -s 5 -r 0

./numerical_analysis -o ${OUTDIR}/Nm1.d2.s2.fig. -n 1 -d 2 -s 2 -r 0
./numerical_analysis -o ${OUTDIR}/Nm1.d2.s4.fig. -n 1 -d 2 -s 4 -r 0
./numerical_analysis -o ${OUTDIR}/Nm1.d2.s8.fig. -n 1 -d 2 -s 8 -r 0
./numerical_analysis -o ${OUTDIR}/Nm1.d5.s5.fig. -n 1 -d 5 -s 5 -r 0
./numerical_analysis -o ${OUTDIR}/Nm1.d5.s10.fig. -n 1 -d 5 -s 10 -r 0
./numerical_analysis -o ${OUTDIR}/Nm1.d5.s15.fig. -n 1 -d 5 -s 15 -r 0