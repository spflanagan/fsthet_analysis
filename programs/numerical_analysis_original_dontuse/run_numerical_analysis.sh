#!/bin/bash

file="parameters.txt"
time=0
while read -r Npop Nsub mRate
do
	if [ $time > 0 ];  then
		filename="n${Npop}_s${Nsub}_m${mRate}"
		echo "running numerical_analysis on $filename"
		./numerical_analysis -n ${Npop} -s ${Nsub} -m mRate -r 10000 -o ${filename}
	fi
	time=$(($time+1))
done < $file
