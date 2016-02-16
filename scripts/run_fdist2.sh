#!/bin/bash
#to compile fdist2: 
#gcc -o fdist2 fdist2.c -lm
#Code to run fdist2 a bunch of times

cd ~/sf_ubuntushare/fdist2
params_file="many_params2.txt"
write_file="fdist_params2.dat"
orig_out="out.dat"

#read file
exec 3<&0
exec 0< $params_file
while read -r rnum ndem	nsub efst nsmp modc nrep
do
	#echo $ndem $nsub $efst $nsmp $modc $nrep
	#write a new fdist_params2.dat file
	echo "$ndem" > $write_file
	echo "$nsub" >> $write_file
	echo "$efst" >> $write_file
	echo "$nsmp" >> $write_file
	echo "$modc" >> $write_file
	echo "$nrep" >> $write_file
	#run fdist2
	echo "Running fdist2, parameter set $rnum"
	./fdist2 
	new_out="out.dat$rnum.txt"
	mv $orig_out $new_out
done
