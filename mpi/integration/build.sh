#!/bin/bash

mpicc -Wall -g main.c integration.c plot.c -lm -o a.out
for n in 1 2 3 4 5 6 7 8 9 10;
do
	mpirun -host tor2 -np $n ./a.out $1 > data.dat
	echo 'Number of processes = '$n
	cat data.dat
done
rm a.out
