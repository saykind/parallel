#!/bin/bash

mpicc -Wall -g -O3 main.c plot.c -lm -o a.out
for n in 1 4;
do
	echo 'Number of processes = '$n
	mpirun -np $n ./a.out $1
done
rm a.out
