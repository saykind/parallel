#!/bin/bash

mpicc -Wall -g main.c -lm -o a.out
for n in 1 2 3 4;
do
	echo 'Number of processes = '$n
	mpirun -np $n ./a.out $1 > data.dat
	cat data.dat
done
rm a.out
