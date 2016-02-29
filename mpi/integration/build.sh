#!/bin/bash

mpicc -Wall -g main.c integration.c plot.c -lm -o a.out
for n in  3;
do
	mpirun -np $n ./a.out $1 >> data.dat
done
rm a.out
