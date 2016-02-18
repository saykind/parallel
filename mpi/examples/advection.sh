#!/bin/bash

mpicc -Wall -g linear.c -lm -o a.out
#for n in {1..18}
for n in 1 2;
do
	mpirun -np $n ./a.out $1 >> data.dat
done
rm a.out
