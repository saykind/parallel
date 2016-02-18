#!/bin/bash

gcc -Wall -g -fopenmp linear.c -lm -o a.out
for n in {1..4};
do
	./a.out $n >> data4.dat
done
rm a.out
