#!/bin/bash

gcc -Wall -g -fopenmp main.c plot.c -lm -o a.out
for n in {1..4};
do
	./a.out $n >> data.dat
done
cat data.dat
rm a.out
rm data.dat
