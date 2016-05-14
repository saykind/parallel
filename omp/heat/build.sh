#!/bin/bash

gcc -Wall -g -fopenmp main.c plot.c -lm -o a.out
for n in {1..4};
do
	./a.out $n
done
rm a.out
