#!/bin/bash

gcc -Wall -g -fopenmp -O3  main.c plot.c -lm -o a.out
time ./a.out
rm a.out
