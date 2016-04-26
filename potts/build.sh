#!/bin/bash

gcc -Wall -g -fopenmp  main.c plot.c -lm -o a.out
time ./a.out 
rm a.out
