#!/bin/bash

gcc -Wall -g -fopenmp -O3  main.c plot.c -lm -o a.out
./a.out
rm a.out
