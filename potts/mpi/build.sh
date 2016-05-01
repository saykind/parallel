#!/bin/bash
mpicc -Wall -g -O3  main.c plot.c -lm -o a.out
mpirun -np 4 ./a.out
rm a.out
