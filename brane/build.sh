#!/bin/bash

gcc -std=gnu11 -Wall -g -fopenmp main.c storage.c simulate.c -lm -o a.out
gcc -std=gnu11 -Wall -g -fopenmp read.c -lm -o x.out

