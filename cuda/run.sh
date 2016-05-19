#!/bin/bash
echo 'Compiling...'
nvcc -O3 -arch=sm_21 main.cu -o a.out
echo 'Running...'
time ./a.out
rm a.out
