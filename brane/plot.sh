#!/bin/bash
for n in 10 15 20 25; do 
    ./x.out "N=$n";
done
cd gnuplot
gnuplot plot.gp
cd ..