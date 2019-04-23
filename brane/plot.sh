#!/bin/bash
for n in 20 25 30; do 
    ./x.out "N=$n";
done
cd gnuplot
gnuplot plot.gp
cd ..
