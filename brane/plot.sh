#!/bin/bash
for n in 30 35 40 45 50 55 65; do 
    ./x.out "N=$n";
done
cd gnuplot
gnuplot plot.gp
cd ..
