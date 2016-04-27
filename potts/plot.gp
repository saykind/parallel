#!/usr/bin/gnuplot -persist

H = .1

# termianl
set terminal postscript eps enhanced color font 'Helvetica, 20' 

# picture
set grid
set xrange [ 0 : 4 ]
set xlabel "T/J" 
set nokey

set output sprintf("E_plot_H=%.2f.eps",H)
set title "Energy"
set ylabel "E/J" 
set yrange [ -3 : -1.4 ]
set label 1 sprintf("H = %.2f", H) at 2.3, -2.35 
plot sprintf("E_data_H=%.2f.dat",H) u 1:2 w p pt 1 ps .5

set output sprintf("M_plot_H=%.2f.eps",H)
set title "Magnetic Moment"
set ylabel "M" 
set yrange [ * : * ]
set label 1 sprintf("H = %.2f", H) at 2.3, 1.45 
plot sprintf("M_data_H=%.2f.dat",H) u 1:($2-1) w p pt 1 ps .5

#    EOF
