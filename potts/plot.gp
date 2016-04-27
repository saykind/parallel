#!/usr/bin/gnuplot -persist

H = .0

# termianl
set terminal postscript eps enhanced color font 'Helvetica, 20' 

# picture
set grid
set xrange [ 0 : 2 ]
set xlabel "T/J" 
set nokey

set output sprintf("E_plot_H=%.2f.eps",H)
set title "Energy"
set ylabel "E/J" 
set yrange [ -4.1 : -2.9 ]
set label 1 sprintf("H = %.2f", H) at .2, -3.55 
plot sprintf("E_data_H=%.2f.dat",H) u 1:2 w p pt 1 ps .5

set output sprintf("M_plot_H=%.2f.eps",H)
set title "Magnetic Moment"
set ylabel "M" 
set yrange [ -.1 : 1.1 ]
set label 1 sprintf("H = %.2f", H) at .2, .25
plot sprintf("M_data_H=%.2f.dat",H) u 1:($2-1) w p pt 1 ps .5

#    EOF
