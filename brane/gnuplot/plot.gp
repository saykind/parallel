#!/usr/bin/gnuplot -persist

# termianl
set terminal postscript eps enhanced color font 'Helvetica, 20' 

# picture
set grid
set xrange [ 2*pi/70 : pi ]
set logscale xy
set xlabel "q" 
set key box notitle right bottom spacing 1.2

set output sprintf("p=30.eps")
set title "Inverse Green (p8=30)"
set ylabel "G^{-1}"
set yrange [ * : * ]

list=system('ls -1B N=*')
plot x**4 dt '.' lw 3, for [file in list] file w lp lw 2  t file, x**3.2 dt '-' lw 4 lt rgb "red"

#    EOF
