#!/usr/bin/gnuplot -persist

# termianl
set terminal pdf enhanced color font 'Helvetica, 20' size 21 cm, 14.8 cm

# picture
set grid
set xrange [ pi/50 : pi ]
set logscale xy
set xlabel "q" 
set key box notitle right bottom spacing 1.2

set output sprintf("p=1.pdf")
set title "Inverse Green (p8=1)"
set ylabel "G^{-1}"
set yrange [ * : * ]


list=system('ls -1B N=*')
plot x**4 dt '.' lw 3, for [file in list] file w lp lw 2  t file, 1.2*x**3.2 dt '-' lw 4 lt rgb "red"

#    EOF
