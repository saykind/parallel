#!/usr/bin/gnuplot -persist

H = .0

# termianl
set terminal postscript eps enhanced color font 'Helvetica, 20' 

# picture
set grid
set xrange [ 0 : 2 ]
set xtics .2
set xlabel "T/J" 
set key box notitle left center spacing 1.2

set output sprintf("energy_plot_H=%.2f.eps", H)
set title "Energy in Potts model (Z=4)"
set ylabel "E/ZJ"
set yrange [ -1.025 : -.625 ]
plot 	sprintf("E_data_q=%d_H=%.2f.dat", 2, H) w line lt 2 lw 2 t "q=2",\
	sprintf("E_data_q=%d_H=%.2f.dat", 3, H) w line lt 3 lw 2 t "q=3",\
	sprintf("E_data_q=%d_H=%.2f.dat", 4, H) w line lt 4 lw 2 t "q=4",\
	sprintf("E_data_q=%d_H=%.2f.dat", 5, H) w line lt 5 lw 2 t "q=5"

set output sprintf("moment_plot_H=%.2f.eps", H)
set title "Magnetic Moment in Potts model (Z=4)"
set ylabel "M_{normalized}" 
set yrange [ -.1 : 1.1 ]
plot 	sprintf("M_data_q=%d_H=%.2f.dat", 2, H) w line lt 2 lw 2 t "q=2",\
	sprintf("M_data_q=%d_H=%.2f.dat", 3, H) w line lt 3 lw 2 t "q=3",\
	sprintf("M_data_q=%d_H=%.2f.dat", 4, H) w line lt 4 lw 2 t "q=4",\
	sprintf("M_data_q=%d_H=%.2f.dat", 5, H) w line lt 5 lw 2 t "q=5"
#plot 	sprintf("M_data_q=%d_H=%.2f.dat", 2, H) w p pt 1 ps .5 t "q=2",\
#	sprintf("M_data_q=%d_H=%.2f.dat", 3, H) w p pt 1 ps .5 t "q=3",\
#	sprintf("M_data_q=%d_H=%.2f.dat", 4, H) w p pt 1 ps .5 t "q=4",\
#	sprintf("M_data_q=%d_H=%.2f.dat", 5, H) w p pt 1 ps .5 t "q=5"

#    EOF
