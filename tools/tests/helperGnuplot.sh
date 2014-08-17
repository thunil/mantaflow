#!/bin/bash

fileIn="$1"
#echo AA $fileIn $2 
gnuplot=/Applications/Maxima/Gnuplot.app/Contents/Resources/bin/gnuplot

#gnuplot -persist <<PLOT
# set output 'runtimes/${fileIn%%.*}.png'

$gnuplot -persist <<PLOT
unset key
set terminal png
set output 'runtimes/$fileIn.png'
set terminal png size 1000, 700
set ylabel "Time [s]"
set xlabel "$fileIn, Date"
plot '$fileIn' using 1:2 with lines
quit
PLOT

