#!/bin/bash

# input file
fileIn="$1"

# update path if necessary:
gnuplot=/Applications/Maxima/Gnuplot.app/Contents/Resources/bin/gnuplot

# generate plot
$gnuplot -persist <<PLOT
unset key
set terminal png
set output '${fileIn%.*}.png'
set terminal png size 1000, 700
set ylabel "Time [s]"
set xlabel "$fileIn, Date"
plot '$fileIn' using 1:2 with lines
quit
PLOT

