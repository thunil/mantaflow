#!/bin/bash

# input file
fileIn="$1"

# update path if necessary:
gnuplot=gnuplot
#gnuplot=/usr/bin/gnuplot
#e.g. for Mac: /Applications/Maxima/Gnuplot.app/Contents/Resources/bin/gnuplot

# check for gnuplot override
[ -z "${MANTA_GNUPLOT+x}" ] && echo "MANTA_GNUPLOT is not set" || echo "MANTA_GNUPLOT found and is set to \"$MANTA_GNUPLOT\"."
[ -z "${MANTA_GNUPLOT+x}" ] && echo "MANTA_GNUPLOT is not set" || gnuplot=$MANTA_GNUPLOT

echo "Using \"$gnuplot\"." >> tmp.txt

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

echo "Plot generated"
