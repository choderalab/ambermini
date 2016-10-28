#!/bin/sh
# Gnuplot script for showing structure energies 

if test -z "$1"  
  then echo "Usage:"
       echo "  plot_energy.x filename"
       exit 1
fi

if ! hash gnuplot
  then echo "This script uses gnuplot, which is not installed."
       echo "Please install gnuplot or use your own plotting program."
       exit 1
fi

gnuplot -persist << EOF

set style line 1 pt 19 linecolor rgb "red"
set style line 2 pt 17 linecolor rgb "blue"
set style line 3 pt 17 linecolor rgb "black" lw 2
set pointsize 0.2

set ylabel "Energy"
set xlabel "Structure"
set title "$1"

plot "$1" using 1:2 title 'Fit Amber' with linespoints linestyle 1, \
"$1" using 1:4 title 'Initial Amber' with linespoints linestyle 2, \
"$1" using 1:3 title 'Quantum' with linespoints linestyle 3

#set term wxt 
#replot

EOF

