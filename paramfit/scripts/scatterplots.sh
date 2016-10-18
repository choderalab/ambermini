#!/bin/bash

# First the bonds
plotcmd="plot "
num=0
if test -n "$(find . -maxdepth 1 -name '*bondeq' -print -quit)"; then
for i in `ls *bondeq`; do
  title=$(head -n 1 $i)
  title=${title//#}
  title=${title//[[:space:]]}
  plotcmd="$plotcmd '$i' title '$title' with points pt 7 ps 1 lc $num,"
  (( num++ ))
done
plotcmd=${plotcmd%,}


gnuplot -persist << EOF
set key autotitle columnhead
set title "Bond Equilibrium Lengths in Sampled Conformations"
set xlabel "Equilibrium Length (Angstrom)"
set xrange [0:]
unset ytics
unset ylabel
$plotcmd

EOF
fi

if test -n "$(find . -maxdepth 1 -name '*angleq' -print -quit)"; then
# Do the angles
plotcmd="plot "
num=0
for i in `ls *angleq`; do
  title=$(head -n 1 $i)
  title=${title//#}
  title=${title//[[:space:]]}
  plotcmd="$plotcmd '$i' title '$title' with points pt 7 ps 1 lc $num,"
  (( num++ ))
done
plotcmd=${plotcmd%,}

gnuplot -persist << EOF
set key autotitle columnhead
set title "Angle Equilibrium Values in Sampled Conformations"
set xlabel "Equilibrium Phase (radians)"
unset ytics
unset ylabel
set xrange [0:3.14]
$plotcmd

EOF
fi 

if test -n "$(find . -maxdepth 1 -name '*diheq' -print -quit)"; then
# Now the dihedrals
plotcmd="plot "
num=0
for i in `ls *diheq`; do
  title=$(head -n 1 $i)
  title=${title//#}
  title=${title//[[:space:]]}
  plotcmd="$plotcmd '$i' title '$title' with points pt 7 ps 1 lc $num,"
  (( num++ ))
done
plotcmd=${plotcmd%,}

gnuplot -persist << EOF
set key autotitle columnhead
set title "Dihedral Equilibrium Values in Sampled Conformations"
set xlabel "Equilibrium Phase (radians)"
set xrange [0:1.57]
unset ytics
unset ylabel
$plotcmd

EOF
fi
