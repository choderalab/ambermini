To obtain AmberTools, go to [the Amber webpage](http://ambermd.org).

**NOTE: This is not an officially-supported AmberTools distribution.**

`ambermini` is an internal dependency for [omnia](http://omnia.md), and contains a stripped-down set of just antechamber, sqm, and tleap.

To install
==========

./configure --prefix <destination>
make
make install


To uninstall
============

make uninstall


Prerequisites
=============
byacc or yacc
flex
Fortran and C compiler
Python
