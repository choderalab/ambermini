`ambermini` is an internal dependency for [omnia](http://omnia.md), and contains a stripped-down set of just antechamber, sqm, and tleap.

**NOTE: This is not an officially-supported AmberTools distribution.**

To obtain AmberTools, go to [the Amber webpage](http://ambermd.org).

This version has been updated to match AmberTools 15 as of 23 Jun 2015.

### To install
```
./configure --prefix <destination>
make
make install
```

### To uninstall
```
make uninstall
```

### Prerequisites
* byacc or yacc
* flex
* Fortran and C compiler
* Python
