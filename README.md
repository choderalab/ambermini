ambermini
=========

A stripped-down set of just antechamber, sqm, and tleap. This package exists
specifically to provide a lightweight subset of the AmberTools installation so
that they can be installed as dependencies of other programs.

Unless you know what you are doing (and do not need support from the Amber
community), you should download and use the full AmberTools suite available at
http://ambermd.org.

Disclaimer
==========

While this package is derived from the AmberTools source code (see
http://ambermd.org), it is *NOT* supported or maintained by the Amber community.

If you need support or plan to use these tools directly, you are recommended
to use the official release of AmberTools.


For Developers
==============

To install
----------

```
./configure --prefix <destination>
make
make install
```


To uninstall
------------

`make uninstall`


Prerequisites
-------------
byacc, yacc, or bison

flex

Fortran and C compiler

Python
