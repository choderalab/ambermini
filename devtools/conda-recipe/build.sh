#!/bin/bash

# Required for cbio cluster for some stupid reason.
#export LD_LIBRARY_PATH="/opt/gnu/gcc/4.8.1/lib64:/opt/gnu/gcc/4.8.1/lib:/opt/gnu/gmp/lib:/opt/gnu/mpc/lib:/opt/gnu/mpfr/lib"

# See https://github.com/omnia-md/conda-recipes/pull/134
#     https://github.com/omnia-md/conda-recipes/issues/132
# This may require a patched version of gfortran to properly produce staticly linked binaries on osx
if [[ "$OSTYPE" == "darwin"* ]]; then
   export CUSTOMBUILDFLAGS="-static-libgfortran /usr/local/gfortran/lib/libquadmath.a -static-libgcc -lgfortran"
fi

export CFLAGS="-I$PREFIX/include $CFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS"

# Configure build
chmod u+x configure
./configure --prefix=$PREFIX

# Build and install.
make
chmod -R u+x bin
make install

