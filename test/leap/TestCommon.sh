# LEaP tests can source this for common functionality.

CLEAN=0

# Set test tleap
if [ -z "$TESTtleap" ] ; then
  if [ -z "$AMBERHOME" ] ; then
    TESTtleap=../../../../bin/tleap
  else
    TESTtleap=$AMBERHOME/bin/tleap
  fi
fi

# Set sander for calculating single point energies
if [ -z "$TESTsander" ] ; then
  if [ -z "$AMBERHOME" ] ; then
    TESTsander=../../../../bin/sander
  else
    TESTsander=$AMBERHOME/bin/sander
  fi
fi

# Set location of dacdif relative to subdirs
DACDIF=../../dacdif

# Clean all files if present
CleanFiles() {
  while [ ! -z "$1" ] ; do
    if [ -f "$1" ] ; then
      rm $1
    fi
    shift
  done
  if [ "$CLEAN" = 1 ] ; then
    exit 0
  fi
}

# RunTleap <expect topology name>
# Run tleap for leap.in. Check for errors.
RunTleap() {
  $TESTtleap -f leap.in > leap.out
  if [ $? -ne 0 -o ! -f "$1" -o ! -s "$1" ] ; then
    echo "Program error"
    exit 1
  fi
}

if [ "$1" = 'clean' ] ; then
  CLEAN=1
  shift
fi
