#!/bin/bash

# Process parameters: process_gaussian input > output filename
USAGE="Usage: `basename $0` <directory with input files> \n \n"

#Parse command line options
if [ $# != 1 ]; then
	echo -e $USAGE
	exit 1
fi

if [ ! -d $1 ]; then
        echo -e " Directory ${1} not found! \n\n"
	echo -e $USAGE
	exit 1
fi


for i in $(ls $1); do
        energy=$(grep "SCF Done" "$1/$i" | awk '{print $5}')
	echo $energy
done
