#!/bin/bash

#
# Usage: ./gro2xyz <name of .gro file> <stride>
#

if [ $# -ne 2 ]; then
    echo 'Usage: ./gro2xyz <name of .gro file> <stride>'
    exit
fi

grofile=$1
stride=$2
natom=`head -2 $grofile|tail -1`
echo Number of atoms: $natom
nline=$(($natom+3))
npoint=`grep 't=' $grofile|wc -l`
echo Number of points: $npoint
nnonhydrogen=`head -$(($nline-1)) $grofile| awk '{if( NR>2 && $2!~/^H/) {print NR}}'|wc -l`
echo Number of non-hydrogen atoms: $nnonhydrogen
dim=$(($nnonhydrogen*3))
echo Number of dimensions: $dim
filename=`echo $grofile|sed 's/.gro/.xyz/'`
echo Write output to file: $filename
npoint_stride=$(($npoint/$stride))
if [ $(($npoint%$stride)) -gt 0 ]; then
    npoint_stride=$(($npoint_stride+1))
fi
echo Write every $stride steps
echo Number of points to write: $npoint_stride
echo $npoint_stride $dim > $filename
awk '{if( NR%'$nline'>2 && int((NR-1)/'$nline')%'$stride'==0 && $2!~/^H/){print $4,$5,$6}}' $grofile >> $filename

