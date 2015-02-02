#!/bin/bash

# this script was automatically created when using DMap Sampling

startgro=input.gro
tmpstartgro=tmp.gro

natoms=$(sed -n '2p' $startgro)
nlines_per_frame=$((natoms+3))

nlines=`wc -l $startgro| cut -d' ' -f1`
nframes=$((nlines/nlines_per_frame))

for idx in `seq 1 $nframes`; do

  start=$(($nlines_per_frame*(idx-1)+1))
  end=$(($nlines_per_frame*idx))
  sed "$start"','"$end"'!d' $startgro > $tmpstartgro

  # gromacs preprocessing & MD
  grompp_jbm -f grompp.mdp -c $tmpstartgro -p topol.top  &> /dev/null
  mdrun_jbm -nt 1 -dms config.ini -s topol.tpr  &> mdrun.log

done

# remove temporary files
rm -f $tmpstartgro
        