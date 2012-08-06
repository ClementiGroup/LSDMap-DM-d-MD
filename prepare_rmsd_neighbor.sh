#!/bin/bash

source run_parameters

eps_id=5

servername=`hostname|cut -d\. -f2`
echo Job will be set up on \" $servername \"

# write rmsd_neighbor.input
echo ${file_traj} > rmsd_neighbor.input
echo '1 ' ${npoints} >> rmsd_neighbor.input
echo ${neighbor_number} >> rmsd_neighbor.input

# write rmsd_neighbor.pbs
echo \#!/bin/bash > rmsd_neighbor.pbs
echo \#PBS -l nodes=${nodes}:ppn=${ppn},walltime=${walltime} >> rmsd_neighbor.pbs
echo \#PBS -q ${queue} >> rmsd_neighbor.pbs
echo \#PBS -N p_rmsd >> rmsd_neighbor.pbs
echo \#PBS -o err/${file_traj}_rmsd.out >> rmsd_neighbor.pbs
echo \#PBS -e err/${file_traj}_rmsd.err >> rmsd_neighbor.pbs
echo \#PBS -r n >> rmsd_neighbor.pbs
echo \#PBS -V >> rmsd_neighbor.pbs
echo cd \$PBS_O_WORKDIR >> rmsd_neighbor.pbs
echo mpiexec -n $(($nodes*$ppn)) p_rmsd_neighbor\<rmsd_neighbor.input\>rmsd_neighbor.log >> rmsd_neighbor.pbs
echo cat neighbor/${file_traj}_neighbor_9* \> neighbor/${file_traj}_neighbor >> rmsd_neighbor.pbs
echo rm neighbor/${file_traj}_neighbor_9* >> rmsd_neighbor.pbs


# write wlsdmap.input
echo $status_lsdmap $status_localscale > wlsdmap.input
echo ${npoints} >> wlsdmap.input
echo ${file_traj}_rmsd >> wlsdmap.input
echo ${file_weight} >> wlsdmap.input
echo ${file_traj}_dif >> wlsdmap.input
echo 0.0 >> wlsdmap.input
if [ $status_localscale -eq 1 ]; then
    echo $file_eps $eps_column_id >> wlsdmap.input
else
    echo ${constant_eps} >> wlsdmap.input
fi

# write wlsdmap.pbs
echo \#!/bin/bash > wlsdmap.pbs
echo \#PBS -l nodes=${nodes}:ppn=${ppn},walltime=${walltime} >> wlsdmap.pbs
echo \#PBS -q ${queue} >> wlsdmap.pbs
echo \#PBS -N p_wdmap >> wlsdmap.pbs
echo \#PBS -o err/${file_traj}_wlsdmap.out >> wlsdmap.pbs
echo \#PBS -e err/${file_traj}_wlsdmap.err >> wlsdmap.pbs
echo \#PBS -r n >> wlsdmap.pbs
echo \#PBS -V >> wlsdmap.pbs
echo cd \$PBS_O_WORKDIR >> wlsdmap.pbs
echo mpiexec -n $(($nodes*$ppn)) p_wlsdmap\<wlsdmap.input\>wlsdmap.log >> wlsdmap.pbs
