#!/bin/bash

source run_parameters

servername=`hostname|cut -d\. -f2`
echo Job will be set up on \" $servername \"

# write wlsdmap_embed.input
echo $status_lsdmap $status_localscale > wlsdmap_embed.input
echo ${npoints} >> wlsdmap_embed.input
echo ${file_traj}_embed_rmsd >> wlsdmap_embed.input
echo ${file_weight} >> wlsdmap_embed.input
echo ${file_traj}_dif >> wlsdmap_embed.input
echo 0.0 >> wlsdmap_embed.input
if [ $status_localscale -eq 1 ]; then
    echo ${file_traj}_eps $eps_column_id >> wlsdmap_embed.input
else
    echo ${constant_eps} >> wlsdmap_embed.input
fi
echo ${file_traj} >> wlsdmap_embed.input
echo ${file_add} >> wlsdmap_embed.input
echo ${nadd} >> wlsdmap_embed.input
echo ${file_weight_add} >> wlsdmap_embed.input

# write wlsdmap_embed.pbs
echo \#!/bin/bash > wlsdmap_embed.pbs
echo \#PBS -l nodes=${nodes}:ppn=${ppn},walltime=${walltime} >> wlsdmap_embed.pbs
echo \#PBS -q ${queue} >> wlsdmap_embed.pbs
echo \#PBS -N p_wdmap >> wlsdmap_embed.pbs
echo \#PBS -o err/${file_traj}_wlsdmap_embed.out >> wlsdmap_embed.pbs
echo \#PBS -e err/${file_traj}_wlsdmap_embed.err >> wlsdmap_embed.pbs
echo \#PBS -r n >> wlsdmap_embed.pbs
echo \#PBS -V >> wlsdmap_embed.pbs
echo cd \$PBS_O_WORKDIR >> wlsdmap_embed.pbs
echo mpiexec -n $(($nodes*$ppn)) p_wlsdmap_embed\<wlsdmap_embed.input\>wlsdmap_embed.log >> wlsdmap_embed.pbs

# write split_rmsd.input
echo ${file_traj}_rmsd > split_rmsd.input
echo ${npoints} >> split_rmsd.input
echo ${nadd} >> split_rmsd.input
echo $(($nodes*$ppn)) >> split_rmsd.input
echo ${file_traj}_embed_rmsd >> split_rmsd.input

# write split_rmsd.pbs
echo \#!/bin/bash > split_rmsd.pbs
echo \#PBS -l nodes=1:ppn=1,walltime=${walltime} >> split_rmsd.pbs
echo \#PBS -q ${queue} >> split_rmsd.pbs
echo \#PBS -N p_rmsd >> split_rmsd.pbs
echo \#PBS -o err/${file_traj}_embed_rmsd.out >> split_rmsd.pbs
echo \#PBS -e err/${file_traj}_embed_rmsd.err >> split_rmsd.pbs
echo \#PBS -r n >> split_rmsd.pbs
echo \#PBS -V >> split_rmsd.pbs
echo cd \$PBS_O_WORKDIR >> split_rmsd.pbs
echo ./split_rmsd\<split_rmsd.input\>split_rmsd.log >> split_rmsd.pbs
