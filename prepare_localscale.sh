#!/bin/bash

source run_parameters

servername=`hostname|cut -d\. -f2`
echo Job will be set up on \" $servername \"

norder=`python -c "from math import ceil;print int(ceil(float(${npoints})/${npoints_core}))"`

rm -f submitlist

for ((i=0;i<$norder;i++));
do 
    nstart=$(($i*${npoints_core}+1))
    if [ $(($i+1)) -eq $norder ] ; then
	nend=${npoints}
    else
	nend=$((($i+1)*${npoints_core}))
    fi
    echo id:$i starts from $nstart to $nend
    # write input
    input=input/${file_traj}_$(($i+1000)).input
    echo id:$i $input
    echo ${file_traj} > $input
    echo ${file_traj}_neighbor >> $input
    echo $nstart $nend >> $input
    echo $(($i+1000)) >> $input
    echo $neighbor_number $mds_dimension $mds_start $mds_step >> $input
    echo $cutoff_start $cutoff_step $cutoff_number >> $input
    echo $ncore >> $input
    # write pbs
    pbs=input/${file_traj}_$(($i+1000)).pbs
    echo id:$i $pbs
    echo \#!/bin/bash > $pbs
    echo \#PBS -l nodes=${nodes}:ppn=${ppn},walltime=${walltime} >> $pbs
    echo \#PBS -q ${queue} >> $pbs
    echo \#PBS -N p_lmds >> $pbs
    echo \#PBS -o err/${file_traj}_lmds_$(($i+1000)).out >> $pbs
    echo \#PBS -e err/${file_traj}_lmds_$(($i+1000)).err >> $pbs
    echo \#PBS -r n >> $pbs
    echo \#PBS -V >> $pbs
    echo cd \$PBS_O_WORKDIR >> $pbs
    echo mpiexec -n $(($nodes*$ppn)) p_local_mds\<${input}\>${file_traj}_lmds_$(($i+1000)).log >> $pbs
    echo qsub $pbs>> submitlist
done
chmod u+x submitlist

echo "cat localscale/${file_traj}_eps_1*" \> ${file_traj}_eps > merge_localscale_results.sh
echo "rm localscale/${file_traj}_eps_1*" >> merge_localscale_results.sh
chmod u+x merge_localscale_results.sh
