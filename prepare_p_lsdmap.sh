#!/bin/bash

source run_parameters

servername=`hostname|cut -d\. -f2`
servername="hector"

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
    input=ExampleInputFile
    echo id:$i $input
    echo ${file_traj} > $input
    echo $nstart $nend >> $input
    echo $neighbor_number >> $input
    echo $(($i+1000)) >> $input
    echo $mds_dimension $mds_start $mds_step >> $input
    echo $cutoff_start $cutoff_step $cutoff_number >> $input
    echo $ncore >> $input
    echo $status_lsdmap $status_localscale >> $input
    echo ${file_weight} >> $input
    echo ${file_traj}_dif >> $input
    echo 0.0 >> $input
    if [ $status_localscale -eq 1 ]; then
	echo $eps_column_id >> $input
    else
	echo ${constant_eps} >> $input
    fi


    # write pbs
    pbs=ExampleBatchSrcipt.pbs
    echo id:$i $pbs
    echo \#!/bin/bash --login > $pbs
    echo \#PBS -l mppnppn=${ppn} >> $pbs
    echo \#PBS -l mppwidth=${nodes} >> $pbs
    echo \#PBS -l walltime=${walltime} >> $pbs
    echo \#PBS -A e290 >> $pbs
    echo \#PBS -N p_lsdmap >> $pbs
    echo \#PBS -o err/${file_traj}_$(($i+1000)).out >> $pbs
    echo \#PBS -e err/${file_traj}_$(($i+1000)).err >> $pbs
    echo \#PBS -r n >> $pbs
    echo \#PBS -V >> $pbs
    echo cd \$PBS_O_WORKDIR >> $pbs
    echo aprun -n $(($nodes)) -N $(($ppn)) main\<${input}\>${file_traj}_$(($i+1000)).log >> $pbs
    echo qsub $pbs>> submitlist
done
