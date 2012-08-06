#!/bin/bash

hostname=`hostname|cut -d\. -f2`

if [ -f LSDMap_$hostname.inc ]; then
    echo Try compiling programs on \" $hostname \"
    cp LSDMap_$hostname.inc LSDMap.inc
    make
else
    echo Warning: Cannot recognize the host name.
    echo Please modify LSDMap.inc according to the comments.
    cp LSDMap_davinci.inc LSDMap.inc
fi






