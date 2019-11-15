#!/bin/bash

npath=`pwd`
[ ! -d src ] && mkdir src
cd src
if [ ! -d ARModule ]; then
    wget "https://github.com/lwang-astro/ARModule/archive/master.zip" 
    unzip master.zip
    mv ARModule-master ARModule
    rm -f master.zip
fi
if [ ! -d PeTar ]; then
    wget "https://github.com/AICS-SC-GROUP/PeTar/archive/master.zip" 
    unzip master.zip
    mv PeTar-master PeTar
    rm -f master.zip
fi
if  [ ! -d FDPS ]; then
    wget "https://github.com/FDPS/FDPS/archive/master.zip"
    unzip master.zip
    mv FDPS-master FDPS
    rm -f master.zip
fi
cd $npath
flist='interface.cc interface.py interface.h test_interface.cc test_interface.py run_one_cluster_movie.py'
for file in $flist
do
    ln -sf src/PeTar/amuse-interface/$file 
done

