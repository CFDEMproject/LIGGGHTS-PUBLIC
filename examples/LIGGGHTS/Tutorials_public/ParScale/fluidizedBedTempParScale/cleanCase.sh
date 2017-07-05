#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
#--------------------------------------------------------------------------------#



#- clean up case
echo "deleting data at: $casePath ?"
#read
rm -r $casePath/*.e*
rm -r $casePath/*.o*
rm -r $casePath/log*
rm -r $casePath/postParticles/*
rm -r $casePath/postGlobal/*
rm -r $casePath/post/*.*
rm -r $casePath/post/restart/*.restart
rm -r $casePath/log.*
rm -r $casePath/*Restart*
rm -r $casePath/*restart*
rm -r $casePath/*.dat

rm -r $casePath/pascal/0.*

echo "done"


