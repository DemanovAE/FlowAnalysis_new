#!/bin/bash

PID_TYPE=$1
WORK_MODE=$2
ENERGY=$3
SYS=$4
IN_REC=$5
IN_FLAT=$6

if [[ $ENERGY == "11" ]]; then
  sbatch runFlow11GeV.sh $PID_TYPE $WORK_MODE 
fi

if [[ $ENERGY == "14" ]]; then
  sbatch runFlow14GeV.sh $PID_TYPE $WORK_MODE 
fi

if [[ $ENERGY == "19" ]]; then
  sbatch runFlow19GeV.sh $PID_TYPE $WORK_MODE 
fi

if [[ $ENERGY == "27" ]]; then
  sbatch runFlow27GeVRun10.sh $PID_TYPE $WORK_MODE $IN_REC $IN_FLAT $SYS
fi

if [[ $ENERGY == "39" ]]; then
  sbatch runFlow39GeV.sh $PID_TYPE $WORK_MODE 
fi

if [[ $ENERGY == "62" ]]; then
  sbatch runFlow62GeV.sh $PID_TYPE $WORK_MODE 
fi

if [[ $ENERGY == "281" ]]; then
  sbatch runFlow28GeVPer1.sh $PID_TYPE $WORK_MODE $IN_REC $IN_FLAT $SYS
fi

if [[ $ENERGY == "282" ]]; then
  sbatch runFlow28GeVPer2.sh $PID_TYPE $WORK_MODE $IN_REC $IN_FLAT $SYS
fi
