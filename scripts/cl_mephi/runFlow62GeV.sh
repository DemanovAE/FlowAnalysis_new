#!/bin/bash

#SBATCH -D /mnt/pool/rhic/1/demanov/cherenkov/TMP
#SBATCH -J flow
#SBATCH -p compute
#SBATCH --time=12:00:00
#SBATCH -a 1-12
#SBATCH -o /mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/OUT_new/sge_out/slurm_%A_%a.out
#SBATCH -e /mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/OUT_new/sge_err/slurm_%A_%a.err

PID_TYPE=$1
WORK_MODE=$2
export ENERGY=62
export MAIN_DIR=/mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/
export INPUT=$MAIN_DIR/lists/lists${ENERGY}GeV/StRuns${SLURM_ARRAY_TASK_ID}.list
export OUTPUT=$MAIN_DIR/OUT_new/${ENERGY}GeV/${ENERGY}GeV_${WORK_MODE}_${SLURM_ARRAY_JOB_ID}_${PID_TYPE}

source /mnt/pool/rhic/4/parfenovpeter/Soft/Basov/ROOT/build/bin/thisroot.sh

mkdir -p ${MAIN_DIR}/OUT_new/${ENERGY}GeV
cd ${MAIN_DIR}/OUT_new/${ENERGY}GeV
 
mkdir -p $OUTPUT/root
mkdir -p $OUTPUT/log

cd ${MAIN_DIR}/build/
./FemtoDstAnalyzer_${PID_TYPE} -i $INPUT -o $OUTPUT/root/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root -m ${WORK_MODE} -g ${ENERGY} &>>$OUTPUT/log/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log
#./combPID_Gaus -i ${MAIN_DIR}/OUT/27GeV/Histo_27GeVRun10_NewXY.root -f ${MAIN_DIR}/OUT/test.root -o $OUTPUT/root/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root -g 27 -m ${WORK_MODE} -p ${SLURM_ARRAY_TASK_ID} &>>$OUTPUT/log/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log
