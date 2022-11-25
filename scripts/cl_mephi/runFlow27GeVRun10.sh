#!/bin/bash

#mkdir -p /mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/OUT_new/27GeV/sge_out/
#mkdir -p /mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/OUT_new/27GeV/sge_err/

#SBATCH -D /mnt/pool/rhic/1/demanov/cherenkov/TMP
#SBATCH -J flow
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=1-2
#SBATCH -o /mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/OUT_new/27GeV/sge_out/slurm_%A_%a.out
#SBATCH -e /mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/OUT_new/27GeV/sge_err/slurm_%A_%a.err

PID_TYPE=$1
WORK_MODE=$2
IN_REC=$3
IN_FLAT=$4
SYF=$5

export ENERGY=27
export MAIN_DIR=/mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/
export INPUT=$MAIN_DIR/lists/lists${ENERGY}GeVRun10/StRuns${SLURM_ARRAY_TASK_ID}.list
export OUTPUT=$MAIN_DIR/OUT_new/${ENERGY}GeV/${ENERGY}GeV_${WORK_MODE}_${PID_TYPE}_${SYF}_${SLURM_ARRAY_JOB_ID}

source /mnt/pool/rhic/4/parfenovpeter/Soft/Basov/ROOT/build/bin/thisroot.sh

mkdir -p ${MAIN_DIR}/OUT_new/${ENERGY}GeV
cd ${MAIN_DIR}/OUT_new/${ENERGY}GeV
 
mkdir -p $OUTPUT/root
mkdir -p $OUTPUT/log

#cd ${MAIN_DIR}/build2/
cd /home/demanov97/STAR_Analysis/build2/
./FemtoDstAnalyzer_${PID_TYPE} -i $INPUT -o $OUTPUT/root/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root -m ${WORK_MODE} -r ${IN_REC} -f ${IN_FLAT} -g ${ENERGY} &>>$OUTPUT/log/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log
#./combPID_Gaus -i ${MAIN_DIR}/OUT/27GeV/Histo_27GeVRun10_NewXY.root -f ${MAIN_DIR}/OUT/test.root -o $OUTPUT/root/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root -g 27 -m ${WORK_MODE} -p ${SLURM_ARRAY_TASK_ID} &>>$OUTPUT/log/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log
