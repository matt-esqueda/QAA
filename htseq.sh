#!/bin/bash

#SBATCH --partition=bgmp                    ### Partition (like a queue in PBS)
#SBATCH --job-name=htseq                    ### Job Name
#SBATCH --nodes=1                           ### Number of nodes needed for the job
#SBATCH --cpus-per-task=4                   ### Number of cpus per task
#SBATCH --account=bgmp                      ### Account used for job submission
#SBATCH --output=reads-%j.out               ### File to store output
#SBATCH --error=reads-%j.error              ### File to store error

conda activate QAA

# stranded
/usr/bin/time -v htseq-count --stranded=yes /projects/bgmp/mesqueda/bioinfo/Bi623/QAA/19_3F_fox_Aligned.out.sam /projects/bgmp/mesqueda/bioinfo/Bi623/QAA/7_2E_fox_Aligned.out.sam /projects/bgmp/mesqueda/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.107.gtf

# reverse
# /usr/bin/time -v htseq-count --stranded=reverse /projects/bgmp/mesqueda/bioinfo/Bi623/QAA/19_3F_fox_Aligned.out.sam /projects/bgmp/mesqueda/bioinfo/Bi623/QAA/7_2E_fox_Aligned.out.sam /projects/bgmp/mesqueda/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.107.gtf