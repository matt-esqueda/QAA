#!/bin/bash

#SBATCH --partition=bgmp                    ### Partition (like a queue in PBS)
#SBATCH --job-name=f                        ### Job Name
#SBATCH --nodes=1                           ### Number of nodes needed for the job
#SBATCH --cpus-per-task=4                   ### Number of cpus per task
#SBATCH --account=bgmp                      ### Account used for job submission
#SBATCH --output=reads-%j.out               ### File to store output
#SBATCH --error=reads-%j.error              ### File to store error

# 19
# fastqc 19_3F_fox_S14_L008_R1_001.trimmed.fastq.gz 19_3F_fox_S14_L008_R2_001.trimmed.fastq.gz -o /projects/bgmp/mesqueda/bioinfo/Bi623/QAA

# 7 
fastqc 7_2E_fox_S6_L008_R1_001.trimmed.fastq.gz 7_2E_fox_S6_L008_R2_001.trimmed.fastq.gz -o /projects/bgmp/mesqueda/bioinfo/Bi623/QAA