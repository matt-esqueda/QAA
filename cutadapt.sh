#!/bin/bash

#SBATCH --partition=bgmp                    ### Partition (like a queue in PBS)
#SBATCH --job-name=c                        ### Job Name
#SBATCH --nodes=1                           ### Number of nodes needed for the job
#SBATCH --cpus-per-task=1                   ### Number of cpus per task
#SBATCH --account=bgmp                      ### Account used for job submission
#SBATCH --output=reads-%j.out               ### File to store output
#SBATCH --error=reads-%j.error              ### File to store error


# 19_R1 adapter sequence
#/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o "19_3F_fox_S14_L008_R1_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz"

# 19_R2 adapter sequence
#/usr/bin/time -v cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o "19_3F_fox_S14_L008_R2_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R2_001.fastq.gz"

# 7_R1 adapter sequence
#/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o "7_2E_fox_S6_L008_R1_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/demultiplexed/7_2E_fox_S6_L008_R1_001.fastq.gz"

# 7_R2 adapter sequence
/usr/bin/time -v cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o "7_2E_fox_S6_L008_R2_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/demultiplexed/7_2E_fox_S6_L008_R2_001.fastq.gz"