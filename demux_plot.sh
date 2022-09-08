#!/bin/bash

#SBATCH --partition=bgmp                    ### Partition (like a queue in PBS)
#SBATCH --job-name=7_R2                     ### Job Name
#SBATCH --nodes=1                           ### Number of nodes needed for the job
#SBATCH --cpus-per-task=1                   ### Number of cpus per task
#SBATCH --account=bgmp                      ### Account used for job submission
#SBATCH --output=reads-%j.out               ### File to store output
#SBATCH --error=reads-%j.error              ### File to store error



# # read1
/usr/bin/time -v python ./demux_plot.py -l 101 -f "/projects/bgmp/shared/2017_sequencing/demultiplexed/7_2E_fox_S6_L008_R2_001.fastq.gz" -o "7_R2.png"


# index1
# python ./part_1.py -l 8 -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" -o "index1.png"

# # index2
# python ./part_1.py -l 8 -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -o "index2.png"

# read2
# python ./part_1.py -l 101 -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -o "read2.png"