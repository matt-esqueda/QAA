Bi621 - Computational Methods in Genomic Analysis

``````````````````````

Bi 621 – Problem Set 6

``````````````````````


In this assignment, we will be using sample data to write a python script that will provide statistical measures of accuracy for whole genome assembly. 
We will then utilize the Velvet program to build assemblies from Illumina reads and use our script to assay those assemblies.



Part 0 – Write unit tests

Created a test file (Unit_test.fa) containing a few records to test whether the script in part 1 correctly extracts information from the header line and performs the calculations accurately.

Calculated the expected results  and stored them in expected_results.txt.



Part 1 – Contig length distributions

Using the sample data contained in contigs.fa (located on Talapas), extracted the length and k-mer coverage from each header line.

/projects/bgmp/shared/Bi621/contigs.fa

contigs.py will extract and “fix” the necessary header line information to calculate the following:

	contigs.py
	/home/mesqueda/bioinfo/Bi621/PS/ps6-matt-esqueda/Part_1/contigs.py

	number of contigs
    maximum contig length
    mean contig length
    total length of the genome assembly across the contigs
    mean depth of coverage for the contigs
	N50 value of assembly

contigs.py will store results Bucket_lists.tsv and produce a distribution histogram.

head Bucket_lists.tsv
    Contig length	Number of contigs in this category
    200	72
    300	47
    400	20
    500	18
    600	17
    700	16
    800	9
    900	11



Part 2 - Velvet

Installed Velvet in fresh environment with these commands

    $ conda create -n bgmp-velvet -c bioconda/label/cf201901 velvet python=3.10
    $ conda activate bgmp-velvet
    $ velvetg


velvetg - de Bruijn graph construction, error removal and repeat resolution
Version 1.2.10
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)


Default compile time options set by Velvet:

CATEGORIES=4
MAXKMERLENGTH=191
OPENMP=1
LONGSEQUENCES=1

Velvet manual: https://www.animalgenome.org/bioinfo/resources/manuals/velvet.pdf


Used quality trimmed Illumina PE100 reads from the following fastq files (located on Talapas)

	/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1
    /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2
    /projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched


In total, 6 velvet runs were performed with the following parameters:

	k31_cov.auto.len_default
	k41_cov.auto.len_default
	k49_cov.auto.len_default
	k31_cov.20x.len_default
	k31_cov.60x.len_default
	k31_cov.auto.len_500bp


Wrote velveth SLURM script 
    /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS6/velveth.sh

Wrote velvetg SLURM scipt
    /velvetg is run using: /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS6/velvetg.sh


*** Had memory allocation issues while running velvet on talapas. Must use #SBATCH --nodelist=n278         
    ### Set node to n=278 when doing velvet runs. Velvetg kept crashing when running on older nodes and would fail to produce contigs.fa. 

Example of failed run log:


Command exited with non-zero status 1
	Command being timed: "/projects/bgmp/mesqueda/miniconda3/envs/bgmp-velvet/bin/velvetg /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS6/results_k31.cov_20x.len_default_n278 -exp_cov 20 -cov_cutoff auto -ins_length auto"
	User time (seconds): 209.71
	System time (seconds): 2.25
	Percent of CPU this job got: 125%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:48.97
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 14875396
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 425088
	Voluntary context switches: 114
	Involuntary context switches: 404
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 1


A successful velveth/g run will produce the following in the results directory:
	
	contigs.fa	* make sure this file is present when velvet has finished 
	Graph2
	Last Graph
	Log
	PreGraph
	Roadmaps
	Sequences
	stats.txt

Results were collected from each run and analyzed by contigs.py. 
This produced descriptive statistics and distribution histograms for each run, 
stored in the appropriate directory contained in the results folder.

    /home/mesqueda/bioinfo/Bi621/PS/ps6-matt-esqueda/results






``````````````````````
Bi 621 – Problem Set 7

``````````````````````


In this assignment, we will use blastp and the Ensembl database to create a blast-based Reciprocal Best Hit (RBH) pipeline. 
This pipeline will be used to call all 1-to-1 RBH alignments between human and zebrafish proteins in Bi623.


Part 1 – Download protein fasta files and filter to retain just the longest protein per gene

Downloaded the following fasta files from the Ensembl database:

https://uswest.ensembl.org/info/data/ftp/index.html

	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Homo_sapiens.GRCh38.pep.all.fa
	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Danio_rerio.GRCz11.pep.all.fa

*** Don’t put in home directory


Downloaded tables of Gene stable IDs, Gene names, and Protein stable ID for human and zebrafish from Ensembl Biomart:

http://uswest.ensembl.org/biomart/martview/3ac6983a2c0184f0f4e7649ef9653d09
	* couldn’t download second file from site, maybe too many people trying to access Biomart??

	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Human_mart.tsv
	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Zfish_mart.tsv


Wrote two python scripts to parse through fasta files and retain the longest protein record per gene

tool1.py
    /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/tool1.py
    -  this will arrange the sequence lines of each record into one line, will be placed in bioinfo.py module as oneline_fasta()

tool2.py
    /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/tool2.py
    - this will extract the longest protein record for each gene in the read file and write it to a new fasta file, formatted for blast
    - there was a bug in this script, causing the peptide sequence to be added several times to each record. I think I found and fixed it.

Result analysis on Talapas:
    - Loaded the bbmap module on Talapas

- Opened an interactive session on Talapas:
    srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=1:00:00 --cpus-per-task=1 --pty bash

- Ran stats.sh on the fasta files and saved the results:
	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Human_LP_stats.txt
	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Zfish_LP_stats.txt



Part 2 – Blastp

Blast manual: https://www.ncbi.nlm.nih.gov/books/NBK279684/

Built human blast database and zebrafish blast database on Talapas with database.sh (the slurm script will build both):

	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/database.sh


Successful run example:
    -always check for Exit status 0!

usr/bin/time -v 

Note that the EASYBUILD tree is NOT supported by RACS.

Please contact user at sydes@uoregon.edu for questions regarding this tree.


Lmod is automatically replacing "gcc/7.3" with "easybuild".



Note that the EASYBUILD tree is NOT supported by RACS.

Please contact user at sydes@uoregon.edu for questions regarding this tree.



Note that the EASYBUILD tree is NOT supported by RACS.

Please contact user at sydes@uoregon.edu for questions regarding this tree.



Building a new DB, current time: 07/15/2022 10:35:58
New DB name:   /gpfs/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Homo_sapiens_database
New DB title:  Homo_sapiens.fa
Sequence type: Protein
Keep Linkouts: T
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 23506 sequences in 1.40466 seconds.
	Command being timed: "makeblastdb -in Homo_sapiens.fa -parse_seqids -dbtype prot -out Homo_sapiens_database"
	User time (seconds): 1.32
	System time (seconds): 0.07
	Percent of CPU this job got: 94%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.48
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 16200
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 4504
	Voluntary context switches: 60
	Involuntary context switches: 2
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0


Building a new DB, current time: 07/15/2022 10:35:59
New DB name:   /gpfs/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/Danio_rerio_database
New DB title:  Danio_rerio_final.fa
Sequence type: Protein
Keep Linkouts: T
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 30313 sequences in 1.06397 seconds.
	Command being timed: "makeblastdb -in Danio_rerio_final.fa -parse_seqids -dbtype prot -out Danio_rerio_database"
	User time (seconds): 0.98
	System time (seconds): 0.06
	Percent of CPU this job got: 90%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.15
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 16884
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 4471
	Voluntary context switches: 47
	Involuntary context switches: 2
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0


* Must load these modules:

module load easybuild
module load eb-hide/1
module load BLAST+/2.2.31

Created blastp commands to align each fasta file to the other’s database (H_to_Z and Z_to_H)

	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/blast_Homo_sapiens_to_D.sh
	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS7/blast_Danio_rerio_to_H.sh

* Important notes
	- used 8 cores on a node for the blast alignment 
    - included -evalue 1e-6 and -use_sw_tback (Compute locally optimal Smith-Waterman alignments)
    - included bitscore in output format for analysis 
	- forgot to load eb-hide/1 module, had to cancel first job and restart


Bash commands to create summary information:
	number of hits in the file:
		wc -l blast_results.dre_hsa

	The 10 hits with the highest bit scores:
		cat blast_results.dre_hsa | cut -f 1-3 | sort -k3 -n -r | head

	All hits with the lowest evalue, sorted by bitscore (highest to lowest):
		cat blast_results.dre_hsa | grep "9e-99" | sort -k3 -n -r > ZLE.txt






``````````````````````
Bi 621 – Problem Set 8

``````````````````````


In this assignment, we will work with RNA-seq data to map reads to an existing reference genome. This will be performed on Talapas.


Within the projects directory on Talapas, created a new PS8 directory:

	$ mkdir -p /projects/bgmp/<YOU>/bioinfo/Bi621/PS/PS8

The data files were located on Talapas, they will be referenced at these locations in scripts:

	/projects/bgmp/shared/Bi621/dre_WT_ovar12_R1.qtrim.fq.gz
    /projects/bgmp/shared/Bi621/dre_WT_ovar12_R2.qtrim.fq.gz

Created a new directory with PS8, called dre which holds the primary_assembly fasta file:

	$ mkdir -p /projects/bgmp/<YOU>/bioinfo/Bi621/PS/PS8/dre

Downloaded the zebrafish reference genome (wget) by chromosome (FASTA) and gene set (GTF) from the ensemble website : https://uswest.ensembl.org/info/data/ftp/index.html
	* still having trouble downloading files from ensemble on UO secure??
	* downloaded fine at home
	* GTF file placed in /PS8

Installed STAR and samtools

	$ conda activate bgmp_py310
    $ conda install star -c bioconda	
    $ conda install samtools -c bioconda
    $ samtools --version

$ STAR -- version 2.7.10a

$ samtools --version
samtools 1.6
Using htslib 1.6
Copyright (C) 2017 Genome Research Ltd.


NOTE: must conda activate bgmp_py310 to use STAR or samtools with this install


Created directory to contain the STAR database:

    /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/Danio_rerio.GRCz11.dna.ens107.STAR_2.7.10a

Built STAR database with STAR_database.sh 
    /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/STAR_database.sh
	* remember to conda activate bgmp_py310 in the SLURM scripts
    * Must request 8 cores on 1 mode (--runThreadN 8 included in Slurm script)
	* Forgot to due this and first run failed: terminated by signal7 (watch for this!!)


STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/Danio_rerio.GRCz11.dna.ens107.STAR_2.7.10a --genomeFastaFiles /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/dre/Danio_rerio.GRCz11.dna.primary_assembly.fa --sjdbGTFfile /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/Danio_rerio.GRCz11.107.chr.gtf
	STAR version: 2.7.10a   compiled: 2022-01-14T18:50:00-05:00 :/home/dobin/data/STAR/STARcode/STAR.master/source
Jul 13 10:47:28 ..... started STAR run
Jul 13 10:47:28 ... starting to generate Genome files
Jul 13 10:47:44 ..... processing annotations GTF
Jul 13 10:47:50 ... starting to sort Suffix Array. This may take a long time...
Jul 13 10:47:57 ... sorting Suffix Array chunks and saving them to disk...
Command terminated by signal 7
	Command being timed: "STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/Danio_rerio.GRCz11.dna.ens107.STAR_2.7.10a --genomeFastaFiles /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/dre/Danio_rerio.GRCz11.dna.primary_assembly.fa --sjdbGTFfile /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/Danio_rerio.GRCz11.107.chr.gtf"
	User time (seconds): 35.09
	System time (seconds): 3.69
	Percent of CPU this job got: 78%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:49.13
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 4710032
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 15771
	Voluntary context switches: 119
	Involuntary context switches: 769
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0



Ran STAR to align the reads to the reference genome using STAR_alignReads.sh
    /projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/STAR_alignReads.sh
	* Must request 8 cores on 1 mode (--runThreadN 8 included in Slurm script)

-STAR align worked and produced correct SAM file. However, accidently overwrote SAM file while working on script and had to run STAR again. 
Followed the same procedure and everything worked well


Used samtools to convert SAM file to BAM file. Also extracted all reads from chromosome 1 into a new SAM file using samtools	
    
    Note: enter samtools on the command line with no options to get help screen. 

	samtools commands: 
		converts SAM file to BAM file:
		samtools view -S -b Aligned.out.sam > Aligned.out.bam

		creates the index file (BAI) for the BAM file:
		samtools index Aligned.out.bam
		
		extract the data for chromosome 1 and store is SAM file:
		samtools view -h Aligned.out.bam chr1 > Aligned.out.chr1.sam
		
		counts the number of alignments on chromosome 1:
		wc -l Aligned.out.chr1.sam



Part 2: Python and SAM

Wrote SAM.mapped.py to check if each of the reads in the SAM file were properly mapped
    /home/mesqueda/bioinfo/Bi621/PS/ps8-matt-esqueda/SAM.mapped.py
        - reads may be aligned more than once, SAM.mapped.py includes check to ensure that mapped reads are not counted more than once. 
        It will print a count of mapped reads and a count of unmapped reads in the SAM file.

    - still a bit confused about:
    if((flag & 4) != 4):
        mapped = True

    -talked with Jason about flag, it is making more sense now
	
	/projects/bgmp/mesqueda/bioinfo/Bi621/PS/PS8/SAM.mapped.py

	-worked properly on test file
	-worked properly on SAM file






######################################################################################################################################################################

######################################################################################################################################################################

Bi622 -  Genomic Topics


Demultiplexing and Index Swapping Assignment the first

	Part 1: Quality Score Distribution-per-nucleotide

	1. initial data exploration
	
	i.	bash commands used:
			zcat 1294_S1_L008_R1_001.fastq.gz | head
				- used on all four files (R1, R2, R3, R4)
				- determined that R1 = read 1, R2 = read 2
				- determined that R3 = index 1, R4 = index 2
			

	ii. 	zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
				- used on all four files 
				- determined that the read length for the read files is 101 nt
				- determined that the read length for the index files is 8 


	iii.	zcat 1294_S1_L008_R2_001.fastq.gz | head -1000000 | tail -20 > /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/test_index.fq.gz 
				- used to create test index fastq files
			
			zcat 1294_S1_L008_R1_001.fastq.gz | head -1000000 | tail -20 > /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/test_read.fq.gz 
				- used to create test read fastq files
				

	2. Used code from PS4 part 1 to create part_1.py python code to determine the average quality scores at each postion for all reads and indexes. 
	Generated a per nucleotide mean distribution for each of the four files:

		import argparse
		import gzip
		import bioinfo 

		python script:
		part_1.py = /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/part_1.py

		slurm script:
		part_1.sh = /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/part_1.sh

		1294_S1_L008_R1_001.fastq.gz
		1294_S1_L008_R2_001.fastq.gz
		1294_S1_L008_R3_001.fastq.gz
		1294_S1_L008_R4_001.fastq.gz
		
		- all located in /projects/bgmp/shared/2017_sequencing/

		ii. Good quality scores for indexes and reads...

		Used these commands to find how many indexes have undetermined (N) base calls:

		zcat 1294_S1_L008_R2_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l	
		3976613

		zcat 1294_S1_L008_R3_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l
		3328051	


	Part 2: Develop an algorithm to de-multiplex the samples

		- Wrote psuedocode to develop algorithm for de-multiplexing files and reporting the amount of index hopping. 
			/projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/Part_2/psuedocode.txt



Demultiplexing and Index Swapping Assignment the Second
	
	Reviewed three peers psuedocode and provided feedbach through github.

	Gained some insight regarding how to make my own code more efficient. 


Demultiplexing and Index Swapping Assignment the Third

	Used psuedocode from part 1 to write python script to demultiplex the samples. 

	*Look up itertools, permutations, for creating mismatched indexes dictionary
		https://docs.python.org/3/library/itertools.html (very useful)

	Biggest challenges were opening all four read files simultaneously and writing to all the approiate fastq files.
		-used gzip to open zipped files, good to know for the future

	Added reverse_complement() to bioinfo.py

	Figured out the general logic for separating the indexes appropiately, and they code seems to be working for the test files
		-hopped and unknown put into their own folders

	bash command for test:
	./demultiplex.py -i test_indexes  -t1 R1_test_input.fq  -t2 R2_test_input.fq  -t3 R3_test_input.fq  -t4 R4_test_input.fq -ic 30
		-will hardcode the actual files after finished testing, rather than using argparse. Command becoming way too long

	bash command for interactive node:
	srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash
	

	python ./demultiplex.py -ic 30
		*set index cutoff to 30

	Ran with the read files, failed the first time failed because index_mismatch_dict was not set up properly
		-fixed dictionary, run again
	
	Run failed because incorrect variable used in readfile(), need to be careful to not mix up variables, better naming...

	Run failed becuase tsv variables got mixed up, 'pair' is not defnined. Easy fix, run again..

	Run worked, all files created properly and correct counts with appropiate tables and tsv created

	Used pigz to zip all fastq files (works fast!)
		pigz *fq 
	
	Tried creating graphs for stats output. The graphs were either uniformative or too many variables to be readable. Tried to create heatmap but was too difficult.
		-For the future, look to pandas and seaborn for creating heatmap




######################################################################################################################################################################

######################################################################################################################################################################


Bi623 - Topics in Genomics Analysis



Assignment 1 - PS7 continued

Part 3 - Reciprocol Best Hit

Used provided blastp files

	/projects/bgmp/shared/Bi623/PS7_RBH_Bi623/H_to_zfishdb.blastp 
 	/projects/bgmp/shared/Bi623/PS7_RBH_Bi623/Z_to_homodb.blastp 

Unix command to sort blast files by protein id and e-value and pipe to sorted files
	cat H_to_zfishdb.blastp | sort -k1,1 -k11,11g > sorted_H_to_Zdb.blastp
	cat Z_to_homodb.blastp | sort -k1,1 -k11,11g > sorted_H_to_Zdb.blastp
	
	*** These sorted it a more difficult way, do not use. Use sort commands below


Fixed sort commands, these will sort so that each protein is sorted by top e-values
	sort -k1,1 -k11,11g H_to_Zfishdb.blastp > sorted_H_to_Zdb.blastp
	sort -k1,1 -k11,11g Z_to_homodb.blastp > sorted_Z_to_Hdb.blastp

Created test files
	head -50 sorted_H_to_Zdb.blastp > test_H_to_Z
	head -50 sorted_Z_to_Hdb.blastp > test_Z_to_H 

	/projects/bgmp/mesqueda/bioinfo/Bi623/Assignment_1/test_H_to_Z
	/projects/bgmp/mesqueda/bioinfo/Bi623/Assignment_1/test_Z_to_H

This was used to sort through the e-values incorectly however, the logic works and this can be useful in the future
 for v in human_dict.values():
                e_vals.setdefault(v, 0)
                e_vals[v] += 1
            
        my_dict = [(k, v) for k, v in human_dict.items() if e_vals[v] == 1]   

Realized that I had misread the assignment and was tryign to filter the data incorectly. Fortunately, the correct way is easier 
so it should not be difficult to fix the code as necessary.
	
	General approach

	1) create this dicionary for human and zebrafish blast outputs
		human_dict[h_prot] = (z_prot, e_val)

	2) check for matching best hits
		if h_prot == human_id:
        match[h_prot] = fish_id
	
	3) create dictionaries from biomart files for both
		for line in h_mart:
        column = line.split("\t")
        h_mart_dict[column[1]] = (column[0], column[2].strip())
	
	4) open RBH tsv file and write a table of all RBH relationships 
		for h_prot, z_prot in match.items():
			etc... 
			* remember to reset list at end of each loop

Command (w/ argparse) used to run RBH.py:
	./RBH.py -s sorted_H_to_Zdp.blastp -z sorted_Z_to_Hdb.blastp -hm Human_BioMart.tsv -zm Zebrafish_BioMart.tsv

	the python script will create tsv table with RBH relationships
		/projects/bgmp/mesqueda/bioinfo/Bi623/Assignment_1/Human_Zebrafish_RBH.tsv
		wc -l Human_Zebrafish_RBH	7912




Assignment - QAA

	bash command for interactive node:
	srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash

	Part 1 - Read quality score distributions

		Load and run fastqc commands:
		module spider fastqc
		module load fastqc/0.11.5
		echo $PATH
		fastqc -h (help)
			fastqc seqfile1 seqfile2 -o output diretory (optional)
		* everytime I run a fastqc command, the following perl warning pops up:
			perl: warning: Setting locale failed.
			perl: warning: Please check that your locale settings:
        		LANGUAGE = (unset),
        		LC_ALL = (unset),
        		LANG = "C.UTF-8"
    			are supported and installed on your system.
			perl: warning: Falling back to the standard locale ("C").
		-asked Pete about this and he said to ignore these warnings, still got correct output

		Library Assignments:
			19_3F_fox_S14_L008
			7_2E_fox_S6_L008
			located in:	
				/projects/bgmp/shared/2017_sequencing/demultiplexed/
					19_3F_fox_S14_L008_R1_001.fastq.gz
					19_3F_fox_S14_L008_R2_001.fastq.gz
					7_2E_fox_S6_L008_R1_001.fastq.gz
					7_2E_fox_S6_L008_R2_001.fastq.gz
		
		Ran Fastqc on demultiplexed file pairs 
			-outputs fastqc.html and fastqc.zip files for each gz file
		
		Ran demux_plot.sh on each of the gz files
			usr/bin/time -v python ./demux_plot.py -l 101 -f " " -o " "
			-outputs plots
	

	Part 2 - Adaptor Trimming Comparison

		Created a new conda enviroment QAA
			conda create --name QAA
		Installed cutadapt and trimmomatic packages in the QAA environment
			conda install -n QAA cutadapt 
				-version 4.1
			conda install -n QAA trimmomatic
				-version 0.39

		R1 Adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
		R2 Adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

		zcat 19_3F_fox_S14_L008_R1_001.fastq.gz | grep "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" | wc
		zcat 19_3F_fox_S14_L008_R2_001.fastq.gz | grep "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" | wc
		zcat 7_2E_fox_S6_L008_R1_001.fastq.gz | grep "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" | wc
		zcat 7_2E_fox_S6_L008_R2_001.fastq.gz | grep "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" | wc

		Ran cutadapt on files to remove adapters with cutadapt.sh
			output:
				/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/7_2E_fox_S6_L008_R1_001.fastq. gz
				/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/19_3F_fox_S14_L008_R2_001.fastq.gz

		Ran trimmomatic on intermediate files with trimmomatic.sh
			output:
				/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/7_2E_fox_S6_L008_R1_001.trimmed.fastq.gz
				/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/7_2E_fox_S6_L008_R2_001.trimmed.fastq.gz
				/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/19_3F_fox_S14_L008_R1_001.trimmed.fastq.gz
				/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/19_3F_fox_S14_L008_R2_001.trimmed.fastq.gz

		Ran fastqc on trimmed fastq.gz files for 7 and 19 (can use fastqc.sh)
			fastqc 7_3F_fox_S14_L008_R1_001.trimmed.fastq.gz 7_3F_fox_S14_L008_R2_001.trimmed.fastq.gz -o /projects/bgmp/mesqueda/bioinfo/Bi623/QAA
			fastqc 19_3F_fox_S14_L008_R1_001.trimmed.fastq.gz 19_3F_fox_S14_L008_R2_001.trimmed.fastq.gz -o /projects/bgmp/mesqueda/bioinfo/Bi623/QAA

			


	Part 3 Alignment and Strand Specificity 

		Installed the following:

			STAR
				conda install star -c bioconda
				environment location: /projects/bgmp/mesqueda/miniconda3/envs/QAA

			Numpy
				conda install -c anaconda numpy	
				environment location: /projects/bgmp/mesqueda/miniconda3/envs/QAA
			
			pysam
				conda install -c bioconda pysamST
				environment location: /projects/bgmp/mesqueda/miniconda3/envs/QAA

			matplotlib
				conda install matplotlib
				environment location: /projects/bgmp/mesqueda/miniconda3/envs/QAA

			htseq
				conda install -c bioconda htseq
				environment location: /projects/bgmp/mesqueda/miniconda3/envs/QAA
		
		Downloaded mouse primary_assembly and GTF files from FTP server on ensembl:
			/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/mouse_assembly/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
			/projects/bgmp/mesqueda/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.107.gtf.gz
		
		Generated alignment database from assembly and GTF files using STAR
			STAR-database.sh
		
		Aligned the reads from both files using STAR
			STAR_alignReads.sh
			-seemed to work on 7_2E_fox files but not on 19_3F_fox files, need to run again
			-ran again for the 19_3F_fox files, worked the second attempt
			

			Used SAM.mapped.py script from PS8 to count the number of mapped and unmapped reads
				7_2E_fox:
				mapped_reads = 9424753
				unmapped_reads = 341679
	
				19_3F_fox:
				mapped_reads = 30512182
				unmapped_reads = 1287922

			Ran htseq-count twice to count the reads that map to GFF features
				htseq.sh
					htseq-count --stranded=yes <Aligned SAM File> <GTF file>
					htseq-count --stranded=reverse <Aligned SAM File> <GTF file> 
			

			Stranded htseq counts

			19_3F_fox:
			grep '^ENSM' both_stranded_outputs | awk '{sum+=$2} END {print sum}'
			555,571 matches

			7_2E_fox:
			grep '^ENSM' both_stranded_outputs | awk '{sum+=$3} END {print sum}'
			182,750 matches


			Reversed htseq counts
			
			19_3F_fox:
			grep '^ENSM' both_reversed_outputs | awk '{sum+=$2} END {print sum}'
			12,936,400 matches

			7_2E_fox:
			grep '^ENSM' both_reversed_outputs | awk '{sum+=$3} END {print sum}'
			4,027,379 matches


