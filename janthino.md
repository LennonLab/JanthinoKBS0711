# Lennon Lab Genome Analysis: *Janthinobacterium* #
---

## Project Goal: *De Novo* Genome Assembly and Annotation ##

IU has a computer cluster specifically designed for genomic analysis: [Mason](https://kb.iu.edu/data/bbhh.html)

Log-on to Mason from unix terminal

    > ssh [username]@mason.indiana.edu

Log-on to Mason from Putty (using cmd or PowerShell)

    > putty -ssh [username]@mason.indiana.edu

*OR (if you've saved settings for Mason)*

    > putty -load mason

Some of the following may take >20 minutes. 
Mason automatically kills memory/time intensive jobs.
However, you can start an interactive session using the following:

    > qsub -I -q shared -l nodes=1:ppn=4,vmem=10gb,walltime=4:00:00

Or just write a work flow script and submit with torque (see below)

***Time to get the data and start some fun***  
The original sequence files (fastq.gz) can be found at:

    /N/dc2/projects/Lennon_Sequences/2014_Janthino

1) Wild-Type Ancestor - *KBS0711_WT_R*.fastq*
2) White Mutant - *KBS0711_White_R*.fastq

## Assembly ##
### Sequence Quality Checks ###

Before you assemble and annotate the *Janthino* genome, you first need to assess the quality of the raw data
cd 
**A. Copy files into a working directory** 

    > cd /N/dc2/projects/Lennon_Sequences/  
    > mkdir ../(Your Dir)
    > cp ./KBS0711*.fastq.gz ../(Your Dir)
    > cd ../(Your Dir)
    
Run as a background process (recommended)

    > nohup cp ./KBS0711*.fastq.gz ../(Your Dir) &

Did a nohub process finish?

    > ps -U[username]
    OR
    > top -U[username]

**B. Unzip compressed files** (for more information about gunzip see: [Linux / Unix Command: gzip](http://linux.about.com/od/commands/l/blcmdl1_gzip.htm))

    > cd ../(Your Dir)
    > gunzip -fv ./*.gz

Again this might take some time so you might need to run in background

    > nohup gunzip -fv ./*.gz &

**C. Check raw sequence quality with *FastQC*** (for more information see: [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

    > module load java
    > module load fastqc
    > for file in ./*fastq
    > do
    > fastqc $file >> results.out
    > done

Because the output is in HTML, you will need to move this file to your local machine before you can open it. Moving files works best in Quarry so open another ssh session

    > cd /N/dc2/projects/Lennon_Sequences/(Your Dir)/
    > cp ./*.zip /afs/iu.edu/home/[Path to RFS]

Or you can just look at the summary

    > less ./KBS0711*_fastqc/summary.txt
    You can switch between the two files using 
    > :n
    OR
    > :p
    Exit with
    > q

What you will notice is that there are both warnings and failures for both files. Hopefully, we can deal with these by some pre-processing.

### Sequence Pre-Processing ###


**D. Interleave Paired End Reads**  
Paired end sequencing (from HiSeq or MiSeq) yield two files per sample: R1 and R2.
To assemble the raw reads into larger contigs, aligning software needs paired reads in the same file and in the correct order. 
The process used to do this is called *interleaving*. 
There are a few tools out there for interleaving paired end sequences.
The software package *Velvet* includes a perl script: *sufflesSequences_fastq.pl*.
The software package *Khmer* includes a python script: *interleave-reads.py*.
Both should work, but I haven't actually tested this.

***Velvet Method***

    > module load velvet
    > shuffleSequences_fastq.pl  ./KBS0711_w(t/hite)_R1.fastq ./KBS0711_w(t/hite)_R2.fastq ./KBS0711_w(t/hite).fastq

The Khmer tools do not recognize these as paired reads due to changes to files names (still exploring). Use the Khmer script for interleaving.

***Khmer Methods (Preferred)***

    > module load python
    > module load khmer
    > interleave-reads.py ./KBS0711_w(t/hite)_R1.fastq ./KBS0711_w(t/hite)_R2.fastq -o ./KBS0711_w(t/hite).fastq

*This takes a long time - go get a beer!!!!*

It helps to occasionally check the length of the sequences

    > wc KBS0711_w(t/hite).fastq

When I do this I get: WT - 92972872, White - 4554608


**E. Remove Any Low Quality Reads**

If you look at the FastQC output, you will see that there is quite a bit of variation in the quality (Phred quality scores) for each file. Just as a refresher, let's run *FastQC* again:

    > fastqc ./KBS0711_w(t/hite).fastq >> results.out
    > less ./KBS0711_w(t/hite)_fastqc/summary.txt

***Remove Low Quality Reads***

Low quality reads with a command from the FastX-Toolkit package

    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
    > fastq_quality_filter -Q33 -q 30 -p 50 -i ./KBS0711_wt.fastq > ./KBS0711_wt.trim.fastq
    > > fastq_quality_filter -Q33 -q 30 -p 50 -i ./KBS0711_wt.fastq > ./KBS0711_white.trim.fastq

Options:
- - Q33 = Illumina Quality Scores
- -q = minimum quality score
- -p = minimum # of bases that must have quality score

***Check Quality Again***

**F. Remove Orphaned Reads**

Orphaned reads are sequences in a pair end project that have for what ever reason, lost the corresponding pair. This actually causes issues when assembling the sequences. We are going to use a command from the Khmer package to remove any orphans

    > module load khmer
    > extract-paired-reads.py ./KBS0711_wt.trim.fastq
    > extract-paired-reads.py ./KBS0711_white.trim.fastq

The output will be the following: *KBS0711_wt.trim.fastq.pe*, *KBS0711_white.trim.fastq.pe*
       
***Re-Check Sequence Quality with *FastQC****  
Wondering what the data look like now?

***Check the files again**
    > wc -l ./KBS0711*.trim.fastq*

My Output: 
KBS0711_wt.trim.fastq - 83866040
KBS0711_wt.trim.fastq.pe - 77333632
KBS0711_white.trim.fastq - 3948392
KBS0711_white.trim.fastq.pe - 3542400

**G. Digital Normalization**  
The largest issue in genome sequencing is coverage.
We want to make sure that we have good coverage across the entire genome.
However, data shows that we get unequal coverage and though the median coverage may be 50X some areas will have up to 10 times that.
To deal with this tools have been developed to normalize the data.
This benefits the assembly in multiple ways: it equalizes the coverage across the genome, it lowers the error rate, and it greatly reduces the file sizes.
All of these end up benefiting the assembly process.

**Here, we will use a digital normalization process from the *Khmer* package.**
**This method was developed by Titus Brown and Colleagues.**
**They recommend using a three step normalization method (see website), but I've found that the best results are achieved by using only the first step.**

***Normalization***  

The following code will normalize everything to a coverage of 20; keep pairs using ‘-p’:

    > normalize-by-median.py -k 31 -C 25 -N 4 -x 2e9 -p ./KBS0711_wt.trim.fastq.pe
    >  normalize-by-median.py -k 31 -C 25 -N 4 -x 2e9 -p ./KBS0711_white.trim.fastq.pe

Options:
- -K = kmer size, 20 is the recommended setting, I used 31 (yeah, I just picked a number) [reference](http://khmer.readthedocs.org/en/v1.1/choosing-table-sizes.html) 
 - -C = Cutoff (coverage), 20 is the recommended setting, I used 25 to have a bit more
- -N = # tables, just use 4 (recommended)
- -x = table size, just use 2e9
- -p = retains paired reads

This produces a set of ‘.keep’ files for each input

Though the -p command should be retaining paired-end reads, we need to double check just in case

***Check for orphaned reads***  

The process of error trimming could have orphaned reads, so split the PE file into still-interleaved and non-interleaved reads: (This is in Khmer, make sure you load it)

    > extract-paired-reads.py ./KBS0711_wt.trim.fastq.pe.keep
    > extract-paired-reads.py ./KBS0711_white.trim.fastq.pe.keep

**H. Genome Assembly: Using Velvet Optimization**

    > module load bioperl
    > module load velvet
    > module load VelvetOptimiser

    > VelvetOptimiser.pl -s 31 -e 61 -f '-shortPaired -fastq ./KBS0711_wt.trim.fastq.pe.keep.pe' -t 8 -k 'n50*ncon' --p 'wt'\
    > VelvetOptimiser.pl -s 51 -e 111 -f '-shortPaired -fastq ./KBS0711_white.trim.fastq.pe.keep.pe' -t 8 -k 'n50*ncon' --p 'white'

This process may take a while. 
It will go through all hash values from 31 to 61 for WT and 51 to 111 for white (they were sequenced at different lengths) 
On my trials, it will choose 59 for the final assembly of wt, it has also been 57 at times. There will be some variability. 
This might take a while, but at the end you should get less than 100 contigs for wt and ~2000 for white. 

Results that I've gotten:
*Wild-Type*
Velvet hash value: 59
Total number of contigs: 86
n50: 784827
length of longest contig: 1251361
Total bases in contigs: 6093070
Number of contigs > 1k: 27

*White*
Velvet hash value: 97
Total number of contigs: 2077
n50: 161269
length of longest contig: 432110
Total bases in contigs: 6590465
Number of contigs > 1k: 70

## Assembly work flow as a torque script: ## <-- Needs Update

    > #!/bin/bash
    > #PBS -k o
    > #PBS -l nodes=1:ppn=8,vmem=50gb,walltime=10:00:00
    > #PBS -M mmuscare@indiana.edu
    > #PBS -m abe
    > #PBS -j oe
    > module load java
    > module load fastqc
    > module load bioperl
    > module load python
    > module load khmer
    > module load velvet
    > module load VelvetOptimiser
    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
    > cd /N/dc2/projects/Lennon_Sequences/2014_Janthino
    > interleave-reads.py ./KBS0711_wt_R1.fastq ./KBS0711_white_R2.fastq -o ./KBS0711_wt.fastq	
    > interleave-reads.py ./KBS0711_white_R1.fastq ./KBS0711_white_R2.fastq -o ./KBS0711_white.fastq	
    > fastq_quality_filter -Q33 -q 30 -p 50 -i ./KBS0711_wt.fastq > ./KBS0711_wt.trim.fastq
    > fastq_quality_filter -Q33 -q 30 -p 50 -i ./KBS0711_white.fastq > ./KBS0711_white.trim.fastq
    > extract-paired-reads.py ./KBS0711_wt.trim.fastq
    > extract-paired-reads.py ./KBS0711_white.trim.fastq
    > normalize-by-median.py -k 31 -C 25 -N 4 -x 2e9 -p ./KBS0711_wt.trim.fastq.pe
    > normalize-by-median.py -k 31 -C 25 -N 4 -x 2e9 -p ./KBS0711_white.trim.fastq.pe
    > extract-paired-reads.py ./KBS0711_wt.trim.fastq.pe.keep
    > extract-paired-reads.py ./KBS0711_white.trim.fastq.pe.keep
    > VelvetOptimiser.pl -s 31 -e 61 -f '-shortPaired -fastq ./KBS0711_wt.trim.fastq.pe.keep.pe' -t 8 -k 'n50*ncon' --p 'wt'\
    > VelvetOptimiser.pl -s 51 -e 111 -f '-shortPaired -fastq ./KBS0711_white.trim.fastq.pe.keep.pe' -t 8 -k 'n50*ncon' --p 'white'

## Gene Annotation ##

We are currently using the software package PROKKA to annotate our genome. 
[PROKKA](http://www.vicbioinformatics.com/prokka-manual.html))

This package integrate numerous other software into a single workflow that predicts genes and annotates both coding and noncoding genes.

This is the base command I have been using:
    > prokka 'input.fasta' --force --centre IU
    > egrep -v '^(ACCESSION|VERSION)' 'input'.gbk > 'new.gbk'

Wild-Type
    > prokka wt.contigs.59.fasta --force --centre IU
White Mutant 
    > prokka white.contigs.97.fasta --force --centre IU
    
###NCBI compliant PROKKA annotation

compliant = Force Genbank/ENA/DDJB compliance: --genes --mincontiglen 200 --centre XXX **(default OFF)**

centre = Sequencing center

outdir = NCBI assigned BioProject ID

locustag = NCBI assigned locus tag

prefix = Locus tag with "Chr1" to specify chromosome counting

Genus, species and strain are what they seem.
	
	> prokka -compliant --centre MSUGC --outdir PRJNA277721 --locustag VM94 --prefix VM94-Chr1 --genus Janthinobacterium --species sp. --strain KBS0711 contigs.fa

Thoughts: commands to add/modify
1. --mincontigslen 1000
2. --rnammer (uses RNAmmer for rRNA prediction)
3. --centre IU (sequencing senter, IU)
4. --genus Janthinobacter
5. --strain KBS0711 (KBS0711_wt, KBS0711_white)
6. --addgenes

This is still the area of exploration.





## Other details not to worry about yet ##
## Gene Predication ##

###A: Using prodigal ###

###B: Using RAST ###


## Submitting to NCBI/GenBank##

### Using tbl2asn

First, download the command-line program from the [NCBI website](http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/). For the basic submission of an annotated genome that contains more than one contig (referred to by NCBI as a [Whole Genome Shotgun Submission](http://www.ncbi.nlm.nih.gov/genbank/wgs)) you want three files: your **suffix.sbt**, **suffix.fasta**, and **suffix.tbl** files. You generate the Fasta file during assembly, the suffix.tbl during PROKKA annotation, and you can create the suffix.sbt (which contains your submission template) [here](http://www.ncbi.nlm.nih.gov/WebSub/template.cgi). This generates your **suffix.sqn** file, which is what you actually submit to GenBank. **Never manually edit your suffix.sqn file.** tbl2asn automatically finds the suffix.tbl and suffix.fasta in your current directory, so keep that in mind.

	> tbl2asn -p. -t template.sbt -M n -Z discrep -j "[organism=Janthinobacterium sp. KBS0711][strain=KBS0711]"
	
-p = Path to directory. Put a period after the "p" you want to execute the command in the current directory. 

-t = Specifies the template file, also works for full path.

-M = Master genome flag. Just use the normal flag (n) unless you have an exception to the [rules](http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/)

-Z discrep= Runs the discrepancy report **Very useful**

 -j = Allows for source qualifiers

---

## Other code that may be of interest: ##

**Remove Adapters with with *FastX-toolkit* (This is my preferred method)**

Add the FastX-toolkit to path

    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13

Remove R1 Adapters

    > fastx_clipper -v -Q 33 -a GCTCTTCCGATCT -i ./711_ATTCCT_L007_R1_001.fastq -o ./trimmed/janthino.trim.R1.fastq

Remove R2 Adapters

    > fastx_clipper -v -Q 33 -a AGATCGGAAGAGC -i ./711_ATTCCT_L007_R2_001.fastq -o ./janthino.trim.R2.fastq

**Remove Adapters with *cutadapt***

    > module load cutadapt
    > cutadapt -a GCTCTTCCGATCT ./711_ATTCCT_L007_R1_001.fastq -o ./trimmed/janthiono.trim2.R1.fastq
    > cutadapt -a AGATCGGAAGAGC ./711_ATTCCT_L007_R2_001.fastq -o ./        trimmed/janthino.trim2.R2.fastq

Remove Adapters with *cutadapt* for **paired-end reads**

*cutadapt* stores reads in temporary files so it can match reads from R1 and R2. Only khmer versions after 1.0 are compatible with the *cutadapt* outputted fastq file.

	> module load cutadapt
	> cutadapt -q 10 -a AGATCGGAAGAGC --minimum-length 20 -o ./tmp.1.fastq -p ./tmp.2.fastq ./711_ATTCCT_L007_R1_001.fastq ./711_ATTCCT_L007_R2_001.fastq
	> cutadapt -q 10 -a AGATCGGAAGAGC --minimum-length 20 -o ./711_ATTCCT_L007_R1_001_no_adapt.fastq -p ./711_ATTCCT_L007_R2_001_no_adapt.fastq ./tmp.1.fastq ./tmp.2.fastq
	> rm tmp.1.fastq tmp.2.fastq

Method | Pros | Cons 
:--------: |:-----: |:------: 
*FastX-toolkit*|Removes adapters & low quality/short seqs|Takes a while <br> Have not tested it on paired-end reads
*cutadapt* |Removes adapters and adapters with minor errors <br> We have tested it on paired-end reads and know that it works |Takes a longer time


###Removing Contamination###
It's always possible to just cut your contamination out of a contig, but that breaks the contig, potentially leaving you with a more fragmented genome. Using read mapping software we can remove whatever contamination we have prior knowledge and a nucleotide sequence of before we begin assembly. **Perform this before running khmer.** There are multiple ways to map reads, so as long as you're keeping track of you're indexed reference file and your paired-end reads (if you have them) you should be good. You can actually perform *de novo* assembly using BAM and SAM files, but if you want to use khmer you'll have to convert your BAM to FASTQ (Don't forget to take into account single- vs. paired-end reads). These are rich sets of tools that contain a lot of tools and a lot of flags, so getting what you want out of them usually invoves some digging. 

First, put the nucleotide sequences that correspond to your contamination in a text file in FASTA format and run [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) mem. This algorithm finds the local maxima, which works better when mapping reads to small indexed files (ask me why after I take my next bioinformatics course). Then you'll convert the SAM file to a BAM format using [samtools](http://samtools.sourceforge.net/samtools.shtml). Last, you want your reads in a FASTQ format for digital normalizaiton and assembly. Run [bedtools](http://bedtools.readthedocs.org/en/latest/content/tools/bamtofastq.html) to get R1 and R2 paired-end reads, or just paired-end reads.

	> bwa index contamination.fasta 
	> bwa mem -t 4 -p contamination.fasta ./KBS0711_wt_no_adapt.trim_pe.fastq > KBS0711_wt_contam.sam
	
	> samtools view -bT contamination.fasta -f 4 ./KBS0711_wt_contam.sam > ./KBS0711_wt_contam.bam
	> bedtools bamtofastq -i ./KBS0711_wt_contam.bam -fq ./KBS0711_wt_no_contam_R1.fastq -fq2 ./KBS0711_wt_no_contam_R2.fastq 

**BWA**

-t = Number of threads

-p = Paired-end mode

**SAMtools**

-b = Output as BAM format

-T = Specifies reference file

-f [*INT*]= Only outputs alignments with all bits in *INT* present in the flag field. *INT* can be in hex format of /^0x[0-9A-F]+/ [0] or just by position. For unmapped reads you want **4**. Read more about it [here](https://www.biostars.org/p/56246/),

**bedtools**

-i = BAM file

-fq =  First FASTQ output

-fq2 = Second FASTQ output (for paired-end data)