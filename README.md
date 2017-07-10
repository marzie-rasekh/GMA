# Genome Mappaility Analyzer (GMA) #
corrected and optimized version of GMA

## The algorithm ##

The idea is to measure the mappability of each locus of a genome. 
To do so, we simulate a read of given length at each position of the genome.
The reads are mapped back to the reference genome and compared against the real position they originated from. 
We count the number of *missing* and *extra* reads:  
**Missing**: read that was created but did not map back to where it originated from (it went missing, loss).  
**Extra**: read that did not belong where it was mapped to (it is extra coverage, gain).

Then we calculate the corrected GMA score for each position *p* as :

``` GMA[i] = (READ_LENGTH - MISSING[p]) / (READ_LENGTH + EXTRA[p]) ```

## Input ##
My code takes a reference genome as input. The fasta files can be created in the fasta directory by:

```awk '{if(substr($1,1,1)==">") {F=substr($1,2);} print > "fasta/" F ".fasta"}' $REFERENCE ```

Note that the reference should be prepared for bwa aligning. See http://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk

## Parameters ##
The algorithm sets multiple parameters. Edit them in the GMA.sh file. 

```
PROC=8	                       # number of threads used for bwa mem mapping
SNP=0.01                       # error rate for the reads (illumina single ended)
COVERAGE=1                     # number of reads to produce per position
READ_LENGTH=101                # %in% O(coverage*|sequence|*|read|) , at read length 100 coverage 1 for human chr1 (250Mbp) makes 24 fastq files of size 2.5G each,
	                       # bam file sizes are neglectible (each ~200-250M)
PADDING=100                    # this is the threshold to consider a read as a mismap, 
                               # any read mapping farther than the PADDING from it's origin is considered a mismap
FASTA_DIR=fasta                # to find $sequence.fasta file, notice the format
FASTQ_DIR=fastq                # to save/read the fastq files, suffix fq
ALIGN_DIR=aln                  # to save/read alignment files
FILE_BREAK=10000000            # setting this line break too low will add over-head, too high will make it slow
	                       # keep in mind that for each read in the fastq, there will be 4 lines
DELIMINATOR="_"                # this is set in the src/GMA_test.cpp code, %TODO fix this

```
Note that the directories will be emptied prior to running.

The code will take each chromosome and in a spearate process simulate reads for it. Each fastq file will be of length FILE_BREAK. 
Setting this value too high will create larger files making the IO slower. 
Setting it too low makes too many files and makes retriving files slower, and also adds overhead to the BWA call.

## Output ##
The output of the code will the produced in the ALIGN_DIR folder. For each sequence header in the reference you will see:
1. sequence.gma
    This file holds all edges of mismaped reads. 
    There are 5 columns: the originating_chromosome, originating_basepair, mapped_chromosome, mapped_basepair, mapped_strand. 
    Since all our reads are from the forward strand, the mapped strand might be inetersting. 
    A star indicates that the read was not mapped.
    This file can be used for graph purposes.
2. sequence.missing.bed
    Consists of 5 columns in bed format: originating_chromosome, originating_basepair, originating_end_position, mapped_strand, number_of_mapped
    Since a read can be mapped to multiple locations. This file is made from the sequence.gma file.
3. sequence.extra.bed
    Consists of 5 columns in bed format: mapped_chromosome, mapped_basepair, mapped_end_position, mapped_strand
    This file is made from the sequence.gma file.
4. Finally for each line in a bed file, one can simply get the ```bedtools coverage -count -f 1.00 ``` to get the number of missing and extra reads passing that locus/region.
    The result can be used to calculate the final GMA score. (see GMA_TR_analysis.sh)


## How to run ##

This code was written for the BU SCC using qsubs. To run on your local cluster remove the qsubs. 
You can use & to call subprocesses and wait on the child process IDs. 

## References ##

This work is based on :

Hayan Lee, Michael C. Schatz; **Genomic dark matter**: the reliability of short read mapping illustrated by the genome mappability score. Bioinformatics 2012; 28 (16): 2097-2105. doi: 10.1093/bioinformatics/bts330

You can download their code at :
https://sourceforge.net/projects/gma-bio/

I could not run their code. The major differences :
- They use a window of read length and slide it, so the resolution of mappability is blocks of read length
- They create reads then map back. I saw they are mapping back to only the same chromosome, which is an upper bound of mappability. I have corrected that to map to the whole reference.
- They use bwa aln, which is a no-gap aligner and takes longer. I am using the standard bwa mem (Li, Heng. "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM." arXiv preprint arXiv:1303.3997 (2013).)
- They do not count the "extra" reads mapping in the location.
- Their code is super slow because they run all steps for each locus at a time. They try to use hadoop to speed things up. But they have too much hard coded settings. 
- Their code takes a week or two according to themselves. Mine takes about 4 hours for a genome.
- They do not count for the extra reads.


