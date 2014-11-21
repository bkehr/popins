popins
======

Population-scale detection of novel sequence insertions.


Prerequisites
-------------

* SeqAn core library, version 1.4.1 (https://github.com/seqan/seqan)
* bwa (https://github.com/lh3/bwa)
* velvet (https://github.com/dzerbino/velvet)
* samtools (https://github.com/samtools/samtools)
* sickle (https://github.com/najoshi/sickle)

PopIns uses the 'bwa mem' alignment algorithm, thus, requires bwa version 0.7.X.
PopIns was tested with SeqAn r14774, bwa 0.7.10-r789, velvet 1.2.10, samtools 0.1.13 (r926:134), and sickle 1.210.


Installation
------------

1) Download the SeqAn core library and the vcf_io module from the SeqAn extras library.
   You do *not* need to follow the SeqAn install instructions.
   You only need the directory core/include/seqan of SeqAn with all its content.
   In addition, copy the header file vcf_io.h and folder vcf_io/ from extras/include/seqan into this directory.
2) Install all other prerequisites (bwa, velvet, samtools, and sickle).
   Compile velvet with a larger k-mer length than the default if desired, e.g. 63 (necessary for default parameters of PopIns).
3) Set the path to the SeqAn core library by editing the file popins.config.
   Also set the paths to bwa, velveth/velvetg, samtools, and sickle if they are not in your PATH variable.
4) Run 'make' in the popins directory.

If everything is setup correclty, this will create the binary 'popins'.


Usage
-----

PopIns consists of five commands: assemble, merge, contigmap, place, and genotype.
For a short description of each command and an overview of arguments and options, run

  ./popins <COMMAND> --help

When analyzing multiple samples simultaneously, the assemble, contigmap, and genotype commands need to be run for each sample separately, whereas the merge and place commands need to be run only once with input from all samples.
PopIns creates and uses a working directory for each sample, which should be specified with the -d option of the assemble and contigmap commands.


### The assemble command

Usage: ./popins assemble [OPTIONS] <BAM FILE>

The assemble command finds the unmapped reads in a bam files and assembles them using velvet.
If a reference fasta file is specified, the unmapped reads will be remapped to this reference before assembly using bwa-mem.
Only reads that remain unmapped in the remapping step are further processed, i.e. they are quality filtered using sickle and passed to assembly with velvet.


### The merge command

Usage: ./popins merge [OPTIONS] <FA FILE 1> ... <FA FILE N>

The merge command merges all sequences given in the fasta files into a single set of supercontigs.
The algorithm first partitions the sequences into sets of similar sequences using the SWIFT filtering approach, and then aligns each set of contigs into a graph of supercontigs.


### The contigmap command

Usage: ./popins contigmap [OPTIONS] <BAM FILE> <FA FILE>

The contigmap command aligns the unmapped reads found in fastq files in the working directory to a set of contigs specified in the fasta file using bwa-mem.
Subsequently, it merges the bwa output file with the file non_ref.bam in the working directory and completes the read mate's information in all bam records.


### The place command

Usage: ./popins place [OPTIONS] <CONTIG FA FILE> <REF FA FILE> <BAM FILE 1> ... <BAM FILE N>

The place command finds the positions of (super-)contigs in the reference genome.
If a file with locations does not already exist, it identifies approximate locations based on anchoring read pairs found in the bam files.
If bam files with all reads of the individuals are specified, it determines exact positions of insertions from split read alignments for each contig end.
Both steps can be run separately or in a single program call.
The split alignment can be done in batches (e.g. 100 locations per batch) if the approximate locations have been computed before.
It outputs a vcf and a fa record for each identified position.


### The genotype command

Usage: ./popins genotype [OPTIONS] <FA FILE> <BAM FILE> <FA FILE ALT> <BAM FILE ALT> <VCF FILE>

The genotype command takes as input a fasta file of the reference genome, a bam file of a single individual, the fasta file with the supercontigs, the bam file of contig mapped and unmapped reads (<WD>/non_ref.bam), and the VCF file with all predicted insertion positions.
It computes genotype likelihoods by aligning all reads from each insertion location and contig to the reference and to the alternative insertion sequence.
It outputs VCF records with the genotype likelihoods in GT:PL format for the individual to std::out.


References
----------

Kehr B., Melsted P., Jónasdóttir A., Jónasdóttir A., Sigurðsson A., Gylfason A., Guðbjartsson D., Halldórsson B. V., Stefánsson K.
Detecting novel sequence insertions in 3000 individuals from short read sequencing data. (Abstract/Program #38).
Presented at the 64th Annual Meeting of The American Society of Human Genetics, October 20, 2014, San Diego, CA.

Kehr B., Melsted P., Halldórsson B. V.
PopIns: population-scale detection of novel sequence insertions.
Submitted.


Contact
-------

For questions and comments contact birte.kehr [ at ] decode.is
