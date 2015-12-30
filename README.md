popins
======

Population-scale detection of novel sequence insertions.


Prerequisites
-------------

* SeqAn core library, version 1.4.2 (https://github.com/seqan/seqan)
* bwa (https://github.com/lh3/bwa)
* velvet (https://github.com/dzerbino/velvet)
* samtools (https://github.com/samtools/samtools)
* sickle (https://github.com/najoshi/sickle)

PopIns uses the 'bwa mem' alignment algorithm, thus, requires bwa version 0.7.X.
PopIns was tested with bwa 0.7.10-r789, velvet 1.2.10, samtools 1.0 and 0.1.13 (r926:134), and sickle 1.210.


Installation
------------

1. Download the SeqAn library. You do *not* need to follow the SeqAn install instructions.
   You only need the directory .../include/seqan of the SeqAn core library with all its content.
   If it is not present, copy the header file vcf_io.h and folder vcf_io/ from extras/include/seqan into this directory.
2. Install all other prerequisites (bwa, velvet, samtools, and sickle).
   Compile velvet with a larger maximum k-mer length than the default if desired, e.g. MAXKMERLENGTH=63.
   A maximum k-mer length of 47 or higher is necessary for default parameters of PopIns (velvet's default is 31).
3. Set the path to the SeqAn core library by editing the file popins.config.
   Also set the paths to bwa, velveth/velvetg, samtools, and sickle if they are not in your PATH variable.
4. Run 'make' in the popins directory.

If everything is setup correctly, this will create the binary 'popins'.


Usage
-----

PopIns consists of five commands: assemble, merge, contigmap, place, and genotype.
For a short description of each command and an overview of arguments and options, run

    ./popins <COMMAND> --help

When analyzing multiple samples simultaneously, the assemble, contigmap, and genotype commands need to be run for each sample separately, whereas the merge and place commands need to be run only once with input from all samples.
PopIns creates and uses a working directory for each sample, which should be specified with the -d option of the assemble and contigmap commands.

### The assemble command

    ./popins assemble [OPTIONS] <BAM FILE>

The assemble command finds the unmapped reads in a bam files and assembles them using velvet.
If a reference fasta file is specified, the unmapped reads will be remapped to this reference before assembly using bwa-mem.
Only reads that remain unmapped in the remapping step are further processed, i.e. they are quality filtered using sickle and passed to assembly with velvet.


### The merge command

    ./popins merge [OPTIONS] <FA FILE 1> ... <FA FILE N>

The merge command merges all sequences given in the fasta files into a single set of supercontigs.
The algorithm first partitions the sequences into sets of similar sequences using the SWIFT filtering approach, and then aligns each set of contigs into a graph of supercontigs.


### The contigmap command

    ./popins contigmap [OPTIONS] <FA FILE>

The contigmap command aligns the unmapped reads found in fastq files in the working directory to a set of contigs specified in the fasta file using bwa-mem.
Subsequently, it merges the bwa output file with the file non_ref.bam in the working directory and completes the read mate's information in all bam records.
Finally, it determines approximate insertion locations for contigs with anchoring read pairs.


### The place command

    ./popins place [OPTIONS] <CONTIG FA FILE> <REF FA FILE>

The place command identifies the positions of (super-)contigs in the reference genome.
If a file with merged locations (-ml option) does not already exist, it requires the -l option to be set and merges locations files.
If bam files with all reads of the individuals are specified, it determines exact positions of insertions from split read alignments for each contig end.
Both steps can be run separately or in a single program call.
The split alignment can be done in batches by genomic region (e.g. only locations in chr3:40000000-41000000) if the locations files have been merged before.
It outputs a vcf record for each identified position, which references contigs and contig positions in the (super-)contigs file.


### The genotype command

    ./popins genotype [OPTIONS] <FA FILE> <BAM FILE> <FA FILE ALT> <BAM FILE ALT> <VCF FILE>

The genotype command takes as input a fasta file of the reference genome, a bam file of a single individual, the fasta file with the supercontigs, the bam file of contig mapped and unmapped reads (&lt;WD&gt;/non_ref.bam), and the VCF file with all predicted insertion positions.
It computes genotype likelihoods by aligning all reads from each insertion location and contig to the reference and to the alternative insertion sequence.
It outputs VCF records with the genotype likelihoods in GT:PL format for the individual to std::out.


References
----------

Kehr B., Melsted P., Halldórsson B. V. (2015).
PopIns: population-scale detection of novel sequence insertions.
Bioinformatics, btv273.

Kehr B., Melsted P., Jónasdóttir A., Jónasdóttir A., Sigurðsson A., Gylfason A., Guðbjartsson D., Halldórsson B. V., Stefánsson K. (2014).
Detecting novel sequence insertions in 3000 individuals from short read sequencing data. (Abstract/Program #38).
Presented at the 64th Annual Meeting of The American Society of Human Genetics, October 20, 2014, San Diego, CA.


Contact
-------

For questions and comments contact birte.kehr [ at ] decode.is
