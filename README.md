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
PopIns was tested with bwa 0.7.10-r789, velvet 1.2.10, samtools 1.3, and sickle 1.210.


Installation
------------

1. Download the SeqAn library. You do *not* need to follow the SeqAn install instructions.
   You only need the directory .../include/seqan of the SeqAn core library with all its content.
   Look for the header file vcf_io.h and folder vcf_io/. If it is not present, copy it from extras/include/seqan into your seqan directory.
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

    ./popins place [OPTIONS] <LOCATION FILE> <OUTPUT FILE>

The place command identifies the positions of (super-)contigs in the reference genome and writes them to a VCF file.
VCF records reference contigs and contig positions in the (super-)contigs file.

The place command consists of four steps that can be run in one program call or in separate calls.

To run all four steps together, the options -l, -b, -c, and -r need to be specified and the <OUTPUT FILE>'s ending needs to be 'vcf'. The <LOCATION FILE> has to list the location files for all samples.

When running the steps separately, the specified parameters determine which step is being run:

1. First the contig locations determined for all samples (files listed in <LOCATION FILE>) need to be merged into one set of locations (written to <OUTPUT FILE>).
2. Then, prefixes/suffixes of all contigs (specify with -c option) are aligned to these locations (specify as <LOCATION FILE>) in the reference genome (specify with -r option) and VCF records are written to <OUTPUF FILE> if the alignment is successful.
   The <OUTPUT FILE>'s ending has to be 'vcf'.
3. Next, contigs (specify -c option) that do not align to the reference genome (specify -r option) are passed on to split-read alignment.
   This step is run by sample.
   The sample's original BAM file needs to be specified.
   The program arguments, the <LOCATION FILE> and <OUTPUT FILE>, need to be the sample's locs\_unplaced.txt file and a locs\_placed.txt file.
4. Finally, the results from split-read alignment (the locs_placed.txt files) of all samples (input files listed in <LOCATION FILE>) are being combined and appended the <OUTPUF FILE> (file ending has to be 'vcf').


### The genotype command

    ./popins genotype [OPTIONS] <FA FILE> <BAM FILE> <FA FILE ALT> <BAM FILE ALT> <VCF FILE>

The genotype command takes as input a fasta file of the reference genome, a bam file of a single individual, the fasta file with the supercontigs, the bam file of contig mapped and unmapped reads (<WD>/non_ref.bam), and the VCF file with all predicted insertion positions.
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
