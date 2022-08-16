__Please note that we have published a more scalable successor of PopIns with a new merging step: [PopIns2](https://github.com/kehrlab/PopIns2).__

popins
======

Population-scale detection of novel-sequence insertions.


Prerequisites
-------------

* GCC version >= 4.9 (supports C++14)
* SeqAn core library, version 2.2.0 (https://github.com/seqan/seqan)
* bwa (https://github.com/lh3/bwa)
* velvet (https://github.com/dzerbino/velvet)
* samtools, version >= 1.3 (https://github.com/samtools/samtools)
* sickle (https://github.com/najoshi/sickle)

PopIns uses the 'bwa mem' alignment algorithm, thus, requires bwa version 0.7.X.
PopIns was tested with bwa 0.7.10-r789, velvet 1.2.10, samtools 1.3, and sickle 1.210.


Installation
------------

1. Download the SeqAn library. You do not need to follow the SeqAn install instructions. You only need the directory .../include/seqan with all its content (the SeqAn core library).
2. If you decide to save the seqan directory not in the popins directory, make a symbolic link to the seqan directory in the popins directory, i.e. type 'ln -s /path/to/seqan seqan' in the popins directory.
3. Install all prerequisites (bwa, velvet, samtools, and sickle).
   Compile velvet with a larger maximum k-mer length than the default if desired, e.g. MAXKMERLENGTH=63.
   A maximum k-mer length of 47 or higher is necessary for default parameters of PopIns (velvet's default is 31).
4. Set the paths to bwa, velveth, velvetg, samtools, and sickle in the file popins.config if they are not in your PATH variable.
5. Run 'make' in the popins directory.

If everything is setup correctly, this will create the binary 'popins'.

*Note*: Please check that you have used the correct version of SeqAn by running `popins genotype -h`. The last line of the help message needs to be `Seqan version: 2.2.0` .


Usage
-----

PopIns consists of seven commands: assemble, merge, contigmap, place-refalign, place-splitalign, place-finish, and genotype.
For a short description of each command and an overview of arguments and options, run

    ./popins <COMMAND> --help
    
The only input needed is a set of BAM files and the reference genome. Additionally, a sequence collection can be specified to filter out contamination.
By default, the reference genome is assumed to be present in the current directory and named `genome.fa` together with its index `genome.fa.fai`, but a different file path can be specified via the options of the place and genotype commands.

When analyzing multiple samples simultaneously, the assemble, contigmap, place-splitalign and genotype commands need to be run for each sample separately, whereas the merge and place commands need to be run only once for all samples together.

PopIns creates and uses a working directory for each sample, which is named by the sample ID (retrieved from BAM file header or user specified).
By default, the sample directories are created in the current directory. Use the --prefix option if you want them to reside in another location.

Once all steps have been run, each sample directory contains the following files:
- `POPINS_SAMPLE_INFO`: Meta information of the sample, e.g. the path to the original BAM file.
- `contigs.fa`: Contigs assembled from the reads without high-quality alignment to the reference genome.
- `insertions.vcf`: **Genotype likelihoods of the sample (GT:PL) for all predicted insertions.**
- `locations.txt`: Candidate insertion locations for the supercontigs based on reads from only this sample.
- `locations_placed.txt`: Split-read alignment results for this sample.
- `locations_unplaced.txt`: Split-read alignment input for this sample.
- `non_ref.bam`: Mates of the reads without high-quality alignments.
- `non_ref_new.bam`: Contig-aligned reads and their mates from `non_ref.bam`.
- `non_ref_new.bam.bai`: BAM index for `non_ref_new.bam`. 
- `paired.1.fastq`: First reads of those read pairs where both reads have no high-quality alignment to the reference genome.
- `paired.2.fastq`: Second reads of those read pairs where both reads have no high-quality alignment to the reference genome.
- `single.fastq`:  Reads without high-quality alignment to the reference genome but whose mates align with high-quality.

In addition to sample-specific files, a number of output files are written (by default in the current directory):
- `insertions.vcf`: Insertion positions without genotype likelihoods of the samples.
- `locations.txt`: Candidate insertion locations for the supercontigs.
- `groups.txt`: Groups of similar contigs for which only a single VCF record is written. For information purposes only.
- `skipped_contigs.fa`: Contigs that are ignored during contig merging. Optional and for information purposes only.
- `supercontigs.fa`: Contigs assembled from unaligned reads and merged from all samples.


### The assemble command

    ./popins assemble [OPTIONS] <BAM FILE>

The assemble command finds reads without high-quality alignment in the input BAM file, quality filters them using SICKLE and assembles them into contigs using VELVET.
If a reference fasta file is specified, the reads are first remapped to this reference using BWA-MEM and only reads that remain without high-quality alignment after remapping are quality-filtered and assembled.
Make sure that the reference fasta file is BWA-indexed, i.e. run `bwa index /path/to/reference.fa` before running the assemble command.


### The merge command

    ./popins merge [OPTIONS]

The merge command merges the contigs in `<prefix>/*/contigs.fa` into a single set of supercontigs.
The input contigs are first partitioned into sets of similar sequences using the SWIFT filtering algorithm, and then each set of sequences is aligned into a graph of supercontigs.


### The contigmap command

    ./popins contigmap [OPTIONS] <SAMPLE ID>

The contigmap command aligns the reads with low-quality alignments of a sample to the set of supercontigs using BWA-MEM.
The BWA output file is merged with the sample's `non_ref.bam` file into a `non_ref_new.bam` file where information about read mates is set.


### The place-refalign command

    ./popins place-refalign [OPTIONS]

This is the first of three place-* commands, which together identify insertion positions of the (super-)contigs in the reference genome and write them to a VCF file.
The place-refalign command merges contig locations in the sample directories into one file of locations and aligns prefixes/suffixes of contigs to the merged locations on the reference genome. VCF records are written if the alignment is successful. Locations of contigs that do not align to the reference genome are written to additional output files `locations_unplaced.txt` in the sample directories.


### The place-splitalign command

    ./popins place-splitalign [OPTIONS] <SAMPLE_ID>

This is the second of the three place-* commands. The place-splitalign command split-read aligns all locations in a sample's `locations_unplaced.txt` and writes the results to a file `locations_placed.txt` in the sample directory.


### The place-finish command

    ./popins place-finish [OPTIONS]
    
This is the third of the three place-* commands. The place-finish command combines the results from split-read alignment (the `locations_placed.txt` files) of all samples and appends them to the VCF output file.


### The genotype command

    ./popins genotype [OPTIONS] <SAMPLE ID>

The genotype command computes genotype likelihoods for a sample for all insertions given in the input VCF file by aligning all reads, which are mapped to the reference genome around the insertion breakpoint or to the contig, to the reference and to the alternative insertion sequence.
VCF records with the genotype likelihoods in GT:PL format for the individual are written to a file `insertions.vcf` in the sample directory.


Example
-------

    mkdir popins_example && cd popins_example/
    ln -s /path/to/hg38.fa genome.fa
    ln -s /path/to/hg38.fa genome.fa.fai
    
    ./popins assemble --sample sample1 /path/to/first_sample.bam
    ./popins assemble --sample sample2 /path/to/second_sample.bam
    ./popins assemble --sample sample3 /path/to/third_sample.bam
    
    ./popins merge
    
    ./popins contigmap sample1
    ./popins contigmap sample2
    ./popins contigmap sample3
    
    ./popins place-refalign
    ./popins place-splitalign sample1
    ./popins place-splitalign sample2
    ./popins place-splitalign sample3
    ./popins place-finish
    
    ./popins genotype sample1
    ./popins genotype sample2
    ./popins genotype sample3
    


References
----------

Kehr B., Helgadóttir A., Melsted P., Jónsson H., Helgason H., Jónasdóttir Að., Jónasdóttir As.,	Sigurðsson Á., Gylfason A., Halldórsson G. H., Kristmundsdóttir S., Þorgeirsson G., Ólafsson Í., Holm H., Þorsteinsdóttir U., Sulem P., Helgason A., Guðbjartsson D. F., Halldórsson B. V., Stefánsson K. (2017).
Diversity in non-repetitive human sequences not found in the reference genome.
[Nature Genetics,](http://rdcu.be/pDbJ) [doi:10.1038/ng.3801](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3801.html).

Kehr B., Melsted P., Halldórsson B. V. (2016).
PopIns: population-scale detection of novel sequence insertions.
[Bioinformatics, 32(7):961-967](http://bioinformatics.oxfordjournals.org/content/32/7/961.abstract).


Contact
-------

For questions and comments contact birte.kehr [ at ] bihealth.de
