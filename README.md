popins
======

Population-scale detection of novel-sequence insertions.


Prerequisites
-------------

* bwa (https://github.com/lh3/bwa)
* velvet (https://github.com/dzerbino/velvet)
* samtools (https://github.com/samtools/samtools)
* sickle (https://github.com/najoshi/sickle)

PopIns uses the 'bwa mem' alignment algorithm, thus, requires bwa version 0.7.X.
PopIns was tested with bwa 0.7.10-r789, velvet 1.2.10, samtools 1.3, and sickle 1.210.


Installation
------------

1. Install all prerequisites (bwa, velvet, samtools, and sickle).
   Compile velvet with a larger maximum k-mer length than the default if desired, e.g. MAXKMERLENGTH=63.
   A maximum k-mer length of 47 or higher is necessary for default parameters of PopIns (velvet's default is 31).
2. Set the paths to bwa, velveth, velvetg, samtools, and sickle in the file popins.config if they are not in your PATH variable.
3. Run 'make' in the popins directory.

If everything is setup correctly, this will create the binary 'popins'.


Usage
-----

PopIns consists of five commands: assemble, merge, contigmap, place, and genotype.
For a short description of each command and an overview of arguments and options, run

    ./popins <COMMAND> --help
    
The only input needed is a set of BAM files (PopIns was only tested on BAM files created with BWA-MEM) and the reference genome.
By default, the reference genome is assumed to be present in the current directory and named `genome.fa`, but a different file path can be specified in those commands that need it.

When analyzing multiple samples simultaneously, the assemble, contigmap, and genotype commands as well as a substep of the place command need to be run for each sample separately, whereas the merge command and the remaining steps of the place commands need to be run only once for all samples together.

PopIns creates and uses a working directory for each sample, which is named by the sample ID (retrived from BAM file header or specify).
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


### The merge command

    ./popins merge [OPTIONS]

The merge command merges the contigs in `<prefix>/*/contigs.fa` into a single set of supercontigs.
The input contigs are first partitioned into sets of similar sequences using the SWIFT filtering algorithm, and then each set of sequences is aligned into a graph of supercontigs.


### The contigmap command

    ./popins contigmap [OPTIONS] <SAMPLE ID>

The contigmap command aligns the reads with low-quality alignments of a sample to the set of supercontigs using BWA-MEM.
The BWA output file is merged with the sample's `non_ref.bam` file into a `non_ref_new.bam` file where information about read mates is set.


### The place command

    ./popins place [OPTIONS]

The place command identifies insertion positions of the (super-)contigs in the reference genome and writes them to a VCF file.
The placing consists of four steps.
Only the third step needs to be run per sample unsing the option `--sample <SAMPLE ID>`.

Step 1: The contig locations in the sample directories are merged into one file of locations.

Step 2: Prefixes/suffixes of contigs are aligned to the merged locations on the reference genome and VCF records are written if the alignment is successful.
Locations of contigs that do not align to the reference genome are written to additional output files `locations_unplaced.txt` in the sample directories.

Step 3: All locations in a sample's `locations_unplaced.txt` are split-read aligned and the results are written to a file `locations_placed.txt` in the sample directory.

Step 4: The results from split-read alignment (the `locations_placed.txt` files) of all samples are combined and appended to the VCF output file.


### The genotype command

    ./popins genotype [OPTIONS] <SAMPLE ID>

The genotype command computes genotype likelihoods for a sample for all insertions given in the input VCF file by aligning all reads, which are mapped to the reference genome around the insertion breakpoint or to the contig, to the reference and to the alternative insertion sequence.
VCF records with the genotype likelihoods in GT:PL format for the individual are written to a file `insertions.vcf` in the sample directory.


Example
-------

    mkdir popins_example && cd popins_example/
    ln -s /path/to/hg38.fa genome.fa
    
    ./popins assemble --sample sample1 /path/to/first_sample.bam
    ./popins assemble --sample sample2 /path/to/second_sample.bam
    ./popins assemble --sample sample3 /path/to/third_sample.bam
    
    ./popins merge
    
    ./popins contigmap sample1
    ./popins contigmap sample2
    ./popins contigmap sample3
    
    ./popins place              # locations.txt does not exist --> runs substeps 1 and 2
    ./popins place -s sample1   # a sample ID is given --> runs substep 3
    ./popins place -s sample2
    ./popins place -s sample3
    ./popins place              # locations.txt now exists --> runs substeps 4
    
    ./popins genotype sample1
    ./popins genotype sample2
    ./popins genotype sample3
    


References
----------

Kehr B., Melsted P., Halldórsson B. V. (2016).
PopIns: population-scale detection of novel sequence insertions.
[Bioinformatics, 32(7):961-967.](http://bioinformatics.oxfordjournals.org/content/32/7/961.abstract)

Kehr B., Melsted P., Jónasdóttir A., Jónasdóttir A., Sigurðsson A., Gylfason A., Guðbjartsson D., Halldórsson B. V., Stefánsson K. (2014).
Detecting novel sequence insertions in 3000 individuals from short read sequencing data. (Abstract/Program #38).
Presented at the 64th Annual Meeting of The American Society of Human Genetics, October 20, 2014, San Diego, CA.


Contact
-------

For questions and comments contact birte.kehr [ at ] decode.is
