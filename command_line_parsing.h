#ifndef POPINS_CLP_H_
#define POPINS_CLP_H_

#include <string>
#include <seqan/arg_parse.h>

#include "popins_utils.h"

using namespace seqan;

// ==========================================================================
// struct <Command>Options
// ==========================================================================

struct AssemblyOptions {
    CharString mappingFile;
    CharString matepairFile;
    CharString referenceFile;

    CharString prefix;
    CharString sampleID;

    unsigned kmerLength;
    CharString adapters;
    int humanSeqs;

    unsigned threads;
    CharString memory;

    AssemblyOptions () :
        matepairFile(""), referenceFile(""), prefix("."), sampleID(""),
      kmerLength(47), humanSeqs(maxValue<int>()), threads(1), memory("768M")
    {}
};

struct MergingOptions {
    CharString prefix;
    CharString outputFile;
    CharString skippedFile;
    std::fstream outputStream;
    std::fstream skippedStream;
    bool verbose;

    double errorRate;
    int minimalLength;
    unsigned qgramLength;
    int matchScore;
    int errorPenalty;
    int minScore;
    int minTipScore;

    double minEntropy;

    MergingOptions() :
        prefix("."), outputFile("supercontigs.fa"), skippedFile(""), verbose(false),
        errorRate(0.01), minimalLength(60), qgramLength(47), matchScore(1), errorPenalty(-5), minScore(90), minTipScore(30), minEntropy(0.75)
    {}
};

struct ContigMapOptions {
    CharString prefix;
    CharString sampleID;
    CharString contigFile;

    bool bestAlignment;
    int maxInsertSize;
    bool deleteNonRefNew;

    unsigned threads;
    CharString memory;

    ContigMapOptions() :
       prefix("."), sampleID(""), contigFile("supercontigs.fa"),
      bestAlignment(false), maxInsertSize(800), deleteNonRefNew(false), threads(1), memory("768M")
    {}
};

struct RefAlign_;
typedef Tag<RefAlign_> RefAlign;
struct SplitAlign_;
typedef Tag<SplitAlign_> SplitAlign;
struct SplitCombine_;
typedef Tag<SplitCombine_> SplitCombine;

template<typename TTag>
struct PlacingOptions {
    CharString prefix;
    CharString sampleID;
    CharString outFile;

    CharString locationsFile;
    CharString groupsFile;
    CharString supercontigFile;
    CharString referenceFile;

    double minLocScore;
    unsigned minAnchorReads;
    unsigned readLength;
    unsigned maxInsertSize;
    unsigned groupDist;

    PlacingOptions() :
        prefix("."), sampleID(""), outFile("insertions.vcf"), locationsFile("locations.txt"), groupsFile("groups.txt"),
        supercontigFile("supercontigs.fa"), referenceFile("genome.fa"),
        minLocScore(0.3), minAnchorReads(2), readLength(100), maxInsertSize(800), groupDist(100)
    {}
};

struct GenotypingOptions {
    CharString prefix;
    CharString sampleID;

    CharString referenceFile;
    CharString supercontigFile;
    CharString vcfFile;

    CharString genotypingModel;
    int regionWindowSize;
    bool addReadGroup;

    int maxInsertSize;
    int bpQclip;
    int minSeqLen;
    double minReadProb;
    int maxBARcount;

    int match;
    int mismatch;
    int gapOpen;
    int gapExtend;
    int minAlignScore;

    // hidden options
    bool verbose;
    bool callBoth;
    bool useReadCounts;
    bool fullOverlap;

    GenotypingOptions() :
        prefix("."), sampleID(""), referenceFile("genome.fa"), supercontigFile("supercontigs.fa"), vcfFile("insertions.vcf"),
      genotypingModel("RANDOM"), regionWindowSize(50), addReadGroup(false),
        maxInsertSize(500), bpQclip(0), minSeqLen(10), minReadProb(0.00001), maxBARcount(200),
        match(1), mismatch(-2), gapOpen(-4), gapExtend(-1), minAlignScore(55),
      verbose(false), callBoth(false), useReadCounts(false), fullOverlap(false)
    {}
};

// ==========================================================================
// Function hideOptions()
// ==========================================================================

void
setHiddenOptions(ArgumentParser & parser, bool hide, AssemblyOptions &)
{
   hideOption(parser, "matePair", hide);
   hideOption(parser, "kmerLength", hide);
}

void
setHiddenOptions(ArgumentParser & parser, bool hide, MergingOptions &)
{
   hideOption(parser, "e", hide);
   hideOption(parser, "l", hide);
   hideOption(parser, "k", hide);
   hideOption(parser, "m", hide);
   hideOption(parser, "mm", hide);
   hideOption(parser, "a", hide);
   hideOption(parser, "t", hide);
   hideOption(parser, "v", hide);
}

void
setHiddenOptions(ArgumentParser & parser, bool hide, ContigMapOptions &)
{
   hideOption(parser, "b", hide);
   hideOption(parser, "e", hide);
   hideOption(parser, "d", hide);
}

void
setHiddenOptions(ArgumentParser & parser, bool hide, PlacingOptions<RefAlign> &)
{
   hideOption(parser, "groupDist", hide);
}

void
setHiddenOptions(ArgumentParser & /*parser*/, bool /*hide*/, PlacingOptions<SplitAlign> &)
{
	// Nothing to be done.
}

void
setHiddenOptions(ArgumentParser & /*parser*/, bool /*hide*/, PlacingOptions<SplitCombine> &)
{
	// Nothing to be done.
}

void
setHiddenOptions(ArgumentParser & parser, bool hide, GenotypingOptions &)
{
   hideOption(parser, "addReadGroup", hide);

   hideOption(parser, "maxInsertSize", hide);
   hideOption(parser, "qual", hide);
   hideOption(parser, "minSeqLen", hide);
   hideOption(parser, "minReadProb", hide);
   hideOption(parser, "maxReadCount", hide);

   hideOption(parser, "match", hide);
   hideOption(parser, "mismatch", hide);
   hideOption(parser, "gapOpen", hide);
   hideOption(parser, "gapExtend", hide);
   hideOption(parser, "minScore", hide);
}

// ==========================================================================
// Function setupParser()
// ==========================================================================

void
setupParser(ArgumentParser & parser, AssemblyOptions & options)
{
    setShortDescription(parser, "Assembly of unmapped reads.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM_FILE\\fP");
    addDescription(parser, "Finds reads without high-quality alignment in the \\fIBAM FILE\\fP, quality filters them "
          "using SICKLE and assembles them into contigs using VELVET. If the option \'--reference \\fIFASTA FILE\\fP\' "
          "is set, the reads are first remapped to this reference using BwA-MEM and only reads that remain without "
          "high-quality alignment after remapping are quality-filtered and assembled.");

    // Require a bam file as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("s", "sample", "An ID for the sample.", ArgParseArgument::STRING, "SAMPLE_ID"));
    addOption(parser, ArgParseOption("mp", "matePair", "", ArgParseArgument::INPUT_FILE, "BAM FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("a", "adapters", "Enable adapter removal for Illumina reads. Default: \\fIno adapter removal\\fP.", ArgParseArgument::STRING, "STR"));
    addOption(parser, ArgParseOption("r", "reference", "Remap reads to this reference before assembly. Default: \\fIno remapping\\fP.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("f", "filter", "Treat reads aligned to all but the first INT reference sequences after remapping as high-quality aligned even if their alignment quality is low. "
          "Recommended for non-human reference sequences.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("k", "kmerLength", "The k-mer size for velvet assembly.", ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "Compute resource options");
    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for BWA and samtools sort.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "memory", "Maximum memory per thread for samtools sort; suffix K/M/G recognized.", ArgParseArgument::STRING, "STR"));

    // Set valid and default values.
    setValidValues(parser, "adapters", "HiSeq HiSeqX");
    setValidValues(parser, "reference", "fa fna fasta");
    setMinValue(parser, "threads", "1");

    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "sample", "retrieval from BAM file header");
    setDefaultValue(parser, "kmerLength", options.kmerLength);
    setDefaultValue(parser, "threads", options.threads);
    setDefaultValue(parser, "memory", options.memory);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}

void
setupParser(ArgumentParser & parser, MergingOptions & options)
{
    setShortDescription(parser, "Merging of insertion contigs into supercontigs.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP]");
    addDescription(parser, "Merges the contigs in \'<prefix>/*/contigs.fa\' into a single set of supercontigs. The "
          "input contigs are first partitioned into sets of similar sequences using the SWIFT filtering algorithm, "
          "and then each set of sequences is aligned into a graph of supercontigs.");

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of supercontigs output file.", ArgParseArgument::OUTPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("s", "skipped", "Write skipped contigs to a file. Default: \\fIdo not write skipped contigs\\fP", ArgParseArgument::OUTPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output of components."));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("y", "minEntropy", "Ignore low-complexity contigs with entropy below FLOAT. Use 0 to disable.", ArgParseArgument::DOUBLE, "FLOAT"));

    addOption(parser, ArgParseOption("e", "errRate", "Maximal error rate for SWIFT filtering.", ArgParseArgument::DOUBLE, "FLOAT"));
    addOption(parser, ArgParseOption("l", "minLength", "Minimal length for SWIFT filtering.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("k", "kmerLength", "Length of k-mers for SWIFT filtering.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "match", "Match score for Smith-Waterman alignment.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("mm", "penalty", "Error penalty for Smith-Waterman alignment.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("a", "minScore", "Minimal score for Smith-Waterman alignment.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("t", "minTipScore", "Minimal score for tips in supercontig graph.", ArgParseArgument::INTEGER, "INT"));

    // Set valid values.
    setValidValues(parser, "c", "fa fna fasta");
    setValidValues(parser, "s", "fa fna fasta");
    setMinValue(parser, "y", "0");
    setMaxValue(parser, "y", "1");
    setMinValue(parser, "e", "0");
    setMaxValue(parser, "e", "0.25");
    setMinValue(parser, "l", "3");
    setMinValue(parser, "k", "3");
    setMinValue(parser, "t", "0");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "c", options.outputFile);
    setDefaultValue(parser, "y", options.minEntropy);

    setDefaultValue(parser, "e", options.errorRate);
    setDefaultValue(parser, "l", options.minimalLength);
    setDefaultValue(parser, "k", options.qgramLength);
    setDefaultValue(parser, "m", options.matchScore);
    setDefaultValue(parser, "mm", options.errorPenalty);
    setDefaultValue(parser, "a", options.minScore);
    setDefaultValue(parser, "t", options.minTipScore);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}

void
setupParser(ArgumentParser & parser, ContigMapOptions & options)
{
    setShortDescription(parser, "Alignment of unmapped reads to assembled contigs.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fISAMPLE_ID\\fP");
    addDescription(parser, "Aligns the reads with low-quality alignments of a sample to the set of supercontigs using "
            "BWA-MEM. Merges the BWA output file with the sample's non_ref.bam file into a non_ref_new.bam file where "
    		"information about read mates is set.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "SAMPLE_ID"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of (super-)contigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("b", "best", "Do not use BWA-mem's -a option to output all alignments of a read."));
    addOption(parser, ArgParseOption("e", "maxInsertSize", "The maximum expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("d", "noNonRefNew", "Delete the non_ref_new.bam file after writing locations."));

    addSection(parser, "Compute resource options");
    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for BWA and samtools sort.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "memory", "Maximum memory per thread for samtools sort; suffix K/M/G recognized.", ArgParseArgument::STRING, "STR"));

    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "c", options.contigFile);
    setDefaultValue(parser, "best", "false");
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);
    setDefaultValue(parser, "noNonRefNew", "false");
    setDefaultValue(parser, "threads", options.threads);
    setDefaultValue(parser, "memory", options.memory);
    setMinValue(parser, "threads", "1");

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}

void
setupParser(ArgumentParser & parser, PlacingOptions<RefAlign> & options)
{
    setShortDescription(parser, "Contig placing by alignment of contig ends to reference genome.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP]");

    addDescription(parser, "This is step 1/3 of contig placing. The contig locations in the sample directories are "
    		"merged into one file of locations. Next, prefixes/suffixes of contigs are aligned to the merged locations "
    		"on the reference genome and VCF records are written if the alignment is successful. Locations of contigs "
    		"that do not align to the reference genome are written to additional output files \\fIlocations_unplaced.txt\\fP "
    		"in the sample directories. Further, groups of contigs that can be placed at the same position and whose "
    		"prefixes/suffixes align to each other are written to another output file; only a single VCF record is "
    		"written per group.");

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("l", "locations", "Name of merged locations file.", ArgParseArgument::OUTPUT_FILE, "FILE"));
    addOption(parser, ArgParseOption("i", "insertions", "Name of VCF output file.", ArgParseArgument::OUTPUT_FILE, "VCF_FILE"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of supercontigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("g", "groups", "Name of file containing groups of contigs that can be placed at the same position and whose prefixes/suffixes align.", ArgParseArgument::OUTPUT_FILE, "FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("", "minScore", "Minimal anchoring score for a location.", ArgParseArgument::DOUBLE, "FLOAT"));
    addOption(parser, ArgParseOption("", "minReads", "Minimal number of anchoring read pairs for a location.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "maxInsertSize", "The maximum expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "readLength", "The length of the reads.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "groupDist", "Minimal distance between groups of locations.", ArgParseArgument::INTEGER, "INT"));

    // Set valid values.
    setMinValue(parser, "minScore", "0");
    setMaxValue(parser, "minScore", "1");
    setValidValues(parser, "contigs", "fa fna fasta");
    setValidValues(parser, "reference", "fa fna fasta");
    setValidValues(parser, "insertions", "vcf");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "locations", options.locationsFile);
    setDefaultValue(parser, "insertions", options.outFile);
    setDefaultValue(parser, "contigs", options.supercontigFile);
    setDefaultValue(parser, "reference", options.referenceFile);
    setDefaultValue(parser, "groups", options.groupsFile);

    setDefaultValue(parser, "minScore", options.minLocScore);
    setDefaultValue(parser, "minReads", options.minAnchorReads);
    setDefaultValue(parser, "groupDist", options.groupDist);
    setDefaultValue(parser, "readLength", options.readLength);
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}

void
setupParser(ArgumentParser & parser, PlacingOptions<SplitAlign> & options)
{
    setShortDescription(parser, "Contig placing by split-read alignment.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fISAMPLE_ID\\fP");

    addDescription(parser, "This is step 2/3 of contig placing. All locations in a sample's "
    		"\\fIlocations_unplaced.txt\\fP are split-read aligned and the results are written to a file "
    		"\\fIlocations_placed.txt\\fP in the sample directory.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "SAMPLE_ID"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of supercontigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("", "maxInsertSize", "The maximum expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "readLength", "The length of the reads.", ArgParseArgument::INTEGER, "INT"));

    // Set valid values.
    setValidValues(parser, "contigs", "fa fna fasta");
    setValidValues(parser, "reference", "fa fna fasta");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "contigs", options.supercontigFile);
    setDefaultValue(parser, "reference", options.referenceFile);

    setDefaultValue(parser, "readLength", options.readLength);
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}

void
setupParser(ArgumentParser & parser, PlacingOptions<SplitCombine> & options)
{
    setShortDescription(parser, "Combining breakpoint positions from split-read alignment.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP]");

    addDescription(parser, "This is step 3/3 of contig placing. The results from split-read alignment (the "
    		"\\fIlocations_placed.txt\\fP files) of all samples are combined and appended to the VCF output file.");

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("i", "insertions", "Name of VCF output file.", ArgParseArgument::OUTPUT_FILE, "VCF_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    // Set valid values.
    setValidValues(parser, "reference", "fa fna fasta");
    setValidValues(parser, "insertions", "vcf");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "insertions", options.outFile);
    setDefaultValue(parser, "reference", options.referenceFile);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}

void
setupParser(ArgumentParser & parser, GenotypingOptions & options)
{
    setShortDescription(parser, "Genotyping a sample for insertions.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fISAMPLE_ID\\fP");
    addDescription(parser, "Computes genotype likelihoods for a sample for all insertions given in the input VCF file "
    		"by aligning all reads, which are mapped to the reference genome around the insertion breakpoint or to the "
    		"contig, to the reference and to the alternative insertion sequence. VCF records with the genotype "
    		"likelihoods in GT:PL format for the individual are written to a file \\fIinsertions.vcf\\fP in the sample "
            "directory.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "SAMPLE_ID"));

    // Setup the options.
    addSection(parser, "Input/output options");
    addOption(parser, ArgParseOption("p", "prefix", "Path to the sample directories.", ArgParseArgument::STRING, "PATH"));
    addOption(parser, ArgParseOption("i", "insertions", "Name of VCF input file.", ArgParseArgument::INPUT_FILE, "VCF_FILE"));
    addOption(parser, ArgParseOption("c", "contigs", "Name of supercontigs file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));
    addOption(parser, ArgParseOption("r", "reference", "Name of reference genome file.", ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

    addSection(parser, "Algorithm options");
    addOption(parser, ArgParseOption("m", "model", "Model used for genotyping", ArgParseArgument::STRING, "GENOTYPING_MODEL"));
    addOption(parser, ArgParseOption("w", "window", "Region window size.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("rg", "addReadGroup", "Add read group."));

    addSection(parser, "Read(-pair) options");
    addOption(parser, ArgParseOption("", "maxInsertSize", "Maximum read pair insert size.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "qual", "Quality score threshold for read trimming.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "minSeqLen", "Minimum read length after trimming.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "minReadProb", "Minimum read probability.", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("", "maxReadCount", "Maximum number of reads to consider in region window.", ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "Alignment options");
    addOption(parser, ArgParseOption("", "match", "Cost of match.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "mismatch", "Cost of mismatch.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "gapOpen", "Cost of gap open.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "gapExtend", "Cost of gap extend.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("", "minScore", "Minimum alignment score.", ArgParseArgument::INTEGER, "INT"));

    // Misc hidden options.
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
   hideOption(parser, "verbose", true);
    addOption(parser, ArgParseOption("", "callBoth", "Call both models."));
   hideOption(parser, "callBoth", true);
    addOption(parser, ArgParseOption("", "readCounts", "Use read counts."));
   hideOption(parser, "readCounts", true);
    addOption(parser, ArgParseOption("", "fullOverlap", "Full overlap of read."));
   hideOption(parser, "fullOverlap", true);

    // Set valid values.
    setValidValues(parser, "contigs", "fa fna fasta");
    setValidValues(parser, "reference", "fa fna fasta");
    setValidValues(parser, "insertions", "vcf");
    setValidValues(parser, "model", "DUP RANDOM");

    // Set default values.
    setDefaultValue(parser, "prefix", "\'.\'");
    setDefaultValue(parser, "insertions", options.vcfFile);
    setDefaultValue(parser, "contigs", options.supercontigFile);
    setDefaultValue(parser, "reference", options.referenceFile);

    setDefaultValue(parser, "model", options.genotypingModel);
    setDefaultValue(parser, "window", options.regionWindowSize);
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);
    setDefaultValue(parser, "qual", options.bpQclip);
    setDefaultValue(parser, "minReadProb", options.minReadProb);
    setDefaultValue(parser, "maxReadCount", options.maxBARcount);

    setDefaultValue(parser, "match", options.match);
    setDefaultValue(parser, "mismatch", options.mismatch);
    setDefaultValue(parser, "gapOpen", options.gapOpen);
    setDefaultValue(parser, "gapExtend", options.gapExtend);
    setDefaultValue(parser, "minScore", options.minAlignScore);

    // Hide some options from default help.
    setHiddenOptions(parser, true, options);
}

// ==========================================================================
// Function getOptionValues()
// ==========================================================================

void
getOptionValues(AssemblyOptions & options, ArgumentParser const & parser)
{
    getArgumentValue(options.mappingFile, parser, 0);

    if (isSet(parser, "prefix"))
       getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "sample"))
       getOptionValue(options.sampleID, parser, "sample");
    if (isSet(parser, "matePair"))
        getOptionValue(options.matepairFile, parser, "matePair");
    if (isSet(parser, "adapters"))
        getOptionValue(options.adapters, parser, "adapters");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");
    if (isSet(parser, "filter"))
        getOptionValue(options.humanSeqs, parser, "filter");
    if (isSet(parser, "kmerLength"))
        getOptionValue(options.kmerLength, parser, "kmerLength");
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");
}

void
getOptionValues(MergingOptions & options, ArgumentParser & parser)
{
    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "contigs"))
        getOptionValue(options.outputFile, parser, "contigs");
    if (isSet(parser, "skipped"))
        getOptionValue(options.skippedFile, parser, "skipped");
    if (isSet(parser, "verbose"))
        options.verbose = true;

    if (isSet(parser, "minEntropy"))
        getOptionValue(options.minEntropy, parser, "minEntropy");

    if (isSet(parser, "errRate"))
        getOptionValue(options.errorRate, parser, "errRate");
    if (isSet(parser, "minLength"))
        getOptionValue(options.minimalLength, parser, "minLength");
    if (isSet(parser, "kmerLength"))
        getOptionValue(options.qgramLength, parser, "kmerLength");
    if (isSet(parser, "match"))
        getOptionValue(options.matchScore, parser, "match");
    if (isSet(parser, "minScore"))
        getOptionValue(options.minScore, parser, "minScore");
    if (isSet(parser, "penalty"))
        getOptionValue(options.errorPenalty, parser, "penalty");
    if (isSet(parser, "minTipScore"))
        getOptionValue(options.minTipScore, parser, "minTipScore");
}

void
getOptionValues(ContigMapOptions & options, ArgumentParser & parser)
{
    getArgumentValue(options.sampleID, parser, 0);

    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "contigs"))
        getOptionValue(options.contigFile, parser, "contigs");
    if (isSet(parser, "best"))
        options.bestAlignment = true;
    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "noNonRefNew"))
        options.deleteNonRefNew = true;
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");
}

void
getOptionValues(PlacingOptions<RefAlign> & options, ArgumentParser & parser)
{
    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "locations"))
        getOptionValue(options.locationsFile, parser, "locations");
    if (isSet(parser, "contigs"))
        getOptionValue(options.supercontigFile, parser, "contigs");
    if (isSet(parser, "insertions"))
        getOptionValue(options.outFile, parser, "insertions");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");
    if (isSet(parser, "groups"))
        getOptionValue(options.groupsFile, parser, "groups");

    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "readLength"))
        getOptionValue(options.readLength, parser, "readLength");
    if (isSet(parser, "minScore"))
        getOptionValue(options.minLocScore, parser, "minScore");
    if (isSet(parser, "minReads"))
        getOptionValue(options.minAnchorReads, parser, "minReads");
    if (isSet(parser, "groupDist"))
        getOptionValue(options.groupDist, parser, "groupDist");
}

void
getOptionValues(PlacingOptions<SplitAlign> & options, ArgumentParser & parser)
{
    getArgumentValue(options.sampleID, parser, 0);

    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "contigs"))
        getOptionValue(options.supercontigFile, parser, "contigs");
    if (isSet(parser, "insertions"))
        getOptionValue(options.outFile, parser, "insertions");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");

    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "readLength"))
        getOptionValue(options.readLength, parser, "readLength");
}

void
getOptionValues(PlacingOptions<SplitCombine> & options, ArgumentParser & parser)
{
    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "insertions"))
        getOptionValue(options.outFile, parser, "insertions");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");
}

void
getOptionValues(GenotypingOptions & options, ArgumentParser & parser)
{
    getArgumentValue(options.sampleID, parser, 0);

    if (isSet(parser, "prefix"))
        getOptionValue(options.prefix, parser, "prefix");
    if (isSet(parser, "contigs"))
        getOptionValue(options.supercontigFile, parser, "contigs");
    if (isSet(parser, "insertions"))
        getOptionValue(options.vcfFile, parser, "insertions");
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");

    if (isSet(parser, "model"))
        getOptionValue(options.genotypingModel, parser, "model");
    if(isSet(parser, "window"))
        getOptionValue( options.regionWindowSize, parser, "window");
    options.addReadGroup = isSet(parser, "addReadGroup");

    if (isSet(parser, "match"))
        getOptionValue(options.match, parser, "match");
    if (isSet(parser, "mismatch"))
        getOptionValue(options.mismatch, parser, "mismatch");
    if (isSet(parser, "gapOpen"))
        getOptionValue(options.gapOpen, parser, "gapOpen");
    if (isSet(parser, "gapExtend"))
        getOptionValue(options.gapExtend, parser, "gapExtend");
    if (isSet(parser, "minScore"))
        getOptionValue(options.minAlignScore, parser, "minScore");

    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "minReadProb"))
        getOptionValue(options.minReadProb, parser, "minReadProb");
    if (isSet(parser, "maxReadCount"))
        getOptionValue(options.maxBARcount, parser, "maxReadCount");
    if (isSet(parser, "qual"))
        getOptionValue(options.bpQclip, parser, "qual");
    if (isSet(parser, "minSeqLen"))
        getOptionValue(options.minSeqLen, parser, "minSeqLen");

    options.verbose = isSet(parser, "verbose");
    options.callBoth = isSet(parser, "callBoth");
    options.fullOverlap = isSet(parser, "fullOverlap");
    options.useReadCounts = isSet(parser, "readCounts");
}

// ==========================================================================
// Function parseCommandLine()
// ==========================================================================

template<typename TOptions>
ArgumentParser::ParseResult
parseCommandLine(TOptions & options, int argc, char const ** argv)
{
    // Concatenate program name from name and command.
    CharString prog_name = argv[0];
    prog_name += " ";
    prog_name += argv[1];

    ++argv;
    --argc;

    // Setup the parser.
    ArgumentParser parser(toCString(prog_name));
    addOption(parser, ArgParseOption("H", "fullHelp", "Display the help message with the full list of options."));
    setupParser(parser, options);

    // Parse the command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Check if full help option is set.
    if (isSet(parser, "fullHelp"))
    {
       setHiddenOptions(parser, false, options);
       printHelp(parser);
        return ArgumentParser::PARSE_HELP;
    }

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Collect the option values.
    getOptionValues(options, parser);

    return res;
}

// ==========================================================================

#endif // #ifndef POPINS_CLP_H_
