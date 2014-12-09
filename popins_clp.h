#include <seqan/arg_parse.h>

#ifndef POPINS_CLP_H_
#define POPINS_CLP_H_

using namespace seqan;

// ==========================================================================

// Returns true if file exists, otherwise false.
inline bool exists(CharString const & filename)
{
  struct stat buffer;
  return (stat(toCString(filename), &buffer) == 0);
}

// ==========================================================================
// struct <Command>Options
// ==========================================================================

struct AssemblyOptions {
    CharString mappingFile;
    CharString referenceFile;
    CharString workingDirectory;
    CharString tmpDir;
    
    unsigned kmerLength;
    CharString adapters;
    int humanSeqs;
    unsigned threads;
    CharString memory;
    
    AssemblyOptions () :
        kmerLength(47), humanSeqs(maxValue<int>()), threads(1), memory("500000000")
    {}
};

struct MergingOptions {
    String<CharString> contigFiles;

    CharString outputFile;
    std::fstream outputStream;
    bool verbose;
    bool veryVerbose;
    
    double errorRate;
    int minimalLength;
    unsigned qgramLength;
    int matchScore;
    int errorPenalty;
    int minScore;
    int minTipScore;

    MergingOptions() :
        outputFile("supercontigs.fa"), verbose(false), veryVerbose(false),
        errorRate(0.01), minimalLength(60), qgramLength(47), matchScore(1), errorPenalty(-5), minScore(90), minTipScore(30)
    {} 
};

struct ContigMapOptions {
    CharString contigFile;
    CharString remappedFile;
    CharString workingDirectory;
    unsigned threads;
    CharString memory;
    CharString tmpDir;
    bool allAlignment;
    
    ContigMapOptions() :
        threads(1), memory("500000000"), allAlignment(false)
    {}
};

struct PlacingOptions {
    CharString supercontigFile;
    CharString referenceFile;
    String<CharString> nonRefFiles;
    
    CharString bamFilesFile;

    CharString locationsFile;
    CharString vcfInsertionsFile;
    CharString faInsertionsFile;
    
    unsigned batchIndex;
    unsigned batchSize;
    
    double minLocScore;
    
    unsigned readLength;
    unsigned maxInsertSize;
    
    bool verbose;
    
    PlacingOptions() :
        locationsFile("locations.txt"), vcfInsertionsFile("insertions.vcf"), faInsertionsFile("insertions.fa"),
        batchIndex(0), batchSize(maxValue<unsigned>()), minLocScore(0.3), readLength(100), maxInsertSize(800), verbose(false)
    {}
};

struct GenotypingOptions {
    CharString referenceFile;
    CharString bamFile;
    CharString altFastaFile;
    CharString altBamFile;
    CharString vcfFile;
    CharString sampleName;

    int match;
    int mismatch;
    int gapOpen;
    int gapExtend;
    int minAlignScore;

    int maxInsertSize;
    int bpQclip;
    int minSeqLen;
    double minReadProb;

    int regionWindowSize;
    bool addReadGroup;
    bool verbose;

    // hidden options
    bool callBoth;
    bool useReadCounts;
    bool fullOverlap;

    GenotypingOptions() : 
        sampleName("sample"), match(1), mismatch(-4), gapOpen(-10), gapExtend(-1), minAlignScore(55),
        maxInsertSize( 500 ), bpQclip(0), minSeqLen(10), minReadProb(0.0001),
        regionWindowSize(50), addReadGroup(false), verbose(false),
        callBoth(false), useReadCounts(false), fullOverlap(false)
    {}
};

// ==========================================================================

inline CharString
getFileName(CharString const & path, CharString const & name)
{
    CharString filename = path;
    filename += "/";
    filename += name;
    return filename;
}

inline void
removeFile(CharString const & path, const char * filename)
{
    CharString file = path;
    file += "/";
    file += filename;
    remove(toCString(file));
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
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM FILE\\fP");
    addDescription(parser, "Finds the unmapped reads in a bam files. If a fasta file is specified, the unmapped reads "
                           "will first be remapped to this reference using bwa and only reads that remain unmapped are "
                           "further processed. All unmapped reads are quality filtered using sickle and passed to "
                           "assembly with velvet.");

    // Require a bam file as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "BAMFILE"));
    
    
    // Setup the options.
    addOption(parser, ArgParseOption("d", "directory", "Path to working directory.", ArgParseArgument::STRING, "PATH"));
    setDefaultValue(parser, "directory", "current directory");
    
    addOption(parser, ArgParseOption("tmp", "tmpdir", "Path to a temporary directory ending with XXXXXX.", ArgParseArgument::STRING, "PATH"));
    setDefaultValue(parser, "tmpdir", "same as working directory");
    
    addOption(parser, ArgParseOption("k", "kmerLength", "The k-mer size for velvet assembly.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "kmerLength", options.kmerLength);
    
    addOption(parser, ArgParseOption("a", "adapters", "Enable adapter removal for Illumina reads. Default: \\fIno adapter removal\\fP.", ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "adapters", "HiSeq HiSeqX");

    addOption(parser, ArgParseOption("r", "reference", "Fasta file with reference sequences for remapping. Default: \\fIno remapping\\fP.", ArgParseArgument::INPUTFILE, "FILE"));
    setValidValues(parser, "reference", "fa fna fasta");
    
    addOption(parser, ArgParseOption("f", "filter", "Consider reads with low quality alignments as unmapped only for first INT sequences in the reference file. Requires reference file for remapping to be set.", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for bwa.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", options.threads);
    setMinValue(parser, "threads", "1");
    
    addOption(parser, ArgParseOption("m", "memory", "Maximum memory for samtools sort.", ArgParseArgument::STRING, "STR"));
    setDefaultValue(parser, "memory", options.memory);
}

void
setupParser(ArgumentParser & parser, MergingOptions & options)
{
    setShortDescription(parser, "Merging of insertion contigs into supercontigs.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE); 
    
    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFA FILE 1\\fP ... \\fIFA FILE N\\fP");
    addDescription(parser, "Merges the sequences given in fasta files into a single set of supercontigs. The algorithm "
                           "first partitions the sequences into sets of similar sequences using the SWIFT filtering "
                           "approach, and then aligns each set of contigs into a graph of supercontigs.");
                           
    // Require a list of fasta files as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FAFILE", true));
    
    // Program parameters
    addSection(parser, "Parameters");
    addOption(parser, ArgParseOption("e", "errRate", "Maximal error rate for SWIFT filtering.", ArgParseArgument::DOUBLE, "FLOAT"));
    addOption(parser, ArgParseOption("l", "minLength", "Minimal length for SWIFT filtering.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("k", "kmerLength", "Length of k-mers for SWIFT filtering.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "matchScore", "Match score for Smith-Waterman alignment.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("p", "penalty", "Error penalty for Smith-Waterman alignment.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("s", "minScore", "Minimal score for Smith-Waterman alignment.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("t", "minTipScore", "Minimal score for tips in supercontig graph.", ArgParseArgument::INTEGER, "INT"));
    
    // Output file options.
    addSection(parser, "Output options");
    addOption(parser, ArgParseOption("o", "outFile", "Name of output fasta file for the supercontigs.", ArgParseArgument::OUTPUTFILE, "CONTIGFILE"));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose screen output."));
    addOption(parser, ArgParseOption("vv", "veryVerbose", "Enable very verbose screen output."));

    // Set minimal/maximal/lists of valid values.
    setMinValue(parser, "e", "0");
    setMaxValue(parser, "e", "0.25");
    setMinValue(parser, "l", "3");
    setMinValue(parser, "k", "3");
    setMinValue(parser, "t", "0");
    setValidValues(parser, "o", "fa fna fasta");

    // Set default values.
    setDefaultValue(parser, "e", options.errorRate);
    setDefaultValue(parser, "l", options.minimalLength);
    setDefaultValue(parser, "k", options.qgramLength);
    setDefaultValue(parser, "m", options.matchScore);
    setDefaultValue(parser, "p", options.errorPenalty);
    setDefaultValue(parser, "s", options.minScore);
    setDefaultValue(parser, "t", options.minTipScore);
    setDefaultValue(parser, "o", options.outputFile);
}

void
setupParser(ArgumentParser & parser, ContigMapOptions & options)
{
    setShortDescription(parser, "Alignment of unmapped reads to assembled contigs.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE); 
    
    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFA FILE\\fP");
    addDescription(parser, "Aligns unmapped reads from fastq files in working directory to a set of contigs specified "
                           "in the fasta file using bwa-mem. Merges the bwa output file with the file non_ref.bam in "
                           "the working directory and sets the read mate's information in all bam records. Note that "
                           "the fasta file needs to be indexed for alignment with bwa-mem.");

    // Require a fasta file as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FAFILE"));

    // Setup the option.
    addOption(parser, ArgParseOption("d", "directory", "Path to working directory.", ArgParseArgument::STRING, "PATH"));
    setDefaultValue(parser, "directory", "current directory"); 

    addOption(parser, ArgParseOption("tmp", "tmpdir", "Path to a temporary directory ending with XXXXXX.", ArgParseArgument::STRING, "PATH"));
    setDefaultValue(parser, "tmpdir", "same as working directory");

    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for bwa.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", options.threads);
    setMinValue(parser, "threads", "1");
    
    addOption(parser, ArgParseOption("m", "memory", "Maximum memory for samtools sort.", ArgParseArgument::STRING, "STR"));
    setDefaultValue(parser, "memory", options.memory);

    addOption(parser, ArgParseOption("a", "all", "Use bwa-mem's -a option to output all alignments of a read."));
    setDefaultValue(parser, "all", "false");
}

void
setupParser(ArgumentParser & parser, PlacingOptions & options)
{
    setShortDescription(parser, "Finding positions of contigs in reference genome.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE); 
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fICONTIGFILE\\fP \\fIREFFILE\\fP \\fIBAM FILE 1\\fP ... \\fIBAM FILE N\\fP");
    addDescription(parser, "Finds the positions of (super-)contigs in the reference genome. Identifies approximate "
                           "locations based on anchoring read pairs found in the bam files if a file with locations "
                           "does not already exist. Determines exact positions of insertions using SplazerS for split "
                           "read alignment for each contig end if bam files with all reads of the individuals are "
                           "specified. Outputs a vcf and fa record for each identified position. The split alignment "
                           "can be done in batches (e.g. 100 locations per batch) if the approximate locations have "
                           "been computed before.");
    
    // Require fasta file with merged contigs as arguments.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "CONTIGFILE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "REFFILE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "BAMFILE", true));
    
    // Setup (input) options.
    addSection(parser, "Main options");
    addOption(parser, ArgParseOption("l", "locations", "Name of file with approximate insertion locations. Computed if not exists.", ArgParseArgument::STRING, "LOCATIONFILE"));
    addOption(parser, ArgParseOption("m", "minScore", "Minimal score of a location to be passed to split mapping.", ArgParseArgument::DOUBLE, "FLOAT"));
    addOption(parser, ArgParseOption("b", "bamFiles", "File listing original, full bam files of individuals, one per line. Specify to determine exact insertion positions from split reads.", ArgParseArgument::INPUTFILE, "FILE"));
    addOption(parser, ArgParseOption("s", "batchSize", "Number of locations per batch. Specify to split computation into smaller batches. Requires locations file to exist, bam files, and batch number.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("i", "batchIndex", "Number of batch. Specify to split computation into smaller batches. Requires locations file to exist, bam files, and batch size.", ArgParseArgument::INTEGER, "INT"));
    
    addOption(parser, ArgParseOption("r", "readLength", "The length of the reads.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("e", "maxInsertSize", "The maximal expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    
    // Output file options.
    addSection(parser, "Output options");
    addOption(parser, ArgParseOption("ov", "outVcf", "Name of output file for vcf records.", ArgParseArgument::OUTPUTFILE, "VCFFILE"));
    addOption(parser, ArgParseOption("of", "outFa", "Name of output file for insertion sequences.", ArgParseArgument::OUTPUTFILE, "FAFILE"));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));

    setValidValues(parser, "of", "fa fna fasta");
    setValidValues(parser, "ov", "vcf");
    setMinValue(parser, "m", "0");
    setMaxValue(parser, "m", "1");
    
    // Set default values.
    setDefaultValue(parser, "l", options.locationsFile);
    setDefaultValue(parser, "m", options.minLocScore);
    setDefaultValue(parser, "r", options.readLength);
    setDefaultValue(parser, "e", options.maxInsertSize);
    setDefaultValue(parser, "ov", options.vcfInsertionsFile);
    setDefaultValue(parser, "of", options.faInsertionsFile);
}

void
setupParser(ArgumentParser & parser, GenotypingOptions & options)
{
    setShortDescription(parser, "Genotyping an individual for the insertions.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE); 
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFASTAFILE\\fP \\fIBAMFILE\\fP \\fIFASTAFILEALT\\fP \\fIBAMFILEALT\\fP \\fIVCFFILE\\fP");
    addDescription(parser, "The genotype command takes as input a fasta file of the reference genome, a bam file of a "
                           "single individual, the fasta file with the supercontigs, the bam file of contig mapped and "
                           "unmapped reads (<WD>/non_ref.bam), and the VCF file with all predicted insertion "
                           "positions. It computes genotype likelihoods by aligning all reads from each insertion "
                           "location and contig to the reference and to the alternative insertion sequence. It outputs "
                           "VCF records with the genotype likelihoods in GT:PL format for the individual to std::out.");
    
    // Require five files as arguments.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "fastafile"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "bamfile"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "fastafilealt"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "bamfilealt"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "vcffile"));
    
    // Options.
    addSection(parser, "Alignment options");
    addOption(parser, ArgParseOption("m", "match", "Cost of match.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("mm", "mismatch", "Cost of mismatch.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("go", "gapopen", "Cost of gap open.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ge", "gapextend", "Cost of gap extend.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("s", "minscore", "Minimum alignment score.", ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "Read(-pair) options");
    addOption(parser, ArgParseOption("i", "insertSize", "Maximum read pair insert size.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("q", "qual", "Quality score threshold for read trimming.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("l", "minSeqLen", "Minimum read length after trimming.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("pm", "minreadprob", "Minimum read probability.", ArgParseArgument::DOUBLE, "DOUBLE"));

    addSection(parser, "Misc options");
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("w", "window", "Region window size.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("rg", "addreadgroup", "Add read group."));
    addOption(parser, ArgParseOption("sa", "samplename", "Name of sample for vcf output", ArgParseArgument::STRING, "SAMPLENAME"));

    // Hidden options.
    addOption(parser, ArgParseOption("b", "callboth", "callboth"));
    hideOption(parser, "b");
    addOption(parser, ArgParseOption("r", "readcounts", "use read counts"));
    hideOption(parser, "r");    
    addOption(parser, ArgParseOption("f", "fulloverlap", "full overlap of read"));
    hideOption(parser, "f");

    // Defualt values.
    setDefaultValue(parser, "match", options.match);
    setDefaultValue(parser, "mismatch", options.mismatch);
    setDefaultValue(parser, "gapopen", options.gapOpen);
    setDefaultValue(parser, "gapextend", options.gapExtend);
    setDefaultValue(parser, "minscore", options.minAlignScore);
    setDefaultValue(parser, "insertSize", options.maxInsertSize);
    setDefaultValue(parser, "qual", options.bpQclip);
    setDefaultValue(parser, "minreadprob", options.minReadProb);
    setDefaultValue(parser, "window", options.regionWindowSize);
    setDefaultValue(parser, "samplename", options.sampleName);
}

// ==========================================================================
// Function getOptionValues()
// ==========================================================================

int
getOptionValues(AssemblyOptions & options, ArgumentParser const & parser)
{
    getArgumentValue(options.mappingFile, parser, 0);
    
    if (isSet(parser, "reference"))
        getOptionValue(options.referenceFile, parser, "reference");
    if (isSet(parser, "kmerLength"))
        getOptionValue(options.kmerLength, parser, "kmerLengh");
    if (isSet(parser, "directory"))
        getOptionValue(options.workingDirectory, parser, "directory");
    if (isSet(parser, "tmpdir"))
        getOptionValue(options.tmpDir, parser, "tmpdir");
    if (isSet(parser, "adapters"))
        getOptionValue(options.adapters, parser, "adapters");
    if (isSet(parser, "filter"))
        getOptionValue(options.humanSeqs, parser, "filter");
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");

    return 0;
}

int
getOptionValues(MergingOptions & options, ArgumentParser & parser)
{

    options.contigFiles = getArgumentValues(parser, 0);
    if (length(options.contigFiles) < 2)
    {
        std::cerr << "ERROR: Too few arguments. Please specify at least two fasta files." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    
    // Get parameter values.
    if (isSet(parser, "errRate"))
        getOptionValue(options.errorRate, parser, "errRate");
    if (isSet(parser, "minLength"))
        getOptionValue(options.minimalLength, parser, "minLength");
    if (isSet(parser, "kmerLength"))
        getOptionValue(options.qgramLength, parser, "kmerLength");
    if (isSet(parser, "matchScore"))
        getOptionValue(options.matchScore, parser, "matchScore");
    if (isSet(parser, "minScore"))
        getOptionValue(options.minScore, parser, "minScore");
    if (isSet(parser, "penalty"))
        getOptionValue(options.errorPenalty, parser, "penalty");
    if (isSet(parser, "minTipScore"))
        getOptionValue(options.minTipScore, parser, "minTipScore");

    // Get output options.
    if (isSet(parser, "outFile"))
        getOptionValue(options.outputFile, parser, "outFile");
    if (isSet(parser, "verbose"))
        options.verbose = true;
    if (isSet(parser, "veryVerbose"))
    {
        options.verbose = true;
        options.veryVerbose = true;
    }
    
    return 0;
}

int
getOptionValues(ContigMapOptions & options, ArgumentParser & parser)
{
    getArgumentValue(options.contigFile, parser, 0);
    
    if (isSet(parser, "directory"))
        getOptionValue(options.workingDirectory, parser, "directory");
    if (isSet(parser, "tmpdir"))
        getOptionValue(options.tmpDir, parser, "tmpdir");
    if (isSet(parser, "all"))
        options.allAlignment = true;
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");

    return 0;
}

int
getOptionValues(PlacingOptions & options, ArgumentParser & parser)
{
    getArgumentValue(options.supercontigFile, parser, 0);
    getArgumentValue(options.referenceFile, parser, 1);
    options.nonRefFiles = getArgumentValues(parser, 2);
    
    if (isSet(parser, "locations"))
        getOptionValue(options.locationsFile, parser, "locations");
    if (!exists(options.locationsFile) && (isSet(parser, "batchIndex") || isSet(parser, "batchSize")))
    {
        std::cerr << "ERROR: Locations file " << options.locationsFile << " does not exist."
                  <<       " Please compute locations before splitting exact positioning into batches." << std::endl;
        return 1;
    }
    if (isSet(parser, "minScore"))
        getOptionValue(options.minLocScore, parser, "minScore");

    if (isSet(parser, "bamFiles"))
        getOptionValue(options.bamFilesFile, parser, "bamFiles");
    if (!isSet(parser, "bamFiles") && (isSet(parser, "batchIndex") || isSet(parser, "batchSize")))
    {
        std::cerr << "ERROR: Bam files with all reads of individuals not specified." << std::endl;
        return 1;
    }

    if (isSet(parser, "batchIndex"))
    {
        if (!isSet(parser, "batchSize"))
        {
            std::cerr << "ERROR: Batch size not specified." << std::endl;
            return 1;
        }
        getOptionValue(options.batchIndex, parser, "batchIndex");
    }
    if (isSet(parser, "batchSize"))
    {
        if (!isSet(parser, "batchIndex"))
        {
            std::cerr << "ERROR: Batch index not specified." << std::endl;
            return 1;
        }
        getOptionValue(options.batchSize, parser, "batchSize");
    }

    if (isSet(parser, "readLength"))
        getOptionValue(options.readLength, parser, "readLength");
    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");

    if (isSet(parser, "outVcf"))
        getOptionValue(options.vcfInsertionsFile, parser, "outVcf");
    if (isSet(parser, "outFa"))
        getOptionValue(options.faInsertionsFile, parser, "outFa");
        
    options.verbose = isSet(parser, "verbose");

    return 0;
}

int
getOptionValues(GenotypingOptions & options, ArgumentParser & parser)
{
    getArgumentValue(options.referenceFile, parser, 0);
    getArgumentValue(options.bamFile, parser, 1);
    getArgumentValue(options.altFastaFile, parser, 2);
    getArgumentValue(options.altBamFile, parser, 3);
    getArgumentValue(options.vcfFile, parser, 4);
    
    if (isSet(parser, "match"))
        getOptionValue(options.match, parser, "match");
    if (isSet(parser, "mismatch"))
        getOptionValue(options.mismatch, parser, "mismatch");
    if (isSet(parser, "gapopen"))
        getOptionValue(options.gapOpen, parser, "gapopen");
    if (isSet(parser, "gapextend"))
        getOptionValue(options.gapExtend, parser, "gapextend");
    if (isSet(parser, "minscore"))
        getOptionValue(options.minAlignScore, parser, "minscore");

    if (isSet(parser, "insertSize"))
        getOptionValue(options.maxInsertSize, parser, "insertSize");
    if (isSet(parser, "minreadprob"))
        getOptionValue(options.minReadProb, parser, "minreadprob");
    if (isSet(parser, "qual"))
        getOptionValue(options.bpQclip, parser, "qual");
    if (isSet(parser, "minSeqLen"))
        getOptionValue(options.minSeqLen, parser, "minSeqLen");
    if (isSet(parser, "samplename"))
        getOptionValue(options.sampleName, parser, "samplename");
        
    if(isSet(parser, "window"))
        getOptionValue( options.regionWindowSize, parser, "window");
    options.addReadGroup = isSet(parser, "addreadgroup");
    options.verbose = isSet(parser, "verbose");

    options.callBoth = isSet(parser, "callboth");
    options.fullOverlap = isSet(parser, "fulloverlap");
    options.useReadCounts = isSet(parser, "readcounts");

    return 0;
}

// ==========================================================================
// Function parseCommandLine()
// ==========================================================================

template<typename TOptions>
int parseCommandLine(TOptions & options, int argc, char const ** argv)
{
    // Concatenate program name from name and command.
    CharString prog_name = argv[0];
    prog_name += " ";
    prog_name += argv[1];

    ++argv;
    --argc;

    // Setup the parser.
    ArgumentParser parser(toCString(prog_name));
    setupParser(parser, options);
    
    // Parse the command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return 1;
    
    // Collect the option values.
    getOptionValues(options, parser);

    return 0;
}

// ==========================================================================

#endif // #ifndef POPINS_CLP_H_
