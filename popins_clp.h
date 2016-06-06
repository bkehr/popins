#include <seqan/file.h>
#include <seqan/arg_parse.h>

#ifndef POPINS_CLP_H_
#define POPINS_CLP_H_

#include <string>

using namespace seqan;

// ==========================================================================

// Returns true if file exists, otherwise false.
inline bool exists(CharString const & filename)
{
  struct stat buffer;
  return (stat(toCString(filename), &buffer) == 0);
}

// ==========================================================================
// Function checkFileEnding()
// ==========================================================================

bool
checkFileEnding(CharString & filename, std::string ending)
{
    std::string name = toCString(filename);
    size_t dotPos = name.find_last_of('.');

    if (dotPos == std::string::npos)
        return false;

    return name.substr(dotPos + 1, 3) == ending;
}

// ==========================================================================
// Function readFileNames()
// ==========================================================================

bool
readFileNames(String<CharString> & files, CharString & filenameFile)
{
    if (filenameFile == "") return 0;

    std::fstream stream(toCString(filenameFile), std::fstream::in);
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open file listing files " << filenameFile << std::endl;
        return 1;
    }

    RecordReader<std::fstream, SinglePass<> > reader(stream);

    while (!atEnd(reader))
    {
        CharString file;
        int res = readLine(file, reader);
        if (res != 0)
        {
            std::cerr << "ERROR while reading line from " << filenameFile << std::endl;
            return 1;
        }
        appendValue(files, file);
    }

    return 0;
}

// ==========================================================================

template <typename TValue>
bool readFileNames(String<CharString> & files, String<TValue> & values)
{
    if (length(files) > 1) return 0;
    std::cerr << "ReadFileNames " << length(files) << " " << files << std::endl;
    // Open input file
    CharString filenameFile = files[0];
    std::cerr << filenameFile << " " << files << std::endl;
    std::fstream stream(toCString(filenameFile), std::fstream::in);
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open file listing files " << filenameFile << std::endl;
        return 1;
    }

    RecordReader<std::fstream, SinglePass<> > reader(stream);
    clear(files);

    while (!atEnd(reader))
    {
        // Read the file name
        CharString file;
        int res = readUntilWhitespace(file, reader);
        if (res != 0)
        {
            std::cerr << "ERROR: Failed reading filename from " << filenameFile << std::endl;
            return 1;
        }
        appendValue(files, file);

        skipBlanks(reader);

        // Read the value for this filename
        CharString buffer;
        res = readLine(buffer, reader);
        if (res != 0 || length(buffer) == 0)
        {
            std::cerr << "ERROR: Failed reading second column of " << filenameFile << std::endl;
            return 1;
        }
        TValue val;
        lexicalCast2<TValue>(val, buffer);
        appendValue(values, val);
    }

    return 0;
}

// ==========================================================================

bool
parseInterval(Triple<CharString, unsigned, unsigned> & out, CharString & in)
{
    Iterator<CharString, Rooted>::Type it = begin(in, Rooted());

    unsigned colonPos = 0;
    while (it != end(in))
    {
        if (*it == ':')
        {
            colonPos = position(it);
            break;
        }
        ++it;
    }

    if (colonPos == 0)
    {
        out.i1 = in;
        out.i2 = 0;
        out.i3 = maxValue<unsigned>();
        return 0;
    }

    unsigned dashPos = 0;
    while (it != end(in))
    {
        if (*it == '-')
        {
            dashPos = position(it);
            break;
        }
        ++it;
    }

    if (dashPos == 0)
    {
        std::cerr << "ERROR: Interval is not in format CHR:BEG-END." << std::endl;
        return 1;
    }

    out.i1 = prefix(in, colonPos);
    out.i2 = lexicalCast<unsigned>(infix(in, colonPos + 1, dashPos));
    out.i3 = lexicalCast<unsigned>(suffix(in, dashPos + 1));

    return 0;
}

// ==========================================================================
// struct <Command>Options
// ==========================================================================

struct AssemblyOptions {
    CharString mappingFile;
    CharString matepairFile;
    CharString referenceFile;
    CharString workingDirectory;

    unsigned kmerLength;
    CharString adapters;
    int humanSeqs;
    unsigned threads;
    CharString memory;
    bool matepair;

    AssemblyOptions () :
        kmerLength(47), humanSeqs(maxValue<int>()), threads(1), memory("500000000"), matepair(false)
    {}
};

struct MergingOptions {
    String<CharString> contigFiles;
    String<CharString> componentFiles;
    CharString outputFile;
    CharString skippedFile;
    std::fstream outputStream;
    std::fstream skippedStream;
    bool verbose;
    bool veryVerbose;

    unsigned batchIndex;
    unsigned batches;

    double errorRate;
    int minimalLength;
    unsigned qgramLength;
    int matchScore;
    int errorPenalty;
    int minScore;
    int minTipScore;
    double minEntropy;

    MergingOptions() :
        outputFile("supercontigs.fa"), verbose(false), veryVerbose(false), batchIndex(0), batches(1),
        errorRate(0.01), minimalLength(60), qgramLength(47), matchScore(1), errorPenalty(-5), minScore(90), minTipScore(30), minEntropy(0.75)
    {}
};

struct ContigMapOptions {
    CharString contigFile;
    CharString remappedFile;
    CharString workingDirectory;
    unsigned threads;
    CharString memory;
    bool allAlignment;
    int maxInsertSize;
    bool keepNonRefNew;

    ContigMapOptions() :
        threads(1), memory("500000000"), allAlignment(false), maxInsertSize(800), keepNonRefNew(false)
    {}
};

struct PlacingOptions {
    CharString locationsFile;           // merged from all individuals
    String<CharString> locationsFiles;  // one file per individual

    CharString outFile;
    bool isVcf;

    CharString groupsFile;

    CharString supercontigFile;
    CharString referenceFile;

    CharString bamFile;                     // for one individual
    String<CharString> bamFiles;
    double bamAvgCov;                     // for one individual
    String<double> bamAvgCovs;
    CharString splitAlignFile;             // for one individual

    Triple<CharString, unsigned, unsigned> interval; // chrom, beginPos, endPos

    double minLocScore;
    unsigned minAnchorReads;

    unsigned readLength;
    unsigned maxInsertSize;

    unsigned groupDist;

    bool verbose;

    PlacingOptions() :
        locationsFile("locations.txt"), groupsFile("groups.txt"), bamFile(""),
        minLocScore(0.3), minAnchorReads(2), readLength(100), maxInsertSize(800), groupDist(100),
        verbose(false)
    {}
};

enum genotypingModelType { randomSequenceGenotyping, duplicationGenotyping };

genotypingModelType
genotypingModelEnum( CharString& gmString ){
  if( gmString == "DUP" ){
    return duplicationGenotyping;
  }else{
    return randomSequenceGenotyping;
  }
}

void
genotypingModelName( genotypingModelType gm, CharString& gmString ){
  if( gm == duplicationGenotyping ){
    gmString = CharString( "DUP" );

  }else{
    gmString = CharString( "RANDOM" );
  }
}


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
  genotypingModelType genotypingModel;

    int maxInsertSize;
    int bpQclip;
    int minSeqLen;
    double minReadProb;
    int maxBARcount;

    int regionWindowSize;
    bool addReadGroup;
    bool verbose;

    // hidden options
    bool callBoth;
    bool useReadCounts;
    bool fullOverlap;

    GenotypingOptions() : 
        sampleName("sample"), match(1), mismatch(-2), gapOpen(-4), gapExtend(-1), minAlignScore(55),
	genotypingModel( randomSequenceGenotyping ), maxInsertSize( 500 ), bpQclip(0), minSeqLen(10), 
        minReadProb(0.00001), maxBARcount(200), regionWindowSize(50), addReadGroup(false), verbose(false),
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
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM FILE\\fP [\\fIMP BAM FILE\\fP]");
    addDescription(parser, "Finds the unmapped reads in a bam files. If a fasta file is specified, the unmapped reads "
                           "will first be remapped to this reference using bwa and only reads that remain unmapped are "
                           "further processed. All unmapped reads are quality filtered using sickle and passed to "
                           "assembly with velvet.");

    // Require a bam file as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "BAMFILE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "MPFILE"));

    // Setup the options.
    addOption(parser, ArgParseOption("d", "directory", "Path to working directory.", ArgParseArgument::STRING, "PATH"));
    setDefaultValue(parser, "directory", "current directory");

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
    addOption(parser, ArgParseOption("", "matepair", "Use matepair data for assembly."));
}

void
setupParser(ArgumentParser & parser, MergingOptions & options)
{
    setShortDescription(parser, "Merging of insertion contigs into supercontigs.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFA FILE 1\\fP ... \\fIFA FILE N\\fP");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFILE LIST FILE\\fP");
    addDescription(parser, "Merges the sequences given in fasta files into a single set of supercontigs. The fasta "
                           "files can be listed either on the command line or in a file (FILE LIST FILE). In the "
                           "latter case, a two-column file is expected, the first column giving path/to/contigs.fa and "
                           "the second column specifying the number of contigs in the fasta file. The algorithm first "
                           "partitions the sequences into sets of similar sequences using the SWIFT filtering "
                           "approach, and then aligns each set of sequences into a graph of supercontigs. These two "
                           "steps of the algorithm can be split into several program calls (see Note below).");

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
    addOption(parser, ArgParseOption("y", "minEntropy", "Minimal entropy for discarding low-complexity sequences. Choose as 0.0 to disable filter.", ArgParseArgument::DOUBLE, "FLOAT"));


    // Program mode options
    addSection(parser, "Program mode options");
    addOption(parser, ArgParseOption("b", "batches", "Total number of batches (see Note below).", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("i", "batchIndex", "Batch number (see Note below).", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("c", "componentFiles", "List of files written by the partioning step. Either the name of a single file that lists one component filename per line, or multiple names of component files, i.e. -c <FILE1> -c <FILE2> ...", ArgParseArgument::INPUTFILE, "FILE", true));

    // Output file options.
    addSection(parser, "Output options");
    addOption(parser, ArgParseOption("o", "outFile", "Name of output file. Either in text format for components or fasta format for the supercontigs.", ArgParseArgument::OUTPUTFILE, "OUTPUTFILE"));
    addOption(parser, ArgParseOption("os", "skippedFile", "Name of output file for skipped contigs. Default: no output of skipped contigs.", ArgParseArgument::OUTPUTFILE, "OUTPUTFILE"));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose screen output."));
    addOption(parser, ArgParseOption("vv", "veryVerbose", "Enable very verbose screen output."));

    // Note on usage.
    addSection(parser, "Note");
    addText(parser, "When merging a large number of contigs, the q-gram index in the partitioning step might not fit "
                    "into memory. In this case, you can split the merging into several program calls using different "
                    "program modes: First, generate a set of component files by running 'popins merge' with the -b "
                    "and -i options. Afterwards, combine the component files and construct the supercontigs in a "
                    "separate call to 'popins merge' by specifying the -c option. By specifying the -c option "
                    "together with the -b and -i option, the supercontig construction is split into batches.");

    // Set minimal/maximal/lists of valid values.
    setMinValue(parser, "e", "0");
    setMaxValue(parser, "e", "0.25");
    setMinValue(parser, "l", "3");
    setMinValue(parser, "k", "3");
    setMinValue(parser, "t", "0");
    setMinValue(parser, "y", "0");
    setMaxValue(parser, "y", "1");
    setValidValues(parser, "o", "txt fa fna fasta");
    setValidValues(parser, "os", "fa fna fasta");

    // Set default values.
    setDefaultValue(parser, "e", options.errorRate);
    setDefaultValue(parser, "l", options.minimalLength);
    setDefaultValue(parser, "k", options.qgramLength);
    setDefaultValue(parser, "m", options.matchScore);
    setDefaultValue(parser, "p", options.errorPenalty);
    setDefaultValue(parser, "s", options.minScore);
    setDefaultValue(parser, "t", options.minTipScore);
    setDefaultValue(parser, "y", options.minEntropy);
    setDefaultValue(parser, "b", options.batches);
    setDefaultValue(parser, "i", options.batchIndex);
    setDefaultValue(parser, "o", "supercontigs.fa or components_<BATCH INDEX>.txt");
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

    addOption(parser, ArgParseOption("t", "threads", "Number of threads to use for bwa.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", options.threads);
    setMinValue(parser, "threads", "1");

    addOption(parser, ArgParseOption("m", "memory", "Maximum memory for samtools sort.", ArgParseArgument::STRING, "STR"));
    setDefaultValue(parser, "memory", options.memory);

    addOption(parser, ArgParseOption("a", "all", "Use bwa-mem's -a option to output all alignments of a read."));
    setDefaultValue(parser, "all", "false");

    addOption(parser, ArgParseOption("e", "maxInsertSize", "The maximal expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "maxInsertSize", options.maxInsertSize);

    addOption(parser, ArgParseOption("n", "nonRefNew", "Do not delete the non_ref_new.bam file after writing locations."));
    setDefaultValue(parser, "nonRefNew", "false");
}

void
setupParser(ArgumentParser & parser, PlacingOptions & options)
{
    setShortDescription(parser, "Finding positions of contigs in reference genome.");
    setVersion(parser, VERSION);
    setDate(parser, VERSION_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fILOCFILE\\fP \\fIOUTPUTFILE\\fP");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] -c \\fICONTIGFILE\\fP -r \\fIGENOME\\fP \\fILOCFILE\\fP \\fIOUTPUTFILE\\fP");
    addDescription(parser, "Places (super-)contigs into the reference genome and writes a VCF file with records for "
                           "placed contig ends. The command consists of four steps that can be run in one program "
                           "call or in separate calls. To run all four steps together, the options -l, -b, -c, and -r "
                           "need to be specified and the \\fIOUTPUTFILE\\fP's ending needs to be 'vcf'. The \\fILOCFILE\\fP has "
                           "to list the location files for all samples. When running the steps separately, the "
                           "specified parameters determine which step is being run:");
    addDescription(parser, "1. First, the contig locations determined for all samples (files listed in \\fILOCFILE\\fP) need "
                           "to be merged into one set of locations (written to \\fIOUTPUTFILE\\fP).");
    addDescription(parser, "2. Then, prefixes/suffixes of all contigs (specify with -c option) are aligned to these "
                           "locations (specify as \\fILOCFILE\\fP) in the reference genome (specify with -r option) and VCF "
                           "records are written to \\fIOUTPUTFILE\\fP if the alignment is successful. The \\fIOUTPUTFILE\\fP's "
                           "ending has to be 'vcf'. Additional output files \\fIlocations_unplaced.txt\\fP are being written "
                           "for each sample.");
    addDescription(parser, "3. Next, contigs (specify -c option) that do not align to the reference genome (specify "
                           "-r option) are passed on to split-read alignment. This step is run by sample. The "
                           "sample's original BAM file needs to be specified. The program arguments, the \\fILOCFILE\\fP "
                           "and \\fIOUTPUTFILE\\fP, need to be the sample's \\fIlocations_unplaced.txt\\fP file and a "
                           "\\fIlocations_placed.txt\\fP file for the sample.");
    addDescription(parser, "4. Finally, the results from split-read alignment (the \\fIlocations_placed.txt\\fP files) of all "
                           "samples (input files listed in \\fILOCFILE\\fP) are being combined and appended to the "
                           "\\fIOUTPUTFILE\\fP (file ending has to be 'vcf').");

    // Require a locations file and an output file as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "LOCFILE")); // Name of locations file OR name of file listing locations files, one per line.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "OUTFILE")); // Name of output file. Either VCF or text file listing locations.

    // Setup (input) options.
    addSection(parser, "Main options");

    addOption(parser, ArgParseOption("c", "contigFile", "Name of (super-)contigs file.", ArgParseArgument::INPUTFILE, "FILE"));
    addOption(parser, ArgParseOption("r", "genomeFile", "Name of reference genome file.", ArgParseArgument::INPUTFILE, "FILE"));
    addOption(parser, ArgParseOption("m", "minScore", "Minimal anchoring score for a location.", ArgParseArgument::DOUBLE, "FLOAT"));
    addOption(parser, ArgParseOption("n", "minReads", "Minimal number of anchoring read pairs for a location.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("d", "groupDist", "Minimal distance between groups of locations.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("b", "bamFile", "Full BAM file of an individual. Specify to determine exact insertion positions from split reads.", ArgParseArgument::INPUTFILE, "FILE"));
    addOption(parser, ArgParseOption("a", "bamCov", "Average coverage of the genome in the BAM file. Required if -b option specified.", ArgParseArgument::DOUBLE, "FLOAT"));

    addOption(parser, ArgParseOption("len", "readLength", "The length of the reads.", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("e", "maxInsertSize", "The maximal expected insert size of the read pairs.", ArgParseArgument::INTEGER, "INT"));

    // Output file options.
    addSection(parser, "Output options");
    addOption(parser, ArgParseOption("l", "locations", "Name of merged locations file if placing is run in one program call.", ArgParseArgument::OUTPUTFILE, "FILE"));
    addOption(parser, ArgParseOption("g", "groups", "Name of groups output file.", ArgParseArgument::OUTPUTFILE, "FILE"));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));

    setMinValue(parser, "m", "0");
    setMaxValue(parser, "m", "1");

    setValidValues(parser, "c", "fa fna fasta");
    setValidValues(parser, "r", "fa fna fasta");

    // Set default values.
    setDefaultValue(parser, "m", options.minLocScore);
    setDefaultValue(parser, "n", options.minAnchorReads);
    setDefaultValue(parser, "d", options.groupDist);
    setDefaultValue(parser, "len", options.readLength);
    setDefaultValue(parser, "e", options.maxInsertSize);
    setDefaultValue(parser, "g", options.groupsFile);
    setDefaultValue(parser, "l", options.locationsFile);
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
    addOption(parser, ArgParseOption("mb", "maxBARcount", "Maximum number of reads to consider in region window.", ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "Misc options");
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("w", "window", "Region window size.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("rg", "addreadgroup", "Add read group."));
    addOption(parser, ArgParseOption("sa", "samplename", "Name of sample for vcf output", ArgParseArgument::STRING, "SAMPLENAME"));
    addOption(parser, ArgParseOption("gm", "genotypingmodel", "Model used for genotyping", ArgParseArgument::STRING, "GENOTYPINGMODEL"));

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
    setDefaultValue(parser, "maxBARcount", options.maxBARcount);
    setDefaultValue(parser, "window", options.regionWindowSize);
    setDefaultValue(parser, "samplename", options.sampleName);
    CharString gmName;
    genotypingModelName( options.genotypingModel, gmName ); 
    setDefaultValue(parser, "genotypingmodel", gmName );
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
    if (isSet(parser, "adapters"))
        getOptionValue(options.adapters, parser, "adapters");
    if (isSet(parser, "filter"))
        getOptionValue(options.humanSeqs, parser, "filter");
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");
    if (isSet(parser, "matepair")) {
        options.matepair = true;
        if (!getArgumentValue(options.matepairFile, parser, 1)) {
          return 1;
        }
    }

    return 0;
}

int
getOptionValues(MergingOptions & options, ArgumentParser & parser)
{
    options.contigFiles = getArgumentValues(parser, 0);

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
    if (isSet(parser, "minEntropy"))
        getOptionValue(options.minEntropy, parser, "minEntropy");

    // Get program mode options.
    if (isSet(parser, "batchIndex") && isSet(parser, "batches"))
    {
        getOptionValue(options.batchIndex, parser, "batchIndex");
        getOptionValue(options.batches, parser, "batches");
        if (options.batches <= options.batchIndex)
        {
            std::cerr << "ERROR: Please specify batch index smaller than number of batches." << std::endl;
            return 1;
        }
        if (!isSet(parser, "componentFiles"))
        {
            std::stringstream filename;
            filename << "components_" << options.batchIndex << ".txt";
            options.outputFile = filename.str();
        }
        else
        {
            std::stringstream filename;
            filename << "supercontigs_" << options.batchIndex << ".fa";
            options.outputFile = filename.str();
        }
    }
    else if (isSet(parser, "batchIndex") || isSet(parser, "batches"))
    {
        std::cerr << "ERROR: Please specify both options --batchIndex and --batches." << std::endl;
        return 1;
    }

    if (isSet(parser, "componentFiles"))
    {
        options.componentFiles = getOptionValues(parser, "componentFiles");
        if (length(options.componentFiles) == 1)
        {
            CharString filenamesFile = options.componentFiles[0];
            clear(options.componentFiles);
            if (readFileNames(options.componentFiles, filenamesFile) != 0) return 1;
        }
    }

    // Get output options.
    if (isSet(parser, "outFile"))
        getOptionValue(options.outputFile, parser, "outFile");
    if (isSet(parser, "skippedFile"))
        getOptionValue(options.skippedFile, parser, "skippedFile");
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
    if (isSet(parser, "all"))
        options.allAlignment = true;
    if (isSet(parser, "nonRefNew"))
        options.keepNonRefNew = true;
    if (isSet(parser, "maxInsertSize"))
        getOptionValue(options.maxInsertSize, parser, "maxInsertSize");
    if (isSet(parser, "threads"))
        getOptionValue(options.threads, parser, "threads");
    if (isSet(parser, "memory"))
        getOptionValue(options.memory, parser, "memory");

    return 0;
}

int
getOptionValues(PlacingOptions & options, ArgumentParser & parser)
{
    CharString locFile = options.locationsFile;

    getArgumentValue(options.locationsFile, parser, 0);
    getArgumentValue(options.outFile, parser, 1);

    options.isVcf = checkFileEnding(options.outFile, "vcf");

    if (isSet(parser, "contigFile") && isSet(parser, "genomeFile"))
    {
        // Step 2 or 3 or all
        getOptionValue(options.supercontigFile, parser, "contigFile");
        getOptionValue(options.referenceFile, parser, "genomeFile");

        if (isSet(parser, "maxInsertSize"))
            getOptionValue(options.maxInsertSize, parser, "maxInsertSize");

        if (isSet(parser, "bamFile"))
        {
            // Step 3 or all
            getOptionValue(options.bamFile, parser, "bamFile");

            if (isSet(parser, "bamCov"))
                getOptionValue(options.bamAvgCov, parser, "bamCov");
            if (isSet(parser, "readLength"))
                getOptionValue(options.readLength, parser, "readLength");

            if (options.isVcf)
            {
                // All steps
                if (readFileNames(options.locationsFiles, options.locationsFile) != 0)
                    return 1;

                if (isSet(parser, "locations"))
                    getOptionValue(locFile, parser, "locations");
                options.locationsFile = locFile;

                appendValue(options.bamFiles, options.bamFile);
                if (readFileNames(options.bamFiles, options.bamAvgCovs) != 0)
                    return 1;
            }
        }
        if (isSet(parser, "minScore"))
            getOptionValue(options.minLocScore, parser, "minScore");
        if (isSet(parser, "minReads"))
            getOptionValue(options.minAnchorReads, parser, "minReads");
        if (isSet(parser, "groups"))
            getOptionValue(options.groupsFile, parser, "groups");
        if (isSet(parser, "groupDist"))
            getOptionValue(options.groupDist, parser, "groupDist");
    }
    else
    {
        // Step 1 or 4
        CharString locationsFilesFile;
        if (readFileNames(options.locationsFiles, options.locationsFile) != 0)
            return 1;

        // Step 4
        if (isSet(parser, "genomeFile"))
            getOptionValue(options.referenceFile, parser, "genomeFile");
        else if (options.isVcf)
        {
            std::cerr << "ERROR: Please specify the reference genome file with the -r (--genomeFile) option for combining split-read alignment results." << std::endl;
            return 1;
        }

    }

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
    if (isSet(parser, "maxBARcount"))
        getOptionValue(options.maxBARcount, parser, "maxBARcount");
    if (isSet(parser, "qual"))
        getOptionValue(options.bpQclip, parser, "qual");
    if (isSet(parser, "minSeqLen"))
        getOptionValue(options.minSeqLen, parser, "minSeqLen");
    if (isSet(parser, "samplename"))
        getOptionValue(options.sampleName, parser, "samplename");
    if (isSet(parser, "genotypingmodel")){
      CharString gmString;
      getOptionValue( gmString, parser, "genotypingmodel");
      options.genotypingModel = genotypingModelEnum( gmString );
    }        
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
    if (getOptionValues(options, parser) != 0)
        return 1;

    return 0;
}

// ==========================================================================

#endif // #ifndef POPINS_CLP_H_
