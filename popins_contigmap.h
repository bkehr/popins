#ifndef POPINS_CONTIGMAP_H_
#define POPINS_CONTIGMAP_H_

#include <sstream>

#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include "popins_clp.h"
#include "popins_crop_unmapped.h"
#include "popins_location.h"

using namespace seqan;

// ==========================================================================
// Function write_fastq()
// ==========================================================================

bool
write_fastq(CharString & fastqFirst,
            CharString & fastqSecond,
            CharString & fastqSingle,
            CharString & unmappedBam)
{
    typedef std::map<CharString, Pair<CharString> > TFastqMap;
    
    // Create maps for fastq records (first read in pair and second read in pair).
    TFastqMap firstReads, secondReads;

    // Open bam file.
    BamStream inStream(toCString(unmappedBam), BamStream::READ);
    if (!isGood(inStream)) 
    {
        std::cerr << "ERROR while opening input bam file " << unmappedBam << std::endl;
        return 1;
    }

    // Open the output fastq files.    
    SequenceStream fastqFirstStream, fastqSecondStream, fastqSingleStream;
    if (openFastq(fastqFirstStream, fastqFirst) != 0 || openFastq(fastqSecondStream, fastqSecond) != 0 ||
        openFastq(fastqSingleStream, fastqSingle) != 0) return 1;
    
    // Iterate over bam file and append fastq records.
    BamAlignmentRecord record;
    while(!atEnd(inStream))
    {
        readRecord(record, inStream);
        if (hasFlagUnmapped(record))
            appendFastqRecord(fastqFirstStream, fastqSecondStream, firstReads, secondReads, record);
    }

    // Write the fastq files.
    if (writeFastq(fastqFirstStream, fastqSecondStream, fastqSingleStream, firstReads, secondReads) != 0) return 1;

    return 0;
}

// ==========================================================================
// Function fill_sequences()
// ==========================================================================

bool
fill_sequences(CharString & outFile, CharString & inFile)
{
    typedef Position<Dna5String>::Type TPos;

    BamStream inStream(toCString(inFile));
    if (!isGood(inStream))
    {
        std::cerr << "ERROR: Could bwa output file " << inFile << "" << std::endl;
        return 1;
    }
    BamStream outStream(toCString(outFile), BamStream::WRITE);
    if (!isGood(outStream))
    {
        std::cerr << "ERROR: Could not open output file " << outFile << "" << std::endl;
        return 1;
    }
    outStream.header = inStream.header;

    BamAlignmentRecord firstRecord, nextRecord;
    while (!atEnd(inStream))
    {
        readRecord(nextRecord, inStream);
        
        if (firstRecord.qName != nextRecord.qName || hasFlagFirst(firstRecord) != hasFlagFirst(nextRecord))
        {
            // update first record
            firstRecord = nextRecord;
            if (length(firstRecord.seq) == 0 || length(firstRecord.qual) == 0)
            {
                std::cerr << "ERROR: First record of read " << firstRecord.qName << " has no sequence." << std::endl;
                return 1;
            }
        }
        else
        {
            // fill sequence field and quality string
            if (length(nextRecord.seq) == 0 || length(nextRecord.qual) == 0)
            {
                TPos last = length(nextRecord.cigar)-1;
                if (nextRecord.cigar[0].operation == 'H' || nextRecord.cigar[last].operation == 'H')
                {
                    TPos begin = 0;
                    if (nextRecord.cigar[0].operation == 'H')
                        begin = nextRecord.cigar[0].count;

                    TPos end = length(firstRecord.seq);
                    if (nextRecord.cigar[last].operation == 'H')
                        end -= nextRecord.cigar[last].count;

                    nextRecord.seq = infix(firstRecord.seq, begin, end);
                    nextRecord.qual = infix(firstRecord.qual, begin, end);
                }
                else
                {
                    nextRecord.seq = firstRecord.seq;
                    nextRecord.qual = firstRecord.qual;
                }
            }
        }
        
        writeRecord(outStream, nextRecord);
    }

    close(outStream);
    clear(outStream.header);
    clear(inStream.header);

    return 0;
}

inline int
removeTmpDir(CharString & tmpDir, CharString & workingDir)
{
    if (tmpDir != workingDir)
    {
        if (remove(toCString(tmpDir)) != 0)
        {
            std::cerr << "ERROR: Could not remove temporary directory " << tmpDir << std::endl;
            return 1;
        }
        std::cerr << "[" << time(0) << "] " << "Temporary directory " << tmpDir << " removed." << std::endl;
    }
    return 0;
}


// ==========================================================================
// Function popins_contigmap()
// ==========================================================================

int popins_contigmap(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    ContigMapOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    // Check for input files to exist.
    CharString fastqFirst = getFileName(options.workingDirectory, "paired.1.fastq");
    CharString fastqSecond = getFileName(options.workingDirectory, "paired.2.fastq");
    CharString fastqSingle = getFileName(options.workingDirectory, "single.fastq");
    CharString nonRefBam = getFileName(options.workingDirectory, "non_ref.bam");
    CharString nonRefNew = getFileName(options.workingDirectory, "non_ref_new.bam");
    CharString locationsFile = getFileName(options.workingDirectory, "locations.txt");

    if (!exists(fastqFirst) || !exists(fastqSecond) || !exists(fastqSingle) || !exists(nonRefBam))
    {
        std::cerr << "ERROR: Could not find all input files ";
        std::cerr << fastqFirst << ", " << fastqSecond << ", " << fastqSingle << ", and " << nonRefBam << std::endl;
        return 1;
    }

    // Create temporary directory.
    if (options.tmpDir == "")
    {
        options.tmpDir = options.workingDirectory;
    }
    else
    {
        errno = 0;
        char* tempDir = mkdtemp(toCString(options.tmpDir));
        if (errno != 0)
        {
            if (errno == EINVAL) std::cerr << "ERROR: Temporary directory does not end in XXXXXX: " << options.tmpDir << std::endl;
            else std::cerr << "ERROR: Could not create temporary directory at " << options.tmpDir << std::endl;
            return 1;
        }
        options.tmpDir = tempDir;
        std::cerr << "[" << time(0) << "] Using temporary directory " << tempDir << std::endl;
    }

    // Create names of temporary files.
    CharString mappedSam = getFileName(options.tmpDir, "contig_mapped_unsorted.sam");
    CharString mappedBamUnsorted = getFileName(options.tmpDir, "contig_mapped_unsorted.bam");
    CharString mappedBam = getFileName(options.tmpDir, "contig_mapped.bam");
    CharString mergedBam = getFileName(options.tmpDir, "merged.bam");

    // Remapping to contigs with bwa.
    std::cerr << "[" << time(0) << "] Mapping reads to contigs using " << BWA << std::endl;
    std::stringstream cmd;
    cmd.str("");
    if (options.allAlignment) cmd << BWA << " mem -a ";
    else cmd << BWA << " mem ";
    cmd << "-t " << options.threads << " " << options.contigFile << " " << fastqFirst << " " << fastqSecond << " > " << mappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFirst << " and " << fastqSecond << std::endl;
        removeTmpDir(options.tmpDir, options.workingDirectory);
        return 1;
    }
    //remove(toCString(fastqFirst));
    //remove(toCString(fastqSecond));

    cmd.str("");
    if (options.allAlignment) cmd << BWA << " mem -a ";
    else cmd << BWA << " mem ";
    cmd << "-t " << options.threads << " " << options.contigFile << " " << fastqSingle << " | awk '$1 !~ /@/' >> " << mappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqSingle << std::endl;
        removeTmpDir(options.tmpDir, options.workingDirectory);
        return 1;
    }
    //remove(toCString(fastqSingle));

    // Fill in sequences in bwa output.
    std::cerr << "[" << time(0) << "] Filling in sequences of secondary records in bwa output" << std::endl;
    if (fill_sequences(mappedBamUnsorted, mappedSam) != 0)
    {
        removeTmpDir(options.tmpDir, options.workingDirectory);
        return 1;
    }
    remove(toCString(mappedSam));

    // Sort <WD>/contig_mapped.bam by read name
    std::cerr << "[" << time(0) << "] " << "Sorting " << mappedBamUnsorted << " by read name using " << SAMTOOLS << std::endl;
    cmd.str("");
    cmd << SAMTOOLS << " sort -n -m " << options.memory << " " << mappedBamUnsorted << " " << options.tmpDir << "/contig_mapped";
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << mappedBamUnsorted << " by read name using " << SAMTOOLS << std::endl;
        removeTmpDir(options.tmpDir, options.workingDirectory);
        return 1;
    }
    remove(toCString(mappedBamUnsorted));

    // Merge non_ref.bam with contig_mapped and set the mates.
    if (merge_and_set_mate(mergedBam, nonRefBam, mappedBam) != 0)
    {
        removeTmpDir(options.tmpDir, options.workingDirectory);
        return 1;
    }
    remove(toCString(mappedBam));
    //remove(toCString(nonRefBam));

    // Sort <WD>/merged.bam by beginPos, output is <WD>/non_ref.bam.
    std::cerr << "[" << time(0) << "] " << "Sorting " << mergedBam << " using " << SAMTOOLS << std::endl;
    cmd.str("");
    cmd << SAMTOOLS << " sort -m " << options.memory << " " << mergedBam << " " << options.workingDirectory << "/non_ref_new";
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << mergedBam << " by beginPos using " << SAMTOOLS << std::endl;
        removeTmpDir(options.tmpDir, options.workingDirectory);
        return 1;
    }
    remove(toCString(mergedBam));

    // Index <WD>/non_ref_new.bam.
    std::cerr << "[" << time(0) << "] " << "Indexing " << nonRefNew << " by beginPos using " << SAMTOOLS << std::endl;
    cmd.str("");
    cmd << SAMTOOLS << " index " << nonRefNew;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while indexing " << nonRefNew << " using " << SAMTOOLS << std::endl;
        removeTmpDir(options.tmpDir, options.workingDirectory);
        return 1;
    }
    
    // Remove temporary directory.
    if (removeTmpDir(options.tmpDir, options.workingDirectory) != 0)
        return 1;

    // Find anchoring locations of contigs for this individual.
    std::cerr << "[" << time(0) << "] " << "Computing contig locations from anchoring reads in " << nonRefNew << std::endl;
    String<Location> locations;
    findLocations(locations, nonRefNew);
    scoreLocations(locations);
    if (writeLocations(locationsFile, locations) != 0) return 1;

    return 0;
}

#endif // #ifndef POPINS_CONTIGMAP_H_

