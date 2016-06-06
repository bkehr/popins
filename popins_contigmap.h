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
        std::cerr << "ERROR: Could not open bwa output file " << inFile << "" << std::endl;
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

    // Create names of temporary files.
    CharString mappedSam = getFileName(options.workingDirectory, "contig_mapped_unsorted.sam");
    CharString mappedBamUnsorted = getFileName(options.workingDirectory, "contig_mapped_unsorted.bam");
    CharString mappedBam = getFileName(options.workingDirectory, "contig_mapped.bam");
    CharString mergedBam = getFileName(options.workingDirectory, "merged.bam");

    std::stringstream cmd;

    std::ostringstream msg;
    msg << "Mapping reads to contigs using " << BWA;
    printStatus(msg);

    // Remapping to contigs with bwa.
    cmd.str("");
    if (options.allAlignment) cmd << BWA << " mem -a ";
    else cmd << BWA << " mem ";
    cmd << "-t " << options.threads << " " << options.contigFile << " " << fastqFirst << " " << fastqSecond << " > " << mappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFirst << " and " << fastqSecond << std::endl;
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
        return 1;
    }
    //remove(toCString(fastqSingle));

    printStatus("Filling in sequences of secondary records in bwa output");

    // Fill in sequences in bwa output.
    if (fill_sequences(mappedBamUnsorted, mappedSam) != 0)
    {
        return 1;
    }
    remove(toCString(mappedSam));

    msg.str("");
    msg << "Sorting " << mappedBamUnsorted << " by read name using " << SAMTOOLS;
    printStatus(msg);

    // Sort <WD>/contig_mapped.bam by read name
    cmd.str("");
    cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " -o " << options.workingDirectory << "/contig_mapped.bam" << " " << mappedBamUnsorted;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << mappedBamUnsorted << " by read name using " << SAMTOOLS << std::endl;
        return 1;
    }
    remove(toCString(mappedBamUnsorted));

    // Merge non_ref.bam with contig_mapped and set the mates.
    if (merge_and_set_mate(mergedBam, nonRefBam, mappedBam) != 0)
        return 1;

    remove(toCString(mappedBam));
    //remove(toCString(nonRefBam));

    msg.str("");
    msg << "Sorting " << mergedBam << " using " << SAMTOOLS;
    printStatus(msg);

    // Sort <WD>/merged.bam by beginPos, output is <WD>/non_ref.bam.
    cmd.str("");
    cmd << SAMTOOLS << " sort -@ " << options.threads << " -m " << options.memory << " -o " << options.workingDirectory << "/non_ref_new.bam" << " " << mergedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << mergedBam << " by beginPos using " << SAMTOOLS << std::endl;
        return 1;
    }
    remove(toCString(mergedBam));

    msg.str("");
    msg << "Indexing " << nonRefNew << " by beginPos using " << SAMTOOLS;
    printStatus(msg);

    // Index <WD>/non_ref_new.bam.
    cmd.str("");
    cmd << SAMTOOLS << " index " << nonRefNew;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while indexing " << nonRefNew << " using " << SAMTOOLS << std::endl;
        return 1;
    }

    msg.str("");
    msg << "Computing contig locations from anchoring reads in " << nonRefNew;
    printStatus(msg);

    // Find anchoring locations of contigs for this individual.
    String<Location> locations;
    findLocations(locations, nonRefNew, options.maxInsertSize);
    scoreLocations(locations);
    if (writeLocations(locationsFile, locations) != 0) return 1;

    // Remove the non_ref_new.bam file.
    if (!options.keepNonRefNew)
        remove(toCString(nonRefNew));

    return 0;
}

#endif // #ifndef POPINS_CONTIGMAP_H_

