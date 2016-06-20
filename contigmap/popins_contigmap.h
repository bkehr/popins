#ifndef POPINS_CONTIGMAP_H_
#define POPINS_CONTIGMAP_H_

#include <sstream>

#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include "../popins_utils.h"
#include "../command_line_parsing.h"
#include "../assemble/crop_unmapped.h"
#include "../place/location.h"

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
    BamFileIn inStream(toCString(unmappedBam));

    // Read bam header and clear it since we don't need it.
    BamHeader header;
    readHeader(header, inStream);
    clear(header);

    // Open the output fastq files.    
    SeqFileOut fastqFirstStream(toCString(fastqFirst));
    SeqFileOut fastqSecondStream(toCString(fastqSecond));
    SeqFileOut fastqSingleStream(toCString(fastqSingle));

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

    BamFileIn inStream(toCString(inFile));
    BamFileOut outStream(context(inStream), toCString(outFile));

    BamHeader header;
    readHeader(header, inStream);
    writeHeader(outStream, header);

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

    return 0;
}


// ==========================================================================
// Function popins_contigmap()
// ==========================================================================

int popins_contigmap(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    ContigMapOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    CharString workingDirectory = getFileName(options.prefix, options.sampleID);

    // Check for input files to exist.
    CharString fastqFirst = getFileName(workingDirectory, "paired.1.fastq");
    CharString fastqSecond = getFileName(workingDirectory, "paired.2.fastq");
    CharString fastqSingle = getFileName(workingDirectory, "single.fastq");
    CharString nonRefBam = getFileName(workingDirectory, "non_ref.bam");
    CharString nonRefNew = getFileName(workingDirectory, "non_ref_new.bam");
    CharString locationsFile = getFileName(workingDirectory, "locations.txt");

    if (!exists(fastqFirst) || !exists(fastqSecond) || !exists(fastqSingle) || !exists(nonRefBam))
    {
        std::cerr << "ERROR: Could not find all input files ";
        std::cerr << fastqFirst << ", " << fastqSecond << ", " << fastqSingle << ", and " << nonRefBam << std::endl;
        return 7;
    }

    // Create names of temporary files.
    CharString mappedSam = getFileName(workingDirectory, "contig_mapped_unsorted.sam");
    CharString mappedBamUnsorted = getFileName(workingDirectory, "contig_mapped_unsorted.bam");
    CharString mappedBam = getFileName(workingDirectory, "contig_mapped.bam");
    CharString mergedBam = getFileName(workingDirectory, "merged.bam");

    std::stringstream cmd;

    CharString indexFile = options.contigFile;
    indexFile += ".bwt";
    if (!exists(indexFile))
    {
        std::ostringstream msg;
        msg << "Indexing contigs int \'" << options.contigFile << "\'using " << BWA;
        printStatus(msg);

        cmd.str("");
        cmd << BWA << " index " << options.contigFile;
        if (system(cmd.str().c_str()) != 0)
        {
            std::cerr << "ERROR while indexing \'" << options.contigFile << "\' using " << BWA << std::endl;
            return 7;
        }
    }

    std::ostringstream msg;
    msg << "Mapping reads to contigs using " << BWA;
    printStatus(msg);

    // Remapping to contigs with bwa.
    cmd.str("");
    if (!options.bestAlignment) cmd << BWA << " mem -a ";
    else cmd << BWA << " mem ";
    cmd << "-t " << options.threads << " " << options.contigFile << " " << fastqFirst << " " << fastqSecond << " > " << mappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFirst << " and " << fastqSecond << std::endl;
        return 7;
    }
    //remove(toCString(fastqFirst));
    //remove(toCString(fastqSecond));

    cmd.str("");
    if (!options.bestAlignment) cmd << BWA << " mem -a ";
    else cmd << BWA << " mem ";
    cmd << "-t " << options.threads << " " << options.contigFile << " " << fastqSingle << " | awk '$1 !~ /@/' >> " << mappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqSingle << std::endl;
        return 7;
    }
    //remove(toCString(fastqSingle));

    printStatus("Filling in sequences of secondary records in bwa output");

    // Fill in sequences in bwa output.
    if (fill_sequences(mappedBamUnsorted, mappedSam) != 0)
    {
        return 7;
    }
    remove(toCString(mappedSam));

    msg.str("");
    msg << "Sorting " << mappedBamUnsorted << " by read name using " << SAMTOOLS;
    printStatus(msg);

    // Sort <WD>/contig_mapped.bam by read name
    cmd.str("");
    cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " -o " << mappedBam << " " << mappedBamUnsorted;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << mappedBamUnsorted << " by read name using " << SAMTOOLS << std::endl;
        return 7;
    }
    remove(toCString(mappedBamUnsorted));

    // Merge non_ref.bam with contig_mapped and set the mates.
    if (merge_and_set_mate(mergedBam, nonRefBam, mappedBam) != 0)
        return 7;

    remove(toCString(mappedBam));
    //remove(toCString(nonRefBam));

    msg.str("");
    msg << "Sorting " << mergedBam << " using " << SAMTOOLS;
    printStatus(msg);

    // Sort <WD>/merged.bam by beginPos, output is <WD>/non_ref.bam.
    cmd.str("");
    cmd << SAMTOOLS << " sort -@ " << options.threads << " -m " << options.memory << " -o " << nonRefNew << " " << mergedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << mergedBam << " by beginPos using " << SAMTOOLS << std::endl;
        return 7;
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
        return 7;
    }

    msg.str("");
    msg << "Computing contig locations from anchoring reads in " << nonRefNew;
    printStatus(msg);

    // Find anchoring locations of contigs for this individual.
    String<Location> locations;
    findLocations(locations, nonRefNew, options.maxInsertSize);
    scoreLocations(locations);
    if (writeLocations(locationsFile, locations) != 0) return 7;

    // Remove the non_ref_new.bam file.
    if (options.deleteNonRefNew)
        remove(toCString(nonRefNew));

    return 0;
}

#endif // #ifndef POPINS_CONTIGMAP_H_

