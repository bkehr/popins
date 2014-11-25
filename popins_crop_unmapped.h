#ifndef NOVINS_CROP_UNMAPPED_H_
#define NOVINS_CROP_UNMAPPED_H_

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include "adapter_removal.h"


using namespace seqan;


/**
 * Checks the mapping quality of a read.
 *
 * The mapping quality is considered low if the other read end does not map further
 * than 1000 bp away or in the same orientation and if, in addition, one of the
 * following holds:
 *   - The total number of matches in the cigar string is below 50.
 *   - The read end is soft-clipped by 25 or more bases at both ends.
 *   - The alignment score as indicated by the AS tag is lower than 0.5 * read length.
 *
 * @param record    a read's mapping record from a bam file
 *
 * @returns         true if the read has low mapping quality and otherwise false.
 */
inline bool
hasLowMappingQuality(BamAlignmentRecord & record, int humanSeqs)
{
    typedef Iterator<String<CigarElement<> > >::Type TIter;
    
    // Check for mapping location of other read end. If within 1000 bp and opposite orientation, accept the mapping.
    if (record.rID == record.rNextId && abs(record.beginPos - record.pNext) < 1000 &&
        hasFlagRC(record) != hasFlagNextRC(record))
        return false;

    if (record.rID > humanSeqs) return false;
    
    // Check for less than 50 bp matches ('M') in cigar string.
    unsigned matches = 0;
    TIter itEnd = end(record.cigar);
    for (TIter it = begin(record.cigar); it != itEnd; ++it)
        if ((*it).operation == 'M') matches += (*it).count;
    if (matches < 50) return true;

    // Check for soft-clipping at BOTH ENDS by more than 9 bp.
    if (record.cigar[0].operation == 'S' &&
        record.cigar[0].count > 24 &&
        record.cigar[length(record.cigar)-1].operation == 'S' &&
        record.cigar[length(record.cigar)-1].count > 24)
        return true;
    
    // Check for AS (alignment score) lower than 0.5 * readLength.
    BamTagsDict tagsDict(record.tags);
    unsigned idx;
    if (findTagKey(idx, tagsDict, "AS"))
    {
        unsigned score = 0;
        extractTagValue(score, tagsDict, idx);
        if (score < 0.5*length(record.seq))
            return true;
    }

    return false;
}

inline void
setUnmapped(BamAlignmentRecord & record)
{
    record.flag |= BAM_FLAG_UNMAPPED;
    record.flag &= ~BAM_FLAG_ALL_PROPER;
    record.rID = record.rNextId;
    record.beginPos = record.pNext;
    record.mapQ = 0;
    clear(record.cigar);
    record.tLen = BamAlignmentRecord::INVALID_LEN;
    //clear(reacord.tags);
}

inline void
setMateUnmapped(BamAlignmentRecord & record)
{
    record.flag |= BAM_FLAG_NEXT_UNMAPPED;
    record.flag &= ~BAM_FLAG_ALL_PROPER;
    record.rNextId = record.rID;
    record.pNext = record.beginPos;
    record.tLen = BamAlignmentRecord::INVALID_LEN;
}

// Append a read to map of fastq records.
void
appendFastqRecord(std::map<CharString, Pair<CharString> > & firstReads,
                  std::map<CharString, Pair<CharString> > & secondReads,
                  BamAlignmentRecord const & record)
{
    CharString seq = record.seq;
    CharString qual = record.qual;

    if (hasFlagRC(record))
    {
        reverseComplement(seq);
        reverse(qual);
    }

    if (hasFlagFirst(record))
    {
        if (firstReads.count(record.qName) != 0)
        {
            std::cerr << "[" << time(0) << "] ";
            std::cerr << "WARNING: Multiple records for read " << record.qName << " in bam file." << std::endl;
        }
        firstReads[record.qName] = Pair<CharString>(seq, qual);
    }
    else // hasFlagLast(record)
    {
        if (secondReads.count(record.qName) != 0)
        {
            std::cerr << "[" << time(0) << "] ";
            std::cerr << "WARNING: Multiple records for read " << record.qName << " in bam file." << std::endl;
        }
        secondReads[record.qName] = Pair<CharString>(seq, qual);
    }
}

int
writeFastq(CharString & fastqFirst,
           CharString & fastqSecond,
           CharString & fastqSingle,
           std::map<CharString, Pair<CharString> > const & firstReads,
           std::map<CharString, Pair<CharString> > const & secondReads)
{
    typedef std::map<CharString, Pair<CharString> > TFastqMap;

    // Open the output files.
    SequenceStream fastqFirstStream(toCString(fastqFirst), SequenceStream::WRITE, SequenceStream::FASTQ);
    if (!isGood(fastqFirstStream))
    {
        std::cerr << "ERROR while opening temporary output file " << fastqFirst << std::endl;
        return 1;
    }
    SequenceStream fastqSecondStream(toCString(fastqSecond), SequenceStream::WRITE, SequenceStream::FASTQ);
    if (!isGood(fastqSecondStream))
    {
        std::cerr << "ERROR while opening temporary output file " << fastqFirst << std::endl;
        return 1;
    }
    SequenceStream fastqSingleStream(toCString(fastqSingle), SequenceStream::WRITE, SequenceStream::FASTQ);
    if (!isGood(fastqSingleStream))
    {
        std::cerr << "ERROR while opening temporary output file " << fastqFirst << std::endl;
        return 1;
    }
    
    // Initialize iterators over reads in fastq maps.
    TFastqMap::const_iterator firstIt = firstReads.begin();
    TFastqMap::const_iterator firstEnd = firstReads.end();
    TFastqMap::const_iterator secondIt = secondReads.begin();
    TFastqMap::const_iterator secondEnd = secondReads.end();
    
    // Iterate over reads and output to fastq files (paired.1 and paired.2, or single).
    while (firstIt != firstEnd && secondIt != secondEnd)
    {
        if (firstIt->first < secondIt->first)
        {
            writeRecord(fastqSingleStream, firstIt->first, firstIt->second.i1, firstIt->second.i2);
            ++firstIt;
        }
        else if (firstIt->first == secondIt->first)
        {
            writeRecord(fastqFirstStream, firstIt->first, firstIt->second.i1, firstIt->second.i2);
            writeRecord(fastqSecondStream, secondIt->first, secondIt->second.i1, secondIt->second.i2);
            ++firstIt; ++secondIt;
        }
        else // firstIt->first > secondIt->first
        {
            writeRecord(fastqSingleStream, secondIt->first, secondIt->second.i1, secondIt->second.i2);
            ++secondIt;
        }
    }
    
    // Iterate over remaining reads and output to single.fastq.
    while (firstIt != firstEnd)
    {
        writeRecord(fastqSingleStream, firstIt->first, firstIt->second.i1, firstIt->second.i2);
        ++firstIt;
    }
    while (secondIt != secondEnd)
    {
        writeRecord(fastqSingleStream, secondIt->first, secondIt->second.i1, secondIt->second.i2);
        ++secondIt;
    }

    return 0;
}

template<typename TAdapterTag>
int
crop_unmapped(CharString & fastqFirst,
              CharString & fastqSecond,
              CharString & fastqSingle,
              CharString & unmappedBam,
              CharString const & mappingBam,
              int humanSeqs,
              TAdapterTag tag)
{
    typedef __int32 TPos;
    typedef std::map<CharString, Pair<CharString> > TFastqMap; // Reads to go into fastq files.
    typedef std::map<Pair<TPos>, Pair<CharString, bool> > TOtherMap; // Reads to crop in a second pass of the input file.
    typedef StringSet<Dna5String> TStringSet;

    // Open the input bam file.
    BamStream inStream(toCString(mappingBam));
    BamStream inStream2(toCString(mappingBam));
    if (!isGood(inStream) || !isGood(inStream2))
    {
        std::cerr << "ERROR while opening input bam file " << mappingBam << std::endl;
        return 1;
    }
    
    // Open the bam output file and copy the header from the input file.
    BamStream unmappedBamStream(toCString(unmappedBam), BamStream::WRITE);
    if (!isGood(unmappedBamStream))
    {
        std::cerr << "ERROR while opening output bam file " << unmappedBam << std::endl;
        return 1;
    }
    unmappedBamStream.header = inStream.header;
    
    // Create maps for fastq records (first read in pair and second read in pair) and bam records without mate.
    TFastqMap firstReads, secondReads;
    TOtherMap otherReads;

    // Retrieve the adapter sequences with up to one error.
    TStringSet universal = complementUniversalOneError(tag);
    TStringSet truSeqs = reverseTruSeqsOneError(tag);
    
    // Create suffix index of the adapter sequences.
    Index<TStringSet> indexUniversal(universal);
    Index<TStringSet> indexTruSeqs(truSeqs);
    
    // Iterate over the input file.
    BamAlignmentRecord record;
    while (!atEnd(inStream))
    {
        // Read the next read from input file.
        if (readRecord(record, inStream) != 0)
        {
            std::cerr << "ERROR while reading bam record from " << mappingBam << std::endl;
            return 1;
        }
        
        // Check for flags that indicate 'uninteresting' bam records.
        if (hasFlagDuplicate(record) or hasFlagSecondary(record) or
            hasFlagQCNoPass(record) or hasFlagSupplementary(record)) continue;

        // Check the read's unmapped flag.
        if (hasFlagUnmapped(record))
        {
            // Remove adapter sequences.
            if (removeAdapter(record, indexUniversal, indexTruSeqs, 30, tag) != 2)
            {
                //----- writeRecord(unmappedBamStream, record);
                // Append read to maps of fastq records.
                appendFastqRecord(firstReads, secondReads, record);
            }
        }

        // Check the mate's unmapped flag.
        else if (hasFlagNextUnmapped(record))
        {
            // Remove adapter sequences.
            if (removeAdapter(record, indexUniversal, indexTruSeqs, 30, tag) != 2)
            {
                // Write the read to the output bam file.
                writeRecord(unmappedBamStream, record);
            }
        }
        
        // Check for mapping quality.
        else if (hasLowMappingQuality(record, humanSeqs))
        {
            // Remove adapter sequences.
            if (removeAdapter(record, indexUniversal, indexTruSeqs, 30, tag) != 2)
            {
                // Set read unmapped and append to maps of fastq records.
                setUnmapped(record);
                //----- writeRecord(unmappedBamStream, record);
                appendFastqRecord(firstReads, secondReads, record);
            
                // Append the read to the list of reads that have a mate mapping elsewhere.
                otherReads[Pair<TPos>(record.rNextId, record.pNext)] = Pair<CharString, bool>(record.qName, hasFlagFirst(record));
            }
        }
    }
    
    std::cerr << "[" << time(0) << "] Map of low quality mates has " << otherReads.size() << " records." << std::endl;
    
    // Write the (temporary) fastq files.
    if (writeFastq(fastqFirst, fastqSecond, fastqSingle, firstReads, secondReads) != 0) return 1;
    firstReads.clear();
    secondReads.clear();
    
    std::cerr << "[" << time(0) << "] Unmapped reads written to ";
    std::cerr << fastqFirst << ", " << fastqSecond << ", " << fastqSingle << std::endl;
    
    
    // Find the other read end of the low quality mapping reads and write them to the output bam file. 
    int found = 0;
    if (readRecord(record, inStream2) != 0)
    {
        std::cerr << "ERROR while reading bam record from " << mappingBam << std::endl;
        return 1;
    }
    TOtherMap::const_iterator itEnd = otherReads.end();
    for (TOtherMap::const_iterator it = otherReads.begin(); it != itEnd; ++it)
    {
        // Skip reads not in list.
        while (!atEnd(inStream2) &&
              (record.rID < it->first.i1 ||
               (record.rID == it->first.i1 && record.beginPos < it->first.i2) || 
               (record.rID == it->first.i1 && record.beginPos == it->first.i2 && record.qName != it->second.i1)))
        {
            if (readRecord(record, inStream2) != 0)
            {
                std::cerr << "ERROR while reading bam record from " << mappingBam << std::endl;
                return 1;
            }
        }
        
        // Output record if it matches qName, rID, and beginPos.
        if (!atEnd(inStream2) &&
            record.qName == it->second.i1 && record.rID == it->first.i1 && record.beginPos == it->first.i2)
        {
            if (otherReads.count(Pair<TPos>(record.rNextId, record.pNext)) == 0 &&
                removeAdapter(record, indexUniversal, indexTruSeqs, 30, tag) != 2)
            {
                // Set mate flags and location as unmapped and write record to output bam file.
                setMateUnmapped(record);
                writeRecord(unmappedBamStream, record);
                ++found;
            }
        }
    }
    std::cerr << "[" << time(0) << "] Mapped mates of unmapped reads written to " << unmappedBam << " , " << found << " found in second pass" << std::endl;
    
    return 0;
}

#endif // #ifndef NOVINS_CROP_UNMAPPED_H_
