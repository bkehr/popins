#ifndef NOVINS_CROP_UNMAPPED_H_
#define NOVINS_CROP_UNMAPPED_H_

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include "adapter_removal.h"


using namespace seqan;


// --------------------------------------------------------------------------
// Function hasLowMappingQuality()
// --------------------------------------------------------------------------

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
    if (record.rID != record.rNextId || abs(record.beginPos - record.pNext) >= 1000 || hasFlagRC(record) == hasFlagNextRC(record))
        return true;

    // Check for non-human chromosome IDs - OUT OF DATE SINCE POPINS-v1.0.1?
    if (record.rID > humanSeqs)
        return false;

    // Check for less than 50 bp matches ('M') in cigar string.
    unsigned matches = 0;
    TIter itEnd = end(record.cigar);
    for (TIter it = begin(record.cigar); it != itEnd; ++it)
        if ((*it).operation == 'M') matches += (*it).count;
    if (matches < 50) return true;

    // Check for soft-clipping at BOTH ENDS by more than 24 bp.
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

// --------------------------------------------------------------------------
// Function removeLowQuality()
// --------------------------------------------------------------------------

template<typename TSize_>
inline bool
removeLowQuality(BamAlignmentRecord & record, TSize_ qualThresh)
{
    typedef Iterator<CharString, Rooted>::Type TIter;
    typedef Size<CharString>::Type TSize;

    TSize windowSize = std::max(TSize(5), length(record.qual) / 10);
    TSize windowThresh = qualThresh*windowSize;

    // Initialize windowQual with first windowSize quality values.
    TSize windowQual = 0;
    TIter qualEnd = end(record.qual);
    TIter windowEnd = begin(record.qual) + std::min(windowSize, length(record.qual));
    TIter windowBegin = begin(record.qual);
    for (; windowBegin != windowEnd; ++windowBegin)
        windowQual += *windowBegin - 33;

    // Check quality from the left.
    for (windowBegin = begin(record.qual); windowEnd < qualEnd; ++windowEnd, ++windowBegin)
    {
        if (windowQual >= (TSize)windowThresh)
        {
            while (*windowBegin - 33 < qualThresh) ++windowBegin;
            record.seq = suffix(record.seq, position(windowBegin));
            record.qual = suffix(record.qual, position(windowBegin));
            break;
        }

        windowQual -= *windowBegin - 33;
        windowQual += *windowEnd - 33;
    }
    if (windowEnd == qualEnd) return 1;

    // Initialize windowQual with last windowSize quality values.
    windowQual = 0;
    TIter qualBegin = begin(record.qual);
    windowEnd = end(record.qual) - 1;
    windowBegin = windowEnd - std::min(windowSize, length(record.qual));
    for (; windowEnd != windowBegin; --windowEnd)
        windowQual += *windowEnd - 33;

    // Check quality from the right.
    for (windowEnd = end(record.qual) - 1; windowBegin >= qualBegin; --windowBegin, --windowEnd)
    {
        if (windowQual >= (TSize)windowThresh)
        {
            while (*windowEnd - 33 < qualThresh) --windowEnd;
            record.seq = prefix(record.seq, position(windowEnd) + 1);
            record.qual = prefix(record.qual, position(windowEnd) + 1);
            break;
        }

        windowQual -= *windowEnd - 33;
        windowQual += *windowBegin -33;
    }

    if (length(record.seq) < 30) return 1;
    return 0;
}

// --------------------------------------------------------------------------
// Functions setUnmapped() and setMateUnmapped()
// --------------------------------------------------------------------------

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
    //clear(record.tags);
}

// --------------------------------------------------------------------------

inline void
setMateUnmapped(BamAlignmentRecord & record)
{
    record.flag |= BAM_FLAG_NEXT_UNMAPPED;
    record.flag &= ~BAM_FLAG_ALL_PROPER;
    record.rNextId = record.rID;
    record.pNext = record.beginPos;
    record.tLen = BamAlignmentRecord::INVALID_LEN;
}

// --------------------------------------------------------------------------
// Function appendFastqRecord()
// --------------------------------------------------------------------------

// Append a read to map of fastq records.
bool
appendFastqRecord(SeqFileOut & firstStream,
        SeqFileOut & secondStream,
        std::map<CharString, Pair<CharString> > & firstReads,
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
        if (secondReads.count(record.qName) != 0)
        {
            Pair<CharString> second = secondReads[record.qName];
            writeRecord(firstStream, record.qName, record.seq, record.qual);
            writeRecord(secondStream, record.qName, second.i1, second.i2);
            secondReads.erase(record.qName);
            return 1;
        }
        else
        {
            firstReads[record.qName] = Pair<CharString>(seq, qual);
            return 0;
        }
    }
    else // hasFlagLast(record)
    {
        if (firstReads.count(record.qName) != 0)
        {
            Pair<CharString> first = firstReads[record.qName];
            writeRecord(firstStream, record.qName, first.i1, first.i2);
            writeRecord(secondStream, record.qName, record.seq, record.qual);
            firstReads.erase(record.qName);
            return 1;
        }
        else
        {
            secondReads[record.qName] = Pair<CharString>(seq, qual);
            return 0;
        }
    }
}

// --------------------------------------------------------------------------
// Function writeFastq()
// --------------------------------------------------------------------------

int
writeFastq(SeqFileOut & fastqFirst,
        SeqFileOut & fastqSecond,
        SeqFileOut & fastqSingle,
        std::map<CharString, Pair<CharString> > const & firstReads,
        std::map<CharString, Pair<CharString> > const & secondReads)
{
    typedef std::map<CharString, Pair<CharString> > TFastqMap;

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
            writeRecord(fastqSingle, firstIt->first, firstIt->second.i1, firstIt->second.i2);
            ++firstIt;
        }
        else if (firstIt->first == secondIt->first)
        {
            writeRecord(fastqFirst, firstIt->first, firstIt->second.i1, firstIt->second.i2);
            writeRecord(fastqSecond, secondIt->first, secondIt->second.i1, secondIt->second.i2);
            ++firstIt; ++secondIt;
        }
        else // firstIt->first > secondIt->first
        {
            writeRecord(fastqSingle, secondIt->first, secondIt->second.i1, secondIt->second.i2);
            ++secondIt;
        }
    }

    // Iterate over remaining reads and output to single.fastq.
    while (firstIt != firstEnd)
    {
        writeRecord(fastqSingle, firstIt->first, firstIt->second.i1, firstIt->second.i2);
        ++firstIt;
    }
    while (secondIt != secondEnd)
    {
        writeRecord(fastqSingle, secondIt->first, secondIt->second.i1, secondIt->second.i2);
        ++secondIt;
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function findOtherReads()
// --------------------------------------------------------------------------

template<typename TPos>
int
findOtherReads(BamFileOut & matesStream,
        std::map<Pair<TPos>, Pair<CharString, bool> > & otherReads,
        CharString const & mappingBam)
{
    typedef std::map<Pair<TPos>, Pair<CharString, bool> > TOtherMap;

    int numFound = 0; // Return value.

    // Open input file.
    BamFileIn inStream(toCString(mappingBam));
    BamHeader header;
    readHeader(header, inStream);

    // Load bam index.
    CharString baiFile = mappingBam;
    baiFile += ".bai";
    BamIndex<Bai> bamIndex;
    if (!open(bamIndex, toCString(baiFile)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << baiFile << std::endl;
        return -1;
    }

    __int32 rID = BamAlignmentRecord::INVALID_REFID;
    BamAlignmentRecord record;

    typename TOtherMap::const_iterator itEnd = otherReads.end();
    for (typename TOtherMap::const_iterator it = otherReads.begin(); it != itEnd; ++it)
    {
        if (rID != it->first.i1)
        {
            // Jump to chromosome.
            rID = it->first.i1;
            bool hasAligns;
            jumpToRegion(inStream, hasAligns, rID, it->first.i2, maxValue<TPos>(), bamIndex);
            if (!hasAligns) continue;
            readRecord(record, inStream);
        }

        // Skip reads not in list.
        bool last = false;
        while (record.rID == it->first.i1 && (record.beginPos < it->first.i2 || (record.beginPos == it->first.i2 && record.qName != it->second.i1)))
        {
            if (atEnd(inStream))
            {
                last = true;
                break;
            }
            readRecord(record, inStream);
        }

        // Output record if it matches qName, rID, and beginPos.
        if (!last && record.qName == it->second.i1 && record.rID == it->first.i1 && record.beginPos == it->first.i2)
        {
            // Check if both ends are low-quality mapped and, hence, are already in fastq files.
            if (otherReads.count(Pair<TPos>(record.rNextId, record.pNext)) == 0)
            {
                setMateUnmapped(record);
                writeRecord(matesStream, record);
            }
            ++numFound;
        }
    }

    return numFound;
}

// ==========================================================================
// Function crop_unmapped()
// ==========================================================================

template<typename TAdapterTag>
int
crop_unmapped(double & avgCov,
        Triple<CharString> & fastqFiles,
        CharString & matesBam,
        CharString const & mappingBam,
        int humanSeqs,
        TAdapterTag tag)
{
    typedef __int32 TPos;
    typedef std::map<CharString, Pair<CharString> > TFastqMap; // Reads to go into fastq files.
    typedef std::map<Pair<TPos>, Pair<CharString, bool> > TOtherMap; // Reads to crop in a second pass of the input file.
    typedef StringSet<Dna5String> TStringSet;

    // Open the input and output bam files.
    BamFileIn inStream(toCString(mappingBam));
    BamFileOut matesStream(context(inStream), toCString(matesBam));

    // Copy the header.
    BamHeader header;
    readHeader(header, inStream);
    writeHeader(matesStream, header);

    unsigned long genomeLength = 0;
    for (unsigned i = 0; i < length(header); ++i)
    {
        if (header[i].type != BamHeaderRecordType::BAM_HEADER_REFERENCE)
            continue;

        for (unsigned j = 0; j < length(header[i].tags); ++j)
        {
            if (header[i].tags[j].i1 != "LN")
                continue;

            genomeLength += lexicalCast<unsigned>(header[i].tags[j].i2);
            break;
        }
    }

    // Create maps for fastq records (first read in pair and second read in pair) and bam records without mate.
    TFastqMap firstReads, secondReads;
    TOtherMap otherReads;

    // Open the output fastq files.
    SeqFileOut fastqFirstStream(toCString(fastqFiles.i1));
    SeqFileOut fastqSecondStream(toCString(fastqFiles.i2));
    SeqFileOut fastqSingleStream(toCString(fastqFiles.i3));

    // Retrieve the adapter sequences with up to one error and create indices.
    TStringSet universal = reverseUniversalOneError(tag);
    TStringSet truSeqs = reverseTruSeqsOneError(tag);
    Index<TStringSet> indexUniversal(universal);
    Index<TStringSet> indexTruSeqs(truSeqs);

    // Iterate over the input file.
    BamAlignmentRecord record;
    unsigned long alignedBaseCount = 0;
    while (!atEnd(inStream))
    {
        // Read the next read from input file.
        readRecord(record, inStream);

        // Check for flags that indicate 'uninteresting' bam records.
        if (hasFlagDuplicate(record) or hasFlagSecondary(record) or
                hasFlagQCNoPass(record) or hasFlagSupplementary(record)) continue;

        if (!hasFlagUnmapped(record))
            alignedBaseCount += length(record.seq);

        // Check the read's unmapped flag.
        if (hasFlagUnmapped(record))
        {
            if (removeLowQuality(record, 20) != 1 && removeAdapter(record, indexUniversal, indexTruSeqs, 30, tag) != 2)
                appendFastqRecord(fastqFirstStream, fastqSecondStream, firstReads, secondReads, record);
        }

        // Check the mate's unmapped flag. It's important to check this BEFORE testing hasLowMappingQuality(),
        // since hasLowMappingQuality() assumes both PE reads to be mapped.
        else if (hasFlagNextUnmapped(record))
        {
            writeRecord(matesStream, record);
        }

        // Check for low mapping quality.
        else if (hasLowMappingQuality(record, humanSeqs))
        {
            if (removeLowQuality(record, 20) != 1 && removeAdapter(record, indexUniversal, indexTruSeqs, 30, tag) != 2)
            {
                if (appendFastqRecord(fastqFirstStream, fastqSecondStream, firstReads, secondReads, record) == 0)
                    otherReads[Pair<TPos>(record.rNextId, record.pNext)] = Pair<CharString, bool>(record.qName, hasFlagFirst(record));
            }
        }
    }
    close(inStream);

    avgCov = (double)alignedBaseCount / (double)genomeLength;

    std::ostringstream msg;
    msg << "Map of low quality mates has " << otherReads.size() << " records.";
    printStatus(msg);

    // Write the remaining fastq records.
    if (writeFastq(fastqFirstStream, fastqSecondStream, fastqSingleStream, firstReads, secondReads) != 0) return 1;
    firstReads.clear();
    secondReads.clear();

    msg.str("");
    msg << "Unmapped reads written to " << fastqFiles.i1 << ", " << fastqFiles.i2 << ", " << fastqFiles.i3;
    printStatus(msg);

    // Find the other read end of the low quality mapping reads and write them to the output bam file.
    int found = findOtherReads(matesStream, otherReads, mappingBam);
    if (found == -1) return 1;

    msg.str("");
    msg << "Mapped mates of unmapped reads written to " << matesBam << " , " << found << " found in second pass.";
    printStatus(msg);

    return 0;
}


template<typename TAdapterTag>
int
crop_unmapped(Triple<CharString> & fastqFiles,
        CharString & matesBam,
        CharString const & mappingBam,
        int humanSeqs,
        TAdapterTag tag)
{
    double cov;
    return crop_unmapped(cov, fastqFiles, matesBam, mappingBam, humanSeqs, tag);
}

#endif // #ifndef NOVINS_CROP_UNMAPPED_H_
