#ifndef POPINS_PLACE_SPLIT_ALIGN_H_
#define POPINS_PLACE_SPLIT_ALIGN_H_

#include <seqan/seeds.h>

#include "align_split.h"
#include "popins_location.h"
#include "popins_place_split_align.h"

using namespace seqan;


// ==========================================================================
// Struct BamInfo
// ==========================================================================

struct BamInfo
{
    BamStream stream;
    BamIndex<Bai> index;
    unsigned covThresh;
};

// ==========================================================================
// Struct LocationInfo
// ==========================================================================

template<typename TSeq, typename TPos, typename TSize>
struct LocationInfo
{
    TSeq refSeq;
    TSeq contigSeq;
    TPos refOffset;
    TPos contigOffset;
    TPos concatPos;

    bool highCoverage;
    TSize splitReadCount;
    std::map<Pair<TPos>, TSize> splitPosMap;
    
    LocationInfo() :
        refOffset(0), contigOffset(0), concatPos(0), highCoverage(false), splitReadCount(0)
    {}
};

// ==========================================================================
// Function isCandidateSplitRead()
// ==========================================================================

bool
hasGoodClippedPrefix(BamAlignmentRecord & record, unsigned minQual)
{
    if (length(record.cigar) == 2 && record.cigar[0].operation == 'S')
        return false;

    unsigned windowLength = 10;
    minQual *= windowLength;
    unsigned windowQual = 0;

    unsigned left = 0;
    for (; left < windowLength; ++left)
        windowQual += record.qual[left];

    while (left < length(record.qual) && windowQual < minQual)
    {
        windowQual += record.qual[left] - record.qual[left - windowLength];
        ++left;
    }
    --left; // now points to first good base

    if (left >= windowLength)
    {
        if (left >= length(record.qual) - 30) // 30 or fewer bases are good
            return false;

        unsigned trim = left;
        unsigned i = 0;
        while (i < length(record.cigar) && record.cigar[i].count <= trim)
        {
            trim -= record.cigar[i].count;
            ++i;
        }
        if (length(record.cigar) - i < 2 || (length(record.cigar) - i == 2 && record.cigar[length(record.cigar) - 1].operation == 'S'))
            return false;
    }
    return true;
}

// ==========================================================================

bool
hasGoodClippedSuffix(BamAlignmentRecord & record, unsigned minQual)
{
    if (length(record.cigar) == 2 && record.cigar[1].operation == 'S')
        return false;
    
    unsigned windowLength = 10;
    minQual *= windowLength;
    unsigned windowQual = 0;
        
    int right = length(record.qual) - 1;
    for (; right >= (int)length(record.qual) - (int)windowLength; --right)
        windowQual += record.qual[right];
    
    while (right >= 0 && windowQual < minQual)
    {
        windowQual += record.qual[right] - record.qual[right + windowLength];
         --right;
    }    
    ++right; // now points to last good base

    if ((unsigned)right < length(record.qual) - windowLength)
    {
        if (right < 30) // 30 or fewer bases are good
            return false;
            
        unsigned trim = length(record.qual) - 1 - right;
        int i = length(record.cigar) - 1;
        while (i >= 0 && record.cigar[i].count <= trim)
        {
            trim -= record.cigar[i].count;
            ++i;
        }
        if (i < 2 || (i == 2 && record.cigar[0].operation == 'S'))
            return false;
    }
    return true;
}

// ==========================================================================

bool
isCandidateSplitRead(BamAlignmentRecord & record, bool locOri)
{
    // Check cigar.
    if (length(record.cigar) == 1 || (length(record.cigar) == 3 &&
        record.cigar[1].count < 10 && (record.cigar[1].operation == 'D' || record.cigar[1].operation == 'I')))
        return false; // read matches the reference

    // Check flags.
    if (hasFlagSecondary(record) || hasFlagQCNoPass(record) || hasFlagDuplicate(record) || hasFlagSupplementary(record))
        return false;

    // Check if unmapped and otherwise for quality of interesting read end.
    // Reverse complement the read sequence if necessary.
    if (locOri)
    {
        if (hasFlagUnmapped(record))
        {
            if (!hasFlagNextRC(record))
            {
                if (!hasFlagRC(record))
                {
                    reverseComplement(record.seq);
                    record.flag |= BAM_FLAG_RC;
                }
                return true;
            }
            else
                return false;
        }

        return hasGoodClippedSuffix(record, 53);
    }
    else
    {
        if (hasFlagUnmapped(record))
        {
            if (hasFlagNextRC(record))
            {
                if (hasFlagRC(record))
                {
                    reverseComplement(record.seq);
                    record.flag &= ~BAM_FLAG_RC;
                }
                return true;
            }
            else
                return false;
        }

        return hasGoodClippedPrefix(record, 53);
    }
}

// ==========================================================================
// Function findsplitReads()
// ==========================================================================

template<typename TSeq1, typename TPos, typename TSize, typename TSeq2>
void
alignRead(LocationInfo<TSeq1, TPos, TSize> & locInfo, TSeq2 & read, bool ori)
{
    Gaps<TSeq2> readRowLeft;
    Gaps<TSeq1> refRowLeft;
    Gaps<TSeq2> readRowRight;
    Gaps<TSeq1> refRowRight;
    
    setSource(readRowLeft, read);
    setSource(readRowRight, read);
    
    if (ori)
    {
        setSource(refRowLeft, locInfo.refSeq);
        setSource(refRowRight, locInfo.contigSeq);
    }
    else
    {
        setSource(refRowLeft, locInfo.contigSeq);
        setSource(refRowRight, locInfo.refSeq);
    }

    Score<int, Simple> scoringScheme(1, -2, -5);
    int splitScore = splitAlignment(readRowLeft, refRowLeft, readRowRight, refRowRight, scoringScheme);

    if (splitScore < length(read)*0.7)
        return;

    // Position on the read.
    TPos readPos = toSourcePosition(readRowLeft, clippedEndPosition(readRowLeft));
    unsigned minOverhang = 0.1 * length(read);
    if (readPos < minOverhang || readPos > length(read)-minOverhang)
        return;
    
    // Position on the reference infix and contig.
    TPos refPos;
    TPos contigPos;

    if (ori)
    {
        refPos = toSourcePosition(refRowLeft, clippedEndPosition(refRowLeft));
        contigPos = toSourcePosition(refRowRight, 0);

        // Move the split position to the rightmost possible position.
        while (refPos < length(locInfo.refSeq) && contigPos < length(locInfo.contigSeq) &&
              locInfo.refSeq[refPos] == locInfo.contigSeq[contigPos])
        {
            ++refPos;
            ++contigPos;
        }

        if (refPos > length(locInfo.refSeq)-minOverhang || contigPos < minOverhang)
            return;
    }
    else
    {
        refPos = toSourcePosition(refRowRight, 0);
        contigPos = toSourcePosition(refRowLeft, clippedEndPosition(refRowLeft));

        if (contigPos > length(locInfo.contigSeq)-minOverhang || refPos < minOverhang)
            return;
    }

    // Increase count for this split positions.
    Pair<TPos> posPair(refPos + locInfo.refOffset, contigPos + locInfo.contigOffset);
    if (locInfo.splitPosMap.count(posPair) == 0)
        locInfo.splitPosMap[posPair] = 1;
    else
        ++locInfo.splitPosMap[posPair];

    ++locInfo.splitReadCount;
}

// ==========================================================================

template<typename TSeq1, typename TPos, typename TSize, typename TSeq2>
void
alignRead(LocationInfo<TSeq1, TPos, TSize> & locInfo, TSeq2 & read)
{
    Gaps<TSeq2> readRow;
    Gaps<TSeq1> refRow;
    setSource(readRow, read);
    setSource(refRow, locInfo.refSeq);

    Score<int, Simple> scoringScheme(1, -2, -5);
    int score = globalAlignment(refRow, readRow, scoringScheme, AlignConfig<true, true, true, true>());
    
    if (score < length(read) * 0.7)
        return;
        
    unsigned minOverhang = 0.1 * length(read);
    if (toSourcePosition(refRow, toViewPosition(readRow, 0)) + minOverhang > locInfo.concatPos ||
        toSourcePosition(refRow, toViewPosition(readRow, length(read))) < locInfo.concatPos + minOverhang)
        return;
    
    ++locInfo.splitReadCount;
}

// ==========================================================================

template<typename TSeq, typename TPos, typename TSize>
bool
alignRegion(LocationInfo<TSeq, TPos, TSize> & locInfo,
            BamInfo & bam,
            Location & loc,
            PlacingOptions & options,
            bool isConcat)
{
    // Find the rID in bam file for the locations chromosome.
    int rID = 0;
    BamIOContext<StringSet<CharString> > context = bam.stream.bamIOContext;
    getIdByName(nameStore(context), loc.chr, rID, nameStoreCache(context));

    // Read the first record to find out read length in this bamfile.
    BamAlignmentRecord record;
    if (readRecord(record, bam.stream))
    {
        std::cerr << "ERROR: Failed reading record from bam file." << std::endl;
        return 1;
    }

    // Narrow down the location to relevant read begin positions.
    __int32 locStart = 0, locEnd = 0;
    if (loc.chrOri)
    {
        locStart = loc.chrStart + length(record.seq);
        locEnd = loc.chrEnd + options.maxInsertSize;
    }
    else
    {
        locStart = loc.chrStart - options.maxInsertSize;
        locEnd = loc.chrEnd;
    }
    
    // Set the coverage threshold for this bamfile for this location (to 3 times the average coverage).
    unsigned covThresh = bam.covThresh * (locEnd - locStart) / length(record.seq);

    // Jump to the location in bam file.
    bool hasAlignments;
    jumpToRegion(bam.stream, hasAlignments, rID, locStart, locEnd, bam.index);
    if (!hasAlignments)
        return 0;

    // Iterate reads in region and align candidate split reads.
    unsigned readCount = 0;
    while (!atEnd(bam.stream))
    {
        // Read record from bamfile.
        if (readRecord(record, bam.stream))
        {
            std::cerr << "ERROR: Failed reading record from bam file." << std::endl;
            return 1;
        }

        // Skip records before the region's start.
        if (record.rID == rID && record.beginPos < locStart)
            continue;

        // Check if read's alignment position is still within the location.
        if (record.rID != rID || record.beginPos > locEnd)
            break;
        
        // Check for high coverage.
        ++readCount;
        if (readCount > covThresh)
        {
            locInfo.highCoverage = true;
            return 0;
        }

        // Check quality of record.
        if (!isCandidateSplitRead(record, loc.chrOri))
            continue;

        // (Split-)Align the read to the reference.
        if (isConcat)
            alignRead(locInfo, record.seq);
        else
            alignRead(locInfo, record.seq, loc.chrOri);
    }

    return 0;
}

// ==========================================================================

template<typename TSize, typename TSeq, typename TPos>
bool
findSplitReads(std::map<TSize, LocationInfo<TSeq, TPos, TSize> > & locInfos,
               std::map<TSize, std::set<TSize> > & concatGroups,
               std::map<TSize, std::set<TSize> > & groups,
               String<Location> & locations,
               PlacingOptions & options)
{
    typedef typename std::map<TSize, std::set<TSize> >::iterator TGroupsIter;

    BamInfo bam;

    // Iterate individuals.
    for (unsigned b = 0; b < length(options.bamFiles); ++b)
    {
        if (options.verbose)
            std::cerr << "[" << time(0) << "] " << "Split aligning reads from " << options.bamFiles[b] << std::endl;
        
        // Open the bam file.
        if (open(bam.stream, toCString(options.bamFiles[b])) != 0)
        {
            std::cerr << "ERROR: Could not open bam file " << options.bamFiles[b] << std::endl;
            return 1;
        }

        // Load the bam index.
        CharString baiFile = options.bamFiles[b];
        baiFile += ".bai";
        if (read(bam.index, toCString(baiFile)) != 0)
        {
            std::cerr << "ERROR: Could not read BAI index file " << baiFile << std::endl;
            return 1;
        }

        // Compute max coverage for this bam file.
        bam.covThresh = 3 * options.bamAvgCov[b];

        // Iterate groups of locations with concatenated reference and align reads.
        TGroupsIter cGroupsEnd = concatGroups.end();
        for (TGroupsIter it = concatGroups.begin(); it != cGroupsEnd; ++it)
        {
            TSize g = it->first;
            if (locInfos[g].highCoverage || locInfos[g].splitReadCount >= options.maxSplitReads)
                continue;
            if (alignRegion(locInfos[g], bam, locations[g], options, true) != 0)
                return 1;
        }

        // Iterate groups of locations with separate reference and contig sequence and split-align reads.
        TGroupsIter groupsEnd = groups.end();
        for (TGroupsIter it = groups.begin(); it != groupsEnd; ++it)
        {
            TSize g = it->first;
            if (locInfos[g].highCoverage || locInfos[g].splitReadCount >= options.maxSplitReads)
                continue;
            if (alignRegion(locInfos[g], bam, locations[g], options, false) != 0)
                return 1;
        }
    }

    return 0;
}

// ==========================================================================
// Functions  alignContigPrefixToRef()  and  alignContigSuffixToRef()
// ==========================================================================

template<typename TSeq1, typename TSeq2, typename TScoringScheme>
bool
alignContigPrefixToRef(Pair<size_t> & splitPosPair, TSeq1 & refInfix, TSeq2 & contig, TScoringScheme & scsc)
{
    Gaps<TSeq1> refRow;
    Gaps<typename Prefix<TSeq2>::Type> contigRow;
    
    typename Prefix<TSeq2>::Type contigPrefix = prefix(contig, 50);

    setSource(refRow, refInfix);
    setSource(contigRow, contigPrefix);
    
    int score = localAlignment(refRow, contigRow, scsc);

    if (score < 15)
        return false;

    size_t firstContigPos = toSourcePosition(contigRow, 0);
    
    if (firstContigPos > 5)
        return false;
        
    size_t lastContigPos = toSourcePosition(contigRow, length(refRow));
    size_t lastRefPos = toSourcePosition(refRow, length(refRow));

    splitPosPair = Pair<size_t>(lastRefPos, lastContigPos);

    if (lastContigPos > 45) // 45 = length(contigSuffix) - 5
    {
        size_t firstRefPos = toSourcePosition(refRow, 0);

        Seed<Simple> seed(firstRefPos, firstContigPos, lastRefPos, lastContigPos);
        extendSeed(seed, refInfix, contig, EXTEND_RIGHT, scsc, -2 * scoreGap(scsc), GappedXDrop());
        
        // Re-align seed sequences to get a clean breakpoint position.
        Gaps<typename Infix<TSeq1>::Type> refInfixRow;
        assignSource(refInfixRow, infix(refInfix, beginPositionH(seed), endPositionH(seed)));
        contigPrefix = prefix(contig, endPositionV(seed));
        localAlignment(refInfixRow, contigRow, scsc);
        
        lastContigPos = toSourcePosition(contigRow, length(refInfixRow));
        lastRefPos = beginPositionH(seed) + toSourcePosition(refInfixRow, length(refInfixRow));
    }

    splitPosPair = Pair<size_t>(lastRefPos, lastContigPos);
    return true;
}

// ==========================================================================

template<typename TSeq1, typename TSeq2, typename TScoringScheme>
bool
alignContigSuffixToRef(Pair<size_t> & splitPosPair, TSeq1 & refInfix, TSeq2 & contig, TScoringScheme & scsc)
{
    Gaps<TSeq1> refRow;
    Gaps<typename Suffix<TSeq2>::Type> contigRow;
    
    typename Suffix<TSeq2>::Type contigSuffix = suffix(contig, length(contig) - 50);

    setSource(refRow, refInfix);
    setSource(contigRow, contigSuffix);
    
    int score = localAlignment(refRow, contigRow, scsc);

    if (score < 15)
        return false;

    size_t lastContigPos = toSourcePosition(contigRow, length(refRow));
    
    if (lastContigPos < 45) // 45 = length(contigSuffix) - 5
        return false;
        
    size_t firstContigPos = beginPosition(contigSuffix) + toSourcePosition(contigRow, 0);
    size_t firstRefPos = toSourcePosition(refRow, 0);

    if (firstContigPos < beginPosition(contigSuffix) + 5)
    {
        size_t lastRefPos = toSourcePosition(refRow, length(refRow));
        lastContigPos += beginPosition(contigSuffix);

        Seed<Simple> seed(firstRefPos, firstContigPos, lastRefPos, lastContigPos);
        extendSeed(seed, refInfix, contig, EXTEND_LEFT, scsc, -2 * scoreGap(scsc), GappedXDrop());
        
        // Re-align seed sequences to get a clean breakpoint position.
        Gaps<typename Infix<TSeq1>::Type> refInfixRow;
        assignSource(refInfixRow, infix(refInfix, beginPositionH(seed), toSourcePosition(refRow, length(contigRow))));
        contigSuffix = suffix(contig, beginPositionV(seed));
        localAlignment(refInfixRow, contigRow, scsc);
        
        firstContigPos = beginPosition(contigSuffix) + toSourcePosition(contigRow, 0);
        firstRefPos = beginPositionH(seed) + toSourcePosition(refInfixRow, 0);
    }
    
    splitPosPair = Pair<size_t>(firstRefPos, firstContigPos);
    return true;
}

#endif // #ifndef POPINS_PLACE_SPLIT_ALIGN_H_
