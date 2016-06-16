#ifndef POPINS_PLACE_SPLIT_ALIGN_H_
#define POPINS_PLACE_SPLIT_ALIGN_H_

#include <seqan/bam_io.h>

#include "location.h"
#include <seqan/align_split.h>

using namespace seqan;

// ---------------------------------------------------------------------------------------
// Function hasGoodClippedPrefix()
// ---------------------------------------------------------------------------------------

template<typename TQualString, typename TCigar>
bool
hasGoodClippedPrefix(TQualString & qual, TCigar & cigar, unsigned minQual)
{
    typedef typename Iterator<TQualString, Rooted>::Type TIter;

    if (cigar[0].operation != 'S')
        return false;

    unsigned windowLength = 10;
    minQual *= windowLength;
    unsigned windowQual = 0;

    TIter it = begin(qual);
    TIter itEnd = end(qual);

    for (unsigned i = 0; i < windowLength; ++i, ++it)
        windowQual += *it;

    TIter it2 = begin(qual, Rooted());
    while (it != itEnd && windowQual < minQual)
    {
        windowQual += (*it) - (*it2);
        ++it;
        ++it2;
    }

    // *it2 now points to the first base in the first window with windowQual >= minQual
    unsigned goodPos = position(it2);

    if (goodPos >= length(qual) - 30) // 30 or fewer bases of the read are good
        return false;

    unsigned i = 0;
    while (i < length(cigar) && cigar[i].count <= goodPos)
    {
        goodPos -= cigar[i].count;
        ++i;
    }
    if (length(cigar) - i < 2) // Only matching bases of the read are good.
        return false;

    return true;
}

// ---------------------------------------------------------------------------------------
// Function isCandidateSplitRead()
// ---------------------------------------------------------------------------------------

// Returns true if read is unmapped and it's mate is mapped in correct orientation or
//              if the read is soft-clipped and the clipped prefix/suffix is of good quality.

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

        ModifiedString<CharString, ModReverse> revQual(record.qual);
        ModifiedString<String<CigarElement<> >, ModReverse> revCigar(record.cigar);
        return hasGoodClippedPrefix(revQual, revCigar, 53);
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

        return hasGoodClippedPrefix(record.qual, record.cigar, 53);
    }
}

// ---------------------------------------------------------------------------------------
// Function getSplitPosition()
// ---------------------------------------------------------------------------------------

// Case 1: RIGHT end of insertion and contig in forward orientation (contigSuffix)

std::pair<unsigned, unsigned>
getSplitPosition(Gaps<typename Infix<Dna5String>::Type> & contigGaps,
        Gaps<Dna5String> & refGaps,
        unsigned offset)
        {
    unsigned refPos = offset + toSourcePosition(refGaps, 0);
    unsigned contigPos = length(host(source(contigGaps))) - length(source(contigGaps)) + toSourcePosition(contigGaps, clippedEndPosition(contigGaps));
    return std::pair<unsigned, unsigned>(refPos, contigPos);
        }

// ---------------------------------------------------------------------------------------

// Case 2: LEFT end of insertion and contig in reverse orientation (contigSuffix)  -->  all in RC

std::pair<unsigned, unsigned>
getSplitPosition(Gaps<typename Infix<Dna5String>::Type> & contigGaps,
        Gaps<ModifiedString<ModifiedString<Dna5String, ModComplementDna5>, ModReverse> > & refGaps,
        unsigned offset)
        {
    unsigned refPos = offset + length(source(refGaps)) - toSourcePosition(refGaps, 0) - 1;
    unsigned contigPos = length(source(contigGaps)) - toSourcePosition(contigGaps, clippedEndPosition(contigGaps));
    return std::pair<unsigned, unsigned>(refPos, contigPos);
        }

// ---------------------------------------------------------------------------------------

// Case 3: RIGHT end of insertion and contig in reverse orientation (contigPrefix)

std::pair<unsigned, unsigned>
getSplitPosition(Gaps<ModifiedString<ModifiedString<typename Infix<Dna5String>::Type, ModComplementDna5>, ModReverse> > & contigGaps,
        Gaps<Dna5String> & refGaps,
        unsigned offset)
        {
    unsigned refPos = offset + toSourcePosition(refGaps, 0);
    unsigned contigPos = length(host(host(host(source(contigGaps))))) - length(source(contigGaps)) + toSourcePosition(contigGaps, clippedEndPosition(contigGaps));
    return std::pair<unsigned, unsigned>(refPos, contigPos);
        }

// ---------------------------------------------------------------------------------------

// Case 4: LEFT end of insertion and contig in forward orientation (contigPrefix)  -->  all in RC

std::pair<unsigned, unsigned>
getSplitPosition(Gaps<ModifiedString<ModifiedString<typename Infix<Dna5String>::Type, ModComplementDna5>, ModReverse> > & contigGaps,
        Gaps<ModifiedString<ModifiedString<Dna5String, ModComplementDna5>, ModReverse> > & refGaps,
        unsigned offset)
        {
    unsigned refPos = offset + length(source(refGaps)) - toSourcePosition(refGaps, 0) - 1;
    unsigned contigPos = length(source(contigGaps)) - toSourcePosition(contigGaps, clippedEndPosition(contigGaps));
    return std::pair<unsigned, unsigned>(refPos, contigPos);
        }

// ---------------------------------------------------------------------------------------
// Function alignRead()
// ---------------------------------------------------------------------------------------

template<typename TContigSeq, typename TRefSeq>
bool
alignRead(std::pair<unsigned, unsigned> & insPos, Dna5String readSeq, TContigSeq & contigPrefix, TRefSeq & ref, unsigned refOffset)
{

    Gaps<Dna5String> readRowLeft;
    Gaps<TContigSeq> contigRowLeft;
    Gaps<Dna5String> readRowRight;
    Gaps<TRefSeq> refRowRight;

    setSource(readRowLeft, readSeq);
    setSource(contigRowLeft, contigPrefix);
    setSource(readRowRight, readSeq);
    setSource(refRowRight, ref);

    Score<int, Simple> scoringScheme(1, -3, -4, -5);
    AlignConfig<false, true, true, false> config;
    Pair<int, int> splitScore = splitAlignment(readRowLeft, contigRowLeft, readRowRight, refRowRight, scoringScheme, config);

    //    std::cout << "\nSplit score = " << splitScore.i1 + splitScore.i2 << " (" << splitScore.i1 << " + " << splitScore.i2 << ")" << std::endl;
    //    std::cout << readRowLeft << "\n" << contigRowLeft << "\n\n" << readRowRight << "\n" << refRowRight << std::endl;

    int minOverhang = 0.1 * length(readSeq);
    if (splitScore.i1 < minOverhang || splitScore.i2 < minOverhang || splitScore.i1 + splitScore.i2 < length(readSeq) * 0.5)
        return 1;

    //    // Split position on the read.
    //    unsigned readPos = toSourcePosition(readRowLeft, clippedEndPosition(readRowLeft));
    //    std::cout << "read pos = " << readPos << "   read len =  " << length(readSeq) << std::endl;

    // Find the split position on the reference and contig.
    insPos = getSplitPosition(contigRowLeft, refRowRight, refOffset);
    //    std::cout << "refPos = " << insPos.first << "   contigPos = " << insPos.second << std::endl;

    return 0;
}

// ---------------------------------------------------------------------------------------
// Function splitAlignReads()
// ---------------------------------------------------------------------------------------

template<typename TContigSeq, typename TRefSeq>
bool
splitAlignReads(std::map<std::pair<unsigned, unsigned>, unsigned> & insPos,
        BamFileIn & bamStream,
        BamIndex<Bai> & bai,
        TContigSeq & contigPrefix,
        TRefSeq & ref,
        Location & loc,
        PlacingOptions & options)
{
    // Set the coverage threshold for this BAM file for this location to 3 times the average coverage.
    unsigned covThresh = 3 * options.bamAvgCov * (loc.chrEnd - loc.chrStart) / options.readLength;

    // Find the rID in BAM file for the location's chromosome.
    int rID = 0;
    getIdByName(rID, contigNamesCache(context(bamStream)), loc.chr);

    // Jump to the location in BAM file.
    bool hasAlignments;
    jumpToRegion(bamStream, hasAlignments, rID, loc.chrStart, loc.chrEnd, bai);

    unsigned readCount = 0;

    // Iterate reads in region and align candidate split reads.
    BamAlignmentRecord record;
    std::pair<unsigned, unsigned> posPair;
    while (!atEnd(bamStream))
    {
        // Read record from BAM file.
        readRecord(record, bamStream);

        // Skip records before the region's start.
        if (record.rID == rID && record.beginPos < loc.chrStart)
            continue;

        // Check if read's alignment position is still within the location.
        if (record.rID != rID || record.beginPos > loc.chrEnd)
            return 0;

        // Check for too high coverage.
        ++readCount;
        if (readCount > covThresh)
            return 1;

        // Check quality of record.
        if (!isCandidateSplitRead(record, loc.chrOri))
            continue;

        // Reverse complement record.seq if (loc.chrOri == true).
        Dna5String readSeq = record.seq;
        if (loc.chrOri)
            reverseComplement(readSeq);

        // (Split-)Align the read to the reference.
        if (alignRead(posPair, readSeq, contigPrefix, ref, loc.chrStart) == 0)
        {
            if (insPos.count(posPair) == 0)
                insPos[posPair] = 1;
            else
                ++insPos[posPair];
        }
    }

    return 0;
}

// ---------------------------------------------------------------------------------------
// Function loadContigAndSplitAlign()
// ---------------------------------------------------------------------------------------

template<typename TRefSeq>
bool
loadContigAndSplitAlign(std::map<std::pair<unsigned, unsigned>, unsigned> & insPos,
        BamFileIn & bamStream,
        BamIndex<Bai> & bai,
        TRefSeq & ref,
        Dna5String & contig,
        Location & loc,
        PlacingOptions & options)
{
    typedef Infix<Dna5String>::Type TInfix;
    typedef ModifiedString<TInfix, ModComplementDna5> TComplementInfix;
    typedef ModifiedString<TComplementInfix, ModReverse> TRcInfix;

    unsigned preSufLen = 200;    // TODO: Make this a program parameter.

    if (loc.contigOri)
    {
        unsigned suffixBegPos = length(contig) - _min(preSufLen, length(contig));
        TInfix contigSuffix = infix(contig, suffixBegPos, length(contig));
        return splitAlignReads(insPos, bamStream, bai, contigSuffix, ref, loc, options);
    }
    else
    {
        unsigned prefixEndPos = _min(preSufLen, length(contig));
        TInfix pref = infix(contig, 0, prefixEndPos);
        TRcInfix contigPrefix(pref);
        return splitAlignReads(insPos, bamStream, bai, contigPrefix, ref, loc, options);
    }

}

// ---------------------------------------------------------------------------------------
// Function loadContigAndSplitAlign()
// ---------------------------------------------------------------------------------------

template<typename TStream>
void
writeLocPos(TStream & outStream,
        Location & loc,
        std::map<std::pair<unsigned, unsigned>, unsigned> & insPos,
        bool highCov)
{
    typedef typename std::map<std::pair<unsigned, unsigned>, unsigned>::iterator TIter;
    outStream << loc.chr;
    if (loc.chr != "OTHER")
    {
        outStream << ":";
        outStream << loc.chrStart << "-";
        outStream << loc.chrEnd;
    }
    outStream << "\t" << (loc.chrOri ? "+" : "-");
    outStream << "\t" << loc.contig;
    outStream << "\t" << (loc.contigOri ? "+" : "-");
    outStream << "\t" << loc.numReads;
    outStream << "\t" << loc.score;

    if (highCov)
    {
        outStream << "\t" << "high_coverage" << std::endl;
    }
    else
    {
        TIter it = insPos.begin();
        TIter itEnd = insPos.end();

        if (it != itEnd)
        {
            outStream << "\t" << (it->first).first << "," << (it->first).second << ":" << (it->second);
            ++it;
        }

        while (it != itEnd)
        {
            outStream << ";" << (it->first).first << "," << (it->first).second << ":" << (it->second);
            ++it;
        }

        outStream << std::endl;
    }
}

// =======================================================================================
// Function popins_place_split_read_align()
// =======================================================================================

bool
popins_place_split_read_align(String<LocationInfo> & locs,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions & options)
{
    typedef typename Iterator<String<LocationInfo> >::Type TIter;
    typedef std::pair<CharString, Dna5String> TPair;

    // Open the BAM file.
    BamFileIn bamStream(toCString(options.bamFile));

    // Read the header and clear it since we don't need it.
    BamHeader header;
    readHeader(header, bamStream);
    clear(header);

    // Load the BAI file.
    BamIndex<Bai> bai;
    CharString baiFile = options.bamFile;
    baiFile += ".bai";
    if (!open(bai, toCString(baiFile)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << baiFile << std::endl;
        return 1;
    }

    // Open the output file.
    std::fstream outStream(toCString(options.outFile), std::ios::out);
    if (!outStream.good())
    {
        std::cerr << "ERROR: Could not open output file " << options.outFile << " for writing." << std::endl;
        return 1;
    }

    std::cerr << "0%   10   20   30   40   50   60   70   80   90   100%" << std::endl;
    std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;
    std::cerr << "*" << std::flush;

    double fiftieth = length(locs) / 50.0;
    unsigned progress = 0;

    TIter it = begin(locs);
    TIter itEnd = end(locs);

    unsigned i = 0;
    while (it != itEnd)
    {
        std::map<std::pair<unsigned, unsigned>, unsigned> insPos;
        bool highCov;

        std::vector<TPair>::iterator contigIt = std::lower_bound(contigs.begin(), contigs.end(), TPair((*it).loc.contig, ""));

        if ((*it).loc.chrOri)
        {
            (*it).loc.chrEnd += options.maxInsertSize;

            // Load the genomic region and reverse complement it.
            Dna5String r = loadInterval(fai, (*it).loc.chr, (*it).loc.chrStart, (*it).loc.chrEnd);
            ModifiedString<ModifiedString<Dna5String, ModComplementDna5>, ModReverse> ref(r);

            // Load the contig prefix/suffix and split align.
            highCov = loadContigAndSplitAlign(insPos, bamStream, bai, ref, contigIt->second, (*it).loc, options);

            (*it).loc.chrEnd -= options.maxInsertSize;
        }
        else
        {
            (*it).loc.chrStart -= options.maxInsertSize;

            // Load the genomic region and keep it in forward orientation.
            Dna5String ref = loadInterval(fai, (*it).loc.chr, (*it).loc.chrStart, (*it).loc.chrEnd);

            // Load the contig prefix/suffix and split align.
            highCov = loadContigAndSplitAlign(insPos, bamStream, bai, ref, contigIt->second, (*it).loc, options);

            (*it).loc.chrStart += options.maxInsertSize;
        }

        writeLocPos(outStream, (*it).loc, insPos, highCov);

        while (progress * fiftieth < i)
        {
            std::cerr << "*" << std::flush;
            ++progress;
        }

        ++i;
        ++it;
    }
    std::cerr << std::endl;

    return 0;
}

#endif  // POPINS_PLACE_SPLIT_ALIGN_H_
