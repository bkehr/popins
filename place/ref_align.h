
#ifndef POPINS_PLACE_REF_ALIGN_H_
#define POPINS_PLACE_REF_ALIGN_H_

#include <algorithm>
#include <seqan/align.h>
#include "location.h"
#include "location_info.h"

using namespace seqan;

// =======================================================================================

// ---------------------------------------------------------------------------------------
// Struct SampleLists
// ---------------------------------------------------------------------------------------

struct SampleLists
{
    std::vector<CharString> pns;
    std::vector<std::vector<int> > lists;
};

// ---------------------------------------------------------------------------------------
// Function initSplitAlignLists()
// ---------------------------------------------------------------------------------------

void
initSplitAlignLists(SampleLists & splitAlignLists, String<LocationInfo> & locations)
{
    // Get a set of all samples IDs.

    std::set<CharString> samples;

    Iterator<String<LocationInfo> >::Type locIt = begin(locations);
    Iterator<String<LocationInfo> >::Type locEnd = end(locations);

    while (locIt != locEnd)
    {
        std::map<CharString, unsigned>::iterator sampleIt = (*locIt).loc.bestSamples.begin();
        std::map<CharString, unsigned>::iterator sampleEnd = (*locIt).loc.bestSamples.end();

        while (sampleIt != sampleEnd)
        {
            samples.insert(sampleIt->first);
            ++sampleIt;
        }

        ++locIt;
    }

    // Resize the splitAlignLists members.

    splitAlignLists.pns.resize(samples.size());
    splitAlignLists.lists.resize(samples.size());

    // Copy the samples IDs in set into splitAlignLists.

    std::set<CharString>::iterator sIt = samples.begin();
    std::set<CharString>::iterator sEnd = samples.end();
    unsigned i = 0;

    while (sIt != sEnd)
    {
        splitAlignLists.pns[i] = *sIt;
        ++sIt;
        ++i;
    }
}

// ---------------------------------------------------------------------------------------
// Function setOtherEndBit()
// ---------------------------------------------------------------------------------------

void setOtherEndBit(String<LocationInfo> & locations, bool checkPos = false)
{
    typedef Iterator<String<LocationInfo> >::Type TIter;

    TIter it = begin(locations);
    TIter itSet = begin(locations);
    TIter itEnd = end(locations);

    bool fwd = false;
    bool rev = false;

    while (it != itEnd)
    {
        if ((*it).loc.contig != (*itSet).loc.contig)
        {
            while (itSet != it)
            {
                (*itSet).otherEnd = ((*itSet).loc.contigOri && rev) || (!(*itSet).loc.contigOri && fwd);
                ++itSet;
            }
            fwd = false;
            rev = false;
        }

        if (!checkPos || (*it).insPos != -1)
        {
            if ((*it).loc.contigOri)
                fwd = true;
            else
                rev = true;
        }

        ++it;
    }

    while (itSet != it)
    {
        (*itSet).otherEnd = ((*itSet).loc.contigOri && rev) || (!(*itSet).loc.contigOri && fwd);
        ++itSet;
    }
}

// ---------------------------------------------------------------------------------------
// Function setContigLengths()
// ---------------------------------------------------------------------------------------

bool
setContigLengths(String<LocationInfo> & locations,
        std::vector<std::pair<CharString, Dna5String> > & contigs)
{
    // Assume contigs to be sorted.
    // Assume locations to be sorted by contig.
    // Assume all loc.contig to be present in contigs.

    typedef Iterator<String<LocationInfo> >::Type TLocIter;
    typedef std::vector<std::pair<CharString, Dna5String> >::iterator TContigIter;

    TLocIter locIt = begin(locations);
    TLocIter locEnd = end(locations);

    TContigIter contigIt = contigs.begin();

    while (locIt != locEnd)
    {
        while (contigIt != contigs.end() && contigIt->first < (*locIt).loc.contig)
            ++contigIt;

        if (contigIt == contigs.end() || contigIt->first != (*locIt).loc.contig)
        {
            std::cerr << "ERROR: Location references a contig that was not found in contig file: " << (*locIt).loc.contig << std::endl;
            return 1;
        }
        (*locIt).contigLength = length(contigIt->second);

        ++locIt;
    }

    return 0;
}

Location
otherEnd(Location & loc, unsigned readLength, unsigned maxInsertSize)
{
    Location rc = loc;

    rc.contigOri = !loc.contigOri;
    rc.chrOri = !loc.chrOri;
    rc.numReads = 0;

    if (loc.chrOri)
    {
        if (loc.chrEnd < readLength)
            loc.chrStart = 0;
        else
            rc.chrStart = loc.chrEnd - readLength;
        rc.chrEnd = loc.chrEnd + maxInsertSize;
    }
    else
    {
        if (loc.chrStart < maxInsertSize)
            rc.chrStart = 0;
        else
            rc.chrStart = loc.chrStart - maxInsertSize;
        rc.chrEnd = loc.chrStart + readLength;
    }

    return rc;
}

// ---------------------------------------------------------------------------------------
// Function appendGroups()
// ---------------------------------------------------------------------------------------

void
appendGroups(String<String<unsigned> > & groups, String<String<LocationInfo> > & locGroups)
{
    for (unsigned i = 0; i < length(locGroups); ++i)
    {
        String<unsigned> group;
        resize(group, length(locGroups[i]));

        for (unsigned j = 0; j < length(locGroups[i]); ++j)
            group[j] = locGroups[i][j].idx - 1;

        appendValue(groups, group);
    }
}

// ---------------------------------------------------------------------------------------
// Function findExcludeLocs()
// ---------------------------------------------------------------------------------------

void
findExcludeLocs(std::vector<Pair<CharString, bool> > & exclude, String<LocationInfo> & locations, String<String<unsigned> > & groups)
{
    typedef Iterator<String<String<unsigned> > >::Type TIter;
    typedef Iterator<String<unsigned> >::Type TGroupIter;

    // Sort locations by contig/contigOri.
    std::stable_sort(begin(locations), end(locations), LocationInfoTypeLess());

    // Set otherEndBit to whether other end was reference aligned or not.
    setOtherEndBit(locations, true);

    // Sort locations by genomic position.
    std::stable_sort(begin(locations), end(locations), LocationInfoPosLess());

    // Iterate groups to find all contig ends that can be excluded from split alignment.
    TIter groupsIt = begin(groups);
    TIter groupsEnd = end(groups);

    while (groupsIt != groupsEnd)
    {
        TGroupIter it = begin(*groupsIt);
        TGroupIter itEnd = end(*groupsIt);

        bool placed = false;
        while (it != itEnd)
        {
            placed |= locations[*it].otherEnd;
            ++it;
        }

        if (placed)
        {
            it = begin(*groupsIt);
            while (it != itEnd)
            {
                Location loc = locations[*it].loc;
                exclude.push_back(Pair<CharString, bool>(loc.contig, !loc.contigOri));
                ++it;
            }
        }

        ++groupsIt;
    }

    std::stable_sort(exclude.begin(), exclude.end());
}

// =======================================================================================

// ---------------------------------------------------------------------------------------
// Function writeGroup()
// ---------------------------------------------------------------------------------------

template<typename TStream>
void
writeGroup(TStream & outStream, String<LocationInfo> & group, bool isInsertion)
{
    // Determine the genomic interval of the group.

    CharString chrom = group[0].loc.chr;
    unsigned beginPos = group[0].loc.chrStart;
    unsigned endPos = group[0].loc.chrEnd;

    for (unsigned i = 1; i < length(group); ++i)
    {
        if (group[i].loc.chrStart < beginPos)
            beginPos = group[i].loc.chrStart;
        if (group[i].loc.chrEnd > endPos)
            endPos = group[i].loc.chrEnd;
    }

    // Write the interval and names of all contigs in the group.

    outStream << chrom;
    outStream << "\t" << beginPos;
    outStream << "\t" << endPos;
    outStream << "\t" << (isInsertion?"SV":"REF");
    outStream << "\t" << (group[0].loc.chrOri?"LEFT":"RIGHT");
    outStream << "\t" << group[0].loc.contig << ":" << (group[0].loc.contigOri == group[0].loc.chrOri?"RC":"FW");
    if (group[0].insPos > 0)
        outStream << ":" << group[0].insPos;

    for (unsigned i = 1; i < length(group); ++i)
    {
        outStream << "," << group[i].loc.contig << ":" << (group[i].loc.contigOri == group[i].loc.chrOri?"RC":"FW");
        if (group[i].insPos > 0)
            outStream << ":" << group[i].insPos;
    }

    outStream << std::endl;
}

// ---------------------------------------------------------------------------------------
// Function writeSplitAlignList()
// ---------------------------------------------------------------------------------------

bool
writeSplitAlignList(CharString & filename,
        std::vector<int> & list,
        std::vector<Pair<CharString, bool> > & exclude,
        String<LocationInfo> & locations,
        PlacingOptions<RefAlign> & options)
{
    typedef std::vector<int>::iterator TIter;

    std::fstream outStream(toCString(filename), std::ios::out);
    if (!outStream.good())
    {
        std::cerr << "ERROR: Could not open locations file " << filename << " for writing." << std::endl;
        return 1;
    }

    TIter it = list.begin();
    TIter itEnd = list.end();

    while (it != itEnd)
    {
        Location loc;
        if (*it > 0)
            loc = locations[(*it) - 1].loc;
        else
            loc = otherEnd(locations[-(*it) - 1].loc, options.readLength, options.maxInsertSize);

        Pair<CharString, bool> c(loc.contig, loc.contigOri);
        if (exclude.size() == 0)
        {
            writeLoc(outStream, loc);
        }
        else
        {
            std::vector<Pair<CharString, bool> >::iterator cIt = lower_bound(exclude.begin(), exclude.end(), c);

            if (cIt != exclude.end() && ((*cIt).i1 != c.i1 || (*cIt).i2 != c.i2))
                writeLoc(outStream, loc);
        }

        ++it;
    }

    return 0;
}

// =======================================================================================

// ---------------------------------------------------------------------------------------
// Function isBetterXXX()
// ---------------------------------------------------------------------------------------

inline bool
isBetterRefAligned(LocationInfo & loc, String<LocationInfo> & group)
{
    if (loc.insPos <= 0)
        return false;

    if (group[0].otherEnd && !loc.otherEnd)
        return false;

    if (loc.contigLength <= group[0].contigLength)
        return false;

    return true;
}

inline bool
isBetterUnaligned(LocationInfo & loc, String<LocationInfo> & group)
{
    if (group[0].otherEnd && !loc.otherEnd)
        return false;

    if (loc.loc.numReads > group[0].loc.numReads)
        return true;

    if (loc.loc.numReads < group[0].loc.numReads)
        return false;

    if (loc.contigLength <= group[0].contigLength)
        return false;

    return true;
}

// =======================================================================================

// ---------------------------------------------------------------------------------------
// Function setInsPos()
// ---------------------------------------------------------------------------------------

template<typename TSeqA, typename TSeqB>
void
setInsPos(LocationInfo & a, LocationInfo & b,
        Dna5String & contigA, Dna5String & contigB,
        Gaps<TSeqA> & gapsA, Gaps<TSeqB> & gapsB)
{
    a.refPos = b.refPos;

    if (b.insPos == -1)
    {
        a.insPos = -1;
        return;
    }

    int bPos = b.insPos;
    if (!a.loc.chrOri)
        bPos = length(contigB) - b.insPos;

    if (bPos < (int)beginPosition(gapsB) || bPos >= (int)endPosition(gapsB))
    {
        a.insPos = -1;
        return;
    }

    if (a.loc.chrOri)
        a.insPos = toSourcePosition(gapsA, toViewPosition(gapsB, b.insPos));
    else
        a.insPos = length(contigA) - toSourcePosition(gapsA, toViewPosition(gapsB, bPos));
}

// ---------------------------------------------------------------------------------------
// Function align()
// ---------------------------------------------------------------------------------------

template<typename TSeqA, typename TSeqB>
bool
align(Gaps<TSeqA> & gapsA, Gaps<TSeqB> & gapsB, TSeqA & a, TSeqB & b)
{
    setSource(gapsA, a);
    setSource(gapsB, b);

    Score<int> scoring(1, -3, -4, -5);

    int score = localAlignment(gapsA, gapsB, scoring);

    //    std::cout << gapsA << std::endl << gapsB << std::endl;

    if (score > 25)        // TODO: Make this 25 a program parameter.
        return true;

    return false;
}

// ---------------------------------------------------------------------------------------
// Function alignContigEnds()
// ---------------------------------------------------------------------------------------

bool
contigEndsAlign(LocationInfo & a, LocationInfo & b, std::vector<std::pair<CharString, Dna5String> > & contigs)
{
    unsigned preSufLen = 200;    // TODO: Make this a program parameter.

    typedef ModifiedString<Suffix<Dna5String>::Type, ModComplementDna5> TComplementSuffix;
    typedef ModifiedString<TComplementSuffix, ModReverse> TRCSuffix;
    typedef Prefix<Dna5String>::Type TPrefix;

    typedef std::pair<CharString, Dna5String> TPair;

    std::vector<TPair>::iterator itA = std::lower_bound(contigs.begin(), contigs.end(), TPair(a.loc.contig, ""));
    std::vector<TPair>::iterator itB = std::lower_bound(contigs.begin(), contigs.end(), TPair(b.loc.contig, ""));

    if (b.insPos != -1)
    {
        if (b.loc.chrOri)
            preSufLen = _max((int)preSufLen, b.insPos + 50);
        else
            preSufLen = _max(preSufLen, length(itB->second) - b.insPos + 50);
    }

    SEQAN_ASSERT_EQ(a.loc.chrOri, b.loc.chrOri);

    //    std::cout << itA->first << "  " << length(itA->second) << std::endl;
    //    std::cout << itB->first << "  " << length(itB->second) << std::endl;

    if (!a.loc.contigOri)
    {
        unsigned prefixLengthA = _min(preSufLen, length(itA->second));
        TPrefix prefixA = prefix(itA->second, prefixLengthA);
        Gaps<TPrefix> gapsA;

        if (!b.loc.contigOri)
        {
            // case A: left insertion end and contigs of both a and b in fwd orientation OR
            //         right insertion end and contigs of both a and b in rev orientation
            //         => align contig prefixes of a and b to each other

            unsigned prefixLengthB = _min(preSufLen, length(itB->second));
            TPrefix prefixB = prefix(itB->second, prefixLengthB);

            Gaps<TPrefix> gapsB;
            if (align(gapsA, gapsB, prefixA, prefixB))
            {
                setInsPos(a, b, itA->second, itB->second, gapsA, gapsB);
                return true;
            }
        }
        else
        {
            // case B: left insertion end and contig of a in fwd and contig of b in rev orientation OR
            //         right insertion end and contig of a in rev and contig of b in fwd orientation
            //         => align contig prefix of a to contig suffix of b

            unsigned suffixBeginPosB = length(itB->second) - _min(preSufLen, length(itB->second));
            Suffix<Dna5String>::Type sufB = suffix(itB->second, suffixBeginPosB);
            TRCSuffix suffixB(sufB);

            Gaps<TRCSuffix> gapsB;
            if (align(gapsA, gapsB, prefixA, suffixB))
            {
                setInsPos(a, b, itA->second, itB->second, gapsA, gapsB);
                return true;
            }
        }
    }
    else
    {
        unsigned suffixBeginPosA = length(itA->second) - _min(preSufLen, length(itA->second));
        Suffix<Dna5String>::Type sufA = suffix(itA->second, suffixBeginPosA);
        TRCSuffix suffixA(sufA);
        Gaps<TRCSuffix> gapsA;

        if (!b.loc.contigOri)
        {
            // case C: left insertion end and contig of a in rev and contig of b in fwd orientation OR
            //         right insertion end and contig of a in fwd and contig of b in rev orientation
            //         => align contig suffix of a to contig prefix of b

            unsigned prefixLengthB = _min(preSufLen, length(itB->second));
            TPrefix prefixB = prefix(itB->second, prefixLengthB);

            Gaps<TPrefix> gapsB;
            if (align(gapsA, gapsB, suffixA, prefixB))
            {
                setInsPos(a, b, itA->second, itB->second, gapsA, gapsB);
                return true;
            }
        }
        else
        {
            // case D: left insertion end and contigs of both a and b in rev orientation OR
            //         right insertion end and contigs of both a and b in fwd orientation
            //         => align contig suffixes of a and b to each other

            unsigned suffixBeginPosB = length(itB->second) - _min(preSufLen, length(itB->second));
            Suffix<Dna5String>::Type sufB = suffix(itB->second, suffixBeginPosB);
            TRCSuffix suffixB(sufB);

            Gaps<TRCSuffix> gapsB;
            if (align(gapsA, gapsB, suffixA, suffixB))
            {
                setInsPos(a, b, itA->second, itB->second, gapsA, gapsB);
                return true;
            }
        }
    }

    return false;
}

// ---------------------------------------------------------------------------------------
// Function alignToRefAligned()
// ---------------------------------------------------------------------------------------

Iterator<String<String<LocationInfo> > >::Type
alignToRefAligned(LocationInfo & loc,
        String<String<LocationInfo> > & refAlignedGroups,
        std::vector<std::pair<CharString, Dna5String> > & contigs)
{
    typedef Iterator<String<String<LocationInfo> > >::Type TIter;

    TIter it = begin(refAlignedGroups);
    TIter itEnd = end(refAlignedGroups);

    while (it != itEnd)
    {
        if (contigEndsAlign(loc, (*it)[0], contigs))
            break;
        ++it;
    }

    return it;
}

// ---------------------------------------------------------------------------------------
// Function alignsToRef()
// ---------------------------------------------------------------------------------------

bool
alignsToRef(LocationInfo & loc,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions<RefAlign> & options)
{
    unsigned dist = 10;            // TODO: Make this a program parameter.
    unsigned preSufLen = 100;    // TODO: Make this a program parameter.

    typedef Infix<Dna5String>::Type TInfix;
    typedef ModifiedString<TInfix, ModComplementDna5> TComplementInfix;
    typedef ModifiedString<TComplementInfix, ModReverse> TRcInfix;

    typedef std::pair<CharString, Dna5String> TPair;

    std::vector<TPair>::iterator contigIt = std::lower_bound(contigs.begin(), contigs.end(), TPair(loc.loc.contig, ""));

    if (contigIt == contigs.end())
    {
        std::cerr << "ERROR: Could not find " << loc.loc.contig << " in contig file." << std::endl;
        return 1;
    }

    Dna5String ref;
    Gaps<Dna5String> refGaps;

    if (loc.loc.chrOri)
    {
        unsigned chrStart = loc.loc.chrStart - options.readLength;
        if (loc.loc.chrStart < options.readLength)
            chrStart = 0;
        if (loadInterval(ref, fai, loc.loc.chr, chrStart, loc.loc.chrEnd + options.maxInsertSize) != 0)
            return 1;

        if (loc.loc.contigOri)
        {
            unsigned suffixEndPos = length(contigIt->second);
            unsigned suffixBeginPos = suffixEndPos - _min(preSufLen, suffixEndPos);
            TInfix suf = infix(contigIt->second, suffixBeginPos, suffixEndPos);
            TRcInfix contigSuffix(suf);

            Gaps<TRcInfix> contigGaps;
            while (align(contigGaps, refGaps, contigSuffix, ref))
            {
                loc.insPos = length(contigIt->second) - suffixEndPos + endPosition(contigGaps);
                loc.refPos = chrStart + endPosition(refGaps) - 1;

                if (suffixEndPos - suffixBeginPos - dist > endPosition(contigGaps))
                    return true;

                if (suffixEndPos - suffixBeginPos < preSufLen) // all of the contig aligns to the reference
                        {
                    loc.insPos = -1;
                    return true;
                        }

                suffixEndPos -= preSufLen/2;
                suffixBeginPos = suffixEndPos - _min(preSufLen, suffixEndPos);
                suf = infix(contigIt->second, suffixBeginPos, suffixEndPos);
                contigSuffix = TRcInfix(suf);

                clear(contigGaps);
                clear(refGaps);
            }
        }
        else
        {
            unsigned prefixBeginPos = 0;
            unsigned prefixEndPos = _min(preSufLen, length(contigIt->second));
            TInfix contigPrefix = infix(contigIt->second, prefixBeginPos, prefixEndPos);

            Gaps<TInfix> contigGaps;
            while (align(contigGaps, refGaps, contigPrefix, ref))
            {
                loc.insPos = prefixBeginPos + endPosition(contigGaps);
                loc.refPos = chrStart + endPosition(refGaps) - 1;

                if (prefixEndPos - prefixBeginPos - dist > endPosition(contigGaps))
                    return true;

                if (prefixEndPos - prefixBeginPos < preSufLen) // all of the contig aligns to the reference
                        {
                    loc.insPos = -1;
                    return true;
                        }

                prefixBeginPos += preSufLen/2;
                prefixEndPos = _min(prefixBeginPos + preSufLen, length(contigIt->second));
                contigPrefix = infix(contigIt->second, prefixBeginPos, prefixEndPos);

                clear(contigGaps);
                clear(refGaps);
            }
        }
    }
    else
    {
        unsigned chrStart = loc.loc.chrStart - options.maxInsertSize;
        if (loc.loc.chrStart < options.maxInsertSize)
             chrStart = 0;
        if (loadInterval(ref, fai, loc.loc.chr, chrStart, loc.loc.chrEnd + options.readLength) != 0)
            return 1;

        if (loc.loc.contigOri)
        {
            unsigned suffixEndPos = length(contigIt->second);
            unsigned suffixBeginPos = suffixEndPos - _min(preSufLen, suffixEndPos);
            TInfix contigSuffix = infix(contigIt->second, suffixBeginPos, suffixEndPos);

            Gaps<TInfix> contigGaps;
            while (align(contigGaps, refGaps, contigSuffix, ref))
            {
                loc.insPos = suffixBeginPos + beginPosition(contigGaps);
                loc.refPos = chrStart + beginPosition(refGaps);

                if (beginPosition(contigGaps) > dist)
                    return true;

                if (suffixEndPos - suffixBeginPos < preSufLen)
                {
                    loc.insPos = -1;
                    return true;
                }

                suffixEndPos -= preSufLen/2;
                suffixBeginPos = suffixEndPos - _min(preSufLen, suffixEndPos);
                contigSuffix = infix(contigIt->second, suffixBeginPos, suffixEndPos);

                clear(contigGaps);
                clear(refGaps);
            }
        }
        else
        {
            unsigned prefixBeginPos = 0;
            unsigned prefixEndPos = _min(preSufLen, length(contigIt->second));
            TInfix pref = infix(contigIt->second, prefixBeginPos, prefixEndPos);
            TRcInfix contigPrefix(pref);

            Gaps<TRcInfix> contigGaps;
            while (align(contigGaps, refGaps, contigPrefix, ref))
            {
                loc.insPos = length(contigIt->second) - prefixEndPos + beginPosition(contigGaps);
                loc.refPos = chrStart + beginPosition(refGaps);

                if (beginPosition(contigGaps) > dist)
                    return true;

                if (prefixEndPos - prefixBeginPos < preSufLen)
                {
                    loc.insPos = -1;
                    return true;
                }

                prefixBeginPos += preSufLen/2;
                prefixEndPos = _min(prefixBeginPos + preSufLen, length(contigIt->second));
                pref = infix(contigIt->second, prefixBeginPos, prefixEndPos);
                contigPrefix = TRcInfix(pref);

                clear(contigGaps);
                clear(refGaps);
            }
        }
    }

    return false;
}

// ---------------------------------------------------------------------------------------
// Function addToLists()
// ---------------------------------------------------------------------------------------

void
addToLists(SampleLists  & splitAlignLists,
        LocationInfo & loc)
{
    std::map<CharString, unsigned>::iterator it = loc.loc.bestSamples.begin();
    std::map<CharString, unsigned>::iterator itEnd = loc.loc.bestSamples.end();

    while (it != itEnd)
    {
        std::vector<CharString>::iterator pnIt = std::find(splitAlignLists.pns.begin(), splitAlignLists.pns.end(), it->first);
        std::vector<std::vector<int> >::iterator listsIt = splitAlignLists.lists.begin() + (pnIt - splitAlignLists.pns.begin());

        (*listsIt).push_back(loc.idx);
        ++it;
    }
}

// ---------------------------------------------------------------------------------------
// Function processOtherEnd()
// ---------------------------------------------------------------------------------------

template<typename TStream>
bool
processOtherEnd(TStream & vcfStream,
        SampleLists & splitAlignLists,
        LocationInfo & loc,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions<RefAlign> & options)
{
    if (loc.otherEnd == false && loc.insPos != -1)
    {
        // Reverse complement the location.
        LocationInfo rc(loc);
        rc.loc = otherEnd(loc.loc, options.readLength, options.maxInsertSize);
        rc.insPos = 0;
        rc.idx = -loc.idx;

        // Process the reverse complemented location.
        if (alignsToRef(rc, contigs, fai, options))
        {
            if (loc.insPos != -1 && loc.loc.chrOri && loc.insPos >= rc.insPos)
                return 0;
            if (loc.insPos != -1 && !loc.loc.chrOri && loc.insPos < rc.insPos)
                return 0;
            if (rc.insPos == -1)
                return 0;

            if (writeVcf(vcfStream, rc, 0, fai) != 0)
                return 1;
        }
        else
            addToLists(splitAlignLists, rc);
    }
    return 0;
}

// =======================================================================================

// ---------------------------------------------------------------------------------------
// Function findRefAlignedGroups()
// ---------------------------------------------------------------------------------------

void
findRefAlignedGroups(String<String<LocationInfo> > & refAlignedGroups,
        String<LocationInfo> & unaligned,
        String<LocationInfo> & locations,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions<RefAlign> & options)
{
    Iterator<String<LocationInfo> >::Type it = begin(locations);
    Iterator<String<LocationInfo> >::Type itEnd = end(locations);

    while (it != itEnd)
    {
        Iterator<String<String<LocationInfo> > >::Type refAlignedIt = alignToRefAligned(*it, refAlignedGroups, contigs);
        if (refAlignedIt != end(refAlignedGroups))
        {
            if (isBetterRefAligned(*it, *refAlignedIt))
            {
                appendValue(*refAlignedIt, (*refAlignedIt)[0]);
                (*refAlignedIt)[0] = LocationInfo(*it);
            }
            else
            {
                appendValue(*refAlignedIt, *it);
            }
        }
        else if (alignsToRef(*it, contigs, fai, options))
        {
            //std::cout << "Ins pos: " << (*it).insPos << std::endl;
            String<LocationInfo> newGroup;
            appendValue(newGroup, *it);
            appendValue(refAlignedGroups, newGroup);
        }
        else
        {
            appendValue(unaligned, *it);
        }

        ++it;
    }
}

// ---------------------------------------------------------------------------------------
// Function processRefAlignedGroups()
// ---------------------------------------------------------------------------------------

template<typename TStream1, typename TStream2>
bool
processRefAlignedGroups(TStream1 & vcfStream,
        TStream2 & groupStream,
        String<String<LocationInfo> > & groups,
        SampleLists & splitAlignLists,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions<RefAlign> & options)
{
    typename Iterator<String<String<LocationInfo> > >::Type it = begin(groups);
    typename Iterator<String<String<LocationInfo> > >::Type itEnd = end(groups);

    while (it != itEnd)
    {
        if ((*it)[0].insPos != -1)
        {
            if (writeVcf(vcfStream, (*it)[0], length(*it), fai) != 0)
                return 1;
            writeGroup(groupStream, *it, true);
        }
        else
        {
            writeGroup(groupStream, *it, false);
        }

        if (processOtherEnd(vcfStream, splitAlignLists, (*it)[0], contigs, fai, options) != 0)
            return 1;

        ++it;
    }
    return 0;
}

// ---------------------------------------------------------------------------------------
// Function findUnalignedGroups()
// ---------------------------------------------------------------------------------------

void
findUnalignedGroups(String<String<LocationInfo> > & groups,
        String<LocationInfo> & unaligned,
        std::vector<std::pair<CharString, Dna5String> > & contigs)
{
    Iterator<String<LocationInfo> >::Type it = begin(unaligned);
    Iterator<String<LocationInfo> >::Type itEnd = end(unaligned);

    unsigned largestGroupSize = 0;

    while (it != itEnd)
    {
        bool added = false;
        for (unsigned s = 0; s < largestGroupSize; ++s)
        {
            for (unsigned g = 0; g < length(groups); ++g)
            {
                if (length(groups[g]) < s + 1)
                    continue;

                if (contigEndsAlign(*it, groups[g][s], contigs))
                {
                    if (isBetterUnaligned(*it, groups[g]))
                    {
                        appendValue(groups[g], groups[g][0]);
                        groups[g][0] = *it;
                    }
                    else
                    {
                        appendValue(groups[g], *it);
                    }
                    largestGroupSize = _max(largestGroupSize, length(groups[g]));
                    added = true;
                }
                if (added)
                    break;
            }
            if (added)
                break;
        }

        if (!added)
        {
            String<LocationInfo> newGroup;
            appendValue(newGroup, *it);
            appendValue(groups, newGroup);
            if (largestGroupSize == 0)
                largestGroupSize = 1;
        }

        ++it;
    }
}
// ---------------------------------------------------------------------------------------
// Function processUnalignedGroups()
// ---------------------------------------------------------------------------------------

template<typename TStream1, typename TStream2>
bool
processUnalignedGroups(TStream1 & vcfStream,
        TStream2 & groupStream,
        String<String<LocationInfo> > & groups,
        SampleLists & splitAlignLists,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions<RefAlign> & options)
{
    typename Iterator<String<String<LocationInfo> > >::Type it = begin(groups);
    typename Iterator<String<String<LocationInfo> > >::Type itEnd = end(groups);

    while (it != itEnd)
    {
        addToLists(splitAlignLists, (*it)[0]);
        writeGroup(groupStream, *it, true);
        if (processOtherEnd(vcfStream, splitAlignLists, (*it)[0], contigs, fai, options) != 0)
            return 1;
        ++it;
    }
    return 0;
}

// ---------------------------------------------------------------------------------------
// Function processOverlappingLocations()
// ---------------------------------------------------------------------------------------

template<typename TStream1, typename TStream2>
bool
processOverlappingLocs(TStream1 & vcfStream,
        TStream2 & groupStream,
        String<String<unsigned> > & groups,
        SampleLists & splitAlignLists,
        String<LocationInfo> & locations,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions<RefAlign> & options)
{
    String<String<LocationInfo> > refAlignedGroups;
    String<String<LocationInfo> > unalignedGroups;
    String<LocationInfo> unaligned;

    //    std::cout << "\nProcessing group of " << length(locations) << " overlapping locations." << std::endl;
    //    for (unsigned i = 0; i < length(locations); ++i)
    //        std::cout << locations[i].loc.chr << "\t" << locations[i].loc.chrStart << "\t" << (locations[i].otherEnd?"true":"false") << std::endl;

    // Sort the set of overlapping locations by bit and then by contig length.
    std::stable_sort(begin(locations), end(locations), LocationInfoGreater());

    // Split locations into ref-aligned groups or unaligned.
    findRefAlignedGroups(refAlignedGroups, unaligned, locations, contigs, fai, options);
    clear(locations);

    //    std::cout << "refAlignedGroups: " << length(refAlignedGroups) << std::endl;

    // Handle the ref-aligned groups.
    if (processRefAlignedGroups(vcfStream, groupStream, refAlignedGroups, splitAlignLists, contigs, fai, options) != 0)
        return 1;
    appendGroups(groups, refAlignedGroups);
    clear(refAlignedGroups);

    //    std::cout << "   unaligned: " << length(unaligned) << std::endl;

    // Split the unaligned locations into groups.
    findUnalignedGroups(unalignedGroups, unaligned, contigs);
    clear(unaligned);

    // Handle the unaligned groups.
    if (processUnalignedGroups(vcfStream, groupStream, unalignedGroups, splitAlignLists, contigs, fai, options) != 0)
        return 1;
    appendGroups(groups, unalignedGroups);

    //    std::cout << "   unlignedGroups: " << length(unalignedGroups) << std::endl;

    return 0;
}

// =======================================================================================
// Function popins_place_ref_align()
// =======================================================================================

template<typename TStream>
bool
popins_place_ref_align(TStream & vcfStream,
        String<LocationInfo> & locations,
        std::vector<std::pair<CharString, Dna5String> > & contigs,
        FaiIndex & fai,
        PlacingOptions<RefAlign> & options)
{
    printStatus("Aligning contigs to reference");

    // Initialize splitAlignLists with all sample IDs (using locations).
    SampleLists splitAlignLists;
    initSplitAlignLists(splitAlignLists, locations);

    // Open the groups output file.
    String<String<unsigned> > groups;
    std::fstream outGroups(toCString(options.groupsFile), std::ios_base::out);
    if (!outGroups.is_open())
    {
        std::cerr << "ERROR: Could not open groups output file " << options.groupsFile << std::endl;
        return 1;
    }

    // Sort locations by contig/contigOri.
    std::stable_sort(begin(locations), end(locations), LocationInfoTypeLess());

    // Add bit to each location indicating whether there is a location for the other end or not.
    setOtherEndBit(locations);
    if (setContigLengths(locations, contigs) != 0)
        return 1;

    // Sort locations by genomic position.
    std::stable_sort(begin(locations), end(locations), LocationInfoPosLess());

    // Reset location indices to match sort order.
    for (unsigned i = 0; i < length(locations); ++i)
        locations[i].idx = i + 1;

    // --- Iterate over locations in sets of overlapping genomic positions.

    std::cerr << "0%   10   20   30   40   50   60   70   80   90   100%" << std::endl;
    std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;
    std::cerr << "*" << std::flush;

    double fiftieth = length(locations) / 50.0;
    unsigned progress = 0;

    Iterator<String<LocationInfo> >::Type it = begin(locations);
    Iterator<String<LocationInfo> >::Type itEnd = end(locations);

    String<LocationInfo> fwd;
    String<LocationInfo> rev;

    CharString prevChromFwd = "";
    unsigned prevPosFwd = 0;
    CharString prevChromRev = "";
    unsigned prevPosRev = 0;

    unsigned i = 0;
    while (it != itEnd)
    {
        if ((*it).loc.chrOri)
        {
            if (length(fwd) != 0 && (prevChromFwd != (*it).loc.chr || prevPosFwd + options.groupDist < (*it).loc.chrStart))
            {
                if (processOverlappingLocs(vcfStream, outGroups, groups, splitAlignLists, fwd, contigs, fai, options) != 0)
                    return 1;
                clear(fwd);
            }
            appendValue(fwd, *it);
            prevChromFwd = (*it).loc.chr;
            prevPosFwd = (*it).loc.chrEnd;
        }
        else
        {
            if (length(rev) != 0 && (prevChromRev != (*it).loc.chr || prevPosRev + options.groupDist < (*it).loc.chrStart))
            {
                if (processOverlappingLocs(vcfStream, outGroups, groups, splitAlignLists, rev, contigs, fai, options) != 0)
                    return 1;
                clear(rev);
            }
            appendValue(rev, *it);
            prevChromRev = (*it).loc.chr;
            prevPosRev = (*it).loc.chrEnd;
        }

        while (progress * fiftieth < i)
        {
            std::cerr << "*" << std::flush;
            ++progress;
        }
        ++i;
        ++it;
    }

    if (length(fwd) != 0)
        if (processOverlappingLocs(vcfStream, outGroups, groups, splitAlignLists, fwd, contigs, fai, options) != 0)
            return 1;

    if (length(rev) != 0)
        if (processOverlappingLocs(vcfStream, outGroups, groups, splitAlignLists, rev, contigs, fai, options) != 0)
            return 1;

    while (progress < 50)
    {
    	std::cerr << "*" << std::flush;
    	++progress;
    }
    std::cerr << std::endl;
    printStatus("Writing locations of contigs that do not align to the reference (per sample).");

    // Find locations to exclude from splitAlignLists
    std::vector<Pair<CharString, bool> > exclude;
    findExcludeLocs(exclude, locations, groups);

    // Write splitAlignLists to output files.
    for (unsigned i = 0; i < splitAlignLists.pns.size(); ++i)
    {
        CharString filename = getFileName(options.prefix, splitAlignLists.pns[i]);
        filename += "/locations_unplaced.txt";
        if (writeSplitAlignList(filename, splitAlignLists.lists[i], exclude, locations, options) != 0)
            return 1;
    }

    return 0;
}

#endif /* POPINS_PLACE_REF_ALIGN_H_ */
