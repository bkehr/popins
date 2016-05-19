#ifndef POPINS_LOCATION_H_
#define POPINS_LOCATION_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

using namespace seqan;

// ==========================================================================
// struct GenomicInterval
// ==========================================================================

struct GenomicInterval
{
    typedef Position<CharString>::Type TPos;
    
    CharString chr;
    TPos begin;
    TPos end;
    bool ori;
    
    GenomicInterval(CharString & c, TPos b, TPos e, bool o) :
        chr(c), begin(b), end(e), ori(o)
    {}
    
    GenomicInterval(GenomicInterval const & other) :
        chr(other.chr), begin(other.begin), end(other.end), ori(other.ori)
    {}
};

// ==========================================================================
// struct AnchoringRecord
// ==========================================================================

struct AnchoringRecord
{
    typedef Position<CharString>::Type TPos;

    // TODO replace by GenomicInterval
    CharString chr;
    TPos chrStart;
    TPos chrEnd;
    bool chrOri;

    CharString contig;
    bool contigOri;
};

// ==========================================================================

struct AnchoringRecordLess : public std::binary_function<AnchoringRecord, AnchoringRecord, bool>
{
	AnchoringRecordLess() {}

    inline int compare(AnchoringRecord const & a, AnchoringRecord const & b) const
    {
        if (a.contig > b.contig) return -1;
        if (a.contig < b.contig) return 1;

        if (a.contigOri && !b.contigOri) return -1;
        if (!a.contigOri && b.contigOri) return 1;

        if (a.chr > b.chr) return -1;
        if (a.chr < b.chr) return 1;

        if (a.chrOri && !b.chrOri) return -1;
        if (!a.chrOri && b.chrOri) return 1;

        if (a.chrStart > b.chrStart) return -1;
        if (a.chrStart < b.chrStart) return 1;

        if (a.chrEnd > b.chrEnd) return -1;
        if (a.chrEnd < b.chrEnd) return 1;

        return 0;
    }

    inline bool operator() (AnchoringRecord const & a, AnchoringRecord const & b) const
    {
        return compare(a, b) == 1;
    }
};



// ==========================================================================
// struct Location
// ==========================================================================

struct Location
{
    typedef Position<CharString>::Type TPos;

    // TODO replace by GenomicInterval
    CharString chr;
    TPos chrStart;
    TPos chrEnd;
    bool chrOri;

    CharString contig;
    bool contigOri;

    unsigned numReads;
    double score;
    
    std::map<CharString, unsigned> bestSamples;
    unsigned fileIndex;
    
    Location ()
    {
        chr = "";
        chrStart = 0;
        chrEnd = 0;
        contig = "";
        score = -1;
    }
    
    Location (CharString h, TPos hs, TPos he, bool ho, CharString c, bool co, unsigned n, double s) :
        chr(h), chrStart(hs), chrEnd(he), chrOri(ho), contig(c), contigOri(co), numReads(n), score(s)
    {}
    
    Location (AnchoringRecord const & r) :
        chr(r.chr), chrStart(r.chrStart), chrEnd(r.chrEnd), chrOri(r.chrOri),
        contig(r.contig), contigOri(r.contigOri),
        numReads(1)
    {}
};

// ==========================================================================

struct LocationPosLess : public std::binary_function<Location, Location, bool>
{
    LocationPosLess() {}
    
    inline int compare(Location const & a, Location const & b) const
    {
        bool chrADigit = std::isdigit(a.chr[0]);
        bool chrBDigit = std::isdigit(b.chr[0]);
        if (chrADigit && chrBDigit)
        {
            int chrA, chrB;
            lexicalCast2<int>(chrA, a.chr);
            lexicalCast2<int>(chrB, b.chr);
            if (chrA > chrB) return -1;
            if (chrA < chrB) return 1;
        }
        else if (!chrADigit && !chrBDigit)
        {
            if (a.chr > b.chr) return -1;
            if (a.chr < b.chr) return 1;
        }
        else if (chrADigit && !chrBDigit) return 1;
        else if (!chrADigit && chrBDigit) return -1;
        
        if (a.chrStart > b.chrStart) return -1;
        if (a.chrStart < b.chrStart) return 1;
        
        if (a.contig > b.contig) return -1;
        if (a.contig < b.contig) return 1;
        
        if (a.chrOri && !b.chrOri) return -1;
        if (!a.chrOri && b.chrOri) return 1;
        
        if (a.contigOri && !b.contigOri) return -1;
        if (!a.contigOri && b.contigOri) return 1;
        
        if (a.chrEnd > b.chrEnd) return -1;
        if (a.chrEnd < b.chrEnd) return 1;
        
        return 0;
    }
    
    inline bool operator() (Location const & a, Location const & b) const
    {
        return compare(a, b) == 1;
    }
};

// ==========================================================================

struct LocationTypeLess : public std::binary_function<Location, Location, bool> 
{
    LocationTypeLess() {}
    
    inline int compare(Location const & a, Location const & b) const
    {
        if (a.contig > b.contig) return -1;
        if (a.contig < b.contig) return 1;
        
        if (a.contigOri && !b.contigOri) return -1;
        if (!a.contigOri && b.contigOri) return 1;
        
        bool chrADigit = std::isdigit(a.chr[0]);
        bool chrBDigit = std::isdigit(b.chr[0]);
        if (chrADigit && chrBDigit)
        {
            int chrA, chrB;
            lexicalCast2<int>(chrA, a.chr);
            lexicalCast2<int>(chrB, b.chr);
            if (chrA > chrB) return -1;
            if (chrA < chrB) return 1;
        }
        else if (!chrADigit && !chrBDigit)
        {
            if (a.chr > b.chr) return -1;
            if (a.chr < b.chr) return 1;
        }
        else if (chrADigit && !chrBDigit) return 1;
        else if (!chrADigit && chrBDigit) return -1;
        
        if (a.chrStart > b.chrStart) return -1;
        if (a.chrStart < b.chrStart) return 1;
        
        if (a.chrOri && !b.chrOri) return -1;
        if (!a.chrOri && b.chrOri) return 1;
        
        if (a.chrEnd > b.chrEnd) return -1;
        if (a.chrEnd < b.chrEnd) return 1;
        
        return 0;
    }
    
    inline bool operator() (Location const & a, Location const & b) const
    {
        return compare(a, b) == 1;
    }
};

// ==========================================================================

struct LocationTypeGreater : public std::binary_function <Location, Location, bool> 
{
    LocationTypeGreater() {}
    
    inline int compare(Location const & a, Location const & b) const
    {
        if (a.contig > b.contig) return -1;
        if (a.contig < b.contig) return 1;
        
        if (a.contigOri && !b.contigOri) return -1;
        if (!a.contigOri && b.contigOri) return 1;
        
        bool chrADigit = std::isdigit(a.chr[0]);
        bool chrBDigit = std::isdigit(b.chr[0]);
        if (chrADigit && chrBDigit)
        {
            int chrA, chrB;
            lexicalCast2<int>(chrA, a.chr);
            lexicalCast2<int>(chrB, b.chr);
            if (chrA > chrB) return -1;
            if (chrA < chrB) return 1;
        }
        else if (!chrADigit && !chrBDigit)
        {
            if (a.chr > b.chr) return -1;
            if (a.chr < b.chr) return 1;
        }
        else if (chrADigit && !chrBDigit) return 1;
        else if (!chrADigit && chrBDigit) return -1;
        
        if (a.chrStart > b.chrStart) return -1;
        if (a.chrStart < b.chrStart) return 1;
        
        if (a.chrOri && !b.chrOri) return -1;
        if (!a.chrOri && b.chrOri) return 1;
        
        if (a.chrEnd > b.chrEnd) return -1;
        if (a.chrEnd < b.chrEnd) return 1;
        
        return 0;
    }
    
    inline bool operator() (Location const & a, Location const & b) const
    {
        return compare(a, b) == -1;
    }
};

// ==========================================================================
// Struct LocationsFilter
// ==========================================================================

struct LocationsFilter
{
    bool other;
    unsigned minReads;
    double minScore;
    unsigned maxLength;

    LocationsFilter() :
        other(true), minReads(1), minScore(0), maxLength(-1)
    {}

    LocationsFilter(unsigned r, double s, unsigned l) :
        other(false), minReads(r), minScore(s), maxLength(l)
    {}

};

// ==========================================================================
// Function passesFilter()
// ==========================================================================

bool
passesFilter(Location & loc, LocationsFilter & filter)
{
    if (loc.numReads < filter.minReads || loc.score < filter.minScore)
        return false;

    if (!filter.other && loc.chr == "OTHER")
        return false;

    if (loc.chrEnd - loc.chrStart > filter.maxLength)
        return false;

    return true;
}

// ==========================================================================

inline bool
isComponentOrNode(CharString & name)
{
    typedef Position<CharString>::Type TPos;

    TPos i = 0;
    while (i < length(name) && isdigit(name[i]))
    {
        ++i;
        if (i < length(name) && name[i] == '.') ++i;
    }

    if (i < length(name)-4 && (infix(name, i, i+4) == "NODE" || infix(name, i, i+4) == "COMP"))
        return true;

    return false;
}

inline bool
isChromosome(CharString & name) // TODO: Allow user to specify sequence names.
{
    typedef Position<CharString>::Type TPos;
    TPos i = 0;
    if (length(name) > 3 && prefix(name, 3) == "chr") i = 3;
    
    if ((length(name) == i+1 && (isdigit(name[i]) || name[i] == 'X' || name[i] == 'Y')) ||
        (length(name) == i+2 && isdigit(name[i]) && isdigit(name[i+1])))
        return true;

    return false;
}

// ==========================================================================

inline Pair<CigarElement<>::TCount>
mappedInterval(String<CigarElement<> > & cigar)
{
    typedef CigarElement<>::TCount TSize;
    
    TSize len = 0;
    TSize beginPos = 0;
    TSize endPos = 0;
    
    Iterator<String<CigarElement<> > >::Type itEnd = end(cigar);
    for (Iterator<String<CigarElement<> > >::Type it = begin(cigar); it < itEnd; ++it)
    {
        len += (*it).count;
        
        switch ((*it).operation)
        {
        case 'S':
            if (it == begin(cigar))
                beginPos += (*it).count;
            break;
        case 'D': case 'H':
            len -= (*it).count;
            break;
        case 'M': case 'I':
            endPos = len;
            break;
        }
    }
    
    return Pair<TSize>(beginPos, endPos);
}

// ==========================================================================

inline double
avgQuality(CharString & qual, Pair<CigarElement<>::TCount> & interval)
{
    if (interval.i1 >= interval.i2) return 0;
    if (length(qual) == 0) return 50; // Accept undefined quality strings ('*' in sam format).
    
    Iterator<CharString>::Type it = begin(qual);
    Iterator<CharString>::Type itEnd = begin(qual);
    it += interval.i1;
    itEnd += interval.i2;
    
    unsigned totalQual = 0;
    while (it != itEnd)
    {
        totalQual += *it;
        ++it;
    }
    return totalQual / (interval.i2 - interval.i1);
}

// ==========================================================================

unsigned
alignmentScore(BamAlignmentRecord & record)
{
    BamTagsDict tagsDict(record.tags);
    unsigned idx;
    if (findTagKey(idx, tagsDict, "AS"))
    {
        unsigned score = 0;
        extractTagValue(score, tagsDict, idx);
        
        return score;
    }
    return length(record.seq);
}

// ==========================================================================

bool
isGoodQuality(BamAlignmentRecord & record, Pair<CigarElement<>::TCount> & interval)
{
    if (interval.i2 - interval.i1 < 50)
        return false;

    if (interval.i2 - interval.i1 < length(record.seq) / 2)
        return false;

    if (avgQuality(record.qual, interval) <= 20)
        return false;

    if (alignmentScore(record) < 0.8 * (interval.i2 - interval.i1))
        return false;

    return true;
}

// ==========================================================================

unsigned
distanceToContigEnd(BamAlignmentRecord & record,
        Pair<CigarElement<>::TCount> & interval,
        StringSet<CharString> & names,
        std::map<CharString, unsigned> & contigLengths)
{
    if (hasFlagRC(record))
    {
        return record.beginPos;
    }
    else
    {
        unsigned endPos = record.beginPos + interval.i2 - interval.i1;
        return contigLengths[names[record.rID]] - endPos;
    }
}

// ==========================================================================

inline int
readAnchoringRecord(AnchoringRecord & record,
        std::map<Triple<CharString, CharString, unsigned>, unsigned> & goodReads,
        BamStream & stream,
        StringSet<CharString> & names,
        std::map<CharString, unsigned> & contigLengths)
{
    BamAlignmentRecord r;
    while (!atEnd(stream))
    {
        if (readRecord(r, stream) != 0)
        {
            std::cerr << "ERROR: Failed to read BAM alignment record." << std::endl;
            return 1;
        }

        if (r.rID == r.rNextId || r.rNextId == -1 || r.rID == -1)
            continue;

        Pair<CigarElement<>::TCount> interval = mappedInterval(r.cigar);
        if (!isGoodQuality(r, interval))
            continue;

        bool isContig = isComponentOrNode(names[r.rID]);

        if (!isContig && r.mapQ == 0)
            continue;

        if (isContig && distanceToContigEnd(r, interval, names, contigLengths) > 500)
            continue;

        Triple<CharString, CharString, unsigned> nameChrPos = Triple<CharString, CharString, unsigned>(r.qName, names[r.rNextId], r.pNext);
        if (goodReads.count(nameChrPos) == 0)
        {
            goodReads[Triple<CharString, CharString, unsigned>(r.qName, names[r.rID], r.beginPos)] = r.beginPos + interval.i2 - interval.i1;
        }
        else if (isContig)
        {
            record.chr = names[r.rNextId];
            record.chrStart = r.pNext;
            record.chrEnd = goodReads[nameChrPos];
            record.chrOri = !hasFlagNextRC(r);
            record.contig = names[r.rID];
            record.contigOri = !hasFlagRC(r);

            return 0;
        }
        else
        {
            record.chr = names[r.rID];
            record.chrStart = r.beginPos;
            record.chrEnd = r.beginPos + interval.i2 - interval.i1;
            record.chrOri = !hasFlagRC(r);
            record.contig = names[r.rNextId];
            record.contigOri = !hasFlagNextRC(r);

            return 0;
        }
    }
    return 2;
}

// ==========================================================================

void
listToLocs(String<Location> & locs, String<AnchoringRecord> & list, unsigned maxInsertSize)
{
    typedef Iterator<String<AnchoringRecord> >::Type TIter;

    TIter it = begin(list);
    TIter itEnd = end(list);

    if (it == itEnd)
        return;

    Location loc(*it);
    ++it;

    if (it == itEnd)
        return;

    while (it != itEnd)
    {
        if (loc.contig == (*it).contig && loc.chr == (*it).chr && loc.chrEnd + maxInsertSize >= (*it).chrStart)
        {
            loc.chrEnd = std::max(loc.chrEnd, (*it).chrEnd);
            ++loc.numReads;
        }
        else
        {
            appendValue(locs, loc);
            loc = Location(*it);
        }
        ++it;
    }
    appendValue(locs, loc);
}

// ==========================================================================
// Function findLocations()
// ==========================================================================

int
findLocations(String<Location> & locations, CharString & nonRefFile, unsigned maxInsertSize)
{
    typedef Pair<CharString, unsigned> TContigEnd;
    typedef std::map<TContigEnd, unsigned> TMap;
    typedef TMap::iterator TMapIter;

    BamStream inStream(toCString(nonRefFile), BamStream::READ);
    if (!isGood(inStream))
    {
        std::cerr << "ERROR: Could not open input bam file " << nonRefFile << std::endl;
        return 1;
    }

    std::map<CharString, unsigned> contigLengths;
    for (unsigned  i = 0; i < length(inStream.header.sequenceInfos); ++i)
    contigLengths[inStream.header.sequenceInfos[i].i1] = inStream.header.sequenceInfos[i].i2;

    StringSet<CharString> names = nameStore(inStream.bamIOContext);

    String<String<AnchoringRecord> > lists;
    resize(lists, 4);
    TMap anchorsToOther;

    unsigned i = 0;
    AnchoringRecord record;
    std::map<Triple<CharString, CharString, unsigned>, unsigned> goodReads; // Triple(qName, chrom, beginPos) -> alignEndPos
    while (!atEnd(inStream))
    {
        if (readAnchoringRecord(record, goodReads, inStream, names, contigLengths) == 1)
            return 1;

        if (record.chrOri)
        {
            if (record.contigOri) i = 0;
            else i = 1;
        }
        else
        {
            if (record.contigOri) i = 2;
            else i = 3;
        }

        if (isChromosome(record.chr))
            appendValue(lists[i], record);
        else
            ++anchorsToOther[TContigEnd(record.contig, i%2)];
    }
    
    for (unsigned i = 0; i < length(lists); ++i)
    {
        std::stable_sort(begin(lists[i]), end(lists[i]), AnchoringRecordLess());

        String<Location> locs;
        listToLocs(locs, lists[i], maxInsertSize);
        clear(lists[i]);

        append(locations, locs);
    }
    TMapIter endMap = anchorsToOther.end();
    for (TMapIter it = anchorsToOther.begin(); it != endMap; ++it)
        append(locations, Location("OTHER", 0, 0, true,
                                   (it->first).i1, ((it->first).i2 == 0 ? true : false), it->second, 0));
                                  
    // Sort locations by contig, contigOri, chr, chrStart, chrOri.
    LocationTypeLess less;
    std::stable_sort(begin(locations, Standard()), end(locations, Standard()), less);
    return 0;
}

// ==========================================================================
// Function scoreLocations()
// ==========================================================================

void
scoreLocations(String<Location> & locations)
{
    typedef Iterator<String<Location> >::Type TIterator;
    TIterator itEnd = end(locations);

    std::map<Pair<CharString, bool>, unsigned> readsPerContig;

    // Count total number of reads per contig.
    for (TIterator it = begin(locations); it != itEnd; ++it)
    {
        Pair<CharString, bool> c((*it).contig, (*it).contigOri);
        if (readsPerContig.count(c) == 0)
            readsPerContig[c] = (*it).numReads;
        else
            readsPerContig[c] += (*it).numReads;
    }

    // Compute the score for each location.
    for (TIterator it = begin(locations); it != itEnd; ++it)
    {
        Pair<CharString, bool> c((*it).contig, (*it).contigOri);
        (*it).score = (*it).numReads/(double)readsPerContig[c];
    }
}

// --------------------------------------------------------------------------
// Function readLocation()
// --------------------------------------------------------------------------

bool
readLocation(Location & loc, RecordReader<std::fstream, SinglePass<> > & reader, CharString & locationsFile)
{
    CharString buffer;

    if (readUntilOneOf(loc.chr, reader, ':', '\t', ' ') != 0)
    {
        std::cerr << "ERROR: Reading CHR from locations file " << locationsFile << " failed." << std::endl;
        return 1;
    }

    if (value(reader) == ':')
    {
        skipChar(reader, ':');

        if (readUntilChar(buffer, reader, '-') != 0)
        {
            std::cerr << "ERROR: Reading CHR_START from locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<Location::TPos>(loc.chrStart, buffer);
        skipChar(reader, '-');

        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading CHR_END from locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<Location::TPos>(loc.chrEnd, buffer);
    }
    skipWhitespaces(reader);

    if (value(reader) == '+') loc.chrOri = true;
    else if (value(reader) == '-') loc.chrOri = false;
    else
    {
        std::cerr << "ERROR: Reading CHR_ORI from locations file " << locationsFile << " failed." << std::endl;
        return 1;
    }
    skipNChars(reader, 1);
    skipWhitespaces(reader);

    if (readUntilWhitespace(loc.contig, reader) != 0)
    {
        std::cerr << "ERROR: Reading CONTIG from locations file " << locationsFile << " failed." << std::endl;
        return 1;
    }
    skipWhitespaces(reader);

    if (value(reader) == '+') loc.contigOri = true;
    else if (value(reader) == '-') loc.contigOri = false;
    else
    {
        std::cerr << "ERROR: Reading CONTIG_ORI from locations file " << locationsFile << " failed." << std::endl;
        return 1;
    }
    skipNChars(reader, 1);
    skipWhitespaces(reader);

    clear(buffer);
    if (readDigits(buffer, reader) != 0)
    {
        std::cerr << "ERROR: Reading NUM_READS from locations file " << locationsFile << " failed." << std::endl;
        return 1;
    }
    lexicalCast2<unsigned>(loc.numReads, buffer);
    skipWhitespaces(reader);

    clear(buffer);
    if (readUntilWhitespace(buffer, reader) != 0)
    {
        std::cerr << "ERROR: Reading SCORE from line locations file " << locationsFile << " failed." << std::endl;
        return 1;
    }
    lexicalCast2<double>(loc.score, buffer);

    if (value(reader) != '\r' && value(reader) != '\n') // if not at end of line
    {
        skipWhitespaces(reader);

        CharString sampleName;
        if (readUntilChar(sampleName, reader, ':') != 0)
        {
            std::cerr << "ERROR: Reading SAMPLE_ID from locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        if (loc.bestSamples.count(sampleName) != 0)
        {
            std::cerr << "ERROR: Sample " << sampleName << " listed twice in " << locationsFile << " for " << loc.chr << ":" << loc.chrStart << "-" << loc.chrEnd << "." << std::endl;
            return 1;
        }
        skipChar(reader, ':');

        CharString numBuffer;
        if (readDigits(numBuffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading NUM_READS for sample " << sampleName << " from locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        unsigned readCount;
        lexicalCast2<unsigned>(readCount, numBuffer);
        loc.bestSamples[sampleName] = readCount;

        while (value(reader) != '\r' && value(reader) != '\n') // while not at end of line
        {
            skipChar(reader, ',');

            clear(sampleName);
            if (readUntilChar(sampleName, reader, ':') != 0)
            {
                std::cerr << "ERROR: Reading next SAMPLE_ID from locations file " << locationsFile << " failed." << std::endl;
                return 1;
            }
            if (loc.bestSamples.count(sampleName) != 0)
            {
                std::cerr << "ERROR: Sample " << sampleName << " listed twice in " << locationsFile << " for " << loc.chr << ":" << loc.chrStart << "-" << loc.chrEnd << "." << std::endl;
                return 1;
            }
            skipChar(reader, ':');

            clear(numBuffer);
            if (readDigits(numBuffer, reader) != 0)
            {
                std::cerr << "ERROR: Reading NUM_READS for sample " << sampleName << " from locations file " << locationsFile << " failed." << std::endl;
                return 1;
            }
            lexicalCast2<unsigned>(readCount, numBuffer);
            loc.bestSamples[sampleName] = readCount;
        }
    }
    else
    {
        CharString sampleName = locationsFile; // TODO: Replace by an actual sample ID.
        if (suffix(sampleName, length(sampleName) - 14) == "/locations.txt")
            sampleName = prefix(sampleName, length(sampleName) - 14);
        if (loc.bestSamples.count(sampleName) != 0)
        {
            std::cerr << "ERROR: Sample " << sampleName << " listed twice in " << locationsFile << " for " << loc.chr << ":" << loc.chrStart << "-" << loc.chrEnd << "." << std::endl;
            return 1;
        }
        loc.bestSamples[sampleName] = loc.numReads;
    }
    skipLine(reader);

    return 0;
}

// ==========================================================================
// Function appendLocation()
// ==========================================================================

void
appendLocation(String<Location> & locs, Location & loc)
{
    appendValue(locs, loc);
}

// ==========================================================================
// Function readLocations()
// ==========================================================================

template<typename TLoc>
int
readLocations(String<TLoc> & locations, CharString & locationsFile, LocationsFilter & filterParams)
{
    std::fstream stream(toCString(locationsFile), std::ios::in);
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open locations file " << locationsFile << std::endl;
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > reader(stream);

    for (unsigned line = 0; !atEnd(reader); ++line)
    {
        Location loc;

        if (readLocation(loc, reader, locationsFile) != 0)
            return 1;

        if (passesFilter(loc, filterParams))
            appendLocation(locations, loc);
    }
    return 0;
}

template<typename TLoc>
int
readLocations(String<TLoc> & locations, CharString & locationsFile, Triple<CharString, unsigned, unsigned> & interval, LocationsFilter & filterParams)
{
    std::fstream stream(toCString(locationsFile), std::ios::in);
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open locations file " << locationsFile << std::endl;
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > reader(stream);

    for (unsigned line = 0; !atEnd(reader); ++line)
    {
        Location loc;

        if (readLocation(loc, reader, locationsFile) != 0)
            return 1;

        if (passesFilter(loc, filterParams) && loc.chr == interval.i1 && loc.chrStart >= interval.i2 && loc.chrStart < interval.i3)
            appendLocation(locations, loc);
    }
    return 0;
}

// ==========================================================================
// Function writeLoc()
// ==========================================================================

void
writeLoc(std::fstream & stream, Location & loc)
{
    stream << loc.chr;
    if (loc.chr != "OTHER")
    {
        stream << ":";
        stream << loc.chrStart << "-";
        stream << loc.chrEnd;
    }
    stream << "\t" << (loc.chrOri ? "+" : "-");
    stream << "\t" << loc.contig;
    stream << "\t" << (loc.contigOri ? "+" : "-");
    stream << "\t" << loc.numReads;
    if (loc.score != -1) stream << "\t" << loc.score;
    if (length(loc.bestSamples) > 0)
    {
        String<Pair<unsigned, CharString> > bestSamples;
        for (std::map<CharString, unsigned>::iterator bsIt = loc.bestSamples.begin(); bsIt != loc.bestSamples.end(); ++bsIt)
            appendValue(bestSamples, Pair<unsigned, CharString>(bsIt->second, bsIt->first));

        std::stable_sort(begin(bestSamples), end(bestSamples), std::greater<Pair<unsigned, CharString> >());
        if (length(bestSamples) > 100)
            resize(bestSamples, 100);

        stream << "\t" << bestSamples[0].i2 << ":" << bestSamples[0].i1;
        for (unsigned i = 1; i < length(bestSamples); ++i)
            stream << "," << bestSamples[i].i2 << ":" << bestSamples[i].i1;
    }
    stream << std::endl;
}

// ==========================================================================
// Function writeLocations()
// ==========================================================================

int
writeLocations(std::fstream & stream, String<Location> & locations)
{
    typedef Iterator<String<Location> >::Type TIterator;
    
    // Iterate over locations to output them one per line.
    TIterator itEnd = end(locations);
    for (TIterator it = begin(locations); it != itEnd; ++it)
        writeLoc(stream, *it);
    
    return 0;
}

int
writeLocations(CharString & filename, String<Location> & locations)
{
    std::fstream stream(toCString(filename), std::ios::out);
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open temporary locations file " << filename << " for writing." << std::endl;
        return 1;
    }
    
    return writeLocations(stream, locations);
}

// --------------------------------------------------------------------------
// addLocation()
// --------------------------------------------------------------------------

void
addLocation(Location & prevLoc, String<Location> & locations, Location & loc, unsigned maxInsertSize)
{
    if (prevLoc.contig == "")
    {
        loc.score = -1;
        prevLoc = loc;
    }
    else if (prevLoc.chr != loc.chr || prevLoc.chrEnd + maxInsertSize < loc.chrStart)
    {
        appendValue(locations, prevLoc);
        loc.score = -1;
        prevLoc = loc;
    }
    else
    {
        prevLoc.chrEnd = std::max(prevLoc.chrEnd, loc.chrEnd);
        prevLoc.numReads += loc.numReads;

        for (std::map<CharString, unsigned>::iterator bsIt = loc.bestSamples.begin(); bsIt != loc.bestSamples.end(); ++bsIt)
        {
            if (prevLoc.bestSamples.count(bsIt->first) == 0)
                prevLoc.bestSamples[bsIt->first] = bsIt->second;
            else
                prevLoc.bestSamples[bsIt->first] += bsIt->second;
        }
    }
}

// --------------------------------------------------------------------------
// mergeLocationsBatch()
// --------------------------------------------------------------------------

int mergeLocationsBatch(std::fstream & stream, String<Location> & locations, String<CharString> locationsFiles, size_t offset, size_t batchSize, unsigned maxInsertSize)
{
    typedef RecordReader<std::fstream, SinglePass<> > TReader;
    typedef Pair<std::fstream *, TReader *> TPtrPair;

    Location forward, reverse;
    unsigned contigCount = 0;

    // define min heap of Locations
    std::priority_queue<Location, std::vector<Location>, LocationTypeGreater> heap;

    unsigned last = std::min(offset+batchSize, length(locationsFiles));
    
    
    // Open files and store String of reader pointers.
    String<TPtrPair> readerPtr;
    resize(readerPtr, length(locationsFiles));
    for (unsigned i = offset; i < last; ++i)
    {
        std::fstream * stream = new std::fstream(toCString(locationsFiles[i]), std::ios::in);
        if (!(*stream).good())
        {
            std::cerr << "ERROR: Could not open locations file " << locationsFiles[i] << std::endl;
            return 1;
        }
        readerPtr[i] = TPtrPair(stream, new TReader(*stream));

        // Read the first location record and push it to min heap.
        Location loc;
        if (readLocation(loc, *readerPtr[i].i2, locationsFiles[i]) != 0) return 1;
        loc.fileIndex = i;
        heap.push(loc);
    }

    // Iterate over all files simultaneously using the min heap.
    while (!heap.empty())
    {
        Location loc = heap.top();

        // Output all the locations for a contig.
        if ((forward.contig != "" && (forward.contig != loc.contig || (forward.contig == loc.contig && forward.contigOri != loc.contigOri))) ||
            (reverse.contig != "" && (reverse.contig != loc.contig || (reverse.contig == loc.contig && reverse.contigOri != loc.contigOri))))
        {
            if (forward.contig != "") appendValue(locations, forward);
            if (reverse.contig != "") appendValue(locations, reverse);

            // Compute the score for each location.
            Iterator<String<Location> >::Type itEnd = end(locations);
            for (Iterator<String<Location> >::Type it = begin(locations); it != itEnd; ++it)
                (*it).score = (*it).numReads/(double)contigCount;

            LocationTypeLess less;
            std::stable_sort(begin(locations, Standard()), end(locations, Standard()), less);
            if (length(locations) > 0) writeLocations(stream, locations);

            clear(locations);
            forward = Location();
            reverse = Location();
            contigCount = 0;
        }

        contigCount += loc.numReads;
        if (loc.chrOri) addLocation(forward, locations, loc, maxInsertSize);
        else addLocation(reverse, locations, loc, maxInsertSize);

        heap.pop();
        unsigned i = loc.fileIndex;
        if (!atEnd(*readerPtr[i].i2))
        {
            Location nextLoc;
            if (readLocation(nextLoc, *readerPtr[i].i2, locationsFiles[i]) != 0) return 1;
            nextLoc.fileIndex = i;
            heap.push(nextLoc);
        }
    }

    // Append the remaining locations.
    Location loc = Location();
    if (forward.contig != "") appendValue(locations, forward);
    if (reverse.contig != "") appendValue(locations, reverse);
    
    // Compute the score for each location.
    Iterator<String<Location> >::Type itEnd = end(locations);
    for (Iterator<String<Location> >::Type it = begin(locations); it != itEnd; ++it)
        (*it).score = (*it).numReads/(double)contigCount;

    LocationTypeLess less;
    std::stable_sort(begin(locations, Standard()), end(locations, Standard()), less);
    if (length(locations) > 0) writeLocations(stream, locations);

    // clean-up
    for (unsigned i = offset; i < last; ++i)
    {
        delete readerPtr[i].i2;
        (*readerPtr[i].i1).close();
        delete readerPtr[i].i1;
    }

    return 0;
}

// ==========================================================================
// Function mergeLocations()
// ==========================================================================

int
mergeLocations(std::fstream & stream, String<Location> & locations, String<CharString> & locationsFiles, CharString & outFile, unsigned maxInsertSize, bool verbose)
{
    String<CharString> tmpFiles;
    unsigned batchSize = 500;
    
    if (length(locationsFiles) > batchSize*batchSize)
    {
        std::cerr << "ERROR: Too many locations files, max: " << batchSize*batchSize << ", given: " << length(locationsFiles) << std::endl;
        return 1;
    }
    else if (length(locationsFiles) <= batchSize)
    {
        mergeLocationsBatch(stream, locations, locationsFiles, 0, batchSize, maxInsertSize);
        return 0;
    }

    // Write a temporary file for each batch of files.
    for (unsigned offset = 0; offset < length(locationsFiles); offset += batchSize)
    {
        if (verbose) std::cerr << "[" << time(0) << "] " << "Merging batch " << (offset/batchSize)+1 << " of location files." << std::endl;

        // Create temporary file name.
        std::stringstream tmpName;
        tmpName << outFile << "." << (offset/batchSize)+1;
        CharString tmp = tmpName.str();
        appendValue(tmpFiles, tmp);

        // Open output file.
        std::fstream tmpStream(toCString(tmp), std::ios::out);
        if (!tmpStream.good())
        {
            std::cerr << "ERROR: Could not open temporary locations file " << tmp << " for writing." << std::endl;
            return 1;
        }

        String<Location> locs;
        if (mergeLocationsBatch(tmpStream, locs, locationsFiles, offset, batchSize, maxInsertSize) != 0) return 1;
    }

    // Merge and remove the temporary files.
    if (verbose) std::cerr << "[" << time(0) << "] " << "Merging temporary location files." << std::endl;
    if (mergeLocationsBatch(stream, locations, tmpFiles, 0, length(tmpFiles), maxInsertSize) != 0) return 1;
    for (unsigned i = 0; i < length(tmpFiles); ++i) remove(toCString(tmpFiles[i]));

    return 0;
}

#endif // #ifndef POPINS_LOCATION_H_
