#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

#ifndef POPINS_LOCATION_H_
#define POPINS_LOCATION_H_

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
            std::cerr << "ERROR: Could not read BAM alignment record." << std::endl;
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
            goodReads[ Triple<CharString, CharString, unsigned>(r.qName, names[r.rID], r.beginPos)] = r.beginPos + interval.i2 - interval.i1;
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

bool
addReadToLocs(String<Location> & locs, AnchoringRecord & record)
{
    typedef Iterator<String<Location> >::Type TIter;

    TIter beginIt = begin(locs);
    TIter it = end(locs);
    if (it != beginIt)
    {
        it--;
        while (it >= beginIt && record.chr == (*it).chr && record.chrStart < (*it).chrEnd + 1000)
        {
            if (record.contig == (*it).contig)
            {
                if (record.chrStart < (*it).chrStart) (*it).chrStart = record.chrStart;
                if (record.chrEnd > (*it).chrEnd) (*it).chrEnd = record.chrEnd;
                ++(*it).numReads;
                return 0;
            }
            --it;
        }
    }
    appendValue(locs, Location(record));
    
    return 0;
}

// ==========================================================================
// Function findLocations()
// ==========================================================================

int
findLocations(String<Location> & locations, CharString & nonRefFile)
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

    String<String<Location> > locs;
    resize(locs, 4);
    TMap anchorsToOther;

    unsigned i = 0;
    AnchoringRecord record;
    std::map<Triple<CharString, CharString, unsigned>, unsigned> goodReads; // Triple(qName, chrom, beginPos) -> alignEndPos
    while (!atEnd(inStream))
    {
        readAnchoringRecord(record, goodReads, inStream, names, contigLengths);
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
        if (!isChromosome(record.chr) || addReadToLocs(locs[i], record) == 1)
            ++anchorsToOther[TContigEnd(record.contig, i%2)];
    }
    
    for (unsigned i = 0; i < length(locs); ++i)
    {
        append(locations, locs[i]);
        clear(locs[i]);
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
// Function readLocations()
// ==========================================================================

int
readLocations(String<Location> & locations, CharString & locationsFile, Triple<CharString, unsigned, unsigned> & interval)
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
        if (readLocation(loc, reader, locationsFile) != 0) return 1;
        if (loc.chr == interval.i1 && loc.chrStart >= interval.i2 && loc.chrStart < interval.i3)
            appendValue(locations, loc);
    }
    return 0;
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
    {
        stream << (*it).chr;
        if ((*it).chr != "OTHER")
        {
            stream << ":";
            stream << (*it).chrStart << "-";
            stream << (*it).chrEnd;
        }
        stream << "\t" << ((*it).chrOri ? "+" : "-");
        stream << "\t" << (*it).contig;
        stream << "\t" << ((*it).contigOri ? "+" : "-");
        stream << "\t" << (*it).numReads;
        if ((*it).score != -1) stream << "\t" << (*it).score;
        if (length((*it).bestSamples) > 0)
        {
            String<Pair<unsigned, CharString> > bestSamples;
            for (std::map<CharString, unsigned>::iterator bsIt = (*it).bestSamples.begin(); bsIt != (*it).bestSamples.end(); ++bsIt)
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
addLocation(Location & prevLoc, String<Location> & locations, Location & loc)
{
    if (prevLoc.contig == "")
    {
        loc.score = -1;
        prevLoc = loc;
    }
    else if (prevLoc.chr != loc.chr || prevLoc.chrEnd + 1000 < loc.chrStart)
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

int mergeLocationsBatch(std::fstream & stream, String<Location> & locations, String<CharString> locationsFiles, size_t offset, size_t batchSize)
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
        if (loc.chrOri) addLocation(forward, locations, loc);
        else addLocation(reverse, locations, loc);

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
mergeLocations(std::fstream & stream, String<Location> & locations, String<CharString> & locationsFiles, CharString & outFile, bool verbose)
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
        mergeLocationsBatch(stream, locations, locationsFiles, 0, batchSize);
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
        if (mergeLocationsBatch(tmpStream, locs, locationsFiles, offset, batchSize) != 0) return 1;
    }

    // Merge and remove the temporary files.
    if (verbose) std::cerr << "[" << time(0) << "] " << "Merging temporary location files." << std::endl;
    if (mergeLocationsBatch(stream, locations, tmpFiles, 0, length(tmpFiles)) != 0) return 1;
    for (unsigned i = 0; i < length(tmpFiles); ++i) remove(toCString(tmpFiles[i]));

    return 0;
}

// ==========================================================================

bool
overlaps(GenomicInterval & i1, GenomicInterval & i2)
{
    if (i1.chr != i2.chr) return false;
    if (i1.ori != i2.ori) return false;
    if (i1.end <= i2.begin) return false;
    if (i2.end <= i1.begin) return false;
    return true;
}

// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
void
verifyCandidateGroup(std::map<TSize, std::set<TSize> > & groups, String<TSize> & candidates, String<Location> & locations, std::map<CharString, TSeq> & contigs)
{
    std::map<TSize, Pair<TSize> > aligned; // which location aligns to which other location with which offset (loc1 <- Pair(loc2, offset))
    
    for (TSize i = 0; i < length(candidates)-1; ++i)
    {
        TSeq contig_i = contigs[locations[candidates[i]].contig];
        if (locations[candidates[i]].contigOri) reverseComplement(contig_i);

        for (TSize j = i+1; j < length(candidates); ++j)
        {
            if (aligned.count(candidates[j]) > 0) continue;
            
            TSeq contig_j = contigs[locations[candidates[j]].contig];
            if (locations[candidates[j]].contigOri) reverseComplement(contig_j);

            Align<TSeq> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), prefix(contig_i, std::min((TSize)200, length(contig_i))));
            assignSource(row(align, 1), prefix(contig_j, std::min((TSize)200, length(contig_j))));

            int score = globalAlignment(align, Score<int, Simple>(1, -2, -1, -4), AlignConfig<true, true, true, true>());

            if (score < 50) continue;

            TSize alignBegin = std::max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
            TSize alignEnd = std::min(toViewPosition(row(align, 0), length(contig_i)), toViewPosition(row(align, 1), length(contig_j)));
//            std::cerr << "Score: " << score << ", Relative score: " << (score / ((double) alignEnd - alignBegin)) << std::endl;
//            std::cerr << align << std::endl;
            
            if (score / ((double) alignEnd - alignBegin) < 0.8) continue;

            // add to map
            if (toViewPosition(row(align, 0), 0) == 0)
                aligned[candidates[j]] = Pair<TSize>(candidates[i], toViewPosition(row(align, 1), 0));
            else
                aligned[candidates[j]] = Pair<TSize>(candidates[i], toViewPosition(row(align, 0), 0));
            break;
        }
    }
    
    
    // Make the groups from aligned
    typename Iterator<String<TSize> >::Type it = begin(candidates);
    typename Iterator<String<TSize> >::Type itEnd = end(candidates);
    while (it != itEnd)
    {
        if (aligned.count(*it) == 0)
        {
            groups[*it];
            ++it;
        }
        else
        {
            TSize i = aligned[*it].i1;
            if (aligned.count(i) > 0)
            {
                aligned[*it].i1 = aligned[i].i1;
                aligned[*it].i2 += aligned[i].i2;
            }
            else // i is group representative
            {
                groups[i].insert(*it);
                ++it;
            }
        }
    }
}

// ==========================================================================
// Function groupLocations()
// ==========================================================================

template<typename TSize, typename TSeq>
void
groupLocations(std::map<TSize, std::set<TSize> > & groups, String<Location> & locations, std::map<CharString, TSeq> & contigs)
{
    String<TSize> candidates;
    appendValue(candidates, 0);
    GenomicInterval interval(locations[0].chr, locations[0].chrStart, locations[0].chrEnd, locations[0].chrOri);
    for (TSize i = 1; i < length(locations); ++i)
    {
        // Find candidate group of locations that overlap on reference
        for (; i < length(locations); ++i)
        {
            GenomicInterval locInterval(locations[i].chr, locations[i].chrStart, locations[i].chrEnd, locations[i].chrOri);
            if (!overlaps(interval, locInterval))
            {
                interval = locInterval;
                break;
            }
            interval.begin = _min(interval.begin, locInterval.begin);
            interval.end = _max(interval.end, locInterval.end);
            appendValue(candidates, i);
        }
        verifyCandidateGroup(groups, candidates, locations, contigs);
        clear(candidates);

        if (i < length(locations))
            appendValue(candidates, i);
    }
    
    if (length(candidates) > 0)
        verifyCandidateGroup(groups, candidates, locations, contigs);

//    // Debug code: Output groups.
//    for(typename std::map<TSize, std::set<TSize> >::iterator it = groups.begin(); it != groups.end(); ++it)
//    {
//        std::cerr << it->first << " <-";
//        for(typename std::set<TSize>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2)
//            std::cerr << " " << *it2;
//        std::cerr << std::endl;
//    }
}

// =======================================================================================
// Function loadLocations()
// =======================================================================================

int
loadLocations(String<Location> & locations, PlacingOptions & options)
{
    if (!exists(options.locationsFile))
    {
        if (length(options.locationsFiles) == 0)
        {
            std::cerr << "ERROR: Locations file " << options.locationsFile << "does not exist. Specify -l option to create it." << std::endl;
            return -1;
        }

        // Open output file.
        std::fstream stream(toCString(options.locationsFile), std::ios::out);
        if (!stream.good())
        {
            std::cerr << "ERROR: Could not open locations file " << options.locationsFile << " for writing." << std::endl;
            return -1;
        }

        // Merge approximate locations and write them to a file.
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Merging locations files." << std::endl;
        if (mergeLocations(stream, locations, options.locationsFiles, options.locationsFile, options.verbose) != 0) return -1;
    }
    else
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Locations file exists." << std::endl;
    }

    if (length(options.bamFiles) == 0)
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "No split mapping. Specify -b option for split mapping." << std::endl;
        return 1;
    }
    
    // Read the locations file.
    if (length(locations) == 0)
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Reading locations in " << options.interval.i1 << ":" << options.interval.i2 << "-" << options.interval.i3
                                                                 << " from " << options.locationsFile << std::endl;
        if (readLocations(locations, options.locationsFile, options.interval) != 0) return -1;
        if (options.verbose) std::cerr << "[" << time(0) << "] " << length(locations) << " locations loaded." << std::endl;
    }
    else
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Sorting locations." << std::endl;
        LocationPosLess less;
        std::stable_sort(begin(locations, Standard()), end(locations, Standard()), less);
    }

    // Discard locations with score below options.minLocScore or OTHER or longer than 2*maxInsertSize // TODO: move this to reading function!
    unsigned i = 0;
    while (i < length(locations))
    {
        if (locations[i].score < options.minLocScore || locations[i].chr == "OTHER" ||
            locations[i].chrEnd - locations[i].chrStart > 2*options.maxInsertSize)
            replace(locations, i, i+1, String<Location>());
        else ++i;
    }

    if (length(locations) == 0)
    {
        std::cerr << "[" << time(0) << "] " << "No locations on genome left after filtering for score >= " << options.minLocScore << std::endl;
        return 1;
    }
    
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Keeping " << length(locations) << " locations with score >= "
                  << options.minLocScore << " and shorter than " << (2*options.maxInsertSize) << std::endl;

    return 0;
}

#endif // #ifndef POPINS_LOCATION_H_
