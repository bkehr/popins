#include <iostream>
#include <fstream>
#include <queue>

#include <seqan/sequence.h>
#include <seqan/stream.h>

#ifndef POPINS_LOCATION_H_
#define POPINS_LOCATION_H_

using namespace seqan;

// ==========================================================================
// struct AnchoringRecord
// ==========================================================================

struct AnchoringRecord
{
    typedef Position<CharString>::Type TPos;

    CharString chr;
    TPos chrStart;
    TPos chrEnd;
    bool chrOri;

    CharString contig;
    bool contigOri;
    
    unsigned bamStreamIndex;
};

// ==========================================================================

struct Greater : public ::std::binary_function <AnchoringRecord, AnchoringRecord, bool> 
{
    Greater() {}
    
    inline int compare (AnchoringRecord const & a, AnchoringRecord const & b) const
    {
        if (std::isdigit(a.chr[0]) && std::isdigit(b.chr[0]))
        {
            int chrA, chrB;
            lexicalCast2<int>(chrA, a.chr);
            lexicalCast2<int>(chrB, b.chr);
            if (chrA > chrB) return -1;
            if (chrA < chrB) return 1;
        }
        else
        {
            if (a.chr > b.chr) return -1;
            if (a.chr < b.chr) return 1;
        }

        if (a.chrStart > b.chrStart) return -1;
        if (a.chrStart < b.chrStart) return 1;
        
        if (a.chrEnd > b.chrEnd) return -1;
        if (a.chrEnd < b.chrEnd) return 1;
        
        if (a.chrOri > b.chrOri) return -1;
        if (a.chrOri < b.chrOri) return 1;
        
        if (a.contig > b.contig) return -1;
        if (a.contig < b.contig) return 1;
        
        if (a.contigOri > b.contigOri) return -1;
        if (a.contigOri < b.contigOri) return 1;
        
        return 0;
    }
    
    inline bool operator() (AnchoringRecord const & a, AnchoringRecord const & b) const
    {
        return compare(a, b) == -1;
    }
};

// ==========================================================================
// struct Location
// ==========================================================================

struct Location
{
    typedef Position<CharString>::Type TPos;

    CharString chr;
    TPos chrStart;
    TPos chrEnd;
    bool chrOri;

    CharString contig;
    bool contigOri;

    unsigned numReads;
    double score;
    
    Location () {}
    
    Location (AnchoringRecord const & r) :
        chr(r.chr), chrStart(r.chrStart), chrEnd(r.chrEnd), chrOri(r.chrOri),
        contig(r.contig), contigOri(r.contigOri),
        numReads(1)
    {}
};

// ==========================================================================

bool
isComponentOrNode(CharString & name)
{
    typedef Position<CharString>::Type TPos;

    if (prefix(name, 4) == "COMP")
        return true;

    TPos i = 0;
    for (; i < length(name); ++i)
        if (name[i] == '.') break;
    ++i;
    if (i < length(name)-4 && infix(name, i, i+4) == "NODE")
        return true;

    return false;
}

// ==========================================================================

Pair<CigarElement<>::TCount>
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
        case 'S': case 'H':
            if (it == begin(cigar))
                beginPos += (*it).count;
            break;
        case 'D':
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

double
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
    return totalQual/(interval.i2-interval.i1);
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

int
readAnchoringRecord(AnchoringRecord & record, BamStream & stream, unsigned bamStreamIndex)
{
    StringSet<CharString> names = nameStore(stream.bamIOContext);

    BamAlignmentRecord r;
    while (!atEnd(stream))
    {
        if (readRecord(r, stream) != 0)
        {
            std::cerr << "ERROR: Could not read bam alignment record." << std::endl;
            return 1;
        }
        
        if (r.rID == r.rNextId || r.rNextId == -1 || r.rID == -1) continue;
        
        if (!isComponentOrNode(names[r.rID]) &&    // rId is not NODE or COMP
            isComponentOrNode(names[r.rNextId]))   // rNextId is NODE or COMP
        {
            Pair<CigarElement<>::TCount> interval = mappedInterval(r.cigar);
            if (interval.i2 - interval.i1 >= 50 && interval.i2 - interval.i1 >= length(r.seq) / 2 &&
                avgQuality(r.qual, interval) > 20 &&
                alignmentScore(r) > 0.8 * (interval.i2-interval.i2))
            {
                record.chr = names[r.rID];
                record.chrStart = r.beginPos + interval.i1;
                record.chrEnd = r.beginPos + interval.i2;
                record.chrOri = !hasFlagRC(r);
                record.contig = names[r.rNextId];
                record.contigOri = !hasFlagNextRC(r);
                record.bamStreamIndex = bamStreamIndex;
                return 0;
            }
        }
    }
    return 2;
}

// ==========================================================================

void
addRead(Location & loc, String<Location> & locations, AnchoringRecord & record)
{
    if (loc.chr == "")
    {
        loc = Location(record);
    }
    else if (record.contig != loc.contig || record.chr != loc.chr || record.chrStart > loc.chrStart + 1000)
    {
        appendValue(locations, loc);
        loc = Location(record);
    }
    else
    {
        ++loc.numReads;
        loc.chrEnd = std::max(loc.chrEnd, record.chrEnd);
    }
}

// ==========================================================================
// Function findLocations()
// ==========================================================================

int
findLocations(String<Location> & locations, String<CharString> & nonRefFiles)
{
    Location forwardForward, forwardReverse, reverseForward, reverseReverse;

    // define min heap of AnchoringRecord
    std::priority_queue<AnchoringRecord, std::vector<AnchoringRecord>, Greater> heap;
    
    AnchoringRecord record;

    // Open non ref files and store String of stream pointers
    String<BamStream *> streamPtr;
    resize(streamPtr, length(nonRefFiles));

    for (unsigned i = 0; i < length(nonRefFiles); ++i)
    {
        streamPtr[i] = new BamStream(toCString(nonRefFiles[i]), BamStream::READ);
        if (!isGood(*streamPtr[i]))
        {
            std::cerr << "ERROR: Could not open input bam file " << nonRefFiles[i] << std::endl;
            delete streamPtr[i]; // TODO: delete for all i?
            return 1;
        }
        
        // Read the first anchoring read record and push it to min heap
        int res = readAnchoringRecord(record, *streamPtr[i], i);
        if (res == 0) heap.push(record);
        else if (res == 1) return 1;
    }

    while (!heap.empty())
    {
        record = heap.top();
        
        // Append and reset location if smallest record doesn't support the current location with corresponding
        // orientation. Otherwise add it to location.
        if (record.chrOri)
        {
            if (record.contigOri) addRead(forwardForward, locations, record);
            else addRead(forwardReverse, locations, record);
        }
        else
        {
            if (record.contigOri) addRead(reverseForward, locations, record);
            else addRead(reverseReverse, locations, record);
        }
        
        heap.pop();
        unsigned i = record.bamStreamIndex;
        int res = readAnchoringRecord(record, *streamPtr[i], i);
        if (res == 0) heap.push(record);
        else if (res == 1) return 1;
    }
    
    for (unsigned i = 0; i < length(nonRefFiles); ++i)
        delete streamPtr[i];

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

// ==========================================================================
// Function readLocations()
// ==========================================================================

int
readLocations(String<Location> & locations, CharString & locationsFile, unsigned batchSize, unsigned batchIndex)
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
        if (batchSize != maxValue<unsigned>() && line < batchIndex*batchSize)
        {
            skipLine(reader);
            continue;
        }

        if (line >= batchIndex*batchSize+batchSize) break;
        
        Location loc;
        
        if (readUntilChar(loc.chr, reader, ':') != 0)
        {
            std::cerr << "ERROR: Reading CHR from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        skipChar(reader, ':');
        
        CharString buffer;
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading CHR_START from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<Location::TPos>(loc.chrStart, buffer);
        skipChar(reader, '-');
        
        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading CHR_END from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<Location::TPos>(loc.chrEnd, buffer);
        skipWhitespaces(reader);
        
        if (value(reader) == '+') loc.chrOri = true;
        else if (value(reader) == '-') loc.chrOri = false;
        else
        {
            std::cerr << "ERROR: Reading CHR_ORI from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        skipNChars(reader, 1);
        skipWhitespaces(reader);
        
        if (readUntilWhitespace(loc.contig, reader) != 0)
        {
            std::cerr << "ERROR: Reading CONTIG from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        skipWhitespaces(reader);

        if (value(reader) == '+') loc.contigOri = true;
        else if (value(reader) == '-') loc.contigOri = false;
        else
        {
            std::cerr << "ERROR: Reading CONTIG_ORI from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        skipNChars(reader, 1);
        skipWhitespaces(reader);
        
        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading NUM_READS from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<unsigned>(loc.numReads, buffer);
        skipWhitespaces(reader);
        
        clear(buffer);
        if (readUntilWhitespace(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading SCORE from line " << line << " of locations file " << locationsFile << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<double>(loc.score, buffer);
        skipLine(reader);
        
        appendValue(locations, loc);
    }
    return 0;
}

// ==========================================================================
// Function writeLocations()
// ==========================================================================

int
writeLocations(CharString & locationsFile, String<Location> & locations)
{
    typedef Iterator<String<Location> >::Type TIterator;
    
    // Open output file.
    std::fstream stream(toCString(locationsFile), std::ios::out);
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open locations file " << locationsFile << " for writing." << std::endl;
        return 1;
    }
    
    // Iterate over locations to output them one per line.
    TIterator itEnd = end(locations);
    for (TIterator it = begin(locations); it != itEnd; ++it)
    {
        stream << (*it).chr << ":";
        stream << (*it).chrStart << "-";
        stream << (*it).chrEnd << "\t";
        stream << ((*it).chrOri ? "+" : "-") << "\t";
        stream << (*it).contig << "\t";
        stream << ((*it).contigOri ? "+" : "-") << "\t";
        stream << (*it).numReads << "\t";
        stream << (*it).score << std::endl;
    }
    
    return 0;
}

#endif // #ifndef POPINS_LOCATION_H_
