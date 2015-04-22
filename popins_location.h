#include <iostream>
#include <sstream>
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
    
    unsigned fileIndex;
    
    Location ()
    {
        chr = "";
        chrStart = 0;
        chrEnd = 0;
        contig = "";
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

inline int
readAnchoringRecord(AnchoringRecord & record, BamStream & stream, StringSet<CharString> & names)//, unsigned bamStreamIndex)
{
    BamAlignmentRecord r;
    while (!atEnd(stream))
    {
        if (readRecord(r, stream) != 0)
        {
            std::cerr << "ERROR: Could not read bam alignment record." << std::endl;
            return 1;
        }
        
        if (r.rID == r.rNextId || r.rNextId == -1 || r.rID == -1) continue;
        
        if (isComponentOrNode(names[r.rNextId]))   // rNextId is NODE or COMP
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
                //record.bamStreamIndex = bamStreamIndex;
                return 0;
            }
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
    typedef Iterator<String<CharString> >::Type TIter;
    typedef Pair<CharString, unsigned> TContigEnd;
    typedef std::map<TContigEnd, unsigned> TMap;
    typedef TMap::iterator TMapIter;

    BamStream inStream(toCString(nonRefFile), BamStream::READ);
    if (!isGood(inStream))
    {
        std::cerr << "ERROR: Could not open input bam file " << nonRefFile << std::endl;
        return 1;
    }
    StringSet<CharString> names = nameStore(inStream.bamIOContext);

    String<String<Location> > locs;
    resize(locs, 4);
    TMap anchorsToOther;

    unsigned i = 0;
    AnchoringRecord record;
    while (!atEnd(inStream))
    {
        readAnchoringRecord(record, inStream, names);
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
    skipLine(reader);

    return 0;
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
        if (readLocation(loc, reader, locationsFile) != 0) return 1;
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
        stream << (*it).chr;
        if ((*it).chr != "OTHER")
        {
            stream << ":";
            stream << (*it).chrStart << "-";
            stream << (*it).chrEnd;
        }
        stream << "\t";
        stream << ((*it).chrOri ? "+" : "-") << "\t";
        stream << (*it).contig << "\t";
        stream << ((*it).contigOri ? "+" : "-") << "\t";
        stream << (*it).numReads << "\t";
        stream << (*it).score << std::endl;
    }
    
    return 0;
}

// --------------------------------------------------------------------------
// addLocation()
// --------------------------------------------------------------------------

void
addLocation(Location & prevLoc, String<Location> & locations, Location & loc)
{
    if (prevLoc.contig == "")
    {
        prevLoc = loc;
    }
    else if (prevLoc.contig != loc.contig || prevLoc.chr != loc.chr || prevLoc.chrEnd + 1000 < loc.chrStart)
    {
        appendValue(locations, prevLoc);
        prevLoc = loc;
    }
    else
    {
        prevLoc.chrEnd = std::max(prevLoc.chrEnd, loc.chrEnd);
        prevLoc.numReads += loc.numReads;
    }
}

// --------------------------------------------------------------------------
// mergeLocationsBatch()
// --------------------------------------------------------------------------

int mergeLocationsBatch(String<Location> & locations, String<CharString> locationsFiles, size_t offset, size_t batchSize)
{
    typedef RecordReader<std::fstream, SinglePass<> > TReader;
    typedef Pair<std::fstream *, TReader *> TPtrPair;

    Location forwardForward, forwardReverse, reverseForward, reverseReverse;

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

        if (loc.chrOri)
        {
            if (loc.contigOri) addLocation(forwardForward, locations, loc);
            else addLocation(forwardReverse, locations, loc);
        }
        else
        {
            if (loc.contigOri) addLocation(reverseForward, locations, loc);
            else addLocation(reverseReverse, locations, loc);
        }

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

    // Add the remaining locations.
    if (forwardForward.contig != "") appendValue(locations, forwardForward);
    if (forwardReverse.contig != "") appendValue(locations, forwardReverse);
    if (reverseForward.contig != "") appendValue(locations, reverseForward);
    if (reverseReverse.contig != "") appendValue(locations, reverseReverse);

    return 0;
}

// ==========================================================================
// Function mergeLocations()
// ==========================================================================

int
mergeLocations(String<Location> & locations, String<CharString> & locationsFiles, CharString & outFile, bool verbose)
{
    String<CharString> tmpFiles;
    unsigned batchSize = 100;
    
    if (length(locationsFiles) > batchSize*batchSize)
    {
        std::cerr << "ERROR: Too many locations files, max: " << batchSize*batchSize << ", given: " << length(locationsFiles) << std::endl;
        return 1;
    }
    else if (length(locationsFiles) <= batchSize)
    {
        mergeLocationsBatch(locations, locationsFiles, 0, batchSize);
        return 0;
    }

    // Write a temporary file for each batch of files.
    for (unsigned offset = 0; offset < length(locationsFiles); offset += batchSize)
    {
        if (verbose) std::cerr << "[" << time(0) << "] " << "Merging batch " << (offset/batchSize)+1 << "of locations files." << std::endl;
        String<Location> locs;
        mergeLocationsBatch(locs, locationsFiles, offset, batchSize);

        if (verbose) std::cerr << "[" << time(0) << "] " << "Sorting batch " << (offset/batchSize)+1 << "." << std::endl;
        LocationTypeLess less;
        std::stable_sort(begin(locs, Standard()), end(locs, Standard()), less);

        std::stringstream tmpName;
        tmpName << outFile << "." << (offset/batchSize)+1;
        CharString tmp = tmpName.str();
        if (verbose) std::cerr << "[" << time(0) << "] " << "Writing batch " << (offset/batchSize)+1 << "of locations to temporary file " << tmp << "." << std::endl;
        writeLocations(tmp, locs);
        appendValue(tmpFiles, tmp);
        
    }

    // Merge and remove the temporary files.
    mergeLocationsBatch(locations, tmpFiles, 0, maxValue<unsigned>()-length(locationsFiles));
    for (unsigned i = 0; i < length(tmpFiles); ++i) remove(toCString(tmpFiles[i]));

    return 0;
}

#endif // #ifndef POPINS_LOCATION_H_
