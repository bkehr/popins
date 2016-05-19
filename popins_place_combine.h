#ifndef POPINS_PLACE_COMBINE_H_
#define POPINS_PLACE_COMBINE_H_

#include <vector>
#include <seqan/seq_io.h>
#include "popins_clp.h"
#include "popins_location_info.h"

// ---------------------------------------------------------------------------------------
// Struct PlacedLocation
// ---------------------------------------------------------------------------------------

struct PlacedLocation
{
    typedef std::map<std::pair<unsigned, unsigned>, unsigned> TPosSupport;
    Location loc;
    TPosSupport insPos;

    PlacedLocation() {}

    PlacedLocation(Location & l, TPosSupport & i) :
        loc(l), insPos(i)
    {}
};

// ---------------------------------------------------------------------------------------
// Struct PlacedLocLess()
// ---------------------------------------------------------------------------------------

struct PlacedLocLess : public std::binary_function<PlacedLocation, PlacedLocation, bool>
{
    PlacedLocLess() {}

    inline bool operator() (PlacedLocation const & a, PlacedLocation const & b) const
    {
        LocationPosLess less;
        return less.compare(a.loc, b.loc) == 1;
    }
};

// ---------------------------------------------------------------------------------------
// Function readPlacedLocation()
// ---------------------------------------------------------------------------------------

bool
readPlacedLocation(PlacedLocation & loc, RecordReader<std::fstream, SinglePass<> > & reader, CharString & filename)
{
    CharString buffer;

    if (readUntilOneOf(loc.loc.chr, reader, ':', '\t', ' ') != 0)
    {
        std::cerr << "ERROR: Reading CHR from " << filename << " failed." << std::endl;
        return 1;
    }

    if (value(reader) == ':')
    {
        skipChar(reader, ':');

        if (readUntilChar(buffer, reader, '-') != 0)
        {
            std::cerr << "ERROR: Reading CHR_START from " << filename << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<Location::TPos>(loc.loc.chrStart, buffer);
        skipChar(reader, '-');

        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading CHR_END from " << filename << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<Location::TPos>(loc.loc.chrEnd, buffer);
    }
    skipWhitespaces(reader);

    if (value(reader) == '+') loc.loc.chrOri = true;
    else if (value(reader) == '-') loc.loc.chrOri = false;
    else
    {
        std::cerr << "ERROR: Reading CHR_ORI from " << filename << " failed." << std::endl;
        return 1;
    }
    skipNChars(reader, 1);
    skipWhitespaces(reader);

    if (readUntilWhitespace(loc.loc.contig, reader) != 0)
    {
        std::cerr << "ERROR: Reading CONTIG from " << filename << " failed." << std::endl;
        return 1;
    }
    skipWhitespaces(reader);

    if (value(reader) == '+') loc.loc.contigOri = true;
    else if (value(reader) == '-') loc.loc.contigOri = false;
    else
    {
        std::cerr << "ERROR: Reading CONTIG_ORI from " << filename << " failed." << std::endl;
        return 1;
    }
    skipNChars(reader, 1);
    skipWhitespaces(reader);

    clear(buffer);
    if (readDigits(buffer, reader) != 0)
    {
        std::cerr << "ERROR: Reading NUM_READS from " << filename << " failed." << std::endl;
        return 1;
    }
    lexicalCast2<unsigned>(loc.loc.numReads, buffer);
    skipWhitespaces(reader);

    clear(buffer);
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
    {
        std::cerr << "ERROR: Reading SCORE from " << filename << " failed." << std::endl;
        return 1;
    }
    lexicalCast2<double>(loc.loc.score, buffer);

    if (value(reader) == '\r' || value(reader) == '\n')
    {
        skipLine(reader);
        return 0;
    }

    skipWhitespaces(reader);

    if (value(reader) == 'h')
    {
        clear(buffer);
        readLine(buffer, reader);
        if (buffer == "high_coverage")
            return 0;

        std::cerr << "ERROR: Unexpected word \'" << buffer << "\' in " << filename << "." << std::endl;
        return 1;
    }

    unsigned refPos, contigPos, posSupport;

    clear(buffer);
    if (readDigits(buffer, reader) != 0)
    {
        std::cerr << "ERROR: Reading refPos from " << filename << " failed." << std::endl;
        return 1;
    }
    lexicalCast2<unsigned>(refPos, buffer);
    skipChar(reader, ',');

    clear(buffer);
    if (readDigits(buffer, reader) != 0)
    {
        std::cerr << "ERROR: Reading contigPos from " << filename << " failed." << std::endl;
        return 1;
    }
    lexicalCast2<unsigned>(contigPos, buffer);
    skipChar(reader, ':');

    clear(buffer);
    if (readDigits(buffer, reader) != 0)
    {
        std::cerr << "ERROR: Reading readSupport from " << filename << " failed." << std::endl;
        return 1;
    }
    lexicalCast2<unsigned>(posSupport, buffer);

    loc.insPos[std::pair<unsigned, unsigned>(refPos, contigPos)] = posSupport;

    while (value(reader) != '\r' && value(reader) != '\n') // while not at end of line
    {
        skipChar(reader, ';');

        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading refPos from " << filename << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<unsigned>(refPos, buffer);
        skipChar(reader, ',');

        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading contigPos from " << filename << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<unsigned>(contigPos, buffer);
        skipChar(reader, ':');

        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: Reading readSupport from " << filename << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<unsigned>(posSupport, buffer);

        loc.insPos[std::pair<unsigned, unsigned>(refPos, contigPos)] = posSupport;
    }

    skipLine(reader);
    return 0;
}

// ---------------------------------------------------------------------------------------
// Function loadPlacedLocations()
// ---------------------------------------------------------------------------------------

bool
loadPlacedLocations(std::vector<PlacedLocation> & locs, CharString & filename)
{
    // Open input file.
    std::fstream stream(toCString(filename), std::ios::in);
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open locations file " << filename << std::endl;
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > reader(stream);

    while (!atEnd(reader))
    {
        PlacedLocation loc;
        if (readPlacedLocation(loc, reader, filename) != 0)
            return 1;
        locs.push_back(loc);
    }

    return 0;
}

// ---------------------------------------------------------------------------------------
// Function combine()
// ---------------------------------------------------------------------------------------

PlacedLocation
combine(std::vector<PlacedLocation>::const_iterator & it,
        std::vector<PlacedLocation>::const_iterator & itEnd)
{
    typedef std::map<std::pair<unsigned, unsigned>, unsigned>::const_iterator TPosIter;

    PlacedLocation loc;
    loc.loc = (*it).loc;

    while (it != itEnd)
    {
        TPosIter posIt = (*it).insPos.begin();
        TPosIter posEnd = (*it).insPos.end();

        while (posIt != posEnd)
        {
            if (loc.insPos.count(posIt->first) == 0)
                loc.insPos[posIt->first] = posIt->second;
            else
                loc.insPos[posIt->first] += posIt->second;

            ++posIt;
        }

        ++it;
    }

    return loc;
}

// ---------------------------------------------------------------------------------------
// Function combineLocations()
// ---------------------------------------------------------------------------------------

void
combineLocations(std::vector<PlacedLocation>  & locs)
{
    typedef std::vector<PlacedLocation>::const_iterator TIter;

    std::vector<PlacedLocation> combined;

    TIter it = locs.begin();
    TIter it2 = locs.begin();
    TIter itEnd = locs.end();

    if (it != itEnd)
        ++it;

    PlacedLocation loc;

    PlacedLocLess less;
    while (it != itEnd)
    {
        if (less(*it2, *it))
        {
            loc = combine(it2, it);
            combined.push_back(loc);
        }
        ++it;
    }

    loc = combine(it2, it);
    combined.push_back(loc);

    locs = combined;
}

// ---------------------------------------------------------------------------------------
// Function chooseBestPlacing()
// ---------------------------------------------------------------------------------------

bool
chooseBestPlacing(unsigned & refPos, unsigned & contigPos, unsigned & support, PlacedLocation & loc)
{
    typedef PlacedLocation::TPosSupport::const_iterator TIter;

    if (loc.insPos.size() == 0)
        return 1;

    unsigned totalCount = 0;
    support = 0;

    TIter it = loc.insPos.begin();
    TIter itEnd = loc.insPos.end();

    while (it != itEnd)
    {
        unsigned cnt = it->second;
        totalCount += cnt;
        if (cnt > support)
        {
            support = cnt;
            refPos = (it->first).first;
            contigPos = (it->first).second;
        }
        ++it;
    }

    if (totalCount == 0 || support < 0.5 * totalCount)
        return 1;

    return 0;
}

// ---------------------------------------------------------------------------------------
// Function writeVcf()
// ---------------------------------------------------------------------------------------

template<typename TStream>
void
writeVcf(TStream & outStream, PlacedLocation & loc, unsigned refPos, unsigned contigPos, unsigned support, FaiIndex & fai)
{
    Dna5String ref = loadInterval(fai, loc.loc.chr, refPos, refPos + 1);

    outStream << loc.loc.chr;
    outStream << "\t" << refPos + 1;
    outStream << "\t" << loc.loc.chr << ":" << refPos + 1 << ":" << "FP";
    outStream << "\t" << ref;

    if (loc.loc.chrOri)
        outStream << "\t" << ref << "[" << loc.loc.contig << (!loc.loc.contigOri?"f":"r") << ":" << contigPos << "[";
    else
        outStream << "\t" << "]" << loc.loc.contig << (loc.loc.contigOri?"f":"r") << ":" << contigPos << "]" << ref;

    outStream << "\t" << ".";
    outStream << "\t" << ".";

    if (loc.loc.numReads != 0)
        outStream << "\t" << "AR=" << loc.loc.numReads << ";AS=" << loc.loc.score;
    else
        outStream << "\t" << "PAIRED";
    outStream << ";" << "SR=" << support;    // TODO Write more info fields.

    outStream << "\t" << ".";
    outStream << std::endl;
}

// ---------------------------------------------------------------------------------------

template<typename TStream>
void
writeVcf(TStream & outStream, PlacedLocation & loc, FaiIndex & fai)
{
    unsigned refPos;
    if (loc.loc.chrOri)
        refPos = loc.loc.chrEnd;
    else
        refPos = loc.loc.chrStart;

    Dna5String ref = loadInterval(fai, loc.loc.chr, refPos, refPos + 1);

    outStream << loc.loc.chr;
    outStream << "\t" << refPos + 1;
    outStream << "\t" << loc.loc.chr << ":" << refPos + 1 << ":" << "FP";
    outStream << "\t" << ref;

    if (loc.loc.chrOri)
        outStream << "\t" << ref << "[" << loc.loc.contig << (!loc.loc.contigOri?"f":"r") << "[";
    else
        outStream << "\t" << "]" << loc.loc.contig << (loc.loc.contigOri?"f":"r") << "]" << ref;

    outStream << "\t" << ".";
    outStream << "\t" << ".";

    if (loc.loc.numReads != 0)
        outStream << "\t" << "AR=" << loc.loc.numReads << ";AS=" << loc.loc.score;
    else
        outStream << "\t" << "NOANCHOR";      // TODO Write more info fields.

    outStream << "\t" << ".";
    outStream << std::endl;
}

// =======================================================================================
// Function popins_place_combine()
// =======================================================================================

template<typename TStream>
bool
popins_place_combine(TStream & vcfStream, PlacingOptions & options)
{
    std::vector<PlacedLocation> locs;

    // Open the FAI file of the reference genome.
    FaiIndex fai;
    if (read(fai, toCString(options.referenceFile)) != 0)
    {
        std::cerr << "ERROR: Could not open FAI index for " << options.referenceFile << std::endl;
        return 1;
    }

    if (options.verbose)
    {
        std::cerr << "[" << time(0) << "] " << "Loading the placed locations from " << length(options.locationsFiles) << " locations files." << std::endl;
        std::cerr << "0%   10   20   30   40   50   60   70   80   90   100%" << std::endl;
        std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;
        std::cerr << "*" << std::flush;
    }
    double fiftieth = length(options.locationsFiles) / 50.0;
    unsigned progress = 0;

    // Read placed location files.
    for (unsigned i = 0; i < length(options.locationsFiles); ++i)
    {
        if (loadPlacedLocations(locs, options.locationsFiles[i]) != 0)
            return 1;

        while (options.verbose && progress * fiftieth < i)
        {
            std::cerr << "*" << std::flush;
            ++progress;
        }
    }

    if (options.verbose)
        std::cerr << std::endl;

    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Sorting " << locs.size() << " placed locations." << std::endl;

    // Sort placed locations.
    std::stable_sort(locs.begin(), locs.end(), PlacedLocLess());

    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Combining the placed locations." << std::endl;

    // Combine placing of the same locations.
    combineLocations(locs);

    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Writing " << locs.size() << " combined locations to output file '" << options.outFile << "'." << std::endl;

    // Choose the best position.
    for (unsigned i = 0; i < locs.size(); ++i)
    {
        unsigned refPos = 0, contigPos = 0, support = 0;
        if (chooseBestPlacing(refPos, contigPos, support, locs[i]) == 0)
            writeVcf(vcfStream, locs[i], refPos, contigPos, support, fai);
        else
            writeVcf(vcfStream, locs[i], fai);
    }

    return 0;
}

#endif /* POPINS_PLACE_COMBINE_H_ */
