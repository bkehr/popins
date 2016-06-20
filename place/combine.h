#ifndef POPINS_PLACE_COMBINE_H_
#define POPINS_PLACE_COMBINE_H_

#include <vector>
#include <seqan/seq_io.h>
#include "../command_line_parsing.h"
#include "location_info.h"

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
readPlacedLocation(PlacedLocation & loc, std::stringstream & stream, CharString & filename)
{
    std::string genomicPos;
    stream >> genomicPos;

    size_t colon;
    if ((colon = genomicPos.find(':')) != std::string::npos)
    {
        loc.loc.chr = genomicPos.substr(0, colon);

        size_t dash = genomicPos.find('-');

        std::string start = genomicPos.substr(colon+1, dash-colon-1);
        if (!lexicalCast<Location::TPos>(loc.loc.chrStart, start))
        {
            std::cerr << "ERROR: Could not parse " << start << " as location start position in \'" << filename << "\'." << std::endl;
            return 1;
        }

        std::string end = genomicPos.substr(dash+1);
        if (!lexicalCast<Location::TPos>(loc.loc.chrEnd, end))
        {
            std::cerr << "ERROR: Could not parse " << end << " as location end position in \'" << filename << "\'." << std::endl;
            return 1;
        }
    }
    else
    {
        loc.loc.chr = genomicPos;
    }

    char ori;
    stream >> ori;

    if (ori == '+') loc.loc.chrOri = true;
    else if (ori == '-') loc.loc.chrOri = false;
    else
    {
        std::cerr << "ERROR: Reading CHR_ORI from locations file " << filename << " failed." << std::endl;
        return 1;
    }

    std::string buffer;
    stream >> buffer;
    loc.loc.contig = buffer;

    stream >> ori;
    if (ori == '+') loc.loc.contigOri = true;
    else if (ori == '-') loc.loc.contigOri = false;
    else
    {
        std::cerr << "ERROR: Reading CONTIG_ORI from locations file " << filename << " failed." << std::endl;
        return 1;
    }

    stream >> loc.loc.numReads;
    stream >> loc.loc.score;

    if (stream.eof())
        return 0;

    stream >> std::ws;

    if (stream.peek() == 'h')
    {
        std::string buffer;
        std::getline(stream, buffer);
        if (buffer.compare("high_coverage") == 0)
            return 0;

        std::cerr << "ERROR: Unexpected word \'" << buffer << "\' in " << filename << "." << std::endl;
        return 1;
    }

    std::string refPosStr, contigPosStr, posSupportStr;
    while (std::getline(stream, refPosStr, ',') && std::getline(stream, contigPosStr, ':') && std::getline(stream, posSupportStr, ';'))
    {
        unsigned refPos = 0, contigPos = 0, posSupport = 0;
        lexicalCast<unsigned>(refPos, refPosStr);
        lexicalCast<unsigned>(contigPos, contigPosStr);
        lexicalCast<unsigned>(posSupport, posSupportStr);

        loc.insPos[std::pair<unsigned, unsigned>(refPos, contigPos)] = posSupport;
    }

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

    std::string line;
    while (std::getline(stream, line))
    {
        std::stringstream ss;
        ss.str(line);

        PlacedLocation loc;
        if (readPlacedLocation(loc, ss, filename) != 0)
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
        outStream << "\t" << "NOANCHOR";
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
popins_place_combine(TStream & vcfStream, CharString & prefix, CharString & referenceFile, CharString & outFile)
{
    std::vector<PlacedLocation> locs;

    // Open the FAI file of the reference genome.
    FaiIndex fai;
    if (!open(fai, toCString(referenceFile)))
    {
        std::cerr << "ERROR: Could not open FAI index for " << referenceFile << std::endl;
        return 1;
    }

    CharString filename = "locations_placed.txt";
    String<Pair<CharString> > locationsFiles = listFiles(prefix, filename);

    std::ostringstream msg;
    msg << "Loading the placed locations from " << length(locationsFiles) << " locations files.";
    printStatus(msg);

    std::cerr << "0%   10   20   30   40   50   60   70   80   90   100%" << std::endl;
    std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;
    std::cerr << "*" << std::flush;

    double fiftieth = length(locationsFiles) / 50.0;
    unsigned progress = 0;

    // Read placed location files.
    for (unsigned i = 0; i < length(locationsFiles); ++i)
    {
        if (loadPlacedLocations(locs, locationsFiles[i].i2) != 0)
            return 1;

        while (progress * fiftieth < i)
        {
            std::cerr << "*" << std::flush;
            ++progress;
        }
    }
    std::cerr << std::endl;

    msg.str("");
    msg << "Sorting " << locs.size() << " placed locations.";
    printStatus(msg);

    // Sort placed locations.
    std::stable_sort(locs.begin(), locs.end(), PlacedLocLess());

    msg.str("");
    msg << "Combining the placed locations.";
    printStatus(msg);

    // Combine placing of the same locations.
    combineLocations(locs);

    msg.str("");
    msg << "Writing " << locs.size() << " combined locations to output file '" << outFile << "'.";
    printStatus(msg);

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
