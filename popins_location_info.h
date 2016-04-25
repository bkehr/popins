#ifndef POPINS_LOCATION_INFO_H_
#define POPINS_LOCATION_INFO_H_

#include "popins_location.h"

using namespace seqan;

// ==========================================================================
// struct LocationInfo
// ==========================================================================

struct LocationInfo
{
    Location loc;
    bool otherEnd;
    unsigned idx;
    unsigned contigLength;
    int insPos;
    unsigned refPos;

    LocationInfo (Location & l, unsigned i, unsigned c) :
        loc(l), otherEnd(false), idx(i), contigLength(c), insPos(0), refPos(0)
    {}
};

// ==========================================================================
// Struct LocationInfoPosLess
// ==========================================================================

struct LocationInfoPosLess : public std::binary_function<LocationInfo, LocationInfo, bool>
{
    LocationInfoPosLess() {}

    inline bool operator() (LocationInfo const & a, LocationInfo const & b) const
    {
        LocationPosLess less;
        return less.compare(a.loc, b.loc) == 1;
    }
};

// ==========================================================================
// Struct LocationInfoTypeLess
// ==========================================================================

struct LocationInfoTypeLess : public std::binary_function<LocationInfo, LocationInfo, bool>
{
    LocationInfoTypeLess() {}

    inline int compare(LocationInfo const & a, LocationInfo const & b) const
    {
        LocationTypeLess less;
        return less.compare(a.loc, b.loc);
    }

    inline bool operator() (LocationInfo const & a, LocationInfo const & b) const
    {
        return compare(a, b) == 1;
    }
};

// ==========================================================================
// Struct LocationInfoGreater
// ==========================================================================

struct LocationInfoGreater : public std::binary_function <LocationInfo, LocationInfo, bool>
{
    LocationInfoGreater() {}

    inline int compare(LocationInfo const & a, LocationInfo const & b) const
    {
        if (a.otherEnd && !b.otherEnd) return -1;
        if (!a.otherEnd && b.otherEnd) return 1;

        if (a.contigLength > b.contigLength) return -1;
        if (a.contigLength < b.contigLength) return 1;

        if (a.loc.score > b.loc.score) return -1;
        if (a.loc.score < b.loc.score) return 1;

        if (a.loc.numReads > b.loc.numReads) return -1;
        if (a.loc.numReads < b.loc.numReads) return 1;

        LocationTypeGreater greater;
        return greater.compare(a.loc, b.loc);
    }

    inline bool operator() (LocationInfo const & a, LocationInfo const & b) const
    {
        return compare(a, b) != 1;
    }
};

// ==========================================================================
// Function appendLocation()
// ==========================================================================

void
appendLocation(String<LocationInfo> & locs, Location & loc)
{
    appendValue(locs, LocationInfo(loc, length(locs), 0));
}

// ==========================================================================
// Function loadInterval()
// ==========================================================================

Dna5String
loadInterval(FaiIndex & fai, CharString & chrom, unsigned beginPos, unsigned endPos)
{
    unsigned idx = 0;
    if (!getIdByName(fai, chrom, idx))
    {
        std::cerr << "ERROR: Could not find " << chrom << " in FAI index." << std::endl;
        return 1;
    }

    Dna5String refInfix;
    readRegion(refInfix, fai, idx, beginPos, endPos);

    return refInfix;
}

// ==========================================================================
// Function writeVcf()
// ==========================================================================

template<typename TStream>
void
writeVcf(TStream & outStream, LocationInfo & loc, FaiIndex & fai)
{
    Dna5String ref = loadInterval(fai, loc.loc.chr, loc.refPos, loc.refPos + 1);

    outStream << loc.loc.chr;
    outStream << "\t" << loc.refPos + 1;
    outStream << "\t" << loc.loc.chr << ":" << loc.refPos + 1 << ":" << "FP";
    outStream << "\t" << ref;

    if (loc.loc.chrOri)
        outStream << "\t" << ref << "[" << loc.loc.contig << (!loc.loc.contigOri?"f":"r") << ":" << loc.insPos << "[";
    else
        outStream << "\t" << "]" << loc.loc.contig << (loc.loc.contigOri?"f":"r") << ":" << loc.insPos << "]" << ref;

    outStream << "\t" << ".";
    outStream << "\t" << ".";
    if (loc.loc.numReads != 0)
        outStream << "\t" << "AR=" << loc.loc.numReads << ";AS=" << loc.loc.score; // TODO Write more info fields.
    else
        outStream << "\t" << "PAIRED";
    outStream << std::endl;
}

#endif /* POPINS_LOCATION_INFO_H_ */
