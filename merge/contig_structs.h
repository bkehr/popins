#ifndef CONTIG_STRUCTS_H_
#define CONTIG_STRUCTS_H_

#include<seqan/sequence.h>

using namespace seqan;

// ============================================================================
// struct ContigId
// ============================================================================

struct ContigId
{
    CharString pn;
    CharString contigId;
    bool orientation;

    ContigId()
    {}

    ContigId(CharString & p, CharString & c, bool o) :
        pn(p), contigId(c), orientation(o)
    {}
};

// --------------------------------------------------------------------------

template<typename TLen>
CharString
formattedIndex(unsigned i, TLen max)
{
    unsigned digits = 0;
    while (max) {
        max /= 10;
        digits++;
    }

    std::stringstream s;
    s << std::setfill('0') << std::setw(digits) << i;

    return s.str();
}

// --------------------------------------------------------------------------
// Function operator<<                                               ContigId
// --------------------------------------------------------------------------

template<typename TStream>
inline TStream &
operator<<(TStream & stream, ContigId & id)
{
    stream << id.pn << "." << id.contigId;
    if (!id.orientation) stream << "_rev";
    return stream;
}

// --------------------------------------------------------------------------
// Function operator!=                                               ContigId
// --------------------------------------------------------------------------

inline bool
operator!=(ContigId const & left, ContigId const & right)
{
    return (left.pn != right.pn) ||
            (left.contigId != right.contigId) ||
            (left.orientation != right.orientation);
}

// --------------------------------------------------------------------------
// Function operator<                                                ContigId
// --------------------------------------------------------------------------

inline bool
operator<(ContigId const & left, ContigId const & right)
{
    return (left.pn < right.pn) ||
            (left.pn == right.pn && left.contigId < right.contigId) ||
            (left.pn == right.pn && left.contigId == right.contigId && left.orientation && !right.orientation);
}

// --------------------------------------------------------------------------
// Function operator>                                                ContigId
// --------------------------------------------------------------------------

inline bool
operator>(ContigId const & left, ContigId const & right)
{
    return (left.pn > right.pn) ||
            (left.pn == right.pn && left.contigId > right.contigId) ||
            (left.pn == right.pn && left.contigId == right.contigId && !left.orientation && right.orientation);
}

// ==========================================================================
// struct Contig
// ==========================================================================

template<typename TSeq>
struct Contig
{
    TSeq seq;
    ContigId id;

    Contig() {}

    Contig(TSeq & s, ContigId & i) :
        seq(s), id(i)
    {}
};

// ==========================================================================
// struct ContigComponent
// ==========================================================================

template<typename TSeq>
struct ContigComponent
{
    typedef typename Size<TSeq>::Type TSize;

    StringSet<Contig<TSeq>, Dependent<> > contigs;
    std::set<Pair<TSize> > alignedPairs;

    ContigComponent()
    {}
};

// --------------------------------------------------------------------------

template<typename TSeq>
void
clear(ContigComponent<TSeq> & c)
{
    clear(c.contigs);
    clear(c.alignedPairs);
}

// --------------------------------------------------------------------------

#endif  // #ifndef CONTIG_STRUCTS_H_
