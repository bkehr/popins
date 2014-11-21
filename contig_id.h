#ifndef CONTIG_ID_H_
#define CONTIG_ID_H_

// ============================================================================
// Classes
// ============================================================================

// --------------------------------------------------------------------------
// struct ContigId
// --------------------------------------------------------------------------

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

// ============================================================================
// Functions
// ============================================================================

template<typename TStream>
inline TStream &
operator<<(TStream & stream, ContigId & id)
{
    stream << id.pn << "." << id.contigId;
    if (!id.orientation) stream << "_rev";
    return stream;
}

// --------------------------------------------------------------------------
// Function operator!=                                                ContigId
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


#endif  // #ifndef CONTIG_ID_H_
