#ifndef CONTIG_COMPONENT_H_
#define CONTIG_COMPONENT_H_

#include<seqan/sequence.h>

#include "contig_id.h"

using namespace seqan;

// --------------------------------------------------------------------------
// struct ContigComponent
// --------------------------------------------------------------------------

template<typename TSeq>
struct ContigComponent
{
    typedef typename Size<TSeq>::Type TSize;

    StringSet<ContigId, Dependent<> > ids;
    StringSet<TSeq, Dependent<> > contigs;
    
    std::set<Pair<TSize> > alignedPairs;

    ContigComponent()
    {}
};

template<typename TSeq>
void
clear(ContigComponent<TSeq> & c)
{
    clear(c.ids);
    clear(c.contigs);
    clear(c.alignedPairs);
}

#endif  // #ifndef CONTIG_COMPONENT_H_
