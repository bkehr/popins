#ifndef POPINS_MERGE_PARTITION_H_
#define POPINS_MERGE_PARTITION_H_

#include <seqan/index.h>
#include <seqan/align.h>

#include "contig_structs.h"

using namespace seqan;

// --------------------------------------------------------------------------
// Function pairwiseAlignment()
// --------------------------------------------------------------------------

template<typename TSeq, typename TValueScore>
inline bool
pairwiseAlignment(TSeq & contig1,
        TSeq & contig2,
        Score<int, Simple> scoringScheme,
        int lowerDiag,
        int upperDiag,
        TValueScore minScore)
{
    // setup alignment object
    Align<TSeq, ArrayGaps> align;
    resize(rows(align), 2);
    setSource(row(align, 0), contig1);
    setSource(row(align, 1), contig2);

    // compute local alignment
    int score = localAlignment(align, scoringScheme, lowerDiag, upperDiag);

    // return true if minimal score is reached
    if (score > minScore)
        return true;
    else
        return false;
}

// ==========================================================================
// Function partitionContigs()
// ==========================================================================

template<typename TSize, typename TSeq>
bool
partitionContigs(UnionFind<int> & uf,
        std::set<Pair<TSize> > & alignedPairs,
        String<Contig<TSeq> > & contigs,
        MergingOptions & options)
{
    typedef typename Iterator<String<Contig<TSeq> > >::Type TContigIter;
    typedef StringSet<TSeq, Dependent<> > TStringSet;
    typedef Index<TStringSet, IndexQGram<SimpleShape, OpenAddressing> > TIndex;
    typedef Finder<TSeq, Swift<SwiftLocal> > TFinder;

    printStatus("Partitioning contigs");
    printStatus("- Indexing contigs");

    TSize numComparisons = 0;
    unsigned fwdContigCount = length(contigs)/2;

    // initialization of SWIFT pattern (q-gram index)
    TStringSet seqs;
    StringSet<TSize> indices;
    TContigIter itEnd = end(contigs);
    for (TContigIter it = begin(contigs); it != itEnd; ++it)
        appendValue(seqs, (*it).seq);
    TIndex qgramIndex(seqs);
    resize(indexShape(qgramIndex), options.qgramLength);
    Pattern<TIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);
    indexRequire(qgramIndex, QGramSADir());

    // define scoring scheme
    Score<int, Simple> scoringScheme(options.matchScore, options.errorPenalty, options.errorPenalty);
    int diagExtension = options.minScore/10;

    // print status bar
    printStatus("- Aligning contigs");
    std::cerr << "0%   10   20   30   40   50   60   70   80   90   100%" << std::endl;
    std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;

    unsigned fiftieth = std::max(fwdContigCount / 50, 1u);

    // Iterate over the forward contigs.
    for (int a = 0; a < fwdContigCount; ++a)
    {
        if (a%fiftieth == 0)
            std::cerr << "*" << std::flush;

        // initialization of swift finder
        TFinder swiftFinder(contigs[a].seq, 1000, 1);

        hash(swiftPattern.data_host.data_value->shape, hostIterator(hostIterator(swiftFinder)));
        while (find(swiftFinder, swiftPattern, options.errorRate, options.minimalLength))
        {
            // get index of pattern sequence
            unsigned b = swiftPattern.curSeqNo;

            // align contigs only of different individuals
            if (contigs[a].id.pn == contigs[b].id.pn) continue;

            // align contigs only if not same component already
            if (findSet(uf, a) == findSet(uf, b)) continue;

            // find the contig sequences
            TSeq contigA = haystack(swiftFinder);
            TSeq contigB = indexText(needle(swiftPattern))[b];

            // compute upper and lower diagonal of band.
            int upperDiag = (*swiftFinder.curHit).hstkPos - (*swiftFinder.curHit).ndlPos;
            int lowerDiag = upperDiag - swiftPattern.bucketParams[b].delta - swiftPattern.bucketParams[b].overlap;
            upperDiag += diagExtension;
            lowerDiag -= diagExtension;

            // verify by banded Smith-Waterman alignment
            ++numComparisons;
            if (!pairwiseAlignment(contigA, contigB, scoringScheme, lowerDiag, upperDiag, options.minScore)) continue;
            alignedPairs.insert(Pair<TSize>(a, b));

            // join sets of the two aligned contigs
            joinSets(uf, findSet(uf, a), findSet(uf, b));

            // join sets for reverse complements of the contigs
            unsigned a1 = a < fwdContigCount ? a + fwdContigCount : a - fwdContigCount;
            unsigned b1 = b < fwdContigCount ? b + fwdContigCount : b - fwdContigCount;
            joinSets(uf, findSet(uf, a1), findSet(uf, b1));

            // stop aligning this contig if it is already in a component with more than 100 other contigs
            if (uf._values[findSet(uf, a)] < -100) break;
        }
    }
    std::cerr << std::endl;

    std::ostringstream msg;
    msg << "Number of pairwise comparisons: " << numComparisons;
    printStatus(msg);

    msg.str("");
    msg << "Number of valid alignments:     " << length(alignedPairs);
    printStatus(msg);

    return 0;
}

// --------------------------------------------------------------------------
// Function unionFindToComponents()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
void
unionFindToComponents(std::map<TSize, ContigComponent<TSeq> > & components,
        UnionFind<int> & uf,
        std::set<Pair<TSize> > & alignedPairs,
      unsigned fwdContigCount)
{
    std::ostringstream msg;

    // Determine components from Union-Find data structure by mapping ids to their representative id.
    for (typename std::set<Pair<TSize> >::iterator it = alignedPairs.begin(); it != alignedPairs.end(); ++it)
    {
        TSize rev1 = (*it).i1 < fwdContigCount ? (*it).i1 + fwdContigCount : (*it).i1 - fwdContigCount;
        TSize rev2 = (*it).i2 < fwdContigCount ? (*it).i2 + fwdContigCount : (*it).i2 - fwdContigCount;

        int set = std::min(findSet(uf, (*it).i1), findSet(uf, rev1));

        components[set].alignedPairs.insert(*it);
        components[set].alignedPairs.insert(Pair<TSize>((*it).i2, (*it).i1));
        components[set].alignedPairs.insert(Pair<TSize>(rev1, rev2));
        components[set].alignedPairs.insert(Pair<TSize>(rev2, rev1));
    }

    msg.str("");
    msg << "There are " << components.size() << " components.";
    printStatus(msg);
}

// --------------------------------------------------------------------------
// Function addSingletons()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
void
addSingletons(std::map<TSize, ContigComponent<TSeq> > & components,
        String<Contig<TSeq> > & contigs,
        UnionFind<int> & uf)
{
    unsigned numSingletons = 0;
    for (int i = 0; i < length(contigs)/2; ++i)
    {
        if (components.count(i) == 0 && i == findSet(uf, i))
        {
            components[i];
            ++numSingletons;
        }
    }

    std::ostringstream msg;
    msg << "Added " << numSingletons << " singletons to components.";
    printStatus(msg);
}

#endif // #ifndef POPINS_MERGE_PARTITION_H_
