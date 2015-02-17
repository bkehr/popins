#ifndef POPINS_MERGE_PARTITION_H_
#define POPINS_MERGE_PARTITION_H_

#include <seqan/index.h>
#include <seqan/align.h>

#include "contig_id.h"
#include "contig_component.h"

using namespace seqan;

// ==========================================================================
// Function readAndMergeComponents()
// ==========================================================================

template<typename TSequence, typename TSpec>
bool
readAndMergeComponents(std::map<ContigId, ContigComponent<TSequence> > & components,
                       String<CharString> & componentFiles,
                       StringSet<ContigId, TSpec> & contigIds,
                       bool verbose)
{
    typedef std::map<ContigId, ContigComponent<TSequence> > TComponents;
    typedef typename ContigComponent<TSequence>::TSize TSize;
    typedef typename std::map<TSize, std::set<TSize> >::iterator TPairsIterator;

    std::cerr << "[" << time(0) << "] " << "Reading and merging components files" << std::endl;

    // Initialize Union-Find data structure.
    UnionFind<int> uf;
    resize(uf, length(contigIds));
    
    // Read the component files and join sets.
    String <ContigComponent<TSequence> > comps;
    for (unsigned i = 0; i < length(componentFiles); ++i)
        if (readComponents(comps, uf, componentFiles[i], length(contigIds)/2, verbose) != 0) return 1;

    // Merge elements of comps that belong to the same set by adding them to the same element in components.
    typename Iterator<String<ContigComponent<TSequence> > >::Type compsIt = begin(comps);
    typename Iterator<String<ContigComponent<TSequence> > >::Type compsEnd = end(comps);
    for (; compsIt != compsEnd; ++compsIt)
    {
        int set_i = findSet(uf, (*compsIt).alignedPairs.begin()->first);
        if (set_i >= (int)length(contigIds)/2) continue;

        ContigId id = contigIds[set_i];

        TPairsIterator pairsIt = (*compsIt).alignedPairs.begin();
        TPairsIterator pairsEnd = (*compsIt).alignedPairs.end();
        for (; pairsIt != pairsEnd; ++pairsIt)
            components[id].alignedPairs[pairsIt->first].insert(pairsIt->second.begin(), pairsIt->second.end());
    }

    // Add singleton contigs to components (= those contigs that don't align to any other contig).
    unsigned numSingletons = 0;
    for (unsigned i = 0; i < length(contigIds)/2; ++i)
    {
        if (components.count(contigIds[i]) == 0 && (int)i == findSet(uf, i))
        {
            components[contigIds[i]].alignedPairs[i];
            ++numSingletons;
        }
    }

    std::cerr << "[" << time(0) << "] " << "Added " << numSingletons << " singletons." << std::endl;
    std::cerr << "[" << time(0) << "] " << "There are " << components.size() << " components after merging." << std::endl;

    return 0;
}

// --------------------------------------------------------------------------
// Function pairwiseAlignment()
// --------------------------------------------------------------------------

template<typename TSeq, typename TValueScore>
inline bool
pairwiseAlignment(TSeq & contig1,
                  TSeq & contig2,
                  Score<int, Simple> scoringScheme,
                  unsigned lowerDiag,
                  unsigned upperDiag,
                  TValueScore minScore)
{
    typedef typename Position<TSeq>::Type TPos;

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

template<typename TSize, typename TSeq, typename TFloat, typename TLen, typename TLength, typename TValueMatch, typename TValueError, typename TValueScore>
void
partitionContigs(UnionFind<int> & uf,
                 std::map<TSize, std::set<TSize> > & alignedPairs,
                 StringSet<TSeq> & contigs,
                 StringSet<ContigId> & contigIds,
                 StringSet<TSeq> & contigSubset,
                 StringSet<ContigId> & contigIdSubset,
                 int offset,
                 unsigned len,
                 TFloat errorRate,
                 TLen minimalLength,
                 TLength qgramLength,
                 TValueMatch matchScore,
                 TValueError errorPenalty,
                 TValueScore minScore,
                 bool verbose)
{
    typedef StringSet<TSeq> TStringSet;

    typedef Index<TStringSet, IndexQGram<SimpleShape, OpenAddressing> > TIndex;
    typedef Finder<TSeq, Swift<SwiftLocal> > TFinder;

    TSize numComparisons = 0;
    TSize numAligns = 0;

    // --- Initialization of SWIFT pattern and finders
    TIndex qgramIndex(contigSubset);
    resize(indexShape(qgramIndex), qgramLength);
    Pattern<TIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);
    indexRequire(qgramIndex, QGramSADir());

    // define scoring scheme
    Score<int, Simple> scoringScheme(matchScore, errorPenalty, errorPenalty);
    int diagExtension = minScore/10;
    
    if (verbose) std::cerr << "0%   10   20   30   40   50   60   70   80   90  100%" << std::endl;
    if (verbose) std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;
    unsigned fiftieth = std::max(len/50, 1u);

    for (unsigned a = 0; a < len; ++a)
    {
        if (verbose && a%fiftieth == 0) std::cerr << "*" << std::flush;
        
        TFinder swiftFinder(contigs[a], 1000, 1);

        // --- SWIFT search

        while (find(swiftFinder, swiftPattern, errorRate, minimalLength)) {
            // get index of pattern sequence (subsetB is index in contigSubset and swiftIndex, b is index in uf)
            unsigned subsetB = swiftPattern.curSeqNo;
            unsigned b = subsetB + offset;
            if (subsetB > length(contigSubset)/2) b += len - length(contigSubset)/2;
            SEQAN_ASSERT_EQ(contigIdSubset[subsetB].pn, contigIds[b].pn);
            SEQAN_ASSERT_EQ(contigIdSubset[subsetB].contigId, contigIds[b].contigId);
            SEQAN_ASSERT_EQ(contigIdSubset[subsetB].orientation, contigIds[b].orientation);

            // align contigs only of different individuals and only if not same component already
            if (contigIds[a].pn == contigIdSubset[subsetB].pn) continue;
            if (findSet(uf, a) == findSet(uf, b)) continue;

            // verify by banded Smith-Waterman alignment
            ++numComparisons;
            unsigned upperDiag = (*swiftFinder.curHit).hstkPos - (*swiftFinder.curHit).ndlPos + diagExtension;
            unsigned lowerDiag = upperDiag - swiftPattern.bucketParams[subsetB].delta - swiftPattern.bucketParams[subsetB].overlap - diagExtension;
            if (!pairwiseAlignment(contigs[a], contigSubset[subsetB], scoringScheme, lowerDiag, upperDiag, minScore)) continue;

            ++numAligns;

            // join sets of the two aligned contigs
            joinSets(uf, findSet(uf, a), findSet(uf, b));
            alignedPairs[a].insert(b);
            alignedPairs[b].insert(a);

            // join sets for reverse complements of the contigs
            unsigned a1 = a;
            if (a1 < len) a1 += len; // a was index of contig in forward orientation, set it to index of contig in reverse orientation
            else a1 -= len;         // a was index of contig in reverse orientation, set it to index of contig in forward orientation
            if (b < len) b += len; // b was index of contig in forward orientation, set it to index of contig in reverse orientation
            else b -= len;         // b was index of contig in reverse orientation, set it to index of contig in forward orientation
            joinSets(uf, findSet(uf, a1), findSet(uf, b));
            alignedPairs[a1].insert(b);
            alignedPairs[b].insert(a1);
        }
    }
    if (verbose) std::cerr << std::endl;

    std::cerr << "[" << time(0) << "] " << "Number of pairwise comparisons: \t" << numComparisons << std::endl;
    std::cerr << "[" << time(0) << "] " << "Number of valid alignments:     \t" << numAligns << std::endl;
}

// --------------------------------------------------------------------------

template<typename TSeq, typename TFloat, typename TLen, typename TLength, typename TValueMatch, typename TValueError, typename TValueScore>
void
partitionContigs(std::map<ContigId, ContigComponent<TSeq> > & components,
                 StringSet<TSeq> & contigs,
                 StringSet<ContigId> & contigIds,
                 TFloat errorRate,
                 TLen minimalLength,
                 TLength qgramLength,
                 TValueMatch matchScore,
                 TValueError errorPenalty,
                 TValueScore minScore,
                 bool verbose)
{
    typedef typename Size<TSeq>::Type TSize;

    // Initialize Union-Find data structure.
    UnionFind<int> uf;
    resize(uf, length(contigs));
    
    std::map<TSize, std::set<TSize> > alignedPairs;

    // Partition the contigs.
    partitionContigs(uf, alignedPairs, contigs, contigIds, contigs, contigIds, 0, length(contigs)/2,
                     errorRate, minimalLength, qgramLength, matchScore, errorPenalty, minScore, verbose);
    
    // Determine components from Union-Find data structure by
    // mapping ids and contigs to their representative ContigId.
    TSize len = length(contigs)/2;
    for (TSize i = 0; i < len; ++i)
    {
        int set_i = findSet(uf, i);
        if (set_i <= findSet(uf, i + len))
        {
            ContigId id = contigIds[set_i];
            components[id].alignedPairs[i].insert(alignedPairs[i].begin(), alignedPairs[i].end());
        }
    }

    for (TSize i = len; i < length(contigs); ++i)
    {
        int set_i = findSet(uf, i);
        if (set_i <= findSet(uf, i - len))
        {
            ContigId id = contigIds[set_i];
            components[id].alignedPairs[i].insert(alignedPairs[i].begin(), alignedPairs[i].end());
        }
    }

    std::cerr << "[" << time(0) << "] " << "There are " << components.size() << " components in total." << std::endl;
}

// --------------------------------------------------------------------------

template<typename TSeq, typename TFloat, typename TLen, typename TLength, typename TValueMatch, typename TValueError, typename TValueScore>
void
partitionContigs(std::map<ContigId, ContigComponent<TSeq> > & components,
                 StringSet<TSeq> & contigs,
                 StringSet<ContigId> & contigIds,
                 StringSet<TSeq> & contigSubset,
                 StringSet<ContigId> & contigIdSubset,
                 int offset,
                 TFloat errorRate,
                 TLen minimalLength,
                 TLength qgramLength,
                 TValueMatch matchScore,
                 TValueError errorPenalty,
                 TValueScore minScore,
                 bool verbose)
{
    typedef typename Size<TSeq>::Type TSize;

    // Initialize Union-Find data structure.
    UnionFind<int> uf;
    resize(uf, 2*length(contigs));

    std::map<TSize, std::set<TSize> > alignedPairs;

    // Partition the contigs.
    partitionContigs(uf, alignedPairs, contigs, contigIds, contigSubset, contigIdSubset, offset, length(contigs),
                     errorRate, minimalLength, qgramLength, matchScore, errorPenalty, minScore, verbose);
    
    // Determine components from Union-Find data structure by
    // mapping ids and contigs to their representative ContigId.
    for (TSize i = 0; i < 2*length(contigs); ++i)
    {
        int set_i = findSet(uf, i);

        ContigId id;
        if (set_i < (int)length(contigs))
        {
            id = contigIds[set_i];
        }
        else
        {
            set_i -= length(contigs);
            id = contigIds[set_i];
            id.orientation = false;
        }

        components[id].alignedPairs[i].insert(alignedPairs[i].begin(), alignedPairs[i].end());
    }

    std::cerr << "[" << time(0) << "] " << "There are " << components.size() << " components in total." << std::endl;
}



#endif // #ifndef POPINS_MERGE_PARTITION_H_
