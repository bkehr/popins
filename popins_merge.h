#ifndef POPINS_MERGE_H_
#define POPINS_MERGE_H_

#include <sstream>
#include <iomanip>

#include <seqan/index.h>
#include <seqan/align.h>

#include "popins_clp.h"
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
    
    std::map<TSize, std::set<TSize> > alignedPairs;

    ContigComponent()
    {}
};

// --------------------------------------------------------------------------
// struct Path
// --------------------------------------------------------------------------

template<typename TSeq, typename TVertexDescriptor>
struct Path
{
    typedef typename Position<TSeq>::Type TPos;

    TSeq seq;
    std::map<TPos, TVertexDescriptor> positionMap;
    
    Path()
    {}
    
    Path(Path<TSeq, TVertexDescriptor> & other) :
        seq(other.seq), positionMap(other.positionMap)
    {}
    
    Path(Path<TSeq, TVertexDescriptor> const & other) :
        seq(other.seq), positionMap(other.positionMap)
    {}
};

// --------------------------------------------------------------------------
// struct ComponentGraph
// --------------------------------------------------------------------------

template<typename TSeq>
struct ComponentGraph
{
    typedef Graph<Directed<> > TGraph_;
    typedef typename VertexDescriptor<TGraph_>::Type TVertexDescriptor;

    TGraph_ graph;
    String<TVertexDescriptor> sources;
    String<TSeq> sequenceMap;

    ComponentGraph()
    {}

    ComponentGraph(TSeq & seq)
    {
        TVertexDescriptor v = addVertex(*this, seq);
        appendValue(sources, v);
    }
};

// --------------------------------------------------------------------------

template<typename TSeq1, typename TSeq2>
typename ComponentGraph<TSeq1>::TVertexDescriptor
addVertex(ComponentGraph<TSeq1> & graph, TSeq2 & seq)
{
    typename ComponentGraph<TSeq1>::TVertexDescriptor v = addVertex(graph.graph);
    if (v != length(graph.sequenceMap)) std::cerr << "WARNING: The vertex map is broken.";
    appendValue(graph.sequenceMap, seq);
    return v;
}

// --------------------------------------------------------------------------

template<typename TSeq, typename TSeq1, typename TSeq2>
typename ComponentGraph<TSeq>::TVertexDescriptor
splitVertex(ComponentGraph<TSeq> & graph,
            typename ComponentGraph<TSeq>::TVertexDescriptor & u,
            TSeq1 & uSeq,
            TSeq2 & vSeq)
{
    typedef typename Iterator<typename ComponentGraph<TSeq>::TGraph_, OutEdgeIterator>::Type TOutEdgeIter;
    typedef typename ComponentGraph<TSeq>::TVertexDescriptor TVertexDescriptor;

    TVertexDescriptor v = addVertex(graph, vSeq);
    for (TOutEdgeIter it(graph.graph, u); !atEnd(it); ++it)
    {
        addEdge(graph.graph, v, targetVertex(it));
    }
    
    removeOutEdges(graph.graph, u);
    graph.sequenceMap[u] = uSeq;
    
    addEdge(graph.graph, u, v);
    
    return v;
}

// --------------------------------------------------------------------------
// Function addReverseComplementContigs()
// --------------------------------------------------------------------------

template<typename TSeq>
void
addReverseComplementContigs(StringSet<TSeq> & contigs,
                      StringSet<ContigId> & ids)
{
    typedef typename Size<TSeq>::Type TSize;
    SEQAN_ASSERT_EQ(length(contigs), length(ids));

    TSize len = length(contigs);

    resize(contigs, 2*len);
    resize(ids, 2*len);

    for (TSize i = 0; i < len; ++i)
    {
        TSeq contig_rev = contigs[i];
        reverseComplement(contig_rev);
        contigs[len + i] = contig_rev;

        SEQAN_ASSERT(ids[i].orientation);
        ContigId id_rev = ids[i];
        id_rev.orientation = false;
        ids[len + i] = id_rev;
    }
}

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
// Function readContigs()
// --------------------------------------------------------------------------

template<typename TSeq>
bool
readContigs(StringSet<TSeq> & contigs,
            StringSet<ContigId> & contigIds,
            String<CharString> & filenames,
            bool verbose)
{
    if (verbose) std::cerr << "[" << time(0) << "] " << "Reading contig files" << std::endl;

    unsigned totalContigCount = 0;

    // open and read the files containing contigs of one individual one by one
    for (unsigned i = 0; i < length(filenames); ++i)
    {
        // open fasta file
        SequenceStream stream(toCString(filenames[i]));
        if (!isGood(stream))
        {
            std::cerr << "ERROR: Could not open " << filenames[i] << std::endl;
            return 1;
        }
        
        // get an index for the file
        CharString index = formattedIndex(i, length(filenames));
        
        // read records from fasta file
        unsigned numContigs = 0;
        unsigned basepairs = 0;
        while (!atEnd(stream))
        {
            CharString id;
            TSeq contig;

            // read record
            if (readRecord(id, contig, stream))
            {
                std::cerr << "ERROR: Could not read fasta record from " << filenames[i] << std::endl;
                return 1;
            }

            // append contig and contigId
            appendValue(contigs, contig);
            appendValue(contigIds, ContigId(index, id, true));
            ++numContigs;
            basepairs += length(contig);
        }

        if (verbose) std::cerr << "[" << time(0) << "] " << "Loaded " << filenames[i] << ": " << basepairs << " bp in " << numContigs << " contigs." << std::endl;
        totalContigCount += numContigs;
    }

    if (verbose) std::cerr << "[" << time(0) << "] " << "Total number of contigs: " << totalContigCount << std::endl;

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
    typedef StringSet<TSeq> TStringSet;

    typedef Index<TStringSet, IndexQGram<SimpleShape, OpenAddressing> > TIndex;
    typedef Finder<TSeq, Swift<SwiftLocal> > TFinder;

    if (verbose) std::cerr << "[" << time(0) << "] " << "Partitioning sets of contigs" << std::endl;

    // --- Initialization for determining components
    
    // Initialize Union-Find data structure.
    UnionFind<int> uf;
    resize(uf, length(contigs));
    
    std::map<TSize, std::set<TSize> > alignedPairs;
    
    TSize len = length(contigs)/2;
    TSize numComparisons = 0;
    TSize numAligns = 0;

    // --- Initialization of SWIFT pattern and finders

    TIndex qgramIndex(contigs);
    resize(indexShape(qgramIndex), qgramLength);
    Pattern<TIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);
    indexRequire(qgramIndex, QGramSADir());

    // define scoring scheme
    Score<int, Simple> scoringScheme(matchScore, errorPenalty, errorPenalty);
    
    if (verbose) std::cerr << "0%   10   20   30   40   50   60   70   80   90  100%" << std::endl;
    if (verbose) std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;
    unsigned fiftieth = std::max(length(contigs)/50, (TSize)1);
   
    for (unsigned a = 0; a < length(contigs); ++a)
    {
        if (verbose && a%fiftieth == 0) std::cerr << "*" << std::flush;
        
        TFinder swiftFinder(contigs[a], 1000, 1);
        
        // --- SWIFT search
        
        while (find(swiftFinder, swiftPattern, errorRate, minimalLength)) {
            // get index of pattern sequence
            unsigned b = swiftPattern.curSeqNo;
            
            // align contigs only of different individuals and only if not same component already
            if (contigIds[a].pn == contigIds[b].pn) continue;
            if (findSet(uf, a) == findSet(uf, b)) continue;

            // verify by SW alignment
            ++numComparisons;
            unsigned upperDiag = (*swiftFinder.curHit).hstkPos - (*swiftFinder.curHit).ndlPos;
            unsigned lowerDiag = upperDiag - swiftPattern.bucketParams[b].delta - swiftPattern.bucketParams[b].overlap;
            if (!pairwiseAlignment(contigs[a], contigs[b], scoringScheme, lowerDiag, upperDiag, minScore)) continue;

            ++numAligns;
            alignedPairs[a].insert(b);
            alignedPairs[b].insert(a);
                
            // join sets of the two aligned contigs
            joinSets(uf, findSet(uf, a), findSet(uf, b));

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
    
    // -- Determine components from Union-Find data structure

    // Map ids and contigs to their representative ContigId.
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

    if (verbose)
    {
        std::cerr << "[" << time(0) << "] " << "Number of pairwise comparisons: \t" << numComparisons << std::endl;
        std::cerr << "[" << time(0) << "] " << "Number of valid alignments:     \t" << numAligns << std::endl;
        
        std::cerr << "[" << time(0) << "] " << "There are " << components.size() << " components in total." << std::endl;
    }
}

// --------------------------------------------------------------------------
// Function getSeqsByAlignOrder()
// --------------------------------------------------------------------------

template<typename TSeq, typename TSpec, typename TId>
void
getSeqsByAlignOrder(ContigComponent<TSeq> & component, StringSet<TSeq, TSpec> & contigs, StringSet<TId, TSpec> & contigIds)
{
    typedef typename Size<TSeq>::Type TSize;

    // --- determine the order ---
    
    String<TSize> order;
    std::set<TSize> ordered;
    appendValue(order, component.alignedPairs.begin()->first);
    ordered.insert(component.alignedPairs.begin()->first);
    
    TSize i = 0;
    while (i < length(order))
    {
        typename std::set<TSize>::iterator neighborsEnd = component.alignedPairs[order[i]].end();
        typename std::set<TSize>::iterator neighbor = component.alignedPairs[order[i]].begin();
        for (; neighbor != neighborsEnd; ++neighbor)
        {
            if (ordered.count(*neighbor) == 0)
            {
                // add the neighbor
                appendValue(order, *neighbor);
                ordered.insert(*neighbor);
            }
            
        }
        ++i;
    }
    clear(ordered);
    
    // --- bring contigs and contig ids into the order ---
    for (TSize i = 0; i < length(order); ++i)
    {
        appendValue(component.ids, contigIds[order[i]]);
        appendValue(component.contigs, contigs[order[i]]);
    }
}

// ==========================================================================

template<typename TSeq>
typename Size<String<TSeq> >::Type
enumeratePathsDfs(String<Path<TSeq, typename ComponentGraph<TSeq>::TVertexDescriptor> > & paths,
                  Path<TSeq, typename ComponentGraph<TSeq>::TVertexDescriptor> & prevPath,
                  ComponentGraph<TSeq> & graph,
                  typename ComponentGraph<TSeq>::TVertexDescriptor & v)
{
    typedef typename ComponentGraph<TSeq>::TVertexDescriptor TVertexDescriptor;
    typedef typename Size<String<Path<TSeq, TVertexDescriptor> > >::Type TSize;
    typedef typename Iterator<typename ComponentGraph<TSeq>::TGraph_, OutEdgeIterator>::Type TOutEdgeIter;
    
    TSize len = length(paths);
    
    append(prevPath.seq, graph.sequenceMap[v]);
    prevPath.positionMap[length(prevPath.seq)] = v;
    
    if (outDegree(graph.graph, v) == 0)
    {
        appendValue(paths, prevPath);
        return 1;
    }
    
    for (TOutEdgeIter it(graph.graph, v); !atEnd(it); ++it)
    {
        Path<TSeq, TVertexDescriptor> path(prevPath);
        TVertexDescriptor u = targetVertex(it);
        enumeratePathsDfs(paths, path, graph, u);
    }    
    
    return length(paths) - len;
}

// --------------------------------------------------------------------------
// Function enumeratePaths()
// --------------------------------------------------------------------------

template<typename TSeq>
typename Size<String<TSeq> >::Type
enumeratePaths(String<Path<TSeq, typename ComponentGraph<TSeq>::TVertexDescriptor> > & paths,
               ComponentGraph<TSeq> & graph)
{
    typedef typename ComponentGraph<TSeq>::TVertexDescriptor TVertexDescriptor;
    typedef typename Size<String<TSeq> >::Type TSize;
    
    for (TSize i = 0; i < length(graph.sources); ++i)
    {
        Path<TSeq, TVertexDescriptor> path;
        enumeratePathsDfs(paths, path, graph, graph.sources[i]);
    }

    return length(paths);
}

// --------------------------------------------------------------------------
// Function bestDiagonal()
// --------------------------------------------------------------------------

template<typename TSeq1, typename TSeq2>
int
bestDiagonal(TSeq1 & seq1, TSeq2 & seq2, unsigned qgramLength)
{
    typedef Index<TSeq1, IndexQGram<SimpleShape, OpenAddressing> > TIndex;
    typedef typename Infix<typename Fibre<TIndex, FibreSA>::Type const>::Type TOccurrences;
    typedef typename Iterator<TOccurrences>::Type TOccIter;

    unsigned len1 = length(seq1);
    unsigned len2 = length(seq2);
    
    if (qgramLength > len1 || qgramLength > len2) return maxValue<int>();

    // Build a k-mer index of seq1
    TIndex qgramIndex(seq1);
    resize(indexShape(qgramIndex), qgramLength);
    indexRequire(qgramIndex, QGramSADir());

    // Init diagonal counters
    String<unsigned> counters;
    resize(counters, len1+len2, 0);

    // Init hash function
    Shape<typename Value<TSeq1>::Type, SimpleShape> myShape(qgramLength);
    hashInit(myShape, begin(seq2));

    // Iterate over seq2 to count k-mer hits per diagonal
    for (unsigned i = 0; i < length(seq2) - length(myShape) + 1; ++i)
    {
        // Compute hash of the k-mer
        hashNext(myShape, begin(seq2) + i);

        // Increase counters of diagonals with hits
        TOccurrences occs = getOccurrences(qgramIndex, myShape);
        TOccIter itEnd = end(occs);
        for (TOccIter it = begin(occs); it != itEnd; ++it)
            ++counters[len1 + i - *it];
    }

    // Return the diagonal with the most k-mer hits
    int diag = maxValue<int>();
    unsigned maxCount = 0;
    for (unsigned i = 0; i < length(counters); ++i)
    {
        if (maxCount < counters[i])
        {
             maxCount = counters[i];
             diag = i - len1;
        }
    }
    
    if (diag == maxValue<int>()) return bestDiagonal(seq1, seq2, qgramLength*2/3);
    
    return diag;
}

// --------------------------------------------------------------------------
// Function mergeSeqWithGraph()
// --------------------------------------------------------------------------

template<typename TSeq1, typename TSeq2, typename TVertexDescriptor, typename TLength>
bool mergeSeqWithGraph(ComponentGraph<TSeq1> & compGraph,
                       Path<TSeq1, TVertexDescriptor> & path,
                       TSeq2 & seq,
                       Gaps<TSeq1> & gapsH,
                       Gaps<TSeq2> & gapsV,
                       TLength minBranchLen)
{
    typedef typename Position<TSeq1>::Type TPos;
    typedef typename Size<TSeq1>::Type TSize;

    // --- handle right end of alignment

    TPos endPosV = toSourcePosition(gapsV, length(gapsV)); // end position of alignment in seq
    TPos endPosH = toSourcePosition(gapsH, length(gapsH)); // end position of alignment in path.seq

    if (endPosV < length(seq))
    {
        TPos vPos = path.positionMap.lower_bound(endPosH)->first;
        TVertexDescriptor v = path.positionMap[vPos];
            
        if (endPosH == length(path.seq))
        {
            append(compGraph.sequenceMap[v], suffix(seq, endPosV));
        }
        else if (length(seq) - endPosV > (TSize)minBranchLen)
        {

            if (vPos > endPosH)
            {
                TPos splitPos = length(compGraph.sequenceMap[v]) - (vPos - endPosH);
                TSeq1 prefixSeq = prefix(compGraph.sequenceMap[v], splitPos);
                TSeq1 suffixSeq = suffix(compGraph.sequenceMap[v], splitPos);
                splitVertex(compGraph, v, prefixSeq, suffixSeq);
            }
            TSeq1 suf = suffix(seq, endPosV);
            TVertexDescriptor vBranch = addVertex(compGraph, suf);
            addEdge(compGraph.graph, v, vBranch);
        }
    }

    // --- handle left end of alignment

    TPos beginPosV = toSourcePosition(gapsV, 0); // begin position of alignment in seq
    TPos beginPosH = toSourcePosition(gapsH, 0); // begin position of alignment in path.seq

    if (beginPosV > 0)
    {
        TPos uPos = path.positionMap.lower_bound(beginPosH)->first;
        TVertexDescriptor u = path.positionMap[uPos];
        
        if (beginPosH == 0)
        {
            replace(compGraph.sequenceMap[u], 0, 0, prefix(seq, beginPosV));
        }
        else if (beginPosV > (TSize)minBranchLen)
        {
    
            TVertexDescriptor uSplit = u;
            if (uPos - length(compGraph.sequenceMap[u]) < beginPosH)
            {
                TPos splitPos = length(compGraph.sequenceMap[u]) - (uPos - beginPosH);
                TSeq1 prefixSeq = prefix(compGraph.sequenceMap[u], splitPos);
                TSeq1 suffixSeq = suffix(compGraph.sequenceMap[u], splitPos);
                uSplit = splitVertex(compGraph, u, prefixSeq, suffixSeq);
            }
            TSeq1 pref = prefix(seq, beginPosV);
            TVertexDescriptor uBranch = addVertex(compGraph, pref);
            appendValue(compGraph.sources, uBranch);
            addEdge(compGraph.graph, uBranch, uSplit);
        }
    }

    return true;
}

// --------------------------------------------------------------------------
// Function addSequencesToGraph()
// --------------------------------------------------------------------------

template<typename TSeq1, typename TSeq2, typename TSpec, typename TLength, typename TValueMatch, typename TValueError>
bool
addSequencesToGraph(ComponentGraph<TSeq1> & compGraph,
                    StringSet<TSeq2, TSpec> & seqs,
                    TLength minBranchLen,
                    TValueMatch matchScore,
                    TValueError errorPenalty,
                    unsigned qgramLength)
{
    typedef int TScoreValue;
    typedef ComponentGraph<TSeq1> TGraph;
    typedef Path<TSeq1, typename TGraph::TVertexDescriptor> TPath;
    typedef typename Size<String<TPath> >::Type TSize;

    Score<TScoreValue, Simple> scoringScheme(matchScore, errorPenalty, errorPenalty);

    for (TSize i = 1; i < length(seqs); ++i)
    {
        String<TPath> paths;
        enumeratePaths(paths, compGraph);

        TScoreValue maxScore = minValue<TScoreValue>();
        TPath bestPath;
        Gaps<TSeq1> bestH;
        Gaps<TSeq2> bestV;

        for (TSize j = 0; j < length(paths); ++j)
        {
            Gaps<TSeq1> gapsH(paths[j].seq);
            Gaps<TSeq2> gapsV(seqs[i]);

            int diag = bestDiagonal(seqs[i], paths[j].seq, qgramLength);

            TScoreValue val = 0;
            if (diag == maxValue<int>()) val = localAlignment(gapsH, gapsV, scoringScheme);
            else val = localAlignment(gapsH, gapsV, scoringScheme, diag-25, diag+25);
            if (val > maxScore)
            {
                maxScore = val;
                bestPath = paths[j];
                bestH = gapsH;
                bestV = gapsV;
            }
        }

        mergeSeqWithGraph(compGraph, bestPath, seqs[i], bestH, bestV, minBranchLen);
    }

    return true;
}

// ==========================================================================
// Function mergeSequences()
// ==========================================================================

template<typename TSeq1, typename TSeq2, typename TSpec, typename TLength, typename TValueMatch, typename TValueError>
bool
mergeSequences(String<TSeq1> & mergedSeqs,
               StringSet<TSeq2, TSpec> & seqs,
               TLength & minBranchLen,
               TValueMatch matchScore,
               TValueError errorPenalty,
               unsigned qgramLength,
               bool verbose)
{
    typedef int TScoreValue;    
    typedef ComponentGraph<TSeq1> TGraph;
    typedef Path<TSeq1, typename TGraph::TVertexDescriptor> TPath;
    typedef typename Size<String<TPath> >::Type TSize;

    TGraph compGraph(seqs[0]);
    bool ret = addSequencesToGraph(compGraph, seqs, minBranchLen, matchScore, errorPenalty, qgramLength);
    
    String<TPath> finalPaths;
    enumeratePaths(finalPaths, compGraph);
    
    if (verbose && numVertices(compGraph.graph) > 1)
    {
        std::cout << compGraph.graph;
        std::cout << "Vertex map:" << std::endl;
        for (TSize i = 0; i < length(compGraph.sequenceMap); ++i)
        {
            std::cout << "Vertex: " << i << ", Length: " << length(compGraph.sequenceMap[i]) << std::endl;
        }
    }

    for (TSize i = 0; i < length(finalPaths); ++i)
        appendValue(mergedSeqs, finalPaths[i].seq);
        
    return ret;
}

// ==========================================================================
// Function popins_merge()
// ==========================================================================

int popins_merge(int argc, char const ** argv)
{
    typedef Dna5String TSequence;
    typedef Size<TSequence>::Type TSize;
    typedef Position<TSequence>::Type TPosition;
    typedef std::map<ContigId, ContigComponent<TSequence> > TComponents;

    // Parse the command line to get option values.
    MergingOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    // Read contigs from file.
    StringSet<TSequence, Owner<> > contigs;
    StringSet<ContigId, Owner<> > contigIds;
    if (readContigs(contigs, contigIds, options.contigFiles, options.verbose) != 0) return 1;
 
    // Append reverse complements of contigs to string sets.
    addReverseComplementContigs(contigs, contigIds);

    // *** Partition contigs into sets of similar contigs (components) by computing pairwise alignments. ***
    TComponents components;
    partitionContigs(components, contigs, contigIds, options.errorRate, options.minimalLength,
                     options.qgramLength, options.matchScore, options.errorPenalty, options.minScore, options.verbose);

    // Prepare the output file.
    options.outputStream.open(toCString(options.outputFile), std::ios_base::out);
    if (!options.outputStream.is_open())
    {
        std::cerr << "ERROR: Could not open output file " << options.outputFile << std::endl;
        return 1;
    }
    
    if (options.verbose) std::cerr << "[" << time(0) << "] " << "Constructing supercontigs" << std::endl;

    unsigned numSingleton = 0;
    unsigned numBranching = 0;
    
    // Iterate over the set of components.
    unsigned pos = 0;
    for (TComponents::iterator it = components.begin(); it != components.end(); ++it)
    {
        ContigComponent<TSequence> component = it->second;

        // Sort the contigs for merging.
        getSeqsByAlignOrder(component, contigs, contigIds);

        // Output component if consisting of a single contig.
        if (length(component.contigs) == 1)
        {
            SEQAN_ASSERT_EQ(length(component.ids), 1u);
            options.outputStream << ">" << component.ids[0] << std::endl;
            options.outputStream << component.contigs[0] << std::endl;

            ++numSingleton;
            continue;
        }

        if (options.verbose) std::cout << "COMPONENT_" << pos << " size:" << length(component.contigs) << std::endl;

        // *** Merge contigs of the component. ***
        String<TSequence> mergedSeqs;
        
        /////// DEBUG CODE ///////
        //if (length(component.contigs) == 10 && pos == 1626)
        //{
        //    std::cout << std::endl;
        //    for (unsigned i = 0; i < length(component.ids); ++i)
        //        std::cout << component.ids[i] << std::endl;
        //}
        
        mergeSequences(mergedSeqs, component.contigs,
                       options.minTipScore, options.matchScore, options.errorPenalty, options.qgramLength,
                       options.verbose);

        if (length(mergedSeqs) > 1) ++numBranching;

        if (length(mergedSeqs) <= 25)
        {
            for (TSize i = 0; i < length(mergedSeqs); ++i)
            {
                options.outputStream << ">COMPONENT_" << pos << "_" << char('a'+i)
                                     << "_length_" << length(mergedSeqs[i])
                                     << "_size_" << length(component.contigs) << std::endl;
                options.outputStream << mergedSeqs[i] << std::endl;
            }
        }
        else
        {
            for (TSize i = 0; i < length(mergedSeqs); ++i)
            {
                options.outputStream << ">COMPONENT_" << pos << "_" << char('a'+i/26) << char('a'+i%26)
                                     << "_length_" << length(mergedSeqs[i])
                                     << "_size_" << length(component.contigs) << std::endl;
                options.outputStream << mergedSeqs[i] << std::endl;
            }
        }

        ++pos;
    }

    options.outputStream.close();

    if (options.verbose) 
    {
        std::cerr << "[" << time(0) << "] " << length(components)-numSingleton << " components are merged from several contigs." << std::endl;
        std::cerr << "[" << time(0) << "] " << numSingleton << " contigs did not align with any other contig." << std::endl;
        std::cerr << "[" << time(0) << "] " << numBranching << " components are branching." << std::endl;
    }

    return 0;
}
#endif // #ifndef POPINS_MERGE_H_
