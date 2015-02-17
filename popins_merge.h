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
    SEQAN_ASSERT_EQ(v, length(graph.sequenceMap));
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
                                                                    
template<typename TSeq>
int
getSubsetWithReverseComplementContigs(StringSet<TSeq> & contigSubset,
                                      StringSet<ContigId> & idSubset,
                                      StringSet<TSeq> & contigs,
                                      StringSet<ContigId> & ids,
                                      int totalBatches,
                                      int batchIndex)
{
    typedef typename Size<TSeq>::Type TSize;
    SEQAN_ASSERT_EQ(length(contigs), length(ids));

    int batchSize = length(contigs) / totalBatches;
    if (length(contigs) % totalBatches != 0) ++batchSize;

    int offset = batchSize * batchIndex;
    TSize actualSize = std::min((int)length(contigs) - offset, batchSize);

//    std::cout << "batchSize: " << batchSize << ", actualSize: " << actualSize << ", offset: " << offset << std::endl;

    resize(contigSubset, actualSize*2);
    resize(idSubset, actualSize*2);

    for (TSize i = 0; i < actualSize; ++i)
    {
        contigSubset[i] = contigs[i + offset];
        idSubset[i] = ids[i + offset];

        // Append reverse complements of contigs to subset.
        TSeq contig_rev = contigs[i + offset];
        reverseComplement(contig_rev);
        contigSubset[actualSize + i] = contig_rev;

        // Append id of reversed contig to subset.
        SEQAN_ASSERT(ids[i].orientation);
        ContigId id_rev = ids[i + offset];
        id_rev.orientation = false;
        idSubset[actualSize + i] = id_rev;
    }

    return offset;
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
    std::cerr << "[" << time(0) << "] " << "Reading contig files" << std::endl;

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

    std::cerr << "[" << time(0) << "] " << "Total number of contigs: " << totalContigCount << std::endl;

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

// --------------------------------------------------------------------------
// Function writeComponents()
// --------------------------------------------------------------------------

template<typename TSeq>
bool
writeComponents(CharString & fileName, std::map<ContigId, ContigComponent<TSeq> > & components)
{
    typedef typename Size<TSeq>::Type TSize;
    typedef std::map<ContigId, ContigComponent<TSeq> > TComponents;

    // Open the output file.
    std::fstream outputStream(toCString(fileName), std::ios::out);
    if (!outputStream.good())
    {
        std::cerr << "ERROR: Could not open component output file " << fileName << std::endl;
        return 1;
    }
    
    // Write the components.
    TSize i = 0;
    for (typename TComponents::iterator it = components.begin(); it != components.end(); ++it, ++i)
    {
        if (it->second.alignedPairs.begin()->second.size() == 0) continue;
        for (typename std::map<TSize, std::set<TSize> >::iterator pairIt = it->second.alignedPairs.begin(); pairIt != it->second.alignedPairs.end(); ++pairIt)
        {
            if (pairIt->second.size() == 0) continue;

            typename std::set<TSize>::iterator alignedTo = pairIt->second.begin();
            outputStream << " [" << pairIt->first << "->" << *alignedTo;
            if (alignedTo != pairIt->second.end()) ++alignedTo;
            for (; alignedTo != pairIt->second.end(); ++alignedTo)
                outputStream << "," << *alignedTo;

            outputStream << "]";
        }
        outputStream << "\n";
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function readComponents()
// --------------------------------------------------------------------------

template<typename TSeq, typename TSize>
bool
readComponents(String<ContigComponent<TSeq> > & components, UnionFind<int> & uf, CharString & fileName, TSize len, bool verbose)
{
    typedef String<ContigComponent<TSeq> > TComponents;
    
    unsigned numComps = length(components);
    
    // Open the input file and initialize record reader.
    std::fstream stream(toCString(fileName), std::ios::in);
    
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open components input file " << fileName << std::endl;
        return 1;
    }
    
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    CharString buffer;
    
    // Read the components line by line.
    while (!atEnd(reader))
    {
        if (skipChar(reader, ' ') != 0)
        {
            std::cerr << "ERROR: File format error. Reading whitespace from " << fileName << " failed." << std::endl;
            return 1;
        }

        ContigComponent<TSeq> component;
        while (value(reader) == '[')
        {
            TSize key, val, key_rev, val_rev;
            skipChar(reader, '[');
            
            clear(buffer);
            if (readDigits(buffer, reader) != 0)
            {
                std::cerr << "ERROR: File format error. Reading key from " << fileName << " failed." << std::endl;
                return 1;
            }
            lexicalCast2<TSize>(key, buffer);
            if (key < len) key_rev = key + len;
            else key_rev = key - len;
            
            if (skipChar(reader, '-') != 0 || skipChar(reader, '>') != 0)
            {
                std::cerr << "ERROR: File format error. Reading '->' from " << fileName << " failed." << std::endl;
                return 1;
            }
            
            clear(buffer);
            if (readDigits(buffer, reader) != 0)
            {
                std::cerr << "ERROR: File format error. Reading value from " << fileName << " failed." << std::endl;
                return 1;
            }
            lexicalCast2<TSize>(val, buffer);
            if (val < len) val_rev = val + len;
            else val_rev = val - len;
            
            // Add the aligned pairs.
            component.alignedPairs[key].insert(val);

            // Join sets of key and value.
            joinSets(uf, findSet(uf, key), findSet(uf, val));
            SEQAN_ASSERT_EQ(findSet(uf, key), findSet(uf, val));
            joinSets(uf, findSet(uf, key_rev), findSet(uf, val_rev));
            SEQAN_ASSERT_EQ(findSet(uf, key_rev), findSet(uf, val_rev));
            
            while (value(reader) == ',')
            {
                skipChar(reader, ',');
                clear(buffer);
                if (readDigits(buffer, reader) != 0)
                {
                    std::cerr << "ERROR: File format error. Reading value from " << fileName << " failed." << std::endl;
                    return 1;
                }
                lexicalCast2<TSize>(val, buffer);
                if (val < len) val_rev = val + len;
                else val_rev = val - len;
                
                // Add the aligned pairs.
                component.alignedPairs[key].insert(val);

                // Join sets of key and value.
                joinSets(uf, findSet(uf, key), findSet(uf, val));
                SEQAN_ASSERT_EQ(findSet(uf, key), findSet(uf, val));
                joinSets(uf, findSet(uf, key_rev), findSet(uf, val_rev));
                SEQAN_ASSERT_EQ(findSet(uf, key_rev), findSet(uf, val_rev));
            }
            
            if (skipChar(reader, ']') != 0)
            {
                std::cerr << "ERROR: File format error. Reading ']' from " << fileName << " failed." << std::endl;
                return 1;
            }
            
            if (value(reader) == ' ') skipChar(reader, ' ');
        }
        appendValue(components, component);
        skipLine(reader);
    }
    
    if (verbose) std::cerr << "[" << time(0) << "] " << "Loaded " << fileName << ": " << (length(components) - numComps) << " components." << std::endl;
    
    return 0;
}

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
                       Gaps<TSeq1> & gapsPath,
                       Gaps<TSeq2> & gapsSeq,
                       TLength minBranchLen)
{
    typedef typename Position<TSeq1>::Type TPos;
    typedef typename Size<TSeq1>::Type TSize;

    // --- handle right end of alignment

    TPos alignEndSeq = toSourcePosition(gapsSeq, length(gapsSeq)); // end position of alignment in seq
    TPos alignEndPath = toSourcePosition(gapsPath, length(gapsPath)); // end position of alignment in path.seq

    if (alignEndSeq < length(seq))
    {
        TPos vPos = path.positionMap.lower_bound(alignEndPath)->first; // end position of vertex label on path
        TVertexDescriptor v = path.positionMap[vPos];

        if (alignEndPath == length(path.seq)) // alignment ends at end of path
        {
            append(compGraph.sequenceMap[v], suffix(seq, alignEndSeq));
        }
        else if (length(seq) - alignEndSeq > (TSize)minBranchLen) // unaligned part of seq is longer than minBranchLen
        {
            if (vPos > alignEndPath) // alignment ends before end of vertex label
            {
                TPos splitPos = length(compGraph.sequenceMap[v]) - (vPos - alignEndPath); // relative to vertex label
                TSeq1 prefixSeq = prefix(compGraph.sequenceMap[v], splitPos);
                TSeq1 suffixSeq = suffix(compGraph.sequenceMap[v], splitPos);
                splitVertex(compGraph, v, prefixSeq, suffixSeq);
            }
            TSeq1 suf = suffix(seq, alignEndSeq);
            TVertexDescriptor vBranch = addVertex(compGraph, suf);
            addEdge(compGraph.graph, v, vBranch);
        }
    }

    // --- handle left end of alignment

    TPos alignBeginSeq = toSourcePosition(gapsSeq, 0); // begin position of alignment in seq
    TPos alignBeginPath = toSourcePosition(gapsPath, 0); // begin position of alignment in path.seq

    if (alignBeginSeq > 0)
    {
        TPos uPos = path.positionMap.upper_bound(alignBeginPath)->first;
        TVertexDescriptor u = path.positionMap[uPos];
        
        if (alignBeginPath == 0)
        {
            replace(compGraph.sequenceMap[u], 0, 0, prefix(seq, alignBeginSeq));
        }
        else if (alignBeginSeq > (TSize)minBranchLen)
        {
            TVertexDescriptor uSplit = u;
            if (uPos - length(compGraph.sequenceMap[u]) < alignBeginPath)
            {
                TPos splitPos = length(compGraph.sequenceMap[u]) - (uPos - alignBeginPath);
                TSeq1 prefixSeq = prefix(compGraph.sequenceMap[u], splitPos);
                TSeq1 suffixSeq = suffix(compGraph.sequenceMap[u], splitPos);
                uSplit = splitVertex(compGraph, u, prefixSeq, suffixSeq);
            }
            TSeq1 pref = prefix(seq, alignBeginSeq);
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
        
        if (length(paths) > 30) return false;

        TScoreValue maxScore = minValue<TScoreValue>();
        TPath bestPath;
        Gaps<TSeq1> bestGapsPath;
        Gaps<TSeq2> bestGapsSeq;

        for (TSize j = 0; j < length(paths); ++j)
        {
            Gaps<TSeq1> gapsPath(paths[j].seq);
            Gaps<TSeq2> gapsSeq(seqs[i]);

            int diag = bestDiagonal(seqs[i], paths[j].seq, qgramLength);

            TScoreValue val = 0;
            if (diag == maxValue<int>()) val = localAlignment(gapsPath, gapsSeq, scoringScheme);
            else val = localAlignment(gapsPath, gapsSeq, scoringScheme, diag-25, diag+25);

            if (val > maxScore)
            {
                maxScore = val;
                bestPath = paths[j];
                bestGapsPath = gapsPath;
                bestGapsSeq = gapsSeq;
            }
        }

        mergeSeqWithGraph(compGraph, bestPath, seqs[i], bestGapsPath, bestGapsSeq, minBranchLen);
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
    if (!addSequencesToGraph(compGraph, seqs, minBranchLen, matchScore, errorPenalty, qgramLength))
        return false;

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

    return true;
}

// ==========================================================================
// Function popins_merge()
// ==========================================================================

// TODO clp: option for batch size for merging components
// TODO clp: option for batch number for merging components

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

    TComponents components;
    if (length(options.componentFiles) == 0)
    {
        if (options.partitioningBatches > 1) {

            // Get subset of contigs and add reverse complements only in this subset
            StringSet<TSequence, Owner<> > contigSubset;
            StringSet<ContigId, Owner<> > contigIdSubset;
            int offset = getSubsetWithReverseComplementContigs(contigSubset, contigIdSubset, contigs, contigIds,
                                                               options.partitioningBatches,
                                                               options.partitioningBatchIndex);

            std::cerr << "[" << time(0) << "] " << "Partitioning sets of contigs, batch "
                      << options.partitioningBatchIndex << "/" << options.partitioningBatches << std::endl;

            // *** Partition contigs into sets of similar contigs (components) by computing pairwise alignments. ***
            partitionContigs(components, contigs, contigIds, contigSubset, contigIdSubset, offset, options.errorRate, options.minimalLength,
                             options.qgramLength, options.matchScore, options.errorPenalty, options.minScore, options.verbose);

            // Write set of components
            if (writeComponents(options.outputFile, components) != 0) return 1;

            return 0;
        }
        else
        {
            // Append reverse complements of contigs to string sets.
            addReverseComplementContigs(contigs, contigIds);

            std::cerr << "[" << time(0) << "] " << "Partitioning sets of contigs" << std::endl;

            // *** Partition contigs into sets of similar contigs (components) by computing pairwise alignments. ***
            partitionContigs(components, contigs, contigIds, options.errorRate, options.minimalLength,
                             options.qgramLength, options.matchScore, options.errorPenalty, options.minScore, options.verbose);
        }
    }    
    else
    {
        // Append reverse complements of contigs to string sets.
        addReverseComplementContigs(contigs, contigIds);

        // Read and merge component sets
        if (readAndMergeComponents(components, options.componentFiles, contigIds, options.verbose) != 0) return 1;
    }

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
    unsigned numVeryBranching = 0;
    unsigned numTooLarge = 0;
    
    // Iterate over the set of components.
    unsigned pos = 0;
    for (TComponents::iterator it = components.begin(); it != components.end(); ++it)
    {
        ContigComponent<TSequence> component = it->second;

        // Sort the contigs for merging.
        getSeqsByAlignOrder(component, contigs, contigIds);
        
        if (length(component.contigs) > 10 * length(options.contigFiles))
        {
            if (options.verbose) std::cout << "COMPONENT_" << pos << " size:" << length(component.contigs) << " skipped." << std::endl;
            ++numTooLarge;
            continue;
        }

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
        
        if (!mergeSequences(mergedSeqs, component.contigs,
                            options.minTipScore, options.matchScore, options.errorPenalty, options.qgramLength,
                            options.verbose))
        {
            if (options.verbose) std::cout << "COMPONENT_" << pos << " size:" << length(component.contigs) << " given up." << std::endl;
            ++numVeryBranching;
            ++numBranching;
            continue;
        }

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
        std::cerr << "[" << time(0) << "] " << numBranching << " components are branching, given up on " << numVeryBranching << " of them." << std::endl;
        std::cerr << "[" << time(0) << "] " << numTooLarge << " components exceeded the maximum number of contigs for merging." << std::endl;
    }

    return 0;
}
#endif // #ifndef POPINS_MERGE_H_
