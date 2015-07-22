#ifndef POPINS_MERGE_PARTITION_H_
#define POPINS_MERGE_PARTITION_H_

#include <seqan/index.h>
#include <seqan/align.h>

#include "contig_id.h"
#include "contig_component.h"

using namespace seqan;

// --------------------------------------------------------------------------
// Function readNextContig()
// --------------------------------------------------------------------------

template<typename TSeq, typename TSize>
int
readNextContig(TSeq & contig, ContigId & contigId, SequenceStream & stream, TSize & i, String<CharString> & filenames)
{
    // Open the next file.
    while ((i < (int)length(filenames) && atEnd(stream)) || i == -1)
    {
        ++i;
        open(stream, toCString(filenames[i]));
        if (!isGood(stream))
        {
            std::cerr << "ERROR: Could not open " << filenames[i] << std::endl;
            return -1;
        }
    }
    
    if (atEnd(stream)) return 1;
    
    // Read the next record.
    contigId.orientation = true;
    contigId.pn = formattedIndex(i, length(filenames));
    if (readRecord(contigId.contigId, contig, stream))
    {
        std::cerr << "ERROR: Could not read fasta record from " << filenames[i] << std::endl;
        return -1;
    }
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
                  int lowerDiag,
                  int upperDiag,
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

template<typename TSize, typename TSeq>
bool
partitionContigs(UnionFind<int> & uf,
                 std::set<Pair<TSize> > & alignedPairs,
                 StringSet<TSeq> & contigs,
                 StringSet<ContigId> & contigIds,
                 int totalContigs,
                 unsigned offset,
                 MergingOptions & options)
{
    typedef Index<StringSet<TSeq> , IndexQGram<SimpleShape, OpenAddressing> > TIndex;
    typedef Finder<TSeq, Swift<SwiftLocal> > TFinder;
    
    if (options.verbose) std::cerr << "[" << time(0) << "] " << "Partitioning contigs" << std::endl;
    if (options.verbose) std::cerr << "[" << time(0) << "] " << "- Indexing batch of contigs" << std::endl;

    TSize numComparisons = 0;

    // initialization of SWIFT pattern (q-gram index)
    TIndex qgramIndex(contigs);
    resize(indexShape(qgramIndex), options.qgramLength);
    Pattern<TIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);
    indexRequire(qgramIndex, QGramSADir());

    // define scoring scheme
    Score<int, Simple> scoringScheme(options.matchScore, options.errorPenalty, options.errorPenalty);
    int diagExtension = options.minScore/10;

    // print status bar
    if (options.verbose) std::cerr << "[" << time(0) << "] " << "- Streaming over all contig files" << std::endl;
    if (options.verbose) std::cerr << "0%   10   20   30   40   50   60   70   80   90  100%" << std::endl;
    if (options.verbose) std::cerr << "|----|----|----|----|----|----|----|----|----|----|" << std::endl;
    unsigned fiftieth = std::max((offset+length(contigs)/2)/50, (TSize)1);

    // stream over the contigs
    SequenceStream contigStream;
    int i = -1;
    for (unsigned a = 0; a < offset+length(contigs)/2; ++a)
    {
        if (options.verbose && a%fiftieth == 0) std::cerr << "*" << std::flush;
        
        // read the next contig
        TSeq contig;
        ContigId contigId;
        int ret = readNextContig(contig, contigId, contigStream, i, options.contigFiles);
        SEQAN_ASSERT_NEQ(ret, 1);
        if (ret == -1) return 1;

        // initialization of swift finder
        TFinder swiftFinder(contig, 1000, 1);
        while (find(swiftFinder, swiftPattern, options.errorRate, options.minimalLength)) {
        
            // get index of pattern sequence
            unsigned bSubset = swiftPattern.curSeqNo;
            
            // align contigs only of different individuals 
            if (contigId.pn == contigIds[bSubset].pn) continue;
            
            // convert index bSubset to the index space over all contigs
            unsigned b = bSubset + offset;
            if (bSubset >= length(contigs)/2) b += (totalContigs - length(contigs)/2);

            // align contigs only if not same component already
            if (findSet(uf, a) == findSet(uf, b)) continue;

            // find the contig sequences
            TSeq contigA = haystack(swiftFinder);
            TSeq contigB = indexText(needle(swiftPattern))[bSubset];
    
            // compute upper and lower diagonal of band.
            int upperDiag = (*swiftFinder.curHit).hstkPos - (*swiftFinder.curHit).ndlPos;
            int lowerDiag = upperDiag - swiftPattern.bucketParams[bSubset].delta - swiftPattern.bucketParams[bSubset].overlap;
            upperDiag += diagExtension;
            lowerDiag -= diagExtension;
    
            // verify by banded Smith-Waterman alignment
            ++numComparisons;
            if (!pairwiseAlignment(contigA, contigB, scoringScheme, lowerDiag, upperDiag, options.minScore)) continue;
            alignedPairs.insert(Pair<TSize>(a, b));

            // join sets of the two aligned contigs
            joinSets(uf, findSet(uf, a), findSet(uf, b));

            // join sets for reverse complements of the contigs
            int a1 = a;
            if (a1 < totalContigs) a1 += totalContigs; // a was index of contig in forward orientation, set it to index of contig in reverse orientation
            else a1 -= totalContigs;                   // a was index of contig in reverse orientation, set it to index of contig in forward orientation
            if (b < (unsigned)totalContigs) b += totalContigs;   // b was index of contig in forward orientation, set it to index of contig in reverse orientation
            else b -= totalContigs;                    // b was index of contig in reverse orientation, set it to index of contig in forward orientation
            joinSets(uf, findSet(uf, a1), findSet(uf, b));
            
            // stop aligning this contig if it is already in a component with more than 100 other contigs
            if (uf._values[findSet(uf, a)] < -100) break;
        }
    }
    if (options.verbose) std::cerr << std::endl;

    std::cerr << "[" << time(0) << "] " << "Number of pairwise comparisons: " << numComparisons << std::endl;
    std::cerr << "[" << time(0) << "] " << "Number of valid alignments:     " << length(alignedPairs) << std::endl;

    return 0;
}

// --------------------------------------------------------------------------
// Function writeAlignedPairs()
// --------------------------------------------------------------------------

template<typename TStream, typename TSize>
void
writeAlignedPairs(TStream & outputStream, std::set<Pair<TSize> > & alignedPairs)
{
    typedef typename std::set<Pair<TSize> >::iterator TIter;

    TIter pairsEnd = alignedPairs.end();
    for (TIter pairsIt = alignedPairs.begin(); pairsIt != pairsEnd; ++pairsIt)
        outputStream << (*pairsIt).i1 << " " << (*pairsIt).i2 << "\n";
}

// --------------------------------------------------------------------------
// Function readAlignedPairs()
// --------------------------------------------------------------------------

template<typename TSize>
bool
readAlignedPairs(UnionFind<int> & uf, std::set<Pair<TSize> > & alignedPairs, CharString & fileName, unsigned len, bool verbose)
{
   // Open the input file and initialize record reader.
    std::fstream stream(toCString(fileName), std::ios::in);
    
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open components input file " << fileName << std::endl;
        return 1;
    }
    
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    CharString buffer;
    
    TSize numPairs = 0;
    
    // Read the components line by line.
    while (!atEnd(reader))
    {
        TSize key, val, key_rev, val_rev;
        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: File format error. Reading key from " << fileName << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<TSize>(key, buffer);
        if (key < len) key_rev = key + len;
        else key_rev = key - len;

        skipWhitespaces(reader);

        clear(buffer);
        if (readDigits(buffer, reader) != 0)
        {
            std::cerr << "ERROR: File format error. Reading value from " << fileName << " failed." << std::endl;
            return 1;
        }
        lexicalCast2<TSize>(val, buffer);
        if (val < len) val_rev = val + len;
        else val_rev = val - len;

        skipLine(reader);

        if (findSet(uf, key) == findSet(uf, val)) continue;
        
        // Add the aligned pairs.
        alignedPairs.insert(Pair<TSize>(key, val));
        ++numPairs;

        // Join sets of key and value.
        joinSets(uf, findSet(uf, key), findSet(uf, val));
        SEQAN_ASSERT_EQ(findSet(uf, key), findSet(uf, val));
        joinSets(uf, findSet(uf, key_rev), findSet(uf, val_rev));
        SEQAN_ASSERT_EQ(findSet(uf, key_rev), findSet(uf, val_rev));

    }
    
    if (verbose) std::cerr << "[" << time(0) << "] " << "Loaded " << fileName << ": " << numPairs << " pairs." << std::endl;
    
    return 0;
}

// --------------------------------------------------------------------------
// Function unionFindToComponents()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
std::set<int>
unionFindToComponents(std::map<TSize, ContigComponent<TSeq> > & components,
                      UnionFind<int> & uf,
                      std::set<Pair<TSize> > & alignedPairs,
                      unsigned samples,
                      int totalContigs,
                      bool verbose)
{
    std::set<int> skipped;

    // Determine components from Union-Find data structure by
    // mapping ids to their representative id.
    for (typename std::set<Pair<TSize> >::iterator it = alignedPairs.begin(); it != alignedPairs.end(); ++it)
    {
        int rev1 = (*it).i1;
        if (rev1 < totalContigs) rev1 += totalContigs;
        else rev1 -= totalContigs;
        
        int rev2 = (*it).i2;
        if (rev2 < totalContigs) rev2 += totalContigs;
        else rev2 -= totalContigs;

        int set = std::min(findSet(uf, (*it).i1), findSet(uf, rev1));

        // skip components that are 10 times larger than number of samples
        if (uf._values[set] < -10 * (int)samples) 
        {
            if (verbose && skipped.count(set) == 0)
            {
                std::cerr << "[" << time(0) << "] " << "WARNING: Skipping component of size " << (-1*uf._values[set]) << std::endl;
            }
            skipped.insert(set);
            skipped.insert((*it).i1);
            skipped.insert((*it).i2);
            skipped.insert(rev1);
            skipped.insert(rev2);
            continue;
        }

        components[set].alignedPairs.insert(*it);
        components[set].alignedPairs.insert(Pair<TSize>((*it).i2, (*it).i1));
        components[set].alignedPairs.insert(Pair<TSize>(rev1, rev2));
        components[set].alignedPairs.insert(Pair<TSize>(rev2, rev1));
    }

    std::cerr << "[" << time(0) << "] " << "There are " << components.size() << " components." << std::endl;
    return skipped;
}

// --------------------------------------------------------------------------
// Function addSingletons()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
void
addSingletons(std::map<TSize, ContigComponent<TSeq> > & components, UnionFind<int> & uf, int totalContigs)
{
    unsigned numSingletons = 0;
    for (int i = 0; i < totalContigs; ++i)
    {
        if (components.count(i) == 0 && i == findSet(uf, i))
        {
            components[i];
            ++numSingletons;
        }
    }

    std::cerr << "[" << time(0) << "] " << "Added " << numSingletons << " singletons to components." << std::endl;
}

// ==========================================================================
// Function readAndMergeComponents()
// ==========================================================================

template<typename TSize, typename TSequence>
bool
readAndMergeComponents(std::map<TSize, ContigComponent<TSequence> > & components,
                       std::set<int> & skipped,
                       String<CharString> & componentFiles,
                       unsigned samples,
                       int numContigs,
                       int batchIndex,
                       int batches,
                       bool verbose)
{
    typedef std::map<TSize, ContigComponent<TSequence> > TComponents;
    typedef typename std::set<Pair<TSize> >::iterator TPairsIterator;
    typedef typename TComponents::iterator TCompIterator;

    std::cerr << "[" << time(0) << "] " << "Reading and merging components files" << std::endl;

    // Initialize Union-Find data structure.
    UnionFind<int> uf;
    resize(uf, numContigs*2);
    std::set<Pair<TSize> > alignedPairs;
    
    // Read the aligned pairs from input files and join sets.
    for (unsigned i = 0; i < length(componentFiles); ++i)
        if (readAlignedPairs(uf, alignedPairs, componentFiles[i], numContigs, verbose) != 0) return 1;

    // Convert union-find data structure to components.
    skipped = unionFindToComponents(components, uf, alignedPairs, samples, numContigs, verbose);

    // Add singleton contigs to components (= those contigs that don't align to any other contig).
    addSingletons(components, uf, numContigs);

    // Keep only batch of components (erase all except every batchSize'th component).
    if (batches != 1)
    {
        TCompIterator it = --components.end();
        TSize i = components.size();
        while (i > 0)
        {
            TCompIterator element = it;
            --it; --i;
            if (i % batches != (TSize)batchIndex)
                components.erase(element);
        }
    }

    return 0;
}

#endif // #ifndef POPINS_MERGE_PARTITION_H_
