#ifndef POPINS_MERGE_PARTITION_H_
#define POPINS_MERGE_PARTITION_H_

#include <seqan/index.h>
#include <seqan/align.h>

#include "contig_id.h"
#include "contig_structs.h"

using namespace seqan;

// --------------------------------------------------------------------------
// Function calculateEntropy()
// --------------------------------------------------------------------------

template<typename TSeq>
double
averageEntropy(TSeq & seq)
{
    typedef typename Size<TSeq>::Type TSize;

    // Count dinucleotide occurrences
    String<TSize> diCounts;
    resize(diCounts, 16, 0);
    int counted = 0;
    for (TSize i = 0; i < length(seq)-1; ++i)
    {
        if (seq[i] != 'N' && seq[i+1] != 'N')
        {
            diCounts[ordValue(seq[i]) + 4*ordValue(seq[i+1])] += 1;
            counted += 1;  
        }
    }

    // Calculate entropy for dinucleotide counts
    double entropy = 0;
    typename Iterator<String<TSize> >::Type countEnd = end(diCounts);
    for (typename Iterator<String<TSize> >::Type count = begin(diCounts); count != countEnd; ++count)
    {
        if (*count == 0) continue;
        double p = double(*count) / counted;
        entropy -= p * log(p) / log(2);
    }

    return entropy / 4;
}

// ==========================================================================
// Function filterByEntropy()
// ==========================================================================

template<typename TSize, typename TSeq>
bool
filterByEntropy(std::map<TSize, Contig<TSeq> > & contigs,
                MergingOptions & options)
{
    typedef typename std::map<TSize, Contig<TSeq> >::iterator TIter;
    String<TSize> lowEntropyContigs;
    
    // Iterate contigs and determine entropy
    TIter itEnd = contigs.end();
    for (TIter it = contigs.begin(); it != itEnd; ++it)
    {
        // Entropy calculation
        double entropy = averageEntropy((it->second).seq);

        if (entropy < options.minEntropy)
        {
            options.skippedStream << ">" << (it->second).id << " (entropy filter, entropy: " << entropy << ")" << std::endl;
            options.skippedStream << (it->second).seq << std::endl;
            appendValue(lowEntropyContigs, it->first);
        }
    }

    typename Iterator<String<TSize> >::Type lowEnd = end(lowEntropyContigs);
    for (typename Iterator<String<TSize> >::Type it = begin(lowEntropyContigs); it != lowEnd; ++it)
        contigs.erase(*it);

    if (length(contigs) == 0)
    {
        std::cerr << "There are no contigs that passed the entropy filter." << std::endl;
        return 1;
    }
    
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Passed entropy filter: " << length(contigs) << std::endl;

    return 0;
}

// --------------------------------------------------------------------------
// Function readNextContig()
// --------------------------------------------------------------------------

template<typename TSeq, typename TSize>
int
readNextContig(Contig<TSeq> & contig, SequenceStream & stream, TSize & i, String<CharString> & filenames)
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
    contig.id.orientation = true;
    contig.id.pn = formattedIndex(i, length(filenames));
    if (readRecord(contig.id.contigId, contig.seq, stream))
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
                 std::map<TSize, Contig<TSeq> > & contigs,
                 ContigBatch & batch,
                 MergingOptions & options)
{
    typedef typename std::map<TSize, Contig<TSeq> >::iterator TContigIter;
    typedef StringSet<TSeq, Dependent<> > TStringSet;
    typedef Index<TStringSet, IndexQGram<SimpleShape, OpenAddressing> > TIndex;
    typedef Finder<TSeq, Swift<SwiftLocal> > TFinder;
    
    if (options.verbose) std::cerr << "[" << time(0) << "] " << "Partitioning contigs" << std::endl;
    if (options.verbose) std::cerr << "[" << time(0) << "] " << "- Indexing batch of contigs" << std::endl;

    TSize numComparisons = 0;

    // initialization of SWIFT pattern (q-gram index)
    TStringSet seqs;
    StringSet<TSize> indices;
    TContigIter itEnd = contigs.end();
    for (TContigIter it = contigs.begin(); it != itEnd; ++it)
    {
        appendValue(seqs, (it->second).seq);
        appendValue(indices, it->first);
    }
    TIndex qgramIndex(seqs);
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
    //unsigned fiftieth = std::max((indexOffset(batch)+length(contigs)/2)/50, 1u);
    unsigned fiftieth = std::max((indexOffset(batch)+batchSize(batch))/50, 1);

    // stream over the contigs
    SequenceStream contigStream;
    int i = -1;
    //for (unsigned a = 0; a < indexOffset(batch)+length(contigs)/2; ++a)
    for (int a = 0; a < indexOffset(batch)+batchSize(batch); ++a)
    {
        if (options.verbose && a%fiftieth == 0) std::cerr << "*" << std::flush;
        
        // read the next contig
        Contig<TSeq> contig;
        int ret = readNextContig(contig, contigStream, i, batch.contigFiles);
        SEQAN_ASSERT_NEQ(ret, 1);
        if (ret == -1) return 1;
        if (contigs.count(a) == 0) continue; // skipped contig

        // initialization of swift finder
        TFinder swiftFinder(contig.seq, 1000, 1);
        
        hash(swiftPattern.shape, hostIterator(hostIterator(swiftFinder)));
        while (find(swiftFinder, swiftPattern, options.errorRate, options.minimalLength))
        {
        
            // get index of pattern sequence
            unsigned bSubset = swiftPattern.curSeqNo;
            unsigned b = indices[bSubset];
            
            // align contigs only of different individuals 
            if (contig.id.pn == contigs[b].id.pn) continue;

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
            unsigned a1 = globalIndexRC(a, batch);
            unsigned b1 = globalIndexRC(b, batch);
            joinSets(uf, findSet(uf, a1), findSet(uf, b1));
            
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
                      ContigBatch & batch,
                      bool /*verbose*/)
{
    std::set<int> skipped;

    // Determine components from Union-Find data structure by
    // mapping ids to their representative id.
    for (typename std::set<Pair<TSize> >::iterator it = alignedPairs.begin(); it != alignedPairs.end(); ++it)
    {
        int rev1 = globalIndexRC((*it).i1, batch);        
        int rev2 = globalIndexRC((*it).i2, batch);

        int set = std::min(findSet(uf, (*it).i1), findSet(uf, rev1));

        /*
        // skip components that are 10 times larger than number of samples
        if (uf._values[set] < -10 * (int)length(batch.contigFiles)) 
        {
            if (verbose && skipped.count(set) == 0)
            {
                std::cerr << "[" << time(0) << "] " << "WARNING: Skipping component of size " << (-1*uf._values[set]) << std::endl;
                skipped.insert(set);
            }
            skipped.insert((*it).i1);
            skipped.insert((*it).i2);
            skipped.insert(rev1);
            skipped.insert(rev2);
        }
        else
        */
        {
            components[set].alignedPairs.insert(*it);
            components[set].alignedPairs.insert(Pair<TSize>((*it).i2, (*it).i1));
            components[set].alignedPairs.insert(Pair<TSize>(rev1, rev2));
            components[set].alignedPairs.insert(Pair<TSize>(rev2, rev1));
        }
    }

    std::cerr << "[" << time(0) << "] " << "There are " << components.size() << " components." << std::endl;
    return skipped;
}

// --------------------------------------------------------------------------
// Function addSingletons()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
void
addSingletons(std::map<TSize, ContigComponent<TSeq> > & components,
              std::set<int> & skipped,
              UnionFind<int> & uf,
              int totalContigs)
{
    unsigned numSingletons = 0;
    for (int i = 0; i < totalContigs; ++i)
    {
        if (skipped.count(i) == 0 && components.count(i) == 0 && i == findSet(uf, i))
        {
            components[i];
            ++numSingletons;
        }
    }

    std::cerr << "[" << time(0) << "] " << "Added " << numSingletons << " singletons to components." << std::endl;
}

template<typename TSize, typename TSeq>
void
addSingletons(std::map<TSize, ContigComponent<TSeq> > & components,
              std::map<TSize, Contig<TSeq> > & contigs,
              std::set<int> & skipped,
              UnionFind<int> & uf,
              int totalContigs)
{
    unsigned numSingletons = 0;
    for (int i = 0; i < totalContigs; ++i)
    {
        if (contigs.count(i) > 0 && skipped.count(i) == 0 && components.count(i) == 0 && i == findSet(uf, i))
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
                       ContigBatch & batch,
                       bool verbose)
{
    typedef std::map<TSize, ContigComponent<TSequence> > TComponents;
    typedef typename std::set<Pair<TSize> >::iterator TPairsIterator;
    typedef typename TComponents::iterator TCompIterator;

    std::cerr << "[" << time(0) << "] " << "Reading and merging components files" << std::endl;

    // Initialize Union-Find data structure.
    UnionFind<int> uf;
    resize(uf, batch.contigsInTotal * 2);
    std::set<Pair<TSize> > alignedPairs;
    
    // Read the aligned pairs from input files and join sets.
    for (unsigned i = 0; i < length(componentFiles); ++i)
        if (readAlignedPairs(uf, alignedPairs, componentFiles[i], batch.contigsInTotal, verbose) != 0) return 1;

    // Convert union-find data structure to components.
    skipped = unionFindToComponents(components, uf, alignedPairs, batch, verbose);

    // Add singleton contigs to components (= those contigs that don't align to any other contig).
    addSingletons(components, skipped, uf, batch.contigsInTotal);

    // Keep only batch of components (erase all except every totalBatches'th component).
    unsigned total = totalBatches(batch);
    if (total != 1)
    {
        TCompIterator it = --components.end();
        TSize i = components.size();
        while (i > 0)
        {
            TCompIterator element = it;
            --it; --i;
            if (i % total != batch.number)
                components.erase(element);
        }
    }

    return 0;
}

#endif // #ifndef POPINS_MERGE_PARTITION_H_
