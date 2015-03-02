#ifndef POPINS_MERGE_H_
#define POPINS_MERGE_H_

#include <sstream>
#include <iomanip>

#include "contig_id.h"
#include "contig_component.h"

#include "popins_clp.h"

#include "popins_merge_partition.h"
#include "popins_merge_seqs.h"


using namespace seqan;

// --------------------------------------------------------------------------

int
countContigs(String<CharString> & filenames, String<unsigned> & contigsPerFile)
{
    unsigned count = 0;
    
    // Use counts specified by user.
    for (unsigned i = 0; i < length(contigsPerFile); ++i)
        count += contigsPerFile[i];

    // Iterate over all files to count fasta records.
    if (count == 0)
    {
        CharString id;
        Dna5String contig;

        unsigned prevCount = 0;
        for (unsigned i = 0; i < length(filenames); ++i)
        {
            // Open file.
            SequenceStream stream(toCString(filenames[i]));
            if (!isGood(stream))
            {
                std::cerr << "ERROR: Could not open " << filenames[i] << " as fasta file." << std::endl;
                return -1;
            }
            
            // Count records in file.
            while (!atEnd(stream))
            {
                if (readRecord(id, contig, stream))
                {
                    std::cerr << "ERROR: Could not read fasta record from " << filenames[i] << std::endl;
                    return -1;
                }
                ++count;
            }
            appendValue(contigsPerFile, count - prevCount); 
            prevCount = count;
        }
    }
    
    return (int)count;
}

// --------------------------------------------------------------------------
// Function readContigs()
// --------------------------------------------------------------------------

template<typename TSeq>
bool
readContigFile(StringSet<TSeq> & contigs,
               StringSet<ContigId> & contigIds,
               CharString & filename,
               int min,
               int max,
               CharString & index,
               bool verbose)
{
    // open fasta file
    SequenceStream stream(toCString(filename));
    if (!isGood(stream))
    {
        std::cerr << "ERROR: Could not open " << filename << " as fasta file." << std::endl;
        return 1;
    }

    unsigned numContigsBefore = length(contigs);

    // read records from fasta file
    int numContigs = 0;
    unsigned basepairs = 0;
    while (!atEnd(stream))
    {
        CharString id;
        TSeq contig;

        // read record
        if (readRecord(id, contig, stream))
        {
            std::cerr << "ERROR: Could not read fasta record from " << filename << std::endl;
            return 1;
        }

        if (numContigs >= min && numContigs < max)
        {
            // append contig and contigId
            appendValue(contigs, contig);
            appendValue(contigIds, ContigId(index, id, true));
            basepairs += length(contig);
        }

        ++numContigs;
    }

    if (verbose) std::cerr << "[" << time(0) << "] " << "Loaded " << filename << ": " << basepairs << " bp in "
                                                     << (length(contigs) - numContigsBefore) << " contigs." << std::endl;

    return 0;
}

// --------------------------------------------------------------------------


template<typename TSeq>
int
readContigs(StringSet<TSeq> & contigs,
            StringSet<ContigId> & contigIds,
            String<CharString> & filenames,
            bool verbose)
{
    std::cerr << "[" << time(0) << "] " << "Reading contig files" << std::endl;

    // open and read the files containing contigs of one individual one by one
    for (unsigned i = 0; i < length(filenames); ++i)
    {
        CharString index = formattedIndex(i, length(filenames)); // get an identifier for the file
        if (readContigFile(contigs, contigIds, filenames[i], 0, maxValue<int>(), index, verbose) != 0) return -1;
    }

    std::cerr << "[" << time(0) << "] " << "Total number of contigs: " << length(contigs) << std::endl;

    return length(contigs);
}

// --------------------------------------------------------------------------

template<typename TSeq>
int
readContigs(StringSet<TSeq> & contigs,
            StringSet<ContigId> & contigIds,
            unsigned & batchOffset,
            String<CharString> & filenames,
            String<unsigned> & contigsPerFile,
            int batchIndex,
            int batches,
            bool verbose)
{
    int totalContigCount = countContigs(filenames, contigsPerFile);
    if (totalContigCount == -1) return -1;

    unsigned batchSize = (totalContigCount + batches - 1) / batches;
    batchOffset = batchSize * batchIndex;
    unsigned batchEnd = std::min(batchOffset + batchSize, (unsigned)totalContigCount);

    std::cerr << "[" << time(0) << "] " << "Reading batch " << batchIndex << "/" << batches << " of "
                                        << totalContigCount << " contigs from contig files" << std::endl;
    
    // open and read the files containing contigs of one individual one by one

    unsigned contigCount = 0;
    unsigned i = 0;

    for (; contigCount + contigsPerFile[i] <= batchOffset; ++i)
    {
        SEQAN_ASSERT_LT(i, length(filenames));
        contigCount += contigsPerFile[i];
    }
    
    for (; contigCount < batchEnd; ++i)
    {
        SEQAN_ASSERT_LT(i, length(filenames));

        int fileMin = std::max((int)batchOffset - contigCount, 0u);
        int fileMax = std::min(contigCount + contigsPerFile[i], batchEnd) - contigCount;

        CharString index = formattedIndex(i, length(filenames)); // get an identifier for the file
        if (readContigFile(contigs, contigIds, filenames[i], fileMin, fileMax, index, verbose) != 0) return -1;

        contigCount += contigsPerFile[i];
    }

    std::cerr << "[" << time(0) << "] " << "Number of contigs loaded: " << length(contigs) << std::endl;

    return totalContigCount;
}

// --------------------------------------------------------------------------

void
insertIndex(std::set<unsigned> & indices, unsigned index, int totalContigs)
{
    if (index >= (unsigned)totalContigs) index -= totalContigs;
    if (indices.count(index) == 0) indices.insert(index);
}

// --------------------------------------------------------------------------

template<typename TSeq, typename TSize>
bool
readContigs(std::map<TSize, TSeq> & contigs,
            std::map<TSize, ContigId> & contigIds,
            std::map<TSize, ContigComponent<TSeq> > & components,
            String<CharString> & filenames,
            String<unsigned> & contigsPerFile,
            int totalContigs,
            bool verbose)
{
    typedef std::map<TSize, ContigComponent<TSeq> > TComponents;
    typedef typename std::set<Pair<TSize> >::iterator TPairsIter;
    typedef typename std::set<unsigned>::iterator TListIter;

    std::cerr << "[" << time(0) << "] " << "Reading contigs for a batch of " << length(components) << " components" << std::flush;

    // Get a sorted list of all contig indices.
    std::set<unsigned> indices;
    for (typename TComponents::iterator it = components.begin(); it != components.end(); ++it)
    {
        TPairsIter pairsEnd = it->second.alignedPairs.end();
        for (TPairsIter pairsIt = it->second.alignedPairs.begin(); pairsIt != pairsEnd; ++pairsIt)
        {
            insertIndex(indices, (*pairsIt).i1, totalContigs);
            insertIndex(indices, (*pairsIt).i2, totalContigs);
        }
    }
    
    std::cerr << " (" << indices.size() << " contigs)" << std::endl;
    
    // Read the contigs with these indices.
    
    int fileIndex = 0;
    int firstIndexInCurrFile = 0;
    SequenceStream stream(toCString(filenames[fileIndex]));
    if (!isGood(stream))
    {
        std::cerr << "ERROR: Could not open " << filenames[fileIndex] << " as fasta file." << std::endl;
        return 1;
    }

    unsigned index = firstIndexInCurrFile;
    CharString sample = formattedIndex(fileIndex, length(filenames));

    TListIter itEnd = indices.end();
    for (TListIter it = indices.begin(); it != itEnd; ++it)
    {
        // Open next file.
        if (firstIndexInCurrFile + contigsPerFile[fileIndex] <= *it)
        {
            do
            {
                firstIndexInCurrFile += contigsPerFile[fileIndex];
                ++fileIndex;
            }
            while (firstIndexInCurrFile + contigsPerFile[fileIndex] < *it);

            open(stream, toCString(filenames[fileIndex]));
            if (!isGood(stream))
            {
                std::cerr << "ERROR: Could not open " << filenames[fileIndex] << " as fasta file." << std::endl;
                return 1;
            }
            index = firstIndexInCurrFile;
            sample = formattedIndex(fileIndex, length(filenames));

            if (verbose) std::cerr  << "[" << time(0) << "] " << "Opened " << filenames[fileIndex] << std::endl;
        }
    
        TSeq contig, contigRev;
        ContigId id, idRev;
        id.pn = sample;
        idRev.pn = sample;
        id.orientation = true;
        idRev.orientation = false;

        // Advance file to *it
        while (index <= *it)
        {
            if (readRecord(id.contigId, contig, stream))
            {
                std::cerr << "ERROR: Could not read fasta record from " << filenames[fileIndex] << std::endl;
                return 1;
            }
            ++index;
        }

        idRev.contigId = id.contigId;
        contigRev = contig;
        reverseComplement(contigRev);
        
        unsigned i = *it;
        contigIds[i] = id;
        contigs[i] = contig;

        i += totalContigs;
        contigIds[i] = idRev;
        contigs[i] = contigRev;
    }

    SEQAN_ASSERT_EQ(length(contigs), indices.size());
    std::cerr << "[" << time(0) << "] " << "Number of contigs loaded: " << length(contigs) << std::endl;

    return 0;
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
// Function readInputFiles()
// --------------------------------------------------------------------------

template<typename TSeq, typename TSize>
int
readInputFiles(StringSet<TSeq, Owner<> > & contigs,
               StringSet<ContigId, Owner<> > & contigIds,
               std::map<TSize, TSeq> & contigsMap,
               std::map<TSize, ContigId> & contigIdsMap,
               unsigned & batchOffset,
               std::map<TSize, ContigComponent<TSeq> > & components,
               MergingOptions & options)
{
    int numContigs = 0;

    if (length(options.componentFiles) == 0) // -c option is not set
    {
        if (options.batches == 1) // -i and -b options are not set -> read all contigs
            numContigs = readContigs(contigs, contigIds, options.contigFiles, options.verbose);
        else                      // -i and -b options are set -> read subset of contigs
            numContigs = readContigs(contigs, contigIds, batchOffset, options.contigFiles, options.contigsPerFile,
                                     options.batchIndex, options.batches, options.verbose);
    }

    else // -c option is set
    {
        if (options.batches == 1) // -i and -b options are not set -> read all contigs
        {
            numContigs = readContigs(contigs, contigIds, options.contigFiles, options.verbose);
            if (numContigs == -1) return -1;
            if (readAndMergeComponents(components, options.componentFiles, numContigs,
                                       options.batchIndex, options.batches, options.verbose) != 0) return -1;
        }
        else                      // -i and -b options are set -> read contigs for this batch of components
        {
            int numContigs = countContigs(options.contigFiles, options.contigsPerFile);
            if (numContigs == -1) return -1;
            if (readAndMergeComponents(components, options.componentFiles, numContigs,  // --> popins_merge_partition.h
                                       options.batchIndex, options.batches, options.verbose) != 0) return -1;
            if (readContigs(contigsMap, contigIdsMap, components,
                            options.contigFiles, options.contigsPerFile, numContigs, options.verbose) != 0) return -1;
        }
    }

    return numContigs;
}

// ==========================================================================
// Function popins_merge()
// ==========================================================================

int popins_merge(int argc, char const ** argv)
{
    typedef Dna5String TSequence;
    typedef Size<TSequence>::Type TSize;

    // Parse the command line to get option values.
    MergingOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    // Containers for contigs, contig ids, and components.
    StringSet<TSequence, Owner<> > contigs;
    StringSet<ContigId, Owner<> > contigIds;
    std::map<TSize, TSequence> contigsMap; // only needed for supercontig construction in batches
    std::map<TSize, ContigId> contigIdsMap; // only needed for superconitg construction in batches
    std::map<TSize, ContigComponent<TSequence> > components;

    // Reading of input files (files of contigs, and files of components if -c option is set).
    unsigned batchOffset = 0;
    int totalContigs = readInputFiles(contigs, contigIds, contigsMap, contigIdsMap, batchOffset, components, options);
    if (totalContigs == -1) return 1;
    addReverseComplementContigs(contigs, contigIds);

    // Prepare the output file.
    options.outputStream.open(toCString(options.outputFile), std::ios_base::out);
    if (!options.outputStream.is_open())
    {
        std::cerr << "ERROR: Could not open output file " << options.outputFile << std::endl;
        return 1;
    }

    // Initialize Union-Find data structure for partitioning.
    UnionFind<int> uf;
    resize(uf, 2*totalContigs);
    std::set<Pair<TSize> > alignedPairs;

    // PARTITIONING into components (if -c option is not set).
    if (length(options.componentFiles) == 0)
    {
        if (partitionContigs(uf, alignedPairs, contigs, contigIds, totalContigs, batchOffset, options) != 0) // --> popins_merge_partition.h
            return 1;
        
        // Write aligned pairs or convert union find data structure into set of components.
        if (options.batches != 1)
            writeAlignedPairs(options.outputStream, alignedPairs); // --> popins_merge_partition.h
        else
        {
            unionFindToComponents(components, uf, alignedPairs, totalContigs); // --> popins_merge_partition.h
            addSingletons(components, uf, totalContigs); // --> popins_merge_partition.h
        }
    }

    // SUPERCONTIG CONSTRUCTION (if -c option is set OR -i and -b options are not set).
    if (length(options.componentFiles) != 0 || options.batches == 1)
    {
        if (length(contigs) != 0)
            constructSupercontigs(components, contigs, contigIds, options);  // --> popins_merge_seqs.h
        else
            constructSupercontigs(components, contigsMap, contigIdsMap, options);  // --> popins_merge_seqs.h
    }

    return 0;
}
#endif // #ifndef POPINS_MERGE_H_
