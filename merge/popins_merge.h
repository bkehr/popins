#ifndef POPINS_MERGE_H_
#define POPINS_MERGE_H_

#include <sstream>
#include <iomanip>

#include "contig_id.h"
#include "contig_structs.h"

#include "../command_line_parsing.h"

#include "partition.h"
#include "merge_seqs.h"


using namespace seqan;

// --------------------------------------------------------------------------
// Function readContigFile()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
bool
readContigFile(std::map<TSize, Contig<TSeq> > & contigs,
        CharString & filename,
        int min,
        int max,
        CharString & index)
{
    // open fasta file
    SeqFileIn stream(toCString(filename));

    unsigned numContigsBefore = length(contigs);

    // read records from fasta file
    int numContigs = 0;
    unsigned basepairs = 0;
    while (!atEnd(stream))
    {
        CharString id;
        TSeq seq;

        // read record
        readRecord(id, seq, stream);

        if (numContigs >= min && numContigs < max)
        {
            // append contig and contigId
            ContigId contigId(index, id, true);
            contigs[length(contigs)] = Contig<TSeq>(seq, contigId);
            basepairs += length(seq);
        }

        ++numContigs;
    }

    std::ostringstream msg;
    msg << "Loaded " << filename << ": " << basepairs << " bp in " << (length(contigs) - numContigsBefore) << " contigs.";
    printStatus(msg);

    return 0;
}

// --------------------------------------------------------------------------
// Function readContigs()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
bool
readContigs(std::map<TSize, Contig<TSeq> > & contigs,
        ContigBatch & batch)
{
    printStatus("Reading contig files");

    // open and read the files containing contigs of one individual one by one
    for (unsigned i = 0; i < length(batch.contigFiles); ++i)
    {
        unsigned contigsBefore = length(contigs);
        CharString index = formattedIndex(i, length(batch.contigFiles)); // get an identifier for the file
        if (readContigFile(contigs, batch.contigFiles[i], 0, maxValue<int>(), index) != 0) return 1;
        appendValue(batch.contigsPerFile, length(contigs) - contigsBefore);
    }

    batch.contigsInTotal = length(contigs);

    std::ostringstream msg;
    msg << "Total number of contigs: " << batch.contigsInTotal;
    printStatus(msg);

    return 0;
}

// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
bool
readContigSubset(std::map<TSize, Contig<TSeq> > & contigs,
        ContigBatch & batch)
{
    unsigned batchOffset = indexOffset(batch);
    unsigned batchEnd = std::min(batchOffset + getSize(batch), batch.contigsInTotal);

    std::ostringstream msg;
    msg << "Reading batch " << batch.number << "/" << totalBatches(batch) << " of " << batch.contigsInTotal << " contigs from contig files";
    printStatus(msg);

    // open and read the files containing contigs of one individual one by one

    unsigned contigCount = 0;
    unsigned i = 0;

    for (; contigCount + batch.contigsPerFile[i] <= batchOffset; ++i)
    {
        SEQAN_ASSERT_LT(i, length(batch.contigFiles));
        contigCount += batch.contigsPerFile[i];
    }

    for (; contigCount < batchEnd; ++i)
    {
        SEQAN_ASSERT_LT(i, length(batch.contigFiles));

        int fileMin = std::max((int)batchOffset - contigCount, 0u);
        int fileMax = std::min(contigCount + batch.contigsPerFile[i], batchEnd) - contigCount;

        CharString index = formattedIndex(i, length(batch.contigFiles)); // get an identifier for the file
        if (readContigFile(contigs, batch.contigFiles[i], fileMin, fileMax, index) != 0) return 1;

        contigCount += batch.contigsPerFile[i];
    }

    msg.str("");
    msg << "Number of contigs loaded: " << length(contigs);
    printStatus(msg);

    return 0;
}

// --------------------------------------------------------------------------

void
insertIndex(std::set<unsigned> & indices, unsigned index, int totalContigs)
{
    if (index >= (unsigned)totalContigs) index -= totalContigs;
    if (indices.count(index) == 0) indices.insert(index);
}

// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
bool
readContigs(std::map<TSize, Contig<TSeq> > & contigs,
        std::map<TSize, ContigComponent<TSeq> > & components,
        ContigBatch & batch)
{
    typedef std::map<TSize, ContigComponent<TSeq> > TComponents;
    typedef typename std::set<Pair<TSize> >::iterator TPairsIter;
    typedef typename std::set<unsigned>::iterator TListIter;

    std::ostringstream msg;
    msg << "Reading contigs for a batch of " << length(components) << " components";
    printStatus(msg);

    // Get a sorted list of all contig indices.
    std::set<unsigned> indices;
    for (typename TComponents::iterator it = components.begin(); it != components.end(); ++it)
    {
        insertIndex(indices, it->first, batch.contigsInTotal);
        TPairsIter pairsEnd = it->second.alignedPairs.end();
        for (TPairsIter pairsIt = it->second.alignedPairs.begin(); pairsIt != pairsEnd; ++pairsIt)
        {
            insertIndex(indices, (*pairsIt).i1, batch.contigsInTotal);
            insertIndex(indices, (*pairsIt).i2, batch.contigsInTotal);
        }
    }

    msg.str("");
    msg << "  ... " << indices.size() << " contigs";
    printStatus(msg);

    // Read the contigs with these indices.

    int fileIndex = 0;
    int firstIndexInCurrFile = 0;
    SeqFileIn stream(toCString(batch.contigFiles[fileIndex]));

    unsigned index = firstIndexInCurrFile;
    CharString sample = formattedIndex(fileIndex, length(batch.contigFiles));

    TListIter itEnd = indices.end();
    for (TListIter it = indices.begin(); it != itEnd; ++it)
    {
        // Open next file.
        if (firstIndexInCurrFile + batch.contigsPerFile[fileIndex] <= *it)
        {
            do
            {
                firstIndexInCurrFile += batch.contigsPerFile[fileIndex];
                ++fileIndex;
            }
            while (firstIndexInCurrFile + batch.contigsPerFile[fileIndex] <= *it);

            open(stream, toCString(batch.contigFiles[fileIndex]));
            index = firstIndexInCurrFile;
            sample = formattedIndex(fileIndex, length(batch.contigFiles));

            msg.str("");
            msg << "Opened " << batch.contigFiles[fileIndex];
            printStatus(msg);
        }

        Contig<TSeq> contig;
        contig.id.pn = sample;
        contig.id.orientation = true;

        // Advance file to *it
        while (index <= *it)
        {
            readRecord(contig.id.contigId, contig.seq, stream);
            ++index;
        }

        Contig<TSeq> contigRev;
        contigRev.id.pn = sample;
        contigRev.id.orientation = false;
        contigRev.id.contigId = contig.id.contigId;
        contigRev.seq = contig.seq;
        reverseComplement(contigRev.seq);

        // Add contig and its reverse complement to map.
        contigs[*it] = contig;
        contigs[globalIndexRC(*it, batch)] = contigRev;
    }

    SEQAN_ASSERT_EQ(length(contigs), indices.size());

    msg.str("");
    msg << "Number of contigs loaded: " << length(contigs);
    printStatus(msg);

    return 0;
}

// --------------------------------------------------------------------------
// Function readInputFiles()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
bool
readInputFiles(std::map<TSize, Contig<TSeq> > & contigs,
        std::map<TSize, ContigComponent<TSeq> > & components,
        std::set<int> & skipped,
        ContigBatch & batch,
        MergingOptions & options)
{
    if (length(options.componentFiles) == 0) // -c option is not set
    {
        if (batch.batchesInTotal == 1) // -i and -b options are not set -> read all contigs
        {
            if (readContigs(contigs, batch) != 0) return 1;
        }
        else                     // -i and -b options are set -> read subset of contigs
        {               
            if (countContigs(batch)) return 1;
            if (readContigSubset(contigs, batch) != 0) return 1;
        }
    }
    else // -c option is set
    {
        if (batch.batchesInTotal == 1) // -i and -b options are not set -> read all contigs
        {
            if (readContigs(contigs, batch) != 0) return 1;
            if (readAndMergeComponents(components, skipped, options.componentFiles, batch) != 0) return 1;  // --> popins_merge_partition.h
        }
        else                      // -i and -b options are set -> read contigs for this batch of components
        {
            if (countContigs(batch)) return 1;
            if (readAndMergeComponents(components, skipped, options.componentFiles, batch) != 0) return 1;  // --> popins_merge_partition.h
            if (readContigs(contigs, components, batch) != 0) return 1; // TODO: Exclude singletons of lowEntropy.
        }
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function addReverseComplementContigs()
// --------------------------------------------------------------------------

template<typename TSize, typename TSeq>
void
addReverseComplementContigs(std::map<TSize, Contig<TSeq> > & contigs, ContigBatch & batch)
{
    typedef typename std::map<TSize, Contig<TSeq> >::iterator TIter;

    typedef Pair<TSize, Contig<TSeq> > TPair;
    String<TPair> revContigs;
    TIter itEnd = contigs.end();
    for (TIter it = contigs.begin(); it != itEnd; ++it)
    {

        TSeq revSeq = (it->second).seq;
        reverseComplement(revSeq);

        SEQAN_ASSERT((it->second).id.orientation);
        ContigId revId = (it->second).id;
        revId.orientation = false;

        appendValue(revContigs, TPair(globalIndexRC(it->first, batch), Contig<TSeq>(revSeq, revId)));
    }

    for (unsigned i = 0; i < length(revContigs); ++i)
        contigs[revContigs[i].i1] = revContigs[i].i2;
}

// --------------------------------------------------------------------------
// Function writeSkipped()
// --------------------------------------------------------------------------

template<typename TStream, typename TSize, typename TSeq>
void
writeSkipped(TStream & stream, std::map<TSize, Contig<TSeq> > & contigs, std::set<int> skipped)
{
    typename std::set<int>::iterator itEnd = skipped.end();
    for (typename std::set<int>::iterator it = skipped.begin(); it != itEnd; ++it)
    {
        if (!contigs[*it].id.orientation) continue;
        stream << ">" << contigs[*it].id << " " << "(large component)" << std::endl;
        stream << contigs[*it].seq << std::endl;
    }
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
    std::map<TSize, Contig<TSequence> > contigs;
    std::map<TSize, ContigComponent<TSequence> > components;
    std::set<int> skipped;

    // Open the output files.
    options.outputStream.open(toCString(options.outputFile), std::ios_base::out);
    if (!options.outputStream.is_open())
    {
        std::cerr << "ERROR: Could not open output file " << options.outputFile << std::endl;
        return 1;
    }
    if (options.skippedFile != "")
    {
        options.skippedStream.open(toCString(options.skippedFile), std::ios_base::out);
        if (!options.skippedStream.is_open())
        {
            std::cerr << "ERROR: Could not open output file for skipped contigs " << options.skippedFile << std::endl;
            return 1;
        }
    }

    // Reading of input files (files of contigs, and files of components if -c option is set).
    ContigBatch batch(options.contigFiles, options.batches, options.batchIndex);
    if (readFileNames(batch.contigFiles, batch.contigsPerFile) != 0) return 1;
    if (readInputFiles(contigs, components, skipped, batch, options) != 0) return 1;
    if (filterByEntropy(contigs, options) != 0) return 1;
    addReverseComplementContigs(contigs, batch);

    // Initialize Union-Find data structure for partitioning.
    UnionFind<int> uf;
    resize(uf, 2 * batch.contigsInTotal);
    std::set<Pair<TSize> > alignedPairs;

    // PARTITIONING into components (if -c option is not set)                               --> popins_merge_partition.h
    if (length(options.componentFiles) == 0)
    {
        if (partitionContigs(uf, alignedPairs, contigs, batch, options) != 0) 
            return 1;

        // Write aligned pairs or convert union find data structure into set of components.
        if (options.batches != 1)
            writeAlignedPairs(options.outputStream, alignedPairs);
        else
        {
            skipped = unionFindToComponents(components, uf, alignedPairs, batch);
            addSingletons(components, contigs, skipped, uf, batch.contigsInTotal);
        }
    }

    // Write contigs that were skipped because they form too large components.
    if (options.skippedFile != "")
        writeSkipped(options.skippedStream, contigs, skipped);

    // SUPERCONTIG CONSTRUCTION (if -c option is set OR -i and -b options are not set)           --> popins_merge_seqs.h
    if (length(options.componentFiles) != 0 || options.batches == 1)
    {
        if (length(contigs) != 0)
            constructSupercontigs(components, contigs, options);
        else
            constructSupercontigs(components, contigs, options);
    }

    return 0;
}
#endif // #ifndef POPINS_MERGE_H_