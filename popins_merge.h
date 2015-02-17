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

// --------------------------------------------------------------------------
// Function computeOrReadContigComponents()
// --------------------------------------------------------------------------

template<typename TSequence>
bool
computeOrReadContigComponents(std::map<ContigId, ContigComponent<TSequence> > & components,
                              StringSet<TSequence, Owner<> > & contigs,
                              StringSet<ContigId, Owner<> > & contigIds,
                              MergingOptions & options)
{
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

            // --- PARTITION CONTIGS into sets of similar contigs (components) by computing pairwise alignments. ---
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

            // --- PARTITION CONTIGS into sets of similar contigs (components) by computing pairwise alignments. ---
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

    return 0;
}

// --------------------------------------------------------------------------

template<typename TStream, typename TSeq, typename TSpec>
void
writeSupercontigs(TStream & outputStream, String<TSeq> & mergedSeqs, StringSet<TSeq, TSpec> & contigs, unsigned pos)
{
    typedef typename Size<TSeq>::Type TSize;

    if (length(mergedSeqs) <= 25)
    {
        for (TSize i = 0; i < length(mergedSeqs); ++i)
        {
            outputStream << ">COMPONENT_" << pos << "_" << char('a'+i)
                         << "_length_" << length(mergedSeqs[i])
                         << "_size_" << length(contigs) << std::endl;
            outputStream << mergedSeqs[i] << std::endl;
        }
    }
    else
    {
        for (TSize i = 0; i < length(mergedSeqs); ++i)
        {
            outputStream << ">COMPONENT_" << pos << "_" << char('a'+i/26) << char('a'+i%26)
                         << "_length_" << length(mergedSeqs[i])
                         << "_size_" << length(contigs) << std::endl;
            outputStream << mergedSeqs[i] << std::endl;
        }
    }
}

// --------------------------------------------------------------------------
// Function constructSupercontigs()
// --------------------------------------------------------------------------

template<typename TSequence>
void
constructSupercontigs(std::map<ContigId, ContigComponent<TSequence> > & components,
                      StringSet<TSequence, Owner<> > & contigs,
                      StringSet<ContigId, Owner<> > & contigIds,
                      MergingOptions & options)
{
    typedef std::map<ContigId, ContigComponent<TSequence> > TComponents;
    
    if (options.verbose) std::cerr << "[" << time(0) << "] " << "Constructing supercontigs" << std::endl;

    unsigned numSingleton = 0;
    unsigned numBranching = 0;
    unsigned numVeryBranching = 0;
    unsigned numTooLarge = 0;
    
    // Iterate over the set of components.
    unsigned pos = 0;
    for (typename TComponents::iterator it = components.begin(); it != components.end(); ++it)
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

        // --- MERGE CONTIGS OF THE COMPONENT ---
        String<TSequence> mergedSeqs;
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

        // Output the supercontig.
        writeSupercontigs(options.outputStream, mergedSeqs, component.contigs, pos);

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
}

// ==========================================================================
// Function popins_merge()
// ==========================================================================

// TODO clp: option for batch size for merging components
// TODO clp: option for batch number for merging components

int popins_merge(int argc, char const ** argv)
{
    typedef Dna5String TSequence;

    // Parse the command line to get option values.
    MergingOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    // Read contigs from file.
    StringSet<TSequence, Owner<> > contigs;
    StringSet<ContigId, Owner<> > contigIds;
    if (readContigs(contigs, contigIds, options.contigFiles, options.verbose) != 0) return 1;

    // --- 1) PARTITIONING ---
    std::map<ContigId, ContigComponent<TSequence> > components;
    if (computeOrReadContigComponents(components, contigs, contigIds, options) != 0) return 1;

    // Prepare the output file.
    options.outputStream.open(toCString(options.outputFile), std::ios_base::out);
    if (!options.outputStream.is_open())
    {
        std::cerr << "ERROR: Could not open output file " << options.outputFile << std::endl;
        return 1;
    }
    
    // --- 2) SUPERCONTIG CONSTRUCTION ---
    constructSupercontigs(components, contigs, contigIds, options);

    return 0;
}
#endif // #ifndef POPINS_MERGE_H_
