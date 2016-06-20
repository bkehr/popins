#ifndef POPINS_MERGE_H_
#define POPINS_MERGE_H_

#include <sstream>
#include <iomanip>

#include "contig_structs.h"
#include "../popins_utils.h"
#include "../command_line_parsing.h"

#include "partition.h"
#include "merge_seqs.h"


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

// --------------------------------------------------------------------------
// Function readContigFile()
// --------------------------------------------------------------------------

template<typename TSeq>
bool
readContigFile(String<Contig<TSeq> > & contigs,
        CharString & filename,
        CharString & sampleId,
      MergingOptions & options)
{
    // Open the FASTA file.
    SeqFileIn stream(toCString(filename));

    unsigned numContigsBefore = length(contigs);

    // Read the records from FASTA file.
    unsigned basepairs = 0, numFiltered = 0;
    while (!atEnd(stream))
    {
        CharString contigName;
        TSeq seq;
        readRecord(contigName, seq, stream);
        ContigId contigId(sampleId, contigName, true);

        // Entropy calculation
        double entropy = averageEntropy(seq);

        if (entropy >= options.minEntropy)
        {
           // Append contig and contigId.
           appendValue(contigs, Contig<TSeq>(seq, contigId));
           basepairs += length(seq);
        }
        else if (options.skippedFile != "")
        {
           // Output contig as skipped.
            options.skippedStream << ">" << contigId << " entropy: " << entropy << std::endl;
            options.skippedStream << seq << std::endl;
            ++numFiltered;
        }
    }

    std::ostringstream msg;
    msg << "Loaded " << filename << ": " << basepairs << " bp in " << (length(contigs) - numContigsBefore) << " contigs";
    if (numFiltered > 0)
        msg << " (additional " << numFiltered << " contigs failed the entropy filter)";
    printStatus(msg);

    return 0;
}

// --------------------------------------------------------------------------
// Function readInputFiles()
// --------------------------------------------------------------------------

template<typename TSeq>
bool
readInputFiles(String<Contig<TSeq> > & contigs, MergingOptions & options)
{
   // List all files <prefix>/*/contigs.fa
   CharString filename = "contigs.fa";
   String<Pair<CharString> > contigFiles = listFiles(options.prefix, filename);

   // Read the contig files.
   for (unsigned i = 0; i < length(contigFiles); ++i)
      if (readContigFile(contigs, contigFiles[i].i2, contigFiles[i].i1, options) != 0)
         return 1;

    return 0;
}

// --------------------------------------------------------------------------
// Function addReverseComplementContigs()
// --------------------------------------------------------------------------

template<typename TSeq>
void
addReverseComplementContigs(String<Contig<TSeq> > & contigs)
{
    unsigned fwdContigCount = length(contigs);

    for (unsigned i = 0; i < fwdContigCount; ++i)
    {
        TSeq revSeq = contigs[i].seq;
        reverseComplement(revSeq);

        SEQAN_ASSERT(contigs[i].id.orientation);
        ContigId revId = contigs[i].id;
        revId.orientation = false;

        appendValue(contigs, Contig<TSeq>(revSeq, revId));
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
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Containers for contigs, contig ids, and components.
    String<Contig<TSequence> > contigs;
    std::map<TSize, ContigComponent<TSequence> > components;
    std::set<int> skipped;

    // Open the output files.
    options.outputStream.open(toCString(options.outputFile), std::ios_base::out);
    if (!options.outputStream.is_open())
    {
        std::cerr << "ERROR: Could not open output file " << options.outputFile << std::endl;
        return 7;
    }
    if (options.skippedFile != "")
    {
        options.skippedStream.open(toCString(options.skippedFile), std::ios_base::out);
        if (!options.skippedStream.is_open())
        {
            std::cerr << "ERROR: Could not open output file for skipped contigs " << options.skippedFile << std::endl;
            return 7;
        }
    }

    // Read and filter the contigs and add reverse complements.
    if (readInputFiles(contigs, options) != 0)
       return 7;
    addReverseComplementContigs(contigs);

    // PARTITIONING into components      --> partition.h
    UnionFind<int> uf;
    resize(uf, 2 * length(contigs));
    std::set<Pair<TSize> > alignedPairs;
    if (partitionContigs(uf, alignedPairs, contigs, options) != 0)
        return 7;

    unionFindToComponents(components, uf, alignedPairs, length(contigs)/2);
    addSingletons(components, contigs, uf);

    // SUPERCONTIG CONSTRUCTION           --> merge_seqs.h
    constructSupercontigs(components, contigs, options);

    return 0;
}
#endif // #ifndef POPINS_MERGE_H_
