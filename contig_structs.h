#ifndef CONTIG_STRUCTS_H_
#define CONTIG_STRUCTS_H_

#include<seqan/sequence.h>

#include "contig_id.h"

using namespace seqan;

// ==========================================================================
// struct Contig
// ==========================================================================

template<typename TSeq>
struct Contig
{
    TSeq seq;
    ContigId id;

    Contig() {}

    Contig(TSeq & s, ContigId & i) :
        seq(s), id(i)
    {}
};

// ==========================================================================
// struct ContigComponent
// ==========================================================================

template<typename TSeq>
struct ContigComponent
{
    typedef typename Size<TSeq>::Type TSize;

    StringSet<Contig<TSeq>, Dependent<> > contigs;
    std::set<Pair<TSize> > alignedPairs;

    ContigComponent()
    {}
};

// --------------------------------------------------------------------------

template<typename TSeq>
void
clear(ContigComponent<TSeq> & c)
{
    clear(c.contigs);
    clear(c.alignedPairs);
}

// ==========================================================================
// struct ContigBatch
// ==========================================================================

struct ContigBatch
{
    // General info about the set of contigs
    String<CharString> contigFiles;
    String<unsigned> contigsPerFile;
    unsigned contigsInTotal;

    // Specific info for this batch of contigs
    unsigned number;
    unsigned batchesInTotal;
    int size; // of all batches
    int actualSize; // of this specific batch (last batch may be smaller than the other batches)
    int offset; // start index of this batch

    ContigBatch() :
        number(0), batchesInTotal(1), size(-1), actualSize(-1), offset(-1)
    {}

    ContigBatch(String<CharString> & files, unsigned batches, unsigned batchIndex) :
        contigFiles(files), number(batchIndex), batchesInTotal(batches), size(-1), actualSize(-1), offset(-1)
    {}
};

// --------------------------------------------------------------------------

int getSize(ContigBatch & batch)
{
    if (batch.size == -1)
        batch.size = (batch.contigsInTotal + batch.batchesInTotal - 1) / batch.batchesInTotal;
    return batch.size;
}

// --------------------------------------------------------------------------

int indexOffset(ContigBatch & batch)
{
    if (batch.offset == -1)
        batch.offset = getSize(batch) * batch.number;
    return batch.offset;
}

// --------------------------------------------------------------------------

unsigned globalIndexRC(unsigned i, ContigBatch & batch)
{
    if (i < batch.contigsInTotal)
        return i += batch.contigsInTotal;
    else
        return i -= batch.contigsInTotal;
}

// --------------------------------------------------------------------------

int batchSize(ContigBatch & batch)
{
    if (batch.actualSize == -1)
    {
        if (unsigned(indexOffset(batch) + getSize(batch)) < batch.contigsInTotal)
            batch.actualSize = getSize(batch);
        else
            batch.actualSize = batch.contigsInTotal - indexOffset(batch);
    }
    return batch.actualSize;
}

// --------------------------------------------------------------------------

int totalBatches(ContigBatch & batch)
{
    return batch.batchesInTotal;
}

// --------------------------------------------------------------------------
// Function countContigs()
// --------------------------------------------------------------------------

bool countContigs(ContigBatch & batch)
{
    unsigned count = 0;

    // Case A) Use counts specified by user.
    for (unsigned i = 0; i < length(batch.contigsPerFile); ++i)
        count += batch.contigsPerFile[i];

    // Case B) Iterate over all files to count fasta records.
    if (count == 0)
    {
        CharString id;
        Dna5String contig;

        unsigned prevCount = 0;
        for (unsigned i = 0; i < length(batch.contigFiles); ++i)
        {
            // Open file.
            SequenceStream stream(toCString(batch.contigFiles[i]));
            if (!isGood(stream))
            {
                std::cerr << "ERROR: Could not open " << batch.contigFiles[i] << " as fasta file." << std::endl;
                return 1;
            }

            // Count records in file.
            while (!atEnd(stream))
            {
                if (readRecord(id, contig, stream))
                {
                    std::cerr << "ERROR: Could not read fasta record from " << batch.contigFiles[i] << std::endl;
                    return 1;
                }
                ++count;
            }
            appendValue(batch.contigsPerFile, count - prevCount); 
            prevCount = count;
        }
    }

    batch.contigsInTotal = count;

    return 0;
}

// --------------------------------------------------------------------------

#endif  // #ifndef CONTIG_STRUCTS_H_
