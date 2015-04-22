#include <sstream>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/vcf_io.h>

#include "popins_location.h"
#include "align_split.h"

#ifndef POPINS_PLACE_H_
#define POPINS_PLACE_H_

using namespace seqan;

template<typename TSeq>
bool
loadSequences(std::map<CharString, TSeq> & seqs,
              CharString const & filename)
{
    // Open fasta file.
    SequenceStream stream(toCString(filename));
    if (!isGood(stream))
    {
        std::cerr << "ERROR: Could not open " << filename << std::endl;
        return 1;
    }

    // Read records from file and append to seqs.
    while (!atEnd(stream))
    {
        CharString id;
        TSeq seq;

        // read record
        if (readRecord(id, seq, stream) != 0)
        {
            std::cerr << "ERROR: Could not read fasta record from " << filename << std::endl;
            return 1;
        }
        
        unsigned i = 0;
        for (; i < length(id); ++i) if (id[i] == ' ') break;
        id = prefix(id, i);
        
        seqs[id] = seq;
    }

    return 0;
}

// ==========================================================================
// Function concatenateArtificialReferences()
// ==========================================================================

template<typename TSeq>
bool
artificialReferences(String<Pair<TSeq> > & concatRefs,
                                String<Location> & locations,
                                CharString & referenceFile,
                                CharString & contigFile,
                                unsigned readLength,
                                unsigned maxInsertSize)
{
    typedef String<Location> TLocations;
    typedef Iterator<TLocations>::Type TLocsIter;
    
    // Load the genome and the contig file.
    std::map<CharString, TSeq> genome, contigs;
    if (loadSequences(genome, referenceFile) != 0) return 1;
    if (loadSequences(contigs, contigFile) != 0) return 1;

    // Append reference sequences.
    TLocsIter locsEnd = end(locations);
    for (TLocsIter loc = begin(locations); loc != locsEnd; ++loc)
    {
        if (contigs.count((*loc).contig) == 0)
        {
            std::cerr << "ERROR: Could not find record for " << (*loc).contig << " in contigs file " << contigFile << std::endl;
            return 1;
        }
        if (genome.count((*loc).chr) == 0)
        {
            std::cerr << "ERROR: Could not find record for " << (*loc).chr << " in file of reference genome " << referenceFile << std::endl;
            return 1;
        }

        TSeq contig = contigs[(*loc).contig];
        if ((*loc).contigOri == (*loc).chrOri) reverseComplement(contig);

        if ((*loc).chrOri)
        {
            TSeq chrInfix = infix(genome[(*loc).chr], (*loc).chrStart+readLength, (*loc).chrEnd+maxInsertSize);
            appendValue(concatRefs, Pair<TSeq>(chrInfix, contig));
        }
        else
        {
            TSeq chrInfix = infix(genome[(*loc).chr], (*loc).chrStart-maxInsertSize, (*loc).chrEnd);
            appendValue(concatRefs, Pair<TSeq>(contig, chrInfix));
        }
    }

    return 0;
}

// ==========================================================================

bool
openBamLoadBai(BamStream & bamStream, BamIndex<Bai> & bamIndex, CharString & filename)
{
    // Open bam file.
    if (open(bamStream, toCString(filename)) != 0)
    {
        std::cerr << "ERROR: Could not open bam file " << filename << std::endl;
        return 1;
    }

    // Load the bam index.
    CharString baiFile = filename;
    baiFile += ".bai";
    if (read(bamIndex, toCString(baiFile)) != 0)
    {
        std::cerr << "ERROR: Could not read BAI index file " << baiFile << std::endl;
        return 1;
    }
    
    return 0;
}

// ==========================================================================

template<typename TPos>
bool
jumpToLocation(TPos & locStart, TPos & locEnd, BamStream & bamStream, BamIndex<Bai> & bamIndex, Location & loc, unsigned readLength, unsigned maxInsertSize)
{
    int rID = 0;
    BamIOContext<StringSet<CharString> > context = bamStream.bamIOContext;
    getIdByName(nameStore(context), loc.chr, rID, nameStoreCache(context));

    if (loc.chrOri)
    {
        locStart = loc.chrStart + readLength;
        locEnd = loc.chrEnd + maxInsertSize;
    }
    else
    {
        locStart = loc.chrStart - maxInsertSize;
        locEnd = loc.chrEnd;
    }

    bool hasAlignments;
    jumpToRegion(bamStream, hasAlignments, rID, locStart, locEnd, bamIndex);
    return hasAlignments;
}

// ==========================================================================
// Function splitAlign()
// ==========================================================================

template<typename TPos, typename TSequence>
bool
splitAlign(Pair<TPos> & refPos, Pair<TSequence> & ref, BamAlignmentRecord & record, bool chrOri)
{
    TSequence read = record.seq;

    Align<TSequence> left;
    resize(rows(left), 2);
    setSource(row(left, 0), read);
    setSource(row(left, 1), ref.i1);

    Align<TSequence> right;
    resize(rows(right), 2);
    setSource(row(right, 0), read);
    setSource(row(right, 1), ref.i2);

    TSequence revRead = read;
    reverseComplement(revRead);
    
    Align<TSequence> leftRev;
    resize(rows(leftRev), 2);
    setSource(row(leftRev, 0), revRead);
    setSource(row(leftRev, 1), ref.i1);

    Align<TSequence> rightRev;
    resize(rows(rightRev), 2);
    setSource(row(rightRev, 0), revRead);
    setSource(row(rightRev, 1), ref.i2);

    Score<int, Simple> scoringScheme(1, -2, -5);
    int splitScore = splitAlignment(left, right, scoringScheme);
    int splitScoreRev = splitAlignment(leftRev, rightRev, scoringScheme);
    
    if (splitScoreRev > splitScore)
    {
        splitScore = splitScoreRev;
        right = rightRev;
        left = leftRev;
    }

    if (splitScore < length(read)*0.5) return 1;

    // Position on the read.
    TPos readPos = toSourcePosition(row(left, 0), clippedEndPosition(row(left, 0)));
    if (readPos < 5 || readPos > length(read)-5) return 1;

    // Position on the reference infix and contig.
    refPos = Pair<TPos>(toSourcePosition(row(left, 1), clippedEndPosition(row(left, 1))),
                        toSourcePosition(row(right, 1), 0));

    if (chrOri)
    {
        // Move the split position to the rightmost possible position.
        while (length(ref.i1) > refPos.i1 && length(ref.i2) > refPos.i2 &&
              ref.i1[refPos.i1] == ref.i2[refPos.i2])
        {
            ++refPos.i1;
            ++refPos.i2;
        }
    }
    
    if (refPos.i1 > length(ref.i1)-5 || refPos.i2 < 5) return 1;

    return 0;
}

// ==========================================================================
// Function bestSplitPosition()
// ==========================================================================

template<typename TPos, typename TSize>
bool
bestSplitPosition(TPos & splitPos, TSize & maxCount, TSize & totalCount, std::map<TPos, unsigned> const & map)
{
    typedef typename std::map<TPos, unsigned>::const_iterator TIter;
    
    if (map.size() == 0) return 1;

    TIter it = map.begin();
    TIter itEnd = map.end();
    
    while (it != itEnd)
    {
        unsigned cnt = it->second;
        totalCount += cnt;
        if (cnt > maxCount)
        {
            maxCount = cnt;
            splitPos = it->first;
        }
        ++it;
    }
    
    if (maxCount < 0.5*totalCount) return 1;
    else return 0;
}

// ==========================================================================

template<typename TStream>
bool
initVcfStream(TStream & vcfStream, CharString & filename)
{
    vcfStream.open(toCString(filename), std::ios_base::out);
    if (!vcfStream.is_open())
    {
        std::cerr << "ERROR: Could not open VCF output file " << filename << std::endl;
        return 1;
    }

    // TODO write VCF header!
    
    return 0;
}

// ==========================================================================

template<typename TStream, typename TPos, typename TScore, typename TSize, typename TSize1>
void
writeVcf(TStream & vcfStream,
         CharString & chr, CharString & contig, TPos chrPos, TPos contigPos, bool chrOri, bool contigOri,
         TScore & score, TSize numReads, TSize1 splitReads, unsigned splitReadsSamePosition)
{
    VcfRecord record;
    
    StringSet<CharString> seqNames, sampleNames;
    appendValue(seqNames, chr);
    VcfIOContext context(seqNames, sampleNames);
    
    record.rID = 0;
    record.beginPos = chrPos;
    record.ref = 'N';
    
    // Create ID for alternative haplotype.
    std::stringstream altId;
    altId << "alt_" << chr << "_" << chrPos << "_" << contig << (contigOri!=chrOri ? "f" : "r");
    record.id = altId.str();
    
    // Create the ALT field.
    std::stringstream alt;
    if (chrOri)
    {
        alt << "N[" << contig << (contigOri ? "r" : "f");
        if (contigPos != maxValue<TPos>()) alt << ":" << contigPos;
        alt << "[";
    }
    else
    {
        alt << "]" << contig << (contigOri ? "f" : "r");
        if (contigPos != maxValue<TPos>()) alt << ":" << contigPos;
        alt << "]N";
    }
    record.alt = alt.str();
    
    // Create the info field.
    std::stringstream info;
    info << "AS=" << score << ";" << "RP=" << numReads << ";";
    if (splitReads != 0) info << "SR=" << splitReads << ";" << "SP=" << splitReadsSamePosition << ";";
    record.info = info.str();

    writeRecord(vcfStream, record, context, Vcf());
}

// ==========================================================================
// Function popins_place()
// ==========================================================================

int popins_place(int argc, char const ** argv)
{
    typedef String<Location> TLocations;
    typedef Iterator<TLocations>::Type TLocsIter;
    typedef Position<Dna5String>::Type TPos;

    // Parse the command line to get option values.
    PlacingOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;
        
    TLocations locations;

    if (!exists(options.locationsFile))
    {
        if (length(options.locationsFiles) == 0)
        {
            std::cerr << "ERROR: Locations file " << options.locationsFile << "does not exist. Specify -l option to create it." << std::endl;
            return 1;
        }

        // Merge approximate locations and write them to a file.
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Merging locations files." << std::endl;
        if (mergeLocations(locations, options.locationsFiles, options.verbose) != 0) return 1;
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Sorting locations." << std::endl;
        LocationPosLess less;
        std::stable_sort(begin(locations, Standard()), end(locations, Standard()), less);
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Scoring locations." << std::endl;
        scoreLocations(locations);
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Writing locations to " << options.locationsFile << std::endl;
        writeLocations(options.locationsFile, locations);
    }
    else
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Locations file exists." << std::endl;
    }

    if (length(options.bamFiles) == 0)
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "No split mapping. Specify -b option for split mapping." << std::endl;
        return 0;
    }
    
    // Read the locations file.
    if (length(locations) == 0)
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Reading " << options.batchSize << " locations from " << options.locationsFile << std::endl;
        if (readLocations(locations, options.locationsFile, options.batchSize, options.batchIndex) != 0) return 1;
    }

    // Discard locations with score below options.minLocScore or OTHER or longer than 2*maxInsertSize
    unsigned i = 0;
    while (i < length(locations))
    {
        if (locations[i].score < options.minLocScore || locations[i].chr == "OTHER" ||
            locations[i].chrEnd - locations[i].chrStart > 2*options.maxInsertSize)
            replace(locations, i, i+1, TLocations());
        else ++i;
    }
    
    if (length(locations) == 0)
    {
        std::cerr << "[" << time(0) << "] " << "No locations on genome left after filtering for score >= " << options.minLocScore << std::endl;
        return 0;
    }
    
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Keeping " << length(locations) << " locations with score >= "
                  << options.minLocScore << " and shorter than " << (2*options.maxInsertSize) << std::endl;
    
    // Concatenate the artificial reference for each location.
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Collecting contig and location sequences from "
                                            << options.referenceFile << " and " << options.supercontigFile << std::endl;
    String<Pair<Dna5String> > artificialRefs;
    if (artificialReferences(artificialRefs, locations,
                             options.referenceFile, options.supercontigFile,
                             options.readLength, options.maxInsertSize) != 0) return 1;

    // Split alignment per individual.
    BamStream bamStream;
    BamIndex<Bai> bamIndex;

    unsigned discardedLocs = 0;
    
    String<std::map<Pair<TPos>, unsigned> > splitPosMaps;
    resize(splitPosMaps, length(locations));

    Iterator<String<CharString> >::Type filesEnd = end(options.bamFiles);
    for (Iterator<String<CharString> >::Type file = begin(options.bamFiles); file != filesEnd; ++file)
    {
        if (options.verbose) std::cerr << "[" << time(0) << "] " << "Split aligning reads from " << *file << std::endl;
        if (openBamLoadBai(bamStream, bamIndex, *file) != 0) return 1;
        
        for (unsigned i = 0; i < length(locations); ++i)
        {
            // Jump to the location in bam file.
            TPos locStart = 0, locEnd = 0;
            bool hasAlignments = jumpToLocation(locStart, locEnd, bamStream, bamIndex, locations[i],
                                                options.readLength, options.maxInsertSize);
            if (!hasAlignments) continue;
            
            unsigned readCount = 0;
            bool highCoverage = false;
            
//            std::cout << "Jumped to " << locations[i].chr << ":" << locations[i].chrStart << " " << locations[i].contig << std::endl;

            // Read and split align reads.
            BamAlignmentRecord record;
            record.beginPos = minValue<TPos>();
            while (!atEnd(bamStream) && (TPos)record.beginPos < locEnd)
            {
                if (readCount * options.readLength > 100 * (locEnd - locStart + options.readLength))
                {
                    highCoverage = true;
                    break;
                }
                if (readRecord(record, bamStream) != 0)
                {
                    std::cerr << "ERROR while reading bam alignment record from " << *file << std::endl;
                    return 1;
                }
                if ((TPos)record.beginPos < locStart) continue;
                ++readCount;
                if (length(record.cigar) == 1) continue;
                
                Pair<TPos> refPos;
                if (splitAlign(refPos, artificialRefs[i], record, locations[i].chrOri) != 0) continue;

                // Make refPos.i1 the genomic position and refPos.i2 the contig position.
                if (locations[i].chrOri)
                {
                    refPos.i1 += locations[i].chrStart + options.readLength;
                }
                else
                {
                    TPos help = refPos.i1;
                    refPos.i1 = refPos.i2 + locations[i].chrStart - options.maxInsertSize;
                    refPos.i2 = help;
                }

                if (splitPosMaps[i].count(refPos) == 0) splitPosMaps[i][refPos] = 1;
                else ++splitPosMaps[i][refPos];
            }

            // Discard high coverage locations.
            if (highCoverage)
            {
                replace(locations, i, i+1, TLocations());
                replace(artificialRefs, i, i+1, String<Pair<Dna5String> >());
                replace(splitPosMaps, i, i+1, String<std::map<Pair<TPos>, unsigned> >());
                ++discardedLocs;
                --i;
            }
        }
        if (file == begin(options.bamFiles) && options.verbose)
            std::cerr << "[" << time(0) << "] " << "Discarded " << discardedLocs << " locations because of high coverage." << std::endl;
    }
        
    // Find the best split positions in sets and write the output.
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Identifying best split positions and writing output to "
                                            << options.vcfInsertionsFile << std::endl;
    std::fstream vcfStream;
    if (initVcfStream(vcfStream, options.vcfInsertionsFile) != 0) return 1;
    for (unsigned i = 0; i < length(locations); ++i)
    {
        Location loc = locations[i];
        unsigned maxCount = 0;
        unsigned totalCount = 0;
        Pair<TPos> splitPos;
        if (bestSplitPosition(splitPos, maxCount, totalCount, splitPosMaps[i]) == 0)
        {
            writeVcf(vcfStream, loc.chr, loc.contig, splitPos.i1, splitPos.i2, loc.chrOri, loc.contigOri,
                     loc.score, loc.numReads, totalCount, maxCount);
        }
        else
        {
            if (loc.chrOri)
                writeVcf(vcfStream, loc.chr, loc.contig, loc.chrEnd + options.readLength, maxValue<TPos>(), loc.chrOri, loc.contigOri,
                         loc.score, loc.numReads, 0, 0u);
            else
                writeVcf(vcfStream, loc.chr, loc.contig, loc.chrStart, maxValue<TPos>(), loc.chrOri, loc.contigOri,
                         loc.score, loc.numReads, 0, 0u);
        }
    }

    return 0;
}

#endif // #ifndef POPINS_PLACE_H_
