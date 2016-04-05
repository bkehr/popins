#include <sstream>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/vcf_io.h>
#include <seqan/seq_io.h>

#include "popins_location.h"
#include "popins_place_split_align.h"
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

template<typename TSeq, typename TPos, typename TSize>
bool
artificialReferences(std::map<TSize, LocationInfo<TSeq, TPos, TSize> > & locInfos,
                     std::map<TSize, std::set<TSize> > & concatGroups,
                     std::map<TSize, std::set<TSize> > & groups,
                     String<Location> & locations,
                     std::map<CharString, TSeq> & contigs,
                     PlacingOptions & options)
{
    typedef typename std::map<TSize, std::set<TSize> >::iterator TIter;
    
    // Use fai to jump into reference.
    FaiIndex fai;
    if (read(fai, toCString(options.referenceFile)) != 0)
    {
        std::cerr << "ERROR: Could not open fai index for " << options.referenceFile << std::endl;
        return 1;
    }

    Score<int, Simple> scsc(1, -2, -5);

    // Create reference sequences.
    TIter groupsEnd = groups.end();
    for (TIter it = groups.begin(); it != groupsEnd; )
    {
        Location & loc = locations[it->first];
        LocationInfo<TSeq, TPos, TSize> & locInfo = locInfos[it->first];
        if (contigs.count(loc.contig) == 0)
        {
            std::cerr << "ERROR: Could not find record for " << loc.contig << " in contigs file." << std::endl;
            return 1;
        }
        unsigned idx = 0;
        if (!getIdByName(fai, loc.chr, idx))
        {
            std::cerr << "ERROR: Could not find " << loc.chr << " in FAI index." << std::endl;
            return 1;
        }

        TSeq contig = contigs[loc.contig];
        if (loc.contigOri == loc.chrOri)
            reverseComplement(contig);

        if (loc.chrOri)
        {
            TSeq chrInfix;
            readRegion(chrInfix, fai, idx, loc.chrStart + options.readLength, loc.chrEnd+options.maxInsertSize);
            locInfo.refOffset = loc.chrStart + options.readLength;
            
            Pair<size_t> splitPosPair;
            if (alignContigPrefixToRef(splitPosPair, chrInfix, contig, scsc))
            {
                locInfo.refOffset += splitPosPair.i1;
                locInfo.contigOffset += splitPosPair.i2;

                TSeq concatSeq = infix(chrInfix, std::max(0, (int)splitPosPair.i1 - (int)options.readLength), splitPosPair.i1);
                concatSeq += suffix(contig, splitPosPair.i2);
                locInfo.refSeq = concatSeq;
                locInfo.concatPos = splitPosPair.i1 - std::max(0, (int)splitPosPair.i1 - (int)options.readLength);
                
                concatGroups.insert(*it);
                groups.erase(it++);
            }
            else
            {
                locInfo.refSeq = chrInfix;
                locInfo.contigSeq = contig;
                ++it;
            }
        }
        else
        {
            TSeq chrInfix;
            readRegion(chrInfix, fai, idx, loc.chrStart - options.maxInsertSize, loc.chrEnd);
            locInfo.refOffset = loc.chrStart - options.maxInsertSize;
            
            Pair<size_t> splitPosPair;
            if (alignContigSuffixToRef(splitPosPair, chrInfix, contig, scsc))
            {
                locInfo.refOffset += splitPosPair.i1;
                locInfo.contigOffset += splitPosPair.i2;

                TSeq concatSeq = prefix(contig, splitPosPair.i2);
                concatSeq += infix(chrInfix, splitPosPair.i1, std::min(splitPosPair.i1 + options.readLength, length(chrInfix)));
                locInfo.refSeq = concatSeq;
                locInfo.concatPos = splitPosPair.i2;
                
                concatGroups.insert(*it);
                groups.erase(it++);
            }
            else
            {
                locInfo.refSeq = chrInfix;
                locInfo.contigSeq = contig;
                ++it;
            }
        }
    }

    return 0;
}

// ==========================================================================
// Function bestSplitPosition()
// ==========================================================================

template<typename TPos, typename TSize1, typename TSize2>
bool
bestSplitPosition(TPos & splitPos, TSize1 & maxCount, TSize1 & totalCount, std::map<TPos, TSize2> const & map)
{
    typedef typename std::map<TPos, TSize2>::const_iterator TIter;
    
    if (map.size() == 0) return 1;

    totalCount = 0;
    maxCount = 0;

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
         TScore & score, TSize numReads, TSize1 splitReads, unsigned splitReadsSamePosition, bool groupRepresentative)
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
        if (contigPos != maxValue<TPos>() && groupRepresentative) alt << ":" << contigPos;
        alt << "[";
    }
    else
    {
        alt << "]" << contig << (contigOri ? "f" : "r");
        if (contigPos != maxValue<TPos>() && groupRepresentative) alt << ":" << contigPos;
        alt << "]N";
    }
    record.alt = alt.str();
    
    // Create the info field.
    std::stringstream info;
    info << "AS=" << score << ";" << "RP=" << numReads << ";";
    if (splitReads != 0)
    {
        if (splitReadsSamePosition != 0)
            info << "SR=" << splitReads << ";" << "SP=" << splitReadsSamePosition << ";";
        else
            info << "CR=" << splitReads << ";";
    }
    if (!groupRepresentative)
        info << "GROUPED;";
    record.info = info.str();

    writeRecord(vcfStream, record, context, Vcf());
}

// ==========================================================================
// Function findBestSplit()
// ==========================================================================

template<typename TSize, typename TSeq, typename TPos>
bool
findBestSplitAndWriteVcf(std::fstream & vcfStream,
                         String<Location> & locations,
                         std::map<TSize, std::set<TSize> > &  groups,
                         std::map<TSize, LocationInfo<TSeq, TPos, TSize> > & locInfos,
                         //std::set<TSize> & highCoverageLocs,
                         //std::map<TSize, std::map<Pair<TPos>, unsigned> > & splitPosMaps,
                         //std::map<TSize, Pair<TPos> > & refOffsets,
                         PlacingOptions & options)
{
    for (typename std::map<TSize, std::set<TSize> >::iterator it = groups.begin(); it != groups.end(); ++it)
    {
        LocationInfo<TSeq, TPos, TSize> & locInfo = locInfos[it->first];
        if (locInfo.highCoverage)
            continue;

        unsigned maxCount = 0;
        unsigned totalCount = locInfo.splitReadCount;
        Pair<TPos> splitPos = Pair<TPos>(locInfo.refOffset, locInfo.contigOffset);

        Location loc = locations[it->first];
        if (locInfo.concatPos != 0 || bestSplitPosition(splitPos, maxCount, totalCount, locInfo.splitPosMap) == 0)
        {
            // Write record for group representative
            writeVcf(vcfStream, loc.chr, loc.contig, splitPos.i1, splitPos.i2, loc.chrOri, loc.contigOri,
                     loc.score, loc.numReads, totalCount, maxCount, true);

            // Write records for group members
            for (typename std::set<TSize>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2)
            {
                loc = locations[*it2];
                writeVcf(vcfStream, loc.chr, loc.contig, splitPos.i1, splitPos.i2, loc.chrOri, loc.contigOri,
                         loc.score, loc.numReads, totalCount, maxCount, false);
            }
        }
        else
        {
            if (loc.chrOri)
                writeVcf(vcfStream, loc.chr, loc.contig, loc.chrEnd + options.readLength, maxValue<TPos>(), loc.chrOri, loc.contigOri,
                         loc.score, loc.numReads, 0, 0u, true);
            else
                writeVcf(vcfStream, loc.chr, loc.contig, loc.chrStart, maxValue<TPos>(), loc.chrOri, loc.contigOri,
                         loc.score, loc.numReads, 0, 0u, true);
            for (typename std::set<TSize>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2)
            {
                loc = locations[*it2];
                if (loc.chrOri)
                    writeVcf(vcfStream, loc.chr, loc.contig, loc.chrEnd + options.readLength, maxValue<TPos>(), loc.chrOri, loc.contigOri,
                             loc.score, loc.numReads, 0, 0u, false);
                else
                    writeVcf(vcfStream, loc.chr, loc.contig, loc.chrStart, maxValue<TPos>(), loc.chrOri, loc.contigOri,
                             loc.score, loc.numReads, 0, 0u, false);
            }
        }
    }
    return 0;
}

// ==========================================================================
// Function popins_place()
// ==========================================================================

int popins_place(int argc, char const ** argv)
{
    typedef Dna5String TSeq;
    typedef Size<TSeq>::Type TSize;
    typedef Position<TSeq>::Type TPos;

    // Parse the command line to get option values.
    PlacingOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;
        
    // Load and/or merge the locations.
    String<Location> locations;
    int res = loadLocations(locations, options);
    if (res == -1)
        return 1;
    else if (res == 1)
        return 0;

    // Load the contig file.
    std::map<CharString, TSeq> contigs;
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Reading contig sequences from " << options.supercontigFile << std::endl;
    if (loadSequences(contigs, options.supercontigFile) != 0)
        return 1;
    
    // Determine groups of locations that overlap and where the contig prefix is highly similar
    std::map<TSize, std::set<TSize> > groups, concatGroups;
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Grouping locations by reference position and contig sequence... " << std::flush;
    groupLocations(groups, locations, contigs);
    if (options.verbose)
        std::cerr << groups.size() << " groups." << std::endl;

    // Concatenate the artificial reference for each location.
    std::map<TSize, LocationInfo<TSeq, TPos, TSize> > locInfos;
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Collecting reference sequences for locations from " << options.referenceFile << std::endl;
    if (artificialReferences(locInfos, concatGroups, groups, locations, contigs, options) != 0)
        return 1;
    clear(contigs);

    // Split alignment per individual.
    if (findSplitReads(locInfos, concatGroups, groups, locations, options) != 0)
        return 1;

    // Open the output file.
    std::fstream vcfStream;
    if (initVcfStream(vcfStream, options.vcfInsertionsFile) != 0)
        return 1;

    // Find the best split positions in maps and write the output.
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Identifying best split positions and writing output to "
                                            << options.vcfInsertionsFile << std::endl;
    if (findBestSplitAndWriteVcf(vcfStream, locations, groups, locInfos, options) != 0)
        return 1;
    if (findBestSplitAndWriteVcf(vcfStream, locations, concatGroups, locInfos, options) != 0)
        return 1;

    return 0;
}

#endif // #ifndef POPINS_PLACE_H_
