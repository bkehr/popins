#ifndef POPINS_PLACE_H_
#define POPINS_PLACE_H_

#include <ctime>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/vcf_io.h>
#include <seqan/seq_io.h>

#include "popins_location.h"
#include "popins_location_info.h"
#include "popins_place_ref_align.h"
#include "popins_place_split_align.h"
#include "popins_place_combine.h"

using namespace seqan;

// ==========================================================================
// Functions openVcf(), initVcf()
// ==========================================================================

template<typename TStream>
bool
openVcf(TStream & vcfStream, CharString & filename)
{
    vcfStream.open(toCString(filename), std::ios_base::app);
    if (!vcfStream.is_open())
    {
        std::cerr << "ERROR: Could not open VCF output file " << filename << std::endl;
        return 1;
    }
    
    return 0;
}

template<typename TStream>
bool
initVcf(TStream & vcfStream, PlacingOptions & options, FaiIndex & fai)
{
    vcfStream.open(toCString(options.outFile), std::ios_base::out);
    if (!vcfStream.is_open())
    {
        std::cerr << "ERROR: Could not open VCF output file " << options.outFile << std::endl;
        return 1;
    }
    
    // Get today's date, e.g. '2016-Apr-15'.
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%b-%d", &tstruct);

    // Write the header.
    vcfStream << "##fileformat=VCFv4.2" << std::endl;
    vcfStream << "##fileDate=" << buf << std::endl;
    vcfStream << "##source=PopIns_" << VERSION << std::endl;
    vcfStream << "##assembly=" << options.supercontigFile << std::endl;
    vcfStream << "##reference=" << options.referenceFile << std::endl;

    for (unsigned rID = 0; rID < numSeqs(fai); ++rID)
    {
        CharString seqName = sequenceName(fai, rID);
        if (isChromosome(seqName))
            vcfStream << "##contig=<ID=" << seqName << ",length=" << sequenceLength(fai, rID) << ">" << std::endl;
    }

    vcfStream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcfStream << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">" << std::endl;
    
    vcfStream << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT";
    vcfStream << "\t" << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << "\t" << "FORMAT" << std::endl;

    return 0;
}

// =======================================================================================
// Function mergeLocations()
// =======================================================================================

int
mergeLocations(String<Location> & locations, PlacingOptions & options)
{
    // Open output file.
    std::fstream stream(toCString(options.outFile), std::ios::out);
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open locations file " << options.locationsFile << " for writing." << std::endl;
        return 1;
    }

    // Merge approximate locations and write them to a file.
    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Merging locations files." << std::endl;

    if (mergeLocations(stream, locations, options.locationsFiles, options.locationsFile, options.maxInsertSize, options.verbose) != 0)
        return 1;

    return 0;
}

// =======================================================================================
// Function loadLocations()
// =======================================================================================

int
loadLocations(String<LocationInfo> & locations, PlacingOptions & options)
{
    if (length(locations) == 0)
    {
        if (options.bamFile != "")
            options.minAnchorReads = 0;

        LocationsFilter filter(options.minAnchorReads, options.minLocScore, 3*options.maxInsertSize);
        if (options.interval.i1 == "")
        {
            if (options.verbose)
                std::cerr << "[" << time(0) << "] " << "Reading locations from " << options.locationsFile << std::endl;
            if (readLocations(locations, options.locationsFile, filter) != 0)
                return 1;
        }
        else
        {
            if (options.verbose)
                std::cerr << "[" << time(0) << "] " << "Reading locations in " << options.interval.i1 << ":" << options.interval.i2 << "-" << options.interval.i3
                                                    << " from " << options.locationsFile << std::endl;
                if (readLocations(locations, options.locationsFile, options.interval, filter) != 0)
                    return 1;
        }
    }
    else
    {
        if (options.verbose)
            std::cerr << "[" << time(0) << "] " << "Sorting locations." << std::endl;

        LocationInfoPosLess less;
        std::stable_sort(begin(locations), end(locations), less);
    }

    if (options.verbose)
        std::cerr << "[" << time(0) << "] " << "Loaded " << length(locations) << " locations that pass filters." << std::endl;

    return 0;
}

// ==========================================================================
// Function loadContigs()
// ==========================================================================

template<typename TSeq>
bool
loadContigs(std::vector<std::pair<CharString, TSeq> > & contigs,
            String<LocationInfo> & locs,
            CharString & filename,
            bool verbose)
{
    typedef std::pair<CharString, TSeq> TPair;
    
    if (verbose)
        std::cerr << "[" << time(0) << "] " << "Reading contig sequences from " << filename << std::endl;

    // Prepare the contigs vector to contain only contigs that anchor to a listed location.
    std::set<CharString> contigSet;
    typename Iterator<String<LocationInfo> >::Type it = begin(locs);
    typename Iterator<String<LocationInfo> >::Type itEnd = end(locs);
    while (it != itEnd)
    {
        contigSet.insert((*it).loc.contig);
        ++it;
    }

    // Open fasta file.
    SequenceStream stream(toCString(filename));
    if (!isGood(stream))
    {
        std::cerr << "ERROR: Could not open " << filename << std::endl;
        return 1;
    }

    // Read records from file and append to contigs.
    while (!atEnd(stream) && contigSet.size() != 0)
    {
        CharString id;
        TSeq seq;

        // read record
        if (readRecord(id, seq, stream) != 0)
        {
            std::cerr << "ERROR: Could not read FASTA record from " << filename << std::endl;
            return 1;
        }
        
        unsigned i = 0;
        for (; i < length(id); ++i) if (id[i] == ' ') break;
        id = prefix(id, i);
        
        // Append the contig and erase it from set.
        if (contigSet.count(id) != 0)
        {
            contigs.push_back(TPair(id, seq));
            contigSet.erase(id);
        }
    }

    if (verbose)
        std::cerr << "[" << time(0) << "] " << "Loaded " << contigs.size() << " contig sequences." << std::endl;

    std::stable_sort(contigs.begin(), contigs.end());

    return 0;
}

// ==========================================================================
// Function popins_place()
// ==========================================================================

int popins_place(int argc, char const ** argv)
{
    std::fstream vcfStream;
    String<LocationInfo> locs;
    std::vector<std::pair<CharString, Dna5String> > contigs;
    FaiIndex fai;

    // Parse the command line to get option values.
    PlacingOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    // Step 1: MERGE THE LOCATIONS FROM ALL INDIVIDUALS
    if (!options.isVcf && length(options.locationsFiles) > 0)
    {
        String<Location> locs;
        mergeLocations(locs, options);
    }

    if (length(options.locationsFiles) == 0 || (options.isVcf && options.bamFile != ""))
    {
        // Load the locations.
        if (loadLocations(locs, options) != 0)
            return 1;

        // Load the contig file.
        if (loadContigs(contigs, locs, options.supercontigFile, options.verbose) != 0)
            return 1;

        // Open the FAI file of the reference genome.
        if (read(fai, toCString(options.referenceFile)) != 0)
        {
            std::cerr << "ERROR: Could not open FAI index for " << options.referenceFile << std::endl;
            return 1;
        }
    
        if (options.isVcf)
        {
            // Step 2: DO THE REFERENCE ALIGNMENT FOR ALL LOCATIONS
            if (initVcf(vcfStream, options, fai) != 0)
                return 1;
            if (popins_place_ref_align(vcfStream, locs, contigs, fai, options) != 0)
                return 1;
        }

        // Step 3: DO THE SPLIT-READ ALIGNMENT FOR AN INDIVIDUAL's LOCATIONS
        if (options.bamFile != "" && length(options.bamFiles) == 0)
        {
            if (popins_place_split_read_align(locs, contigs, fai, options) != 0)
                return 1;
            return 0;
        }
        else if (length(options.bamFiles) > 0)
        {
            for (unsigned i = 0; i < length(options.bamFiles); ++i)
            {
                clear(locs);
                // TODO load locations for the sample.
                options.bamFile = options.bamFiles[i];
                options.bamAvgCov = options.bamAvgCovs[i];
                // TODO set output file for the sample options.outFile = ;
                if (popins_place_split_read_align(locs, contigs, fai, options) != 0)
                    return 1;
            }
        }
    }

    // Step 4: COMBINE SPLIT-READ ALIGNMENTS FROM ALL INDIVIDUALS
    if (options.isVcf && length(options.locationsFiles) > 0)
    {
        if (openVcf(vcfStream, options.outFile) != 0)
            return 1;
        if (popins_place_combine(vcfStream, options) != 0)
            return 1;
    }

    return 0;
}

#endif  // POPINS_PLACE_H
