#ifndef POPINS_PLACE_H_
#define POPINS_PLACE_H_

#include <ctime>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/vcf_io.h>
#include <seqan/seq_io.h>

#include "../popins_utils.h"
#include "../command_line_parsing.h"
#include "location.h"
#include "location_info.h"
#include "ref_align.h"
#include "split_align.h"
#include "combine.h"

using namespace seqan;

// ==========================================================================
// Function selectPlacingSteps
// ==========================================================================

void
selectPlacingSteps(PlacingOptions & options)
{
    // User specified --all option.
    if (options.all)
    {
       options.merge = true;
       options.refAlign = true;
       options.splitAlign = true;
       options.combine = true;
       printStatus("Executing all placing steps.");
       return;
    }

    // Step was not explicitly specified.
    if (!(options.all || options.merge || options.refAlign || options.splitAlign || options.combine))
    {
       if (!exists(options.locationsFile))
       {
          options.merge = true;
          options.refAlign = true;
           printStatus("Executing the merging and reference alignment steps.");
          return;
       }
       else if (options.sampleID != "")
       {
          options.splitAlign = true;
          std::ostringstream msg;
          msg << "Executing the splitAlignment step for \'" << options.sampleID << "\'.";
          printStatus(msg);
          return;
       }
       else
       {
          options.combine = true;
          printStatus("Executing the combine step.");
          return;
       }
    }
}

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
    char date[80];
    tstruct = *localtime(&now);
    strftime(date, sizeof(date), "%Y-%b-%d", &tstruct);

    // Write the header.
    vcfStream << "##fileformat=VCFv4.2" << std::endl;
    vcfStream << "##fileDate=" << date << std::endl;
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

    printStatus("Listing locations files.");

    CharString filename = "locations.txt";
    String<Pair<CharString> > locationsFiles = listFiles(options.prefix, filename);

    std::ostringstream msg;
    msg << "Merging " << length(locationsFiles) << " locations files.";
    printStatus(msg);

    // Merge approximate locations and write them to a file.
    if (mergeLocations(stream, locations, locationsFiles, options.locationsFile, options.maxInsertSize) != 0)
        return 1;

    return 0;
}

// =======================================================================================
// Function loadLocations()
// =======================================================================================

int
loadLocations(String<LocationInfo> & locations, CharString filename, double minLocScore, unsigned minAnchorReads, unsigned maxInsertSize)
{
    std::ostringstream msg;
    msg << "Reading locations from " << filename;
    printStatus(msg);

    LocationsFilter filter(minAnchorReads, minLocScore, 3*maxInsertSize);
    if (readLocations(locations, filename, filter) != 0)
        return 1;

    msg.str("");
    msg << "Loaded " << length(locations) << " locations that pass filters.";
    printStatus(msg);

    return 0;
}

// ==========================================================================
// Function loadContigs()
// ==========================================================================

template<typename TSeq>
bool
loadContigs(std::vector<std::pair<CharString, TSeq> > & contigs,
        String<LocationInfo> & locs,
        CharString & filename)
{
    typedef std::pair<CharString, TSeq> TPair;

    std::ostringstream msg;
    msg << "Reading contig sequences from " << filename;
    printStatus(msg);

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
    SeqFileIn stream(toCString(filename));

    // Read records from file and append to contigs.
    while (!atEnd(stream) && contigSet.size() != 0)
    {
        CharString id;
        TSeq seq;
        readRecord(id, seq, stream);

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

    msg.str("");
    msg << "Loaded " << contigs.size() << " contig sequences.";
    printStatus(msg);

    std::stable_sort(contigs.begin(), contigs.end());

    return 0;
}

bool
loadInputAndSplitReadAlign(CharString & samplePath, PlacingOptions & options, FaiIndex & fai)
{
   // Load the POPINS_SAMPLE_INFO file.
   SampleInfo sampleInfo;
   CharString sampleInfoFile = getFileName(samplePath, "POPINS_SAMPLE_INFO");
   if (readSampleInfo(sampleInfo, sampleInfoFile) != 0)
      return 1;

    // Load the locations.
    String<LocationInfo> locs;
    CharString locationsFile = getFileName(samplePath, "locations_unplaced.txt");
    if (loadLocations(locs, locationsFile, options.minLocScore, 0u, options.maxInsertSize) != 0)
        return 1;

   // Load the contig file.
   std::vector<std::pair<CharString, Dna5String> > contigs;
   if (loadContigs(contigs, locs, options.supercontigFile) != 0)
      return 1;

   CharString outfile = getFileName(samplePath, "locations_placed.txt");
    if (popins_place_split_read_align(outfile, locs, contigs, fai, sampleInfo, options.maxInsertSize, options.readLength) != 0)
        return 1;

    return 0;
}

// ==========================================================================
// Function popins_place()
// ==========================================================================

int popins_place(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    PlacingOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    selectPlacingSteps(options);

    // Step 1: MERGE THE LOCATIONS FROM ALL INDIVIDUALS
    if (options.merge)
    {
        String<Location> locs;
        mergeLocations(locs, options);
    }

    // Step 2: REFERENCE ALIGNMENT FOR ALL LOCATIONS
    if (options.refAlign)
    {
        // Open the FAI file of the reference genome.
        FaiIndex fai;
        if (!open(fai, toCString(options.referenceFile)))
        {
            std::cerr << "ERROR: Could not open FAI index for " << options.referenceFile << std::endl;
            return 7;
        }

        // Load the locations.
        String<LocationInfo> locs;
        if (loadLocations(locs, options.locationsFile, options.minLocScore, options.minAnchorReads, options.maxInsertSize) != 0)
            return 7;

        // Load the contig file.
        std::vector<std::pair<CharString, Dna5String> > contigs;
        if (loadContigs(contigs, locs, options.supercontigFile) != 0)
            return 7;

        // Open and initialize the output file.
        std::ofstream vcfStream;
        if (initVcf(vcfStream, options, fai) != 0)
            return 7;

        // Compute the reference alignments.
        if (popins_place_ref_align(vcfStream, locs, contigs, fai, options) != 0)
            return 7;
    }

    // Step 3: SPLIT-READ ALIGNMENT PER INDIVIDUAL
    if (options.splitAlign)
    {
        // Open the FAI file of the reference genome.
        FaiIndex fai;
        if (!open(fai, toCString(options.referenceFile)))
        {
            std::cerr << "ERROR: Could not open FAI index for " << options.referenceFile << std::endl;
            return 7;
        }

        // Do the split read alignment only for the specified sample and return.
       if (options.sampleID != "")
       {
          CharString samplePath = getFileName(options.prefix, options.sampleID);
          if (loadInputAndSplitReadAlign(samplePath, options, fai) != 0)
             return 7;

          return 0;
       }

       // Do the split read alignment only for all samples.
       String<CharString> sampleIDs = listSubdirectories(options.prefix);
       for (unsigned i = 0; i < length(sampleIDs); ++i)
       {
          CharString samplePath = getFileName(options.prefix, sampleIDs[i]);
          if (loadInputAndSplitReadAlign(samplePath, options, fai) != 0)
             return 7;
       }
    }

    // Step 4: COMBINE SPLIT-READ ALIGNMENTS FROM ALL INDIVIDUALS
    if (options.combine)
    {
       // Open the output file.
       std::ofstream vcfStream;
        if (openVcf(vcfStream, options.outFile) != 0)
            return 7;

        // Combine the split-read alignments of all individuals and write VCF records.
        if (popins_place_combine(vcfStream, options.prefix, options.referenceFile, options.outFile) != 0)
            return 7;
    }

    return 0;
}

#endif  // POPINS_PLACE_H
