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

template<typename TStream, typename TTag>
bool
initVcf(TStream & vcfStream, PlacingOptions<TTag> & options, FaiIndex & fai)
{
    vcfStream.open(toCString(options.outFile), std::ios_base::out);
    if (!vcfStream.is_open())
    {
        std::cerr << "ERROR: Could not open VCF output file " << options.outFile << std::endl;
        return 1;
    }

    std::set<CharString> chromosomes;
    if (readChromosomes(chromosomes, options.referenceFile) != 0)
        return 1;

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
        if (isChromosome(seqName, chromosomes))
            vcfStream << "##contig=<ID=" << seqName << ",length=" << sequenceLength(fai, rID) << ">" << std::endl;
    }
    vcfStream << "##INFO=<ID=AR,Number=1,Type=Integer,Description=\"Number of anchoring read pairs.\">" << std::endl;
    vcfStream << "##INFO=<ID=AS,Number=1,Type=Float,Description=\"Anchoring score.\">" << std::endl;
    vcfStream << "##INFO=<ID=GS,Number=1,Type=Integer,Description=\"Number of contigs in the contig group.\">" << std::endl;

    vcfStream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcfStream << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">" << std::endl;
    vcfStream << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;

    vcfStream << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT";
    vcfStream << "\t" << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << std::endl;

    return 0;
}

// =======================================================================================
// Function mergeLocations()
// =======================================================================================

template<typename TTag>
int
mergeLocations(String<Location> & locations, PlacingOptions<TTag> & options)
{
    // Open output file.
    std::fstream stream(toCString(options.locationsFile), std::ios::out);
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

    close(stream);

    return 0;
}

// =======================================================================================
// Function loadLocations()
// =======================================================================================

int
loadLocations(String<LocationInfo> & locations, CharString & sampleID, CharString & filename, double minLocScore, unsigned minAnchorReads, unsigned maxInsertSize)
{
    std::ostringstream msg;
    msg << "Reading locations from " << filename;
    printStatus(msg);

    LocationsFilter filter(minAnchorReads, minLocScore, 3*maxInsertSize);
    if (readLocations(locations, sampleID, filename, filter) != 0)
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
loadInputAndSplitReadAlign(CharString & samplePath, PlacingOptions<SplitAlign> & options, FaiIndex & fai)
{
    // Load the POPINS_SAMPLE_INFO file.
    SampleInfo sampleInfo;
    CharString sampleInfoFile = getFileName(samplePath, "POPINS_SAMPLE_INFO");
    if (readSampleInfo(sampleInfo, sampleInfoFile) != 0)
        return 1;

    // Load the locations.
    String<LocationInfo> locs;
    CharString locationsFile = getFileName(samplePath, "locations_unplaced.txt");
    if (!exists(locationsFile))
    {
        std::ostringstream msg;
        msg << "WARNING: No file \'locations_unplaced.txt\' present for sample \'" << options.sampleID<< "\'.";
        printStatus(msg);
        return 0;
    }
    if (loadLocations(locs, sampleInfo.sample_id, locationsFile, 0.0, 0u, options.maxInsertSize) != 0)
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
// Function popins_place_refalign()
// ==========================================================================

int popins_place_refalign(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    PlacingOptions<RefAlign> options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Placing step 1: MERGE THE LOCATIONS FROM ALL INDIVIDUALS
    if (!exists(options.locationsFile))
    {
    	String<Location> locs;
    	mergeLocations(locs, options);
    }
    else
    {
    	std::ostringstream msg;
    	msg << "WARNING: Locations file \'" << options.locationsFile << "\' exists. Skipping merging step.";
    	printStatus(msg);
    }

    // Placing step 2: REFERENCE ALIGNMENT FOR ALL LOCATIONS

    // Open the FAI file of the reference genome.
    FaiIndex fai;
    if (!open(fai, toCString(options.referenceFile)))
    {
        std::cerr << "ERROR: Could not open FAI index for " << options.referenceFile << std::endl;
        return 7;
    }

    // Load the locations.
    String<LocationInfo> locs;
    CharString sID = "";
    if (loadLocations(locs, sID, options.locationsFile, options.minLocScore, options.minAnchorReads, options.maxInsertSize) != 0)
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

    return 0;
}

// ==========================================================================
// Function popins_place_splitalign()
// ==========================================================================

int popins_place_splitalign(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    PlacingOptions<SplitAlign> options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Placing step 3: SPLIT-READ ALIGNMENT PER INDIVIDUAL

    // Open the FAI file of the reference genome.
    FaiIndex fai;
    if (!open(fai, toCString(options.referenceFile)))
    {
        std::cerr << "ERROR: Could not open FAI index for " << options.referenceFile << std::endl;
        return 7;
    }

    // Do the split read alignment for the specified sample.
    CharString samplePath = getFileName(options.prefix, options.sampleID);
    if (loadInputAndSplitReadAlign(samplePath, options, fai) != 0)
       return 7;

    return 0;
}

// ==========================================================================
// Function popins_place_finish()
// ==========================================================================

int popins_place_finish(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    PlacingOptions<SplitCombine> options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

	// Placing step 4: COMBINE SPLIT-READ ALIGNMENTS FROM ALL INDIVIDUALS

    // Open the output file.
    std::ofstream vcfStream;
     if (openVcf(vcfStream, options.outFile) != 0)
         return 7;

     // Combine the split-read alignments of all individuals and write VCF records.
     if (popins_place_combine(vcfStream, options.prefix, options.referenceFile, options.outFile) != 0)
         return 7;

     return 0;
}

#endif  // POPINS_PLACE_H
