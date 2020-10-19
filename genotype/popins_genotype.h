#ifndef POPINS_GENOTYPE_H_
#define POPINS_GENOTYPE_H_

#include "../popins_utils.h"
#include "../command_line_parsing.h"
#include "variant_caller.h"
#include <algorithm>    // std::sort

using namespace seqan;

// ==========================================================================

#define LL_THRESHOLD -25.5

int phredLikelihoodsToGenotypeQuality(const int phr0, const int phr1, const int phr2){
    int pl[] = {phr0, phr1, phr2};
    int n = sizeof(pl)/sizeof(pl[0]);
    std::sort(pl, pl+n);

    const int pq = pl[1] - pl[0];
    return pq;
}

void
probsToGtString(std::vector<double> & probs, std::string & gtString)
{
    int max = 0;
    double maxP = probs[0];
    if( probs[1] > maxP ){
        maxP = probs[1];
        max = 1;
    }
    if( probs[2] > maxP ){
        maxP = probs[2];
        max = 2;
    }
    std::ostringstream buff;
    if( max == 0 ){
        buff << "0/0:";
    }else if( max == 1 ){
        buff << "0/1:";
    }else{
        buff << "1/1:";
    }
    double lp0 = log10( probs[0]/maxP );
    double lp1 = log10( probs[1]/maxP );
    double lp2 = log10( probs[2]/maxP );
    if( lp0 < LL_THRESHOLD ) lp0 = LL_THRESHOLD;
    if( lp1 < LL_THRESHOLD ) lp1 = LL_THRESHOLD;
    if( lp2 < LL_THRESHOLD ) lp2 = LL_THRESHOLD;

    int phr0 = int( -10*lp0 );
    int phr1 = int( -10*lp1 );
    int phr2 = int( -10*lp2 );
    buff  << phr0 << "," << phr1 << "," << phr2 << ":" << phredLikelihoodsToGenotypeQuality(phr0,phr1,phr2);
    gtString = buff.str();
}

// ==========================================================================
// Function popins_genotype()
// ==========================================================================
using namespace std;
int
popins_genotype(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    GenotypingOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    printStatus("Opening input files.");

    CharString samplePath = getFileName(options.prefix, options.sampleID);

    // Load the POPINS_SAMPLE_INFO file.
    SampleInfo sampleInfo;
    CharString sampleInfoFile = getFileName(samplePath, "POPINS_SAMPLE_INFO");
    if (readSampleInfo(sampleInfo, sampleInfoFile) != 0)
       return 1;

    // Open the input VCF file and prepare output VCF stream.
    VcfFileIn vcfIn(toCString(options.vcfFile));
    CharString outfile = getFileName(samplePath, "insertions.vcf");
    std::ofstream vcfStream(toCString(outfile));
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, vcfStream, Vcf());

    appendName(sampleNamesCache(context(vcfOut)), options.sampleID);

    VcfHeader header;
    readHeader(header, vcfIn);
    writeHeader(vcfOut, header);

    /*
    string chr, f1,f2,f3,f4,f5,f6,f7,f8;
    ifstream f;
    // The chromosomes in the output file need to be in the same order as in the input file
    f.open( toCString( options.vcfFile ) );
    map< string, int> chrs;
    f >> chr >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8;
    while( f ){
        f >> chr >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8;
        if( chrs.count( chr ) == 0 ){
            appendValue(vcfOut.header.sequenceNames, chr.c_str());
            chrs[chr] = 1;
        }
    }
    f.close();
    */

    // Build an index of the fasta file (reference genome).
    FaiIndex faIndex;
    if (!open(faIndex, toCString(options.referenceFile)))
    {
        if (!build(faIndex, toCString(options.referenceFile)))
        {
            std::cerr << "ERROR: Could not find nor build the index of " << options.referenceFile << std::endl;
            return 7;
        }
    }

    // Open the bam file. (A bam file needs the bai index and the bam file stream.)
    BamIndex<Bai> bamIndex;
    BamFileIn bamStream;
    if (initializeBam(toCString(sampleInfo.bam_file), bamIndex, bamStream) != 0)
       return 7;

    // Build an index of the insertion sequences' fasta file.
    FaiIndex faIndexAlt;
    if (!open(faIndexAlt, toCString(options.supercontigFile)))
    {
        if (!build(faIndexAlt, toCString(options.supercontigFile)))
        {
            std::cerr << "ERROR: Could not build the index of " << options.supercontigFile << std::endl;
            return 7;
        }
    }

    // Open the bam file. (A bam file needs the bai index and the bam file stream.)
    BamIndex<Bai> bamIndexAlt;
    BamFileIn bamStreamAlt;
    CharString altBamFile = getFileName(samplePath, "non_ref_new.bam");
    if (initializeBam(toCString(altBamFile), bamIndexAlt, bamStreamAlt))
        return 7;

    std::ostringstream msg;
    msg << "Genotyping VCF records in \'" << options.vcfFile << "\' for " << options.sampleID << ".";
    printStatus(msg);

    // Iterate over VCF file and call the variants.
    VcfRecord record;    
    while (!atEnd(vcfIn))
    {
        readRecord(record, vcfIn);

        std::vector<double> vC(3);
        variantCallRegion(record, vcfIn, faIndex, faIndexAlt, bamIndex, bamStream, bamIndexAlt, bamStreamAlt, options, vC);

        std::string gtString;
        probsToGtString(vC, gtString);
        appendValue(record.genotypeInfos, gtString);
        record.format = "GT:PL:GQ";
        writeRecord(vcfOut, record);
    }

    return 0;
}

#endif // #ifndef POPINS_GENOTYPE_H_
