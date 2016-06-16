#ifndef POPINS_GENOTYPE_H_
#define POPINS_GENOTYPE_H_

#include "popins_clp.h"
#include "variant_caller.h"

using namespace seqan;

// ==========================================================================

#define LL_THRESHOLD -25.5

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
    buff  << phr0 << "," << phr1 << "," << phr2;
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
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    printStatus("Opening input files.");

    // Open the input VCF file and prepare output vcf stream.
    VcfFileIn vcfIn(toCString(options.vcfFile));
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, Vcf());

    appendName(sampleNamesCache(context(vcfOut)), options.sampleName);

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
            std::cerr << "ERROR: Could not build the index of " << options.referenceFile << std::endl;
            return 1;
        }
    }

    // Open the bam file. (A bam file needs the bai index and the bam file stream.)
    BamIndex<Bai> bamIndex;
    BamFileIn bamStream;
    if (initializeBam(toCString(options.bamFile), bamIndex, bamStream) != 0) return 1;

    // Build an index of the insertion sequences' fasta file.
    FaiIndex faIndexAlt;
    if (!open(faIndexAlt, toCString(options.altFastaFile)))
    {
        if (!build(faIndexAlt, toCString(options.altFastaFile)))
        {
            std::cerr << "ERROR: Could not build the index of " << options.altFastaFile << std::endl;
            return 1;
        }
    }

    // Open the bam file. (A bam file needs the bai index and the bam file stream.)
    BamIndex<Bai> bamIndexAlt;
    BamFileIn bamStreamAlt;
    if (initializeBam(toCString(options.altBamFile), bamIndexAlt, bamStreamAlt))
        return 1;

    std::ostringstream msg;
    msg << "Genotyping VCF records in \'" << options.vcfFile << "\' for " << options.sampleName << ".";
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
        record.format = "GT:PL";
        writeRecord(vcfOut, record);
    }

    return 0;
}

#endif // #ifndef POPINS_GENOTYPE_H_
