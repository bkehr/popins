#include <sstream>
#include <cerrno>

#include <seqan/file.h>
#include <seqan/sequence.h>

#include "../popins_utils.h"
#include "../command_line_parsing.h"
#include "crop_unmapped.h"

#ifndef POPINS_ASSEMBLE_H_
#define POPINS_ASSEMBLE_H_

using namespace seqan;

// ==========================================================================

inline void
removeAssemblyDirectory(CharString & path)
{
    removeFile(path, "contigs.fa");
    removeFile(path, "Graph2");
    removeFile(path, "LastGraph");
    removeFile(path, "Log");
    removeFile(path, "PreGraph");
    removeFile(path, "Roadmaps");
    removeFile(path, "Sequences");
    removeFile(path, "stats.txt");
    remove(toCString(path));
}

bool
retrieveSampleID(CharString & sampleID, CharString & mappingBam)
{
    BamFileIn inStream(toCString(mappingBam));

    BamHeader header;
    readHeader(header, inStream);

    for (unsigned i = 0; i < length(header); ++i)
    {
        if (header[i].type != BamHeaderRecordType::BAM_HEADER_READ_GROUP)
            continue;

        for (unsigned j = 0; j < length(header[i].tags); ++j)
        {
            if (header[i].tags[j].i1 != "SM")
                continue;

            sampleID = header[i].tags[j].i2;
            return 0;
        }
    }
    std::cerr << "ERROR: Could not find sample ID in BAM file header." << std::endl;
    return 1;
}

// ==========================================================================
// Function remapping()
// ==========================================================================

inline int
remapping(Triple<CharString> & fastqFilesTemp,
        Triple<CharString> & fastqFiles,
        CharString const & referenceFile,
        CharString const & workingDir,
        unsigned humanSeqs,
        unsigned threads,
        CharString & memory,
        CharString & prefix)
{
    std::stringstream cmd;

    CharString f1 = prefix;
    f1 += "remapped.sam";
    CharString remappedSam = getFileName(workingDir, f1);

    CharString f2 = prefix;
    f2 += "remapped.bam";
    CharString remappedBam = getFileName(workingDir, f2);

    CharString f3 = prefix;
    f3 += "remapped.bam.bai";
    CharString remappedBai = getFileName(workingDir, f3);

    CharString f4 = prefix;
    f4 += "remapped_unsorted.bam";
    CharString remappedUnsortedBam = getFileName(workingDir, f4);

    std::ostringstream msg;
    msg << "Remapping unmapped reads using " << BWA;
    printStatus(msg);

    // Run BWA on unmapped reads (pairs).
    cmd.str("");
    cmd << BWA << " mem -t " << threads << " " << referenceFile << " " << fastqFilesTemp.i1 << " " << fastqFilesTemp.i2 << " > " << remappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFilesTemp.i1 << " and " << fastqFilesTemp.i2 << std::endl;
        return 1;
    }

    remove(toCString(fastqFilesTemp.i1));
    remove(toCString(fastqFilesTemp.i2));

    // Run BWA on unmapped reads (single end).
    cmd.str("");
    cmd << BWA << " mem -t " << threads << " " << referenceFile << " " << fastqFilesTemp.i3 << " | awk '$1 !~ /@/' >> " << remappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFilesTemp.i3 << std::endl;
        return 1;
    }

    remove(toCString(fastqFilesTemp.i3));

    msg.str("");
    msg << "Converting BWA output " << remappedSam << " to bam format.";
    printStatus(msg);

    // Convert BWA output to bam.
    cmd.str("");
    cmd << SAMTOOLS << " view -S -h -b " << remappedSam << " > " << remappedUnsortedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while converting BWA output " << remappedSam << " to bam." << std::endl;
        return 1;
    }
    remove(toCString(remappedSam));

    msg.str("");
    msg << "Sorting " << remappedUnsortedBam << " using " << SAMTOOLS;
    printStatus(msg);

    // Sort bam file.
    cmd.str("");
    cmd << SAMTOOLS << " sort -@ " << threads << " -m " << memory << " -o " << remappedBam << " " << remappedUnsortedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting BWA output " << remappedUnsortedBam << std::endl;
        return 1;
    }

    msg.str("");
    msg << "Indexing " << remappedBam << " using " << SAMTOOLS;
    printStatus(msg);

    // Index bam file.
    cmd.str("");
    cmd << SAMTOOLS << " index " << remappedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while indexing BWA output " << remappedBam << std::endl;
        return 1;
    }

    msg.str("");
    msg << "Cropping unmapped reads from " << remappedBam;
    printStatus(msg);

    // Crop unmapped and create bam file of remapping.
    if (crop_unmapped(fastqFiles, remappedUnsortedBam, remappedBam, humanSeqs, NoAdapters()) != 0)
        return 1;
    remove(toCString(remappedBai));

    msg.str("");
    msg << "Sorting " << remappedUnsortedBam << " by read name using " << SAMTOOLS;
    printStatus(msg);

    // Sort <WD>/remapped.bam by read name.
    cmd.str("");
    cmd << SAMTOOLS << " sort -n -@ " << threads << " -m " << memory << " -o " << remappedBam << " " << remappedUnsortedBam << " ";
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << remappedUnsortedBam << std::endl;
        return 1;
    }

    remove(toCString(remappedUnsortedBam));

    return 0;
}

// ==========================================================================

inline void
setMates(BamAlignmentRecord & record1, BamAlignmentRecord & record2)
{
    SEQAN_ASSERT(!hasFlagFirst(record1) || !hasFlagFirst(record2));
    SEQAN_ASSERT(!hasFlagLast(record1) || !hasFlagLast(record2));

    // Set the next ref id.
    record1.rNextId = record2.rID;
    record2.rNextId = record1.rID;

    // Set the next ref pos.
    record1.pNext = record2.beginPos;
    record2.pNext = record1.beginPos;

    // Fix the next unmapped flag.
    if (hasFlagUnmapped(record2)) record1.flag |= BAM_FLAG_NEXT_UNMAPPED;
    else record1.flag &= ~BAM_FLAG_NEXT_UNMAPPED;
    if (hasFlagUnmapped(record1)) record2.flag |= BAM_FLAG_NEXT_UNMAPPED;
    else record2.flag &= ~BAM_FLAG_NEXT_UNMAPPED;

    // Fix the next reversed flag.
    if (hasFlagRC(record2)) record1.flag |= BAM_FLAG_NEXT_RC;
    else record1.flag &= ~BAM_FLAG_NEXT_RC;
    if (hasFlagRC(record1)) record2.flag |= BAM_FLAG_NEXT_RC;
    else record2.flag &= ~BAM_FLAG_NEXT_RC;

    // Fix first/second in pair flags.
    if (hasFlagFirst(record1)) record2.flag |= BAM_FLAG_LAST;
    if (hasFlagFirst(record2)) record1.flag |= BAM_FLAG_LAST;
    if (hasFlagLast(record1)) record2.flag |= BAM_FLAG_FIRST;
    if (hasFlagLast(record2)) record1.flag |= BAM_FLAG_FIRST;

    // Set flag paired.
    record1.flag |= BAM_FLAG_MULTIPLE;
    record2.flag |= BAM_FLAG_MULTIPLE;
}

// ==========================================================================

// Correct the reference ids of a BamAlignmentRecord for the concatenated header.
template<typename TNameStore>
inline void
readRecordAndCorrectRIds(BamAlignmentRecord & record,
        BamFileIn & stream,
        NameStoreCache<TNameStore> & nameStoreCache)
{
    readRecord(record, stream);

    if (record.rID != BamAlignmentRecord::INVALID_REFID)
    {
        CharString rName = contigNames(context(stream))[record.rID];
        getIdByName(record.rID, nameStoreCache, rName);
    }
    if (record.rNextId != BamAlignmentRecord::INVALID_REFID)
    {
        CharString rNextName = contigNames(context(stream))[record.rNextId];
        getIdByName(record.rNextId, nameStoreCache, rNextName);
    }
}

// ==========================================================================

inline void
mergeHeaders(BamHeader & header,
        FormattedFileContext<BamFileOut, Owner<> >::Type & context,
        BamFileIn & stream1,
        BamFileIn & stream2)
{
    StringSet<CharString> contigNames;
    NameStoreCache<StringSet<CharString> > nameStoreCache;
    String<int32_t> contigLengths;

    // Read and append the two headers. Remove duplicate entries.
    readHeader(header, stream1);
    BamHeader header2;
    readHeader(header2, stream2);
    for (unsigned i = 0; i < length(header2); ++i)
    {
        if (header2[i].type != BAM_HEADER_FIRST)
            appendValue(header, header2[i]);
    }
    std::stable_sort(begin(header, Standard()), end(header, Standard()), BamHeaderRecordTypeLess());

    // Fill sequence names into nameStoreCache.
    for (unsigned i = 0; i < length(header); ++i)
    {
        if (header[i].type == BAM_HEADER_REFERENCE)
        {
            CharString name, len;
            for (unsigned j = 0; j < length(header[i].tags); ++j)
            {
                if (header[i].tags[j].i1 == "SN")
                    name = header[i].tags[j].i2;
                else if (header[i].tags[j].i1 == "LN")
                    len = header[i].tags[j].i2;
            }
            appendName(context._contigNamesCache, name);
            int32_t l;
            lexicalCast<int32_t>(l, len);
            appendValue(context._contigLengths, l);
        }
    }
}

// ==========================================================================

// This function is adapted from samtools code (fuction strnum_cmp in bam_sort.c) to ensure the exact same sort order.
int
compare_qName(CharString & nameA, CharString & nameB)
{
    const char * _a = toCString(nameA);
    const char * _b = toCString(nameB);
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

// ==========================================================================
// Function merge_and_set_mate()
// ==========================================================================

bool
merge_and_set_mate(CharString & mergedBam, CharString & nonRefBam, CharString & remappedBam)
{
    std::ostringstream msg;
    msg << "Merging bam files " << nonRefBam << " and " << remappedBam;
    printStatus(msg);

    // Open the two input streams (can read SAM and BAM files).
    BamFileIn nonRefStream(toCString(nonRefBam));
    BamFileIn remappedStream(toCString(remappedBam));

    printStatus(" - merging headers...");

    // Prepare a header for the output file.
    BamHeader outHeader;
    FormattedFileContext<BamFileOut, Owner<> >::Type context;
    mergeHeaders(outHeader, context, nonRefStream, remappedStream);

    printStatus(" - writing header...");

    // Open the output stream and write the header.
    FormattedFileContext<BamFileOut, Dependent<> >::Type contextDep(context);
    BamFileOut outStream(contextDep, toCString(mergedBam));
    writeHeader(outStream, outHeader);

    printStatus(" - merging read records...");

    // Read the first record from each input file. Correct ids in records from remappedStreams for new header.
    BamAlignmentRecord record1, record2;
    if (!atEnd(nonRefStream)) readRecordAndCorrectRIds(record1, nonRefStream, contigNamesCache(contextDep));
    else record1.qName = "*";
    if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, contigNamesCache(contextDep));
    else record2.qName = "*";

    // Iterate both input files, set mate positions in pairs, and write all records to the output file.
    while (record1.qName != "*" || record2.qName != "*")
    {
        while ((compare_qName(record2.qName, record1.qName) < 0 || record1.qName == "*") && record2.qName != "*")
        {
            writeRecord(outStream, record2);
            if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, contigNamesCache(contextDep));
            else record2.qName = "*";
        }

        bool incr1 = false;
        while (record1.qName == record2.qName && record2.qName != "*")
        {
            incr1 = true;
            setMates(record1, record2);
            writeRecord(outStream, record1);
            writeRecord(outStream, record2);
            if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, contigNamesCache(contextDep));
            else record2.qName = "*";
        }
        if (incr1)
        {
            if (!atEnd(nonRefStream)) readRecordAndCorrectRIds(record1, nonRefStream, contigNamesCache(contextDep));
            else record1.qName = "*";
        }

        while ((compare_qName(record1.qName, record2.qName) < 0 || record2.qName == "*") && record1.qName != "*")
        {
            writeRecord(outStream, record1);
            if (!atEnd(nonRefStream)) readRecordAndCorrectRIds(record1, nonRefStream, contigNamesCache(contextDep));
            else record1.qName = "*";
        }
    }

    return 0;
}

// ==========================================================================
// Function sickle_filtering()
// ==========================================================================

inline bool
sickle_filtering(Triple<CharString> & filteredFiles,
        Triple<CharString> & fastqFiles,
        CharString & workingDirectory)
{
    std::stringstream cmd;

    std::ostringstream msg;
    msg << "Filtering fastq files using " << SICKLE ;
    printStatus(msg);

    cmd.str("");
    cmd << SICKLE << " pe -q 20 -l 60 -x -n -t sanger -f" << fastqFiles.i1 << " -r " << fastqFiles.i2;
    cmd << " -o " << filteredFiles.i1 << " -p " << filteredFiles.i2 << " -s " << filteredFiles.i3;
    if (system(cmd.str().c_str()) != 0) // runs sickle on paired end reads
    {
        std::cerr << "ERROR while filtering " << fastqFiles.i1 << " and " << fastqFiles.i2;
        std::cerr << " using " << SICKLE << std::endl;
        return 1;
    }

    CharString singleFiltered2 = getFileName(workingDirectory, "filtered.single2.fastq");

    cmd.str("");
    cmd << SICKLE << " se -q 20 -l 60 -x -n -t sanger -f " << fastqFiles.i3 << " -o " << singleFiltered2;
    if (system(cmd.str().c_str()) != 0) // runs sickle on single end reads
    {
        std::cerr << "ERROR while filtering " << fastqFiles.i3 << " using " << SICKLE << std::endl;
        return 1;
    }

    cmd.str("");
    cmd << "cat " << singleFiltered2 << " >> " << filteredFiles.i3;
    if (system(cmd.str().c_str()) != 0) // merges the fastq files of single end reads.
    {
        std::cerr << "ERROR while concatenating " << singleFiltered2 << " to " << filteredFiles.i3 << std::endl;
        return 1;
    }

    remove(toCString(singleFiltered2));

    return 0;
}

// ==========================================================================
// Function velvet_assembly()
// ==========================================================================

inline bool
velvet_assembly(Triple<CharString> & filteredFiles, Triple<CharString> & filteredMPFiles, CharString & assemblyDirectory, unsigned kmerLength, bool matepair)
{
    std::stringstream cmd;

    std::ostringstream msg;
    msg << "Preparing assembly of unmapped reads from filtered fastq files using " << VELVETH;
    printStatus(msg);

    cmd.str("");
    cmd << VELVETH << " " << assemblyDirectory << " " << kmerLength << " -short -fastq " << filteredFiles.i3;
    cmd << " -shortPaired -fastq -separate " << filteredFiles.i1 << " " << filteredFiles.i2;
    if (matepair) {
        cmd << " -shortPaired2 -fastq -separate " << filteredMPFiles.i1 << " " << filteredMPFiles.i2;
    }

    if (system(cmd.str().c_str()) != 0) // prepares velvet assembly, use k=47 for longer contigs
    {
        std::cerr << "ERROR while preparing assembly with " << VELVETH << " of ";
        std::cerr << filteredFiles.i3 << ", " << filteredFiles.i1 << ", and " << filteredFiles.i2 << std::endl;
        if (matepair) {
            std::cerr << "and matepair files " << filteredMPFiles.i1 << " and " << filteredMPFiles.i2 << std::endl;
        }
        return 1;
    }

    msg.str("");
    msg << "Assembling unmapped reads from filtered fastq files using " << VELVETG;
    printStatus(msg);

    cmd.str("");
    cmd << VELVETG << " " << assemblyDirectory << " -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no";
    if (matepair) {
        cmd << " -shortMatePaired2 yes";
    }

    if (system(cmd.str().c_str()) != 0) // runs the velvet graph part
    {
        std::cerr << "ERROR while assembling " << assemblyDirectory << " with " << VELVETG << std::endl;
        return 1;
    }

    return 0;
}

// ==========================================================================
// Function popins_assemble()
// ==========================================================================

int popins_assemble(int argc, char const ** argv)
{
    std::ostringstream msg;

    // Parse the command line to get option values.
    AssemblyOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Retrieve the sample ID from the first read group listed in BAM file header.
    if (options.sampleID == "" && retrieveSampleID(options.sampleID, options.mappingFile) == 1)
        return 7;

    // Create working directory if it does not exist.
    CharString workingDirectory = getFileName(options.prefix, options.sampleID);
    if (mkdir(toCString(workingDirectory), 0755) == 0)
    {
        msg.str("");
        msg << "Working directory created at " << workingDirectory;
        printStatus(msg);
    }

    SampleInfo info = initSampleInfo(options.mappingFile, options.sampleID, options.adapters);

    CharString matesBam = getFileName(workingDirectory, "mates.bam");
    CharString nonRefBamTemp = getFileName(workingDirectory, "non_ref_tmp.bam");
    CharString nonRefBam = getFileName(workingDirectory, "non_ref.bam");

    CharString fastqFirst = getFileName(workingDirectory, "paired.1.fastq");
    CharString fastqSecond = getFileName(workingDirectory, "paired.2.fastq");
    CharString fastqSingle = getFileName(workingDirectory, "single.fastq");
    Triple<CharString> fastqFiles = Triple<CharString>(fastqFirst, fastqSecond, fastqSingle);

    // check if files already exits
    std::fstream stream(toCString(fastqFirst));
    if (!stream.is_open())
    {
        msg.str("");
        msg << "Cropping unmapped reads from " << options.mappingFile;
        printStatus(msg);

        // Crop unmapped reads and reads with unreliable mappings from the input bam file.
        if (options.adapters == "HiSeqX")
        {
            if (crop_unmapped(info.avg_cov, fastqFiles, matesBam, options.mappingFile, options.humanSeqs, HiSeqXAdapters()) != 0)
                return 7;
        }
        else if (options.adapters == "HiSeq")
        {
            if (crop_unmapped(info.avg_cov, fastqFiles, matesBam, options.mappingFile, options.humanSeqs, HiSeqAdapters()) != 0)
                return 7;
        }
        else
        {
            if (crop_unmapped(info.avg_cov, fastqFiles, matesBam, options.mappingFile, options.humanSeqs, NoAdapters()) != 0)
                return 7;
        }

        CharString sampleInfoFile = getFileName(workingDirectory, "POPINS_SAMPLE_INFO");
        writeSampleInfo(info, sampleInfoFile);

        msg.str("");
        msg << "Sample info written to \'" << sampleInfoFile << "\'.";
        printStatus(msg);

        msg.str("");
        msg << "Sorting " << matesBam << " using " << SAMTOOLS;
        printStatus(msg);

        // Sort <WD>/mates.bam by read name.
        std::stringstream cmd;
        if (options.referenceFile != "")
            cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " -o " << nonRefBamTemp << " " << matesBam;
        else
            cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " -o " << nonRefBam << " " << matesBam;
        if (system(cmd.str().c_str()) != 0)
        {
            std::cerr << "ERROR while sorting " << matesBam << std::endl;
            return 7;
        }

        remove(toCString(matesBam));

        // Remapping of unmapped with bwa if a fasta reference is given.
        if (options.referenceFile != "")
        {
            Triple<CharString> fastqFilesTemp = fastqFiles;
            fastqFiles = Triple<CharString>(fastqFirst, fastqSecond, fastqSingle);

            // Align with bwa, update fastq files of unaligned reads, and sort remaining bam records by read name.
            CharString remappedBam = getFileName(workingDirectory, "remapped.bam");
            CharString prefix = "";
            if (remapping(fastqFilesTemp, fastqFiles, options.referenceFile, workingDirectory,
                    options.humanSeqs, options.threads, options.memory, prefix) != 0)
                return 7;

            // Set the mate's location and merge non_ref.bam and remapped.bam into a single file.
            if (merge_and_set_mate(nonRefBam, nonRefBamTemp, remappedBam) != 0) return 7;
            remove(toCString(remappedBam));
            remove(toCString(nonRefBamTemp));
        }
    }
    else
    {
        printStatus("Found files, skipping cropping step.");
    }

    CharString firstFiltered = getFileName(workingDirectory, "filtered.paired.1.fastq");
    CharString secondFiltered = getFileName(workingDirectory, "filtered.paired.2.fastq");
    CharString singleFiltered = getFileName(workingDirectory, "filtered.single.fastq");
    Triple<CharString> filteredFiles(firstFiltered, secondFiltered, singleFiltered);

    // Quality filtering/trimming with sickle.
    if (sickle_filtering(filteredFiles, fastqFiles, workingDirectory) != 0)
        return 7;

    // MP handling
    CharString matesMPBam = getFileName(workingDirectory, "MP.mates.bam");
    CharString nonRefBamMPTemp = getFileName(workingDirectory, "MP.non_ref_tmp.bam");
    CharString nonRefMPBam = getFileName(workingDirectory, "MP.non_ref.bam");

    CharString fastqMPFirst = getFileName(workingDirectory, "MP.paired.1.fastq");
    CharString fastqMPSecond = getFileName(workingDirectory, "MP.paired.2.fastq");
    CharString fastqMPSingle = getFileName(workingDirectory, "MP.single.fastq");
    Triple<CharString> fastqMPFiles = Triple<CharString>(fastqMPFirst, fastqMPSecond, fastqMPSingle);

    if (options.matepairFile != "")
    {
        // check if MP files already exits
        std::fstream MPstream(toCString(fastqMPFirst));
        if (!MPstream.is_open())
        {
            msg.str("");
            msg << "Cropping unmapped matepair reads from " << options.matepairFile;
            printStatus(msg);

            // Crop unmapped reads and reads with unreliable mappings from the input bam file.
            if (options.adapters == "HiSeqX")
            {
                if (crop_unmapped(fastqMPFiles, matesMPBam, options.matepairFile, options.humanSeqs, HiSeqXAdapters()) != 0)
                    return 7;
            }
            else if (options.adapters == "HiSeq")
            {
                if (crop_unmapped(fastqMPFiles, matesMPBam, options.matepairFile, options.humanSeqs, HiSeqAdapters()) != 0)
                    return 7;
            }
            else
            {
                if (crop_unmapped(fastqMPFiles, matesMPBam, options.matepairFile, options.humanSeqs, NoAdapters()) != 0)
                    return 7;
            }

            msg.str("");
            msg << "Sorting " << matesMPBam << " using " << SAMTOOLS;
            printStatus(msg);

            // Sort <WD>/mates.bam by read name.
            std::stringstream cmd;
            if (options.referenceFile != "")
                cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " -o " << nonRefBamMPTemp << " " << matesMPBam;
            else
                cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " -o " << nonRefMPBam << " " << matesMPBam;
            if (system(cmd.str().c_str()) != 0)
            {
                std::cerr << "ERROR while sorting " << matesMPBam << std::endl;
                return 7;
            }

            remove(toCString(matesMPBam));

            // Remapping of unmapped with bwa if a fasta reference is given.
            if (options.referenceFile != "")
            {
                Triple<CharString> fastqMPFilesTemp = fastqMPFiles;
                fastqMPFiles = Triple<CharString>(fastqMPFirst, fastqMPSecond, fastqMPSingle);

                // Align with bwa, update fastq files of unaligned reads, and sort remaining bam records by read name.
                CharString remappedMPBam = getFileName(workingDirectory, "MP.remapped.bam");
                CharString prefix = "MP.";
                if (remapping(fastqMPFilesTemp, fastqMPFiles, options.referenceFile, workingDirectory,
                        options.humanSeqs, options.threads, options.memory, prefix) != 0)
                    return 7;

                // Set the mate's location and merge non_ref.bam and remapped.bam into a single file.
                if (merge_and_set_mate(nonRefMPBam, nonRefBamMPTemp, remappedMPBam) != 0)
                    return 7;
                remove(toCString(remappedMPBam));
                remove(toCString(nonRefBamMPTemp));
            }
        }
        else
        {
            printStatus("Found matepair files, skipping cropping step");
        }
    }


    CharString firstMPFiltered = getFileName(workingDirectory, "MP.filtered.paired.1.fastq");
    CharString secondMPFiltered = getFileName(workingDirectory, "MP.filtered.paired.2.fastq");
    CharString singleMPFiltered = getFileName(workingDirectory, "MP.filtered.single.fastq");

    Triple<CharString> filteredMPFiles(firstMPFiltered, secondMPFiltered, singleMPFiltered);

    if (options.matepairFile != "")
    {
        if (sickle_filtering(filteredMPFiles,fastqMPFiles, workingDirectory) != 0)
            return 7;
    }

    // Assembly with velvet.
    CharString assemblyDirectory = getFileName(workingDirectory, "assembly");
    if (velvet_assembly(filteredFiles, filteredMPFiles, assemblyDirectory, options.kmerLength, options.matepairFile == "") != 0)
        return 7;

    remove(toCString(firstFiltered));
    remove(toCString(secondFiltered));
    remove(toCString(singleFiltered));

    if (options.matepairFile == "")
    {
        remove(toCString(firstMPFiltered));
        remove(toCString(secondMPFiltered));
        remove(toCString(singleMPFiltered));
    }

    // Copy contigs file to workingDirectory and remove assembly directory.
    CharString contigFileAssembly = getFileName(assemblyDirectory, "contigs.fa");
    CharString contigFile = getFileName(workingDirectory, "contigs.fa");
    std::ifstream src(toCString(contigFileAssembly), std::ios::binary);
    std::ofstream dst(toCString(contigFile), std::ios::binary);
    dst << src.rdbuf();
    src.close();
    dst.close();
    removeAssemblyDirectory(assemblyDirectory);

    return res;
}

#endif // #ifndef POPINS_ASSEMBLE_H_
