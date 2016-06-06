#include <sstream>
#include <cerrno>

#include <seqan/file.h>
#include <seqan/sequence.h>

#include "popins_clp.h"
#include "popins_crop_unmapped.h"

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

// ==========================================================================
// Function remapping()
// ==========================================================================

inline int
remapping(Triple<CharString> & fastqFilesTemp,
        Triple<CharString> & fastqFiles,
        CharString const & referenceFile,
        CharString const & tempDir,
        unsigned humanSeqs,
        unsigned threads,
        CharString & memory,
        CharString & prefix)
{
    std::stringstream cmd;

    CharString f1 = prefix;
    f1 += "remapped.sam";
    CharString remappedSam = getFileName(tempDir, f1);
    CharString f2 = prefix;
    f2 += "remapped.bam";
    CharString remappedBam = getFileName(tempDir, f2);
    CharString f3 = prefix;
    f3 += "remapped.bam.bai";
    CharString remappedBai = getFileName(tempDir, f3);
    CharString f4 = prefix;
    f4 += "remapped_unsorted.bam";
    CharString remappedUnsortedBam = getFileName(tempDir, f4);

    // Run BWA on unmapped reads (pairs).
    std::cerr << "[" << time(0) << "] Remapping unmapped reads using " << BWA << std::endl;
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

    // Convert BWA output to bam.
    std::cerr << "[" << time(0) << "] Converting BWA output " << remappedSam << " to bam format." << std::endl;
    cmd.str("");
    cmd << SAMTOOLS << " view -S -h -b " << remappedSam << " > " << remappedUnsortedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while converting BWA output " << remappedSam << " to bam." << std::endl;
        return 1;
    }
    remove(toCString(remappedSam));

    // Sort bam file.
    std::cerr << "[" << time(0) << "] Sorting " << remappedUnsortedBam << " using " << SAMTOOLS << std::endl;
    cmd.str("");
    cmd << SAMTOOLS << " sort -@ " << threads << " -m " << memory << " " << remappedUnsortedBam << " " << tempDir << "/remapped";
    cmd << SAMTOOLS << " sort -@ " << threads << " -m " << memory << " " << remappedUnsortedBam << " " << tempDir << "/" << prefix << "remapped";
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting BWA output " << remappedUnsortedBam << std::endl;
        return 1;
    }

    // Index bam file.
    std::cerr << "[" << time(0) << "] Indexing " << remappedBam << " using " << SAMTOOLS << std::endl;
    cmd.str("");
    cmd << SAMTOOLS << " index " << remappedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while indexing BWA output " << remappedBam << std::endl;
        return 1;
    }

    // Crop unmapped and create bam file of remapping.
    std::cerr << "[" << time(0) << "] Cropping unmapped reads from " << remappedBam << std::endl;
    if (crop_unmapped(fastqFiles, remappedUnsortedBam, remappedBam, humanSeqs, NoAdapters()) != 0)
        return 1;
    remove(toCString(remappedBai));

    // Sort <WD>/remapped.bam by read name.
    std::cerr << "[" << time(0) << "] " << "Sorting " << remappedUnsortedBam << " by read name using " << SAMTOOLS << std::endl;

    cmd.str("");
    cmd << SAMTOOLS << " sort -n -@ " << threads << " -m " << memory << " " << remappedUnsortedBam << " " << tempDir << "/remapped";
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
        BamStream & stream,
        TNameStore & nameStor,
        NameStoreCache<TNameStore> & nameStoreCache)
{
    readRecord(record, stream);

    if (record.rID != BamAlignmentRecord::INVALID_REFID)
    {
        CharString rName = nameStore(stream.bamIOContext)[record.rID];
        getIdByName(nameStor, rName, record.rID, nameStoreCache);
    }
    if (record.rNextId != BamAlignmentRecord::INVALID_REFID)
    {
        CharString rNextName = stream.header.sequenceInfos[record.rNextId].i1;
        getIdByName(nameStor, rNextName, record.rNextId, nameStoreCache);
    }
}

// ==========================================================================

inline void
mergeHeaders(BamHeader & header,
        StringSet<CharString> & nameStore,
        NameStoreCache<StringSet<CharString> > & nameStoreCache,
        BamHeader const & header1,
        BamHeader const & header2)
{
    // Write sequenceInfos, name store and cache for this header.
    unsigned idx = 0;
    for (unsigned i = 0; i < length(header1.sequenceInfos); ++i)
    {
        if (!getIdByName(nameStore, header1.sequenceInfos[i].i1, idx, nameStoreCache))
        {
            appendName(nameStore, header1.sequenceInfos[i].i1, nameStoreCache);
            appendValue(header.sequenceInfos, header1.sequenceInfos[i]);
        }
    }

    for (unsigned i = 0; i < length(header2.sequenceInfos); ++i)
    {
        if (!getIdByName(nameStore, header2.sequenceInfos[i].i1, idx, nameStoreCache))
        {
            appendName(nameStore, header2.sequenceInfos[i].i1, nameStoreCache);
            appendValue(header.sequenceInfos, header2.sequenceInfos[i]);
        }
    }

    // Copy all other than sequence records from header1 and header2.
    for (unsigned i = 0; i < length(header1.records); ++i)
    {
        if (header1.records[i].type != BamHeaderRecordType::BAM_HEADER_REFERENCE)
            appendValue(header.records, header1.records[i]);
    }
    for (unsigned i = 0; i < length(header2.records); ++i)
    {
        if (header2.records[i].type != BamHeaderRecordType::BAM_HEADER_REFERENCE && header2.records[i].type != BamHeaderRecordType::BAM_HEADER_FIRST)
            appendValue(header.records, header2.records[i]);
    }

    // Add a record for each sequence in sequenceInfos.
    for (unsigned i = 0; i < length(header.sequenceInfos); ++i)
    {
        BamHeaderRecord r;
        r.type = BamHeaderRecordType::BAM_HEADER_REFERENCE;
        setTagValue("SN", header.sequenceInfos[i].i1, r);
        std::stringstream ss;
        ss << header.sequenceInfos[i].i2;
        setTagValue("LN", ss.str(), r);
        appendValue(header.records, r);
    }

    BamHeaderRecordTypeLess less;
    std::stable_sort(begin(header.records, Standard()), end(header.records, Standard()), less);
}

// ==========================================================================

int
compare_qName(CharString & a, CharString & b)
{
    typedef Iterator<CharString, Rooted>::Type TIter;

    TIter itEndA = end(a);
    TIter itEndB = end(b);

    TIter itA = begin(a);
    TIter itB = begin(b);

    while (itA != itEndA && itB != itEndB)
    {
        if (std::isdigit(*itA) && std::isdigit(*itB))
        {
            while (itA != itEndA && *itA == '0') ++itA;
            while (itB != itEndB && *itB == '0') ++itB;
            while (itA != itEndA && std::isdigit(*itA) && itB != itEndB && std::isdigit(*itB) && *itA == *itB) ++itA, ++itB;
            if (std::isdigit(*itA) && std::isdigit(*itB)) // pointing to two different digits
            {
                int dig_a = *itA;
                int dig_b = *itB;
                while (itA != itEndA && std::isdigit(*itA) && itB != itEndB && std::isdigit(*itB)) ++itA, ++itB; // is one number longer
                return std::isdigit(*itA) ? 1 : std::isdigit(*itB) ? -1 : dig_a - dig_b;
            }
            else if (std::isdigit(*itA)) return 1;
            else if (std::isdigit(*itB)) return -1;
            else if (position(itA) != position(itB)) return position(itA) < position(itB) ? 1 : -1;
        }
        else
        {
            if (*itA != *itB) return (int)*itA - (int)*itB;
            ++itA; ++itB;
        }
    }
    return itA != itEndA ? 1 : itB != itEndB ? -1 : 0;
}

// ==========================================================================
// Function merge_and_set_mate()
// ==========================================================================

bool
merge_and_set_mate(CharString & mergedBam, CharString & nonRefBam, CharString & remappedBam)
{
    typedef StringSet<CharString> TNameStore;

    std::cerr << "[" << time(0) << "] " << "Merging bam files " << nonRefBam << " and " << remappedBam << std::endl;

    // Open the two input streams, can read SAM and BAM files.
    BamStream nonRefStream(toCString(nonRefBam));
    if (!isGood(nonRefStream))
    {
        std::cerr << "ERROR: Could not open input bam/sam file " << nonRefBam << "." << std::endl;
        return 1;
    }
    std::cerr << "[" << time(0) << "]  - opened non_ref file " << nonRefBam << std::endl;

    BamStream remappedStream(toCString(remappedBam));
    if (!isGood(remappedStream))
    {
        std::cerr << "ERROR: Could not open input bam/sam file " << remappedBam << "." << std::endl;
        return 1;
    }
    std::cerr << "[" << time(0) << "]  - opened remapped file " << remappedBam << std::endl;

    // Open output file.
    BamStream outStream(toCString(mergedBam), BamStream::WRITE);
    if (!isGood(remappedStream))
    {
        std::cerr << "ERROR: Could not open output bam file " << mergedBam << "." << std::endl;
        return 1;
    }
    std::cerr << "[" << time(0) << "]  - opened output file " << mergedBam << std::endl;

    // Prepare a header for the output files.
    TNameStore nameStor;
    NameStoreCache<TNameStore> nameStoreCache(nameStor);
    mergeHeaders(outStream.header, nameStor, nameStoreCache, nonRefStream.header, remappedStream.header);
    std::cerr << "[" << time(0) << "]  - merged headers" << std::endl;

    // Read the first record from each input file. Correct ids in records from remappedStreams for new header.
    BamAlignmentRecord record1, record2;
    if (!atEnd(nonRefStream)) readRecord(record1, nonRefStream);
    else record1.qName = "*";
    if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, nameStor, nameStoreCache);
    else record2.qName = "*";

    // Iterate both input files, set mate positions in pairs, and write all records to the output file.
    while (record1.qName != "*" || record2.qName != "*")
    {
        while ((compare_qName(record2.qName, record1.qName) < 0 || record1.qName == "*") && record2.qName != "*")
        {
            writeRecord(outStream, record2);
            if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, nameStor, nameStoreCache);
            else record2.qName = "*";
        }

        bool incr1 = false;
        while (record1.qName == record2.qName && record2.qName != "*")
        {
            incr1 = true;
            setMates(record1, record2);
            writeRecord(outStream, record1);
            writeRecord(outStream, record2);
            if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, nameStor, nameStoreCache);
            else record2.qName = "*";
        }
        if (incr1)
        {
            if (!atEnd(nonRefStream)) readRecord(record1, nonRefStream);
            else record1.qName = "*";
        }

        while ((compare_qName(record1.qName, record2.qName) < 0 || record2.qName == "*") && record1.qName != "*")
        {
            writeRecord(outStream, record1);
            if (!atEnd(nonRefStream)) readRecord(record1, nonRefStream);
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

    std::cerr << "[" << time(0) << "] " << "Filtering fastq files using " << SICKLE << std::endl;
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

    std::cerr << "[" << time(0) << "] " << "Preparing assembly of unmapped reads from filtered fastq files using " << VELVETH << std::endl;
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

    std::cerr << "[" << time(0) << "] " << "Assembling unmapped reads from filtered fastq files using " << VELVETG << std::endl;
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
    bool ret = 0;

    // Parse the command line to get option values.
    AssemblyOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    // Create working directory if it does not exist.
    if (mkdir(toCString(options.workingDirectory), 0755) == 0)
    {
        CharString logMsg = "Working directory created at ";
        logMsg += options.workingDirectory;
        std::cerr << "[" << time(0) << "] " << logMsg << std::endl;
    }

    CharString matesBam = getFileName(options.workingDirectory, "mates.bam");

    CharString fastqFirstTemp = getFileName(options.workingDirectory, "paired.1.fastq");
    CharString fastqSecondTemp = getFileName(options.workingDirectory, "paired.2.fastq");
    CharString fastqSingleTemp = getFileName(options.workingDirectory, "single.fastq");
    CharString nonRefBamTemp = getFileName(options.workingDirectory, "non_ref_tmp.bam");

    CharString fastqFirst = getFileName(options.workingDirectory, "paired.1.fastq");
    CharString fastqSecond = getFileName(options.workingDirectory, "paired.2.fastq");
    CharString fastqSingle = getFileName(options.workingDirectory, "single.fastq");
    CharString nonRefBam = getFileName(options.workingDirectory, "non_ref.bam");

    Triple<CharString> fastqFiles;
    if (options.referenceFile != "")
        fastqFiles = Triple<CharString>(fastqFirstTemp, fastqSecondTemp, fastqSingleTemp);
    else
        fastqFiles = Triple<CharString>(fastqFirst, fastqSecond, fastqSingle);


    // check if files already exits
    std::fstream stream(toCString(fastqFirst));
    if (!stream.is_open()) {
        // Crop unmapped reads and reads with unreliable mappings from the input bam file.
        std::cerr << "[" << time(0) << "] Cropping unmapped reads from " << options.mappingFile << std::endl;
        if (options.adapters == "HiSeqX")
        {
            if (crop_unmapped(fastqFiles, matesBam, options.mappingFile, options.humanSeqs, HiSeqXAdapters()) != 0)
                return 1;
        }
        else if (options.adapters == "HiSeq")
        {
            if (crop_unmapped(fastqFiles, matesBam, options.mappingFile, options.humanSeqs, HiSeqAdapters()) != 0)
                return 1;
        }
        else
        {
            if (crop_unmapped(fastqFiles, matesBam, options.mappingFile, options.humanSeqs, NoAdapters()) != 0)
                return 1;
        }

        // Sort <WD>/mates.bam by read name.
        std::cerr << "[" << time(0) << "] " << "Sorting " << matesBam << " using " << SAMTOOLS << std::endl;
        std::stringstream cmd;
        if (options.referenceFile != "")
            cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " " << matesBam << " " << options.workingDirectory << "/non_ref_tmp";
        else
            cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " " << matesBam << " " << options.workingDirectory << "/non_ref";
        if (system(cmd.str().c_str()) != 0)
        {
            std::cerr << "ERROR while sorting " << matesBam << std::endl;
            return 1;
        }

        remove(toCString(matesBam));

        // Remapping of unmapped with bwa if a fasta reference is given.
        if (options.referenceFile != "")
        {
            Triple<CharString> fastqFilesTemp = fastqFiles;
            fastqFiles = Triple<CharString>(fastqFirst, fastqSecond, fastqSingle);

            // Align with bwa, update fastq files of unaligned reads, and sort remaining bam records by read name.
            CharString remappedBam = getFileName(options.workingDirectory, "remapped.bam");
            CharString prefix = "";
            if (remapping(fastqFilesTemp, fastqFiles, options.referenceFile, options.workingDirectory,
                    options.humanSeqs, options.threads, options.memory, prefix) != 0)
                return 1;

            // Set the mate's location and merge non_ref.bam and remapped.bam into a single file.
            if (merge_and_set_mate(nonRefBam, nonRefBamTemp, remappedBam) != 0) return 1;
            remove(toCString(remappedBam));
            remove(toCString(nonRefBamTemp));
        }
    } else {
        std::cerr << "[" << time(0) << "] Found files, skipping cropping step" << std::endl;
    }


    CharString firstFiltered = getFileName(options.workingDirectory, "filtered.paired.1.fastq");
    CharString secondFiltered = getFileName(options.workingDirectory, "filtered.paired.2.fastq");
    CharString singleFiltered = getFileName(options.workingDirectory, "filtered.single.fastq");
    Triple<CharString> filteredFiles(firstFiltered, secondFiltered, singleFiltered);

    // Quality filtering/trimming with sickle.
    if (sickle_filtering(filteredFiles, fastqFiles, options.workingDirectory) != 0)
        return 1;


    // MP handling
    CharString matesMPBam = getFileName(options.workingDirectory, "MP.mates.bam");
    CharString fastqMPFirstTemp = getFileName(options.workingDirectory, "MP.paired.1.fastq");
    CharString fastqMPSecondTemp = getFileName(options.workingDirectory, "MP.paired.2.fastq");
    CharString fastqMPSingleTemp = getFileName(options.workingDirectory, "MP.single.fastq");
    CharString nonRefBamMPTemp = getFileName(options.workingDirectory, "MP.non_ref_tmp.bam");

    CharString fastqMPFirst = getFileName(options.workingDirectory, "MP.paired.1.fastq");
    CharString fastqMPSecond = getFileName(options.workingDirectory, "MP.paired.2.fastq");
    CharString fastqMPSingle = getFileName(options.workingDirectory, "MP.single.fastq");
    CharString nonRefMPBam = getFileName(options.workingDirectory, "MP.non_ref.bam");

    Triple<CharString> fastqMPFiles;
    if (options.referenceFile != "")
        fastqMPFiles = Triple<CharString>(fastqMPFirstTemp, fastqMPSecondTemp, fastqMPSingleTemp);
    else
        fastqMPFiles = Triple<CharString>(fastqMPFirst, fastqMPSecond, fastqMPSingle);


    if (options.matepair) {
        // check if MP files already exits
        std::fstream MPstream(toCString(fastqMPFirst));
        if (!MPstream.is_open()) {
            // Crop unmapped reads and reads with unreliable mappings from the input bam file.
            std::cerr << "[" << time(0) << "] Cropping unmapped matepair reads from " << options.matepairFile << std::endl;
            if (options.adapters == "HiSeqX")
            {
                if (crop_unmapped(fastqMPFiles, matesMPBam, options.matepairFile, options.humanSeqs, HiSeqXAdapters()) != 0)
                    return 1;
            }
            else if (options.adapters == "HiSeq")
            {
                if (crop_unmapped(fastqMPFiles, matesMPBam, options.matepairFile, options.humanSeqs, HiSeqAdapters()) != 0)
                    return 1;
            }
            else
            {
                if (crop_unmapped(fastqMPFiles, matesMPBam, options.matepairFile, options.humanSeqs, NoAdapters()) != 0)
                    return 1;
            }

            // Sort <WD>/mates.bam by read name.
            std::cerr << "[" << time(0) << "] " << "Sorting " << matesMPBam << " using " << SAMTOOLS << std::endl;
            std::stringstream cmd;
            if (options.referenceFile != "")
                cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " " << matesMPBam << " " << options.workingDirectory << "/MP.non_ref_tmp";
            else
                cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " " << matesMPBam << " " << options.workingDirectory << "/MP.non_ref";
            if (system(cmd.str().c_str()) != 0)
            {
                std::cerr << "ERROR while sorting " << matesMPBam << std::endl;
                return 1;
            }

            remove(toCString(matesMPBam));

            // Remapping of unmapped with bwa if a fasta reference is given.
            if (options.referenceFile != "")
            {
                Triple<CharString> fastqMPFilesTemp = fastqMPFiles;
                fastqMPFiles = Triple<CharString>(fastqMPFirst, fastqMPSecond, fastqMPSingle);

                // Align with bwa, update fastq files of unaligned reads, and sort remaining bam records by read name.
                CharString remappedMPBam = getFileName(options.workingDirectory, "MP.remapped.bam");
                CharString prefix = "MP.";
                if (remapping(fastqMPFilesTemp, fastqMPFiles, options.referenceFile, options.workingDirectory,
                        options.humanSeqs, options.threads, options.memory, prefix) != 0)
                    return 1;

                // Set the mate's location and merge non_ref.bam and remapped.bam into a single file.
                if (merge_and_set_mate(nonRefMPBam, nonRefBamMPTemp, remappedMPBam) != 0) return 1;
                remove(toCString(remappedMPBam));
                remove(toCString(nonRefBamMPTemp));
            }
        } else {
            std::cerr << "[" << time(0) << "] Found matepair files, skipping cropping step" << std::endl;
        }
    }


    CharString firstMPFiltered = getFileName(options.workingDirectory, "MP.filtered.paired.1.fastq");
    CharString secondMPFiltered = getFileName(options.workingDirectory, "MP.filtered.paired.2.fastq");
    CharString singleMPFiltered = getFileName(options.workingDirectory, "MP.filtered.single.fastq");

    Triple<CharString> filteredMPFiles(firstMPFiltered, secondMPFiltered, singleMPFiltered);

    if (options.matepair) {
        if (sickle_filtering(filteredMPFiles,fastqMPFiles, options.workingDirectory) != 0)
            return 1;
    }

    // Assembly with velvet.
    CharString assemblyDirectory = getFileName(options.workingDirectory, "assembly");
    if (velvet_assembly(filteredFiles, filteredMPFiles, assemblyDirectory, options.kmerLength, options.matepair) != 0) return 1;

    remove(toCString(firstFiltered));
    remove(toCString(secondFiltered));
    remove(toCString(singleFiltered));

    if (options.matepair) {
        remove(toCString(firstMPFiltered));
        remove(toCString(secondMPFiltered));
        remove(toCString(singleMPFiltered));
    }

    // Copy contigs file to workingDirectory.
    CharString contigFileAssembly = getFileName(assemblyDirectory, "contigs.fa");
    CharString contigFile = getFileName(options.workingDirectory, "contigs.fa");
    std::ifstream src(toCString(contigFileAssembly), std::ios::binary);
    std::ofstream dst(toCString(contigFile), std::ios::binary);
    dst << src.rdbuf();
    src.close();
    dst.close();

    removeAssemblyDirectory(assemblyDirectory);

    return ret;
}

#endif // #ifndef POPINS_ASSEMBLE_H_
