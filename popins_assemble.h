#include <sstream>

#include <seqan/file.h>
#include <seqan/sequence.h>

#include "popins_clp.h"
#include "popins_crop_unmapped.h"

#ifndef POPINS_ASSEMBLE_H_
#define POPINS_ASSEMBLE_H_

using namespace seqan;

// ==========================================================================

inline void
removeAssemblyDirectory(CharString const & path)
{
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
remapping(CharString & fastqFirst,
          CharString & fastqSecond,
          CharString & fastqSingle,
          CharString const & referenceFile,
          CharString const & workingDirectory,
          unsigned humanSeqs,
          unsigned threads,
          CharString & memory)
{
    std::stringstream cmd;
    
    CharString remappedSam = getFileName(workingDirectory, "remapped.sam");
    CharString remappedBam = getFileName(workingDirectory, "remapped.bam");
    CharString remappedUnsortedBam = getFileName(workingDirectory, "remapped_unsorted.bam");
    
    // Run BWA on unmapped reads (pairs).
    std::cerr << "[" << time(0) << "] Remapping unmapped reads using " << BWA << std::endl;
    cmd.str("");
    cmd << BWA << " mem -t " << threads << " " << referenceFile << " " << fastqFirst << " " << fastqSecond << " > " << remappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFirst << " and " << fastqSecond << std::endl;
        return 1;
    }
    
    // Run BWA on unmapped reads (single end).
    cmd.str("");
    cmd << BWA << " mem -t " << threads << " " << referenceFile << " " << fastqSingle << " | awk '$1 !~ /@/' >> " << remappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqSingle << std::endl;
        return 1;
    }
    
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
    std::cerr << "[" << time(0) << "] Sorting " << remappedBam << " using " << SAMTOOLS << std::endl;
    cmd.str("");
    cmd << SAMTOOLS << " sort -m " << memory << " " << remappedUnsortedBam << " " << workingDirectory << "/remapped";
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting BWA output " << remappedUnsortedBam << std::endl;
        return 1;
    }

    // Crop unmapped and create bam file of remapping.
    std::cerr << "[" << time(0) << "] Cropping unmapped reads from " << remappedBam << std::endl;
    if (crop_unmapped(fastqFirst, fastqSecond, fastqSingle, remappedUnsortedBam, remappedBam, humanSeqs, NoAdapters()) != 0)
        return 1;
    
    // Sort <WD>/remapped.bam by read name.
    std::cerr << "[" << time(0) << "] " << "Sorting " << remappedUnsortedBam << " by read name using " << SAMTOOLS << std::endl;
    
    cmd.str("");
    cmd << SAMTOOLS << " sort -n -m " << memory << " " << remappedUnsortedBam << " " << workingDirectory << "/remapped";
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
    record1.flag &= ~BAM_FLAG_NEXT_UNMAPPED;
    record2.flag &= ~BAM_FLAG_NEXT_UNMAPPED;

    // Fix the next reversed flag.
    if (hasFlagNextRC(record2)) record1.flag |= BAM_FLAG_NEXT_RC;
    else record1.flag &= ~BAM_FLAG_NEXT_RC;
    if (hasFlagNextRC(record1)) record2.flag |= BAM_FLAG_NEXT_RC;
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
             StringSet<CharString> & nameStor,
             NameStoreCache<StringSet<CharString> > & nameStoreCache,
             BamHeader const & header1,
             BamHeader const & header2)
{
    header = header1;
    for (unsigned i = 0; i < length(header2.records); ++i)
    {
        if (header1.records[i].type == 1)
            appendValue(header.records, header2.records[i]);
    }
    BamHeaderRecordTypeLess less;
    std::stable_sort(begin(header.records, Standard()), end(header.records, Standard()), less);
    
    // Write sequenceInfos, name store and cache for this header.
    unsigned idx = 0;
    for (unsigned i = 0; i < length(header.sequenceInfos); ++i)
        if (!getIdByName(nameStor, header.sequenceInfos[i].i1, idx, nameStoreCache))
            appendName(nameStor, header.sequenceInfos[i].i1, nameStoreCache);
    for (unsigned i = 0; i < length(header2.sequenceInfos); ++i)
    {
        if (!getIdByName(nameStor, header2.sequenceInfos[i].i1, idx, nameStoreCache))
        {
            appendName(nameStor, header2.sequenceInfos[i].i1, nameStoreCache);
            appendValue(header.sequenceInfos, header2.sequenceInfos[i]);
        }
    }
}

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
    BamStream remappedStream(toCString(remappedBam));
    if (!isGood(remappedStream))
    {
        std::cerr << "ERROR: Could not open input bam/sam file " << remappedBam << "." << std::endl;
        return 1;
    }

    // Prepare a header for the output files.
    TNameStore nameStor;
    NameStoreCache<TNameStore> nameStoreCache(nameStor);
    BamHeader header;
    mergeHeaders(header, nameStor, nameStoreCache, nonRefStream.header, remappedStream.header);

    // Open output file and set the header.
    BamStream outStream(toCString(mergedBam), BamStream::WRITE);
    if (!isGood(remappedStream))
    {
        std::cerr << "ERROR: Could not open output bam file " << mergedBam << "." << std::endl;
        return 1;
    }
    outStream.header = header;

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
sickle_filtering(CharString & firstFiltered,
                 CharString & secondFiltered,
                 CharString & singleFiltered,
                 CharString & fastqFirst,
                 CharString & fastqSecond,
                 CharString & fastqSingle,
                 CharString & workingDirectory)
{
    std::stringstream cmd;

    std::cerr << "[" << time(0) << "] " << "Filtering fastq files using " << SICKLE << std::endl;
    cmd.str("");
    cmd << SICKLE << " pe -q 20 -l 60 -x -n -t sanger -f" << fastqFirst << " -r " << fastqSecond;
    cmd << " -o " << firstFiltered << " -p " << secondFiltered << " -s " << singleFiltered;
    if (system(cmd.str().c_str()) != 0) // runs sickle on paired end reads
    {
        std::cerr << "ERROR while filtering " << fastqFirst << " and " << fastqSecond;
        std::cerr << " using " << SICKLE << std::endl;
        return 1;
    }
    
    CharString singleFiltered2 = getFileName(workingDirectory, "filtered.single2.fastq");
    
    cmd.str("");
    cmd << SICKLE << " se -q 20 -l 60 -x -n -t sanger -f " << fastqSingle << " -o " << singleFiltered2;
    if (system(cmd.str().c_str()) != 0) // runs sickle on single end reads
    {
        std::cerr << "ERROR while filtering " << fastqSingle << " using " << SICKLE << std::endl;
        return 1;
    }
    
    cmd.str("");
    cmd << "cat " << singleFiltered2 << " >> " << singleFiltered;
    if (system(cmd.str().c_str()) != 0) // merges the fastq files of single end reads.
    {
        std::cerr << "ERROR while concatenating " << singleFiltered2 << " to " << singleFiltered << std::endl;
        return 1;
    }
    
    remove(toCString(singleFiltered2));

    return 0;
}

// ==========================================================================
// Function velvet_assembly()
// ==========================================================================

inline bool
velvet_assembly(CharString & firstFiltered, CharString & secondFiltered, CharString & singleFiltered, CharString & workingDirectory, unsigned kmerLength)
{
    std::stringstream cmd;
    
    CharString assemblyDirectory = getFileName(workingDirectory, "assembly");
    
    std::cerr << "[" << time(0) << "] " << "Preparing assembly of unmapped reads from filtered fastq files using " << VELVETH << std::endl;
    cmd.str("");
    cmd << VELVETH << " " << assemblyDirectory << " " << kmerLength << " -short -fastq " << singleFiltered;
    cmd << " -shortPaired -fastq -separate " << firstFiltered << " " << secondFiltered;
    if (system(cmd.str().c_str()) != 0) // prepares velvet assembly, use k=47 for longer contigs
    {
        std::cerr << "ERROR while preparing assembly with " << VELVETH << " of ";
        std::cerr << singleFiltered << ", " << firstFiltered << ", and " << secondFiltered << std::endl;
        return 1;
    }
    
    std::cerr << "[" << time(0) << "] " << "Assembling unmapped reads from filtered fastq files using " << VELVETG << std::endl;
    cmd.str("");
    cmd << VELVETG << " " << assemblyDirectory << " -exp_cov auto -cov_cutoff 2 -max_coverage 100 -scaffolding no";
    if (system(cmd.str().c_str()) != 0) // runs the velvet graph part
    {
        std::cerr << "ERROR while assembling " << assemblyDirectory << " with " << VELVETG << std::endl;
        return 1;
    }
    
    CharString contigFileAssembly = getFileName(workingDirectory, "assembly/contigs.fa");
    CharString contigFile = getFileName(workingDirectory, "contigs.fa");

    rename(toCString(contigFileAssembly), toCString(contigFile));
    removeAssemblyDirectory(assemblyDirectory);
    
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
    
    CharString fastqFirst = getFileName(options.workingDirectory, "paired.1.fastq");
    CharString fastqSecond = getFileName(options.workingDirectory, "paired.2.fastq");
    CharString fastqSingle = getFileName(options.workingDirectory, "single.fastq");
    CharString matesBam = getFileName(options.workingDirectory, "mates.bam");
    CharString nonRefBam = getFileName(options.workingDirectory, "non_ref.bam");

    // Crop unmapped reads and reads with unreliable mappings from the input bam file.
    std::cerr << "[" << time(0) << "] Cropping unmapped reads from " << options.mappingFile << std::endl;
    if (options.adapters == "HiSeqX")
    {
        if (crop_unmapped(fastqFirst, fastqSecond, fastqSingle, matesBam, options.mappingFile, options.humanSeqs, HiSeqXAdapters()) != 0)
            return 1;
    }
    else if (options.adapters == "HiSeq")
    {
        if (crop_unmapped(fastqFirst, fastqSecond, fastqSingle, matesBam, options.mappingFile, options.humanSeqs, HiSeqAdapters()) != 0)
            return 1;
    }
    else
    {
        if (crop_unmapped(fastqFirst, fastqSecond, fastqSingle, matesBam, options.mappingFile, options.humanSeqs, NoAdapters()) != 0)
            return 1;
    }

    // Sort <WD>/mates.bam by read name.
    std::cerr << "[" << time(0) << "] " << "Sorting " << matesBam << " using " << SAMTOOLS << std::endl;
    std::stringstream cmd;
    cmd << SAMTOOLS << " sort -n -m " << options.memory << " " << matesBam << " " << options.workingDirectory << "/non_ref";
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << matesBam << std::endl;
        return 1;
    }

    remove(toCString(matesBam));
    
    // Remapping of unmapped with bwa if a fasta reference is given.
    if (options.referenceFile != "")
    {
        // Align with bwa, update fastq files of unaligned reads, and sort remaining bam records by read name.
        CharString remappedBam = getFileName(options.workingDirectory, "remapped.bam");
        if (remapping(fastqFirst, fastqSecond, fastqSingle,
                      options.referenceFile, options.workingDirectory,
                      options.humanSeqs, options.threads, options.memory) != 0)
            return 1;
        
        // Set the mate's location and merge non_ref.bam and remapped.bam into a single file.
        CharString mergedBam = getFileName(options.workingDirectory, "merged.bam");
        merge_and_set_mate(mergedBam, nonRefBam, remappedBam);
        remove(toCString(remappedBam));
        remove(toCString(nonRefBam));
        rename(toCString(mergedBam), toCString(nonRefBam));
    }
    
    CharString firstFiltered = getFileName(options.workingDirectory, "filtered.paired.1.fastq");
    CharString secondFiltered = getFileName(options.workingDirectory, "filtered.paired.2.fastq");
    CharString singleFiltered = getFileName(options.workingDirectory, "filtered.single.fastq");
    
    // Quality filtering/trimming with sickle.
    if (sickle_filtering(firstFiltered, secondFiltered, singleFiltered,
                     fastqFirst, fastqSecond, fastqSingle,
                     options.workingDirectory) != 0) return 1;

    // Assembly with velvet.
    if (velvet_assembly(firstFiltered, secondFiltered, singleFiltered, options.workingDirectory, options.kmerLength) != 0) return 1;
    
    remove(toCString(firstFiltered));
    remove(toCString(secondFiltered));
    remove(toCString(singleFiltered));
    
    return ret;
}

#endif // #ifndef POPINS_ASSEMBLE_H_

