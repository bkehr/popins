#ifndef SAMPLE_INFO_H_
#define SAMPLE_INFO_H_

#include <seqan/bam_io.h>

using namespace seqan;

// ==========================================================================
// Struct SampleInfo
// ==========================================================================

struct SampleInfo
{
    CharString sample_id;
    CharString bam_file;
    double avg_cov;
    unsigned read_len;
    CharString adapter_type;
    
    SampleInfo() {}
};

// ==========================================================================
// Function initSampleInfo()
// ==========================================================================

SampleInfo
initSampleInfo(CharString & filename, CharString & adapter_type)
{
	SampleInfo info;
	info.bam_file = filename;

    BamFileIn bamFile(toCString(filename));
    BamHeader header;
    readHeader(header, bamFile);
    
    for (unsigned i = 0; i < length(header); ++i)
    {
        if (header[i].type != BamHeaderRecordType::BAM_HEADER_READ_GROUP)
            continue;
        
        for (unsigned j = 0; j < length(header[i].tags); ++j)
        {
            if (header[i].tags[j].i1 != "SM")
                continue;

            info.sample_id = header[i].tags[j].i2;
            break;
        }
    }
    
    BamAlignmentRecord record;
    readRecord(record, bamFile);
    info.read_len = length(record.seq);
    
    info.adapter_type = adapter_type;

    return info;
}

// ==========================================================================
// Function readSampleInfo()
// ==========================================================================

bool
readSampleInfo(SampleInfo & info, CharString & filename)
{
    std::ifstream stream(toCString(filename));
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open sample info file \'" << filename << "\' for reading." << std::endl;
        return 1;
    }
    
    std::string field, value;
    while (stream >> field >> value)
    {
        if (field.compare("SAMPLE_ID") != 0)
            info.sample_id = value;
        else if (field.compare("BAM_FILE") != 0)
            info.bam_file = value;
        else if (field.compare("AVG_COV") != 0)
            info.avg_cov = lexicalCast<double>(value);
        else if (field.compare("READ_LEN") != 0)
            lexicalCast<unsigned>(info.read_len, value);
        else if (field.compare("ADAPTER_TYPE") != 0)
            info.adapter_type = value;
        else
            std::cerr << "WARNING: Ignoring field \'" << field << "\' in sample info file \'" << filename << "\'." << std::endl;
    }
    
    return 0;
}

// ==========================================================================
// Function write()
// ==========================================================================

bool
writeSampleInfo(SampleInfo & info, CharString & filename)
{
    // open the fileDate
    std::ofstream stream(toCString(filename));
    if (!stream.good())
    {
        std::cerr << "ERROR: Could not open sample info file \'" << filename << "\' for writing." << std::endl;
        return 1;
    }
    
    stream << "SAMPLE_ID" << "\t" << info.sample_id << "\n";
    stream << "BAM_FILE" << "\t" << info.bam_file << "\n";
    stream << "AVG_COV" << "\t" << info.avg_cov << "\n";
    stream << "READ_LEN" << "\t" << info.read_len << "\n";
    stream << "ADAPTER_TYPE" << "\t" << info.adapter_type << "\n";
    
    stream.close();
    return 0;
}

#endif  // SAMPLE_INFO_H_
