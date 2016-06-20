#ifndef POPINS_UILS_H_
#define POPINS_UILS_H_

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
initSampleInfo(CharString & filename, CharString sample_id, CharString & adapter_type)
{
    SampleInfo info;
    info.bam_file = filename;
    info.sample_id = sample_id;

    BamFileIn bamFile(toCString(filename));
    BamHeader header;
    readHeader(header, bamFile);
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

// ==========================================================================

void printStatus(const char * message)
{
        // Get the current date and time.
        char timestamp[80];
        time_t now = time(0);
        struct tm tstruct;
        tstruct = *localtime(&now);
        strftime(timestamp, sizeof(timestamp), "[PopIns %Y-%m-%d %X] ", &tstruct);

        // Print time and message.
        std::cerr << timestamp << message << std::endl;
}

void printStatus(std::ostringstream & message)
{
        std::string msg = message.str();
        printStatus(toCString(msg));
}

// ==========================================================================

// Returns true if file exists, otherwise false.
inline bool exists(CharString const & filename)
{
    struct stat buffer;
    return (stat(toCString(filename), &buffer) == 0);
}

// ==========================================================================

// Lists all files <prefix>/*/<filename>.
String<Pair<CharString> >
listFiles(CharString & prefix, CharString & filename)
{
   String<Pair<CharString> > paths;

    DIR *dir = opendir(toCString(prefix));

    struct dirent *entry = readdir(dir);
    while (entry != NULL)
    {
        if (entry->d_type == DT_DIR)
        {
           CharString sampleID = entry->d_name;
           std::stringstream path;
           path << prefix << "/" << sampleID << "/" << filename;
           CharString pathStr = path.str();
           if (exists(pathStr))
                appendValue(paths, Pair<CharString>(sampleID, pathStr));
           else
              std::cerr << "WARNING: \'" << pathStr << "\' does not exist." << std::endl;
        }
        entry = readdir(dir);
    }

    closedir(dir);

    return paths;
}

// ==========================================================================

// List all directories <prefix>/*/ and return only the basename.
String<CharString>
listSubdirectories(CharString & prefix)
{
   String<CharString> subdirs;

   DIR *dir = opendir(toCString(prefix));

   struct dirent *entry = readdir(dir);
   while (entry != NULL)
   {
      if (entry->d_type == DT_DIR)
         appendValue(subdirs, entry->d_name);
      entry = readdir(dir);
   }

   closedir(dir);

   return subdirs;
}

// ==========================================================================

inline CharString
getFileName(CharString const & path, CharString const & name)
{
    CharString filename = path;
    filename += "/";
    filename += name;
    return filename;
}

// ==========================================================================

inline void
removeFile(CharString const & path, const char * filename)
{
    CharString file = path;
    file += "/";
    file += filename;
    remove(toCString(file));
}

// ==========================================================================
// Function checkFileEnding()
// ==========================================================================

bool
checkFileEnding(CharString & filename, std::string ending)
{
    std::string name = toCString(filename);
    size_t dotPos = name.find_last_of('.');

    if (dotPos == std::string::npos)
        return false;

    return name.substr(dotPos + 1, 3) == ending;
}

// ==========================================================================
// Function readFileNames()
// ==========================================================================

bool
readFileNames(String<CharString> & files, CharString & filenameFile)
{
    if (filenameFile == "") return 0;

    std::fstream stream(toCString(filenameFile), std::fstream::in);
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open file listing files " << filenameFile << std::endl;
        return 1;
    }

    std::string file;
    while (stream >> file)
        appendValue(files, CharString(file));

    return 0;
}

// ==========================================================================

template <typename TValue>
bool readFileNames(String<CharString> & files, String<TValue> & values)
{
    if (length(files) > 1) return 0;
    std::cerr << "ReadFileNames " << length(files) << " " << files << std::endl;
    // Open input file
    CharString filenameFile = files[0];
    std::cerr << filenameFile << " " << files << std::endl;
    std::fstream stream(toCString(filenameFile), std::fstream::in);
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open file listing files " << filenameFile << std::endl;
        return 1;
    }

    clear(files);
    clear(values);

    std::string file;
    TValue val;
    while (stream >> file >> val)
    {
        appendValue(files, file);
        appendValue(values, val);
    }

    return 0;
}

// ==========================================================================

bool
parseInterval(Triple<CharString, unsigned, unsigned> & out, CharString & in)
{
    Iterator<CharString, Rooted>::Type it = begin(in, Rooted());

    unsigned colonPos = 0;
    while (it != end(in))
    {
        if (*it == ':')
        {
            colonPos = position(it);
            break;
        }
        ++it;
    }

    if (colonPos == 0)
    {
        out.i1 = in;
        out.i2 = 0;
        out.i3 = maxValue<unsigned>();
        return 0;
    }

    unsigned dashPos = 0;
    while (it != end(in))
    {
        if (*it == '-')
        {
            dashPos = position(it);
            break;
        }
        ++it;
    }

    if (dashPos == 0)
    {
        std::cerr << "ERROR: Interval is not in format CHR:BEG-END." << std::endl;
        return 1;
    }

    out.i1 = prefix(in, colonPos);
    out.i2 = lexicalCast<unsigned>(infix(in, colonPos + 1, dashPos));
    out.i3 = lexicalCast<unsigned>(suffix(in, dashPos + 1));

    return 0;
}

#endif  // POPINS_UILS_H_
