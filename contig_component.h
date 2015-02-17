#ifndef CONTIG_COMPONENT_H_
#define CONTIG_COMPONENT_H_

#include<seqan/sequence.h>

#include "contig_id.h"

using namespace seqan;

// --------------------------------------------------------------------------
// struct ContigComponent
// --------------------------------------------------------------------------

template<typename TSeq>
struct ContigComponent
{
    typedef typename Size<TSeq>::Type TSize;

    StringSet<ContigId, Dependent<> > ids;
    StringSet<TSeq, Dependent<> > contigs;
    
    std::map<TSize, std::set<TSize> > alignedPairs;

    ContigComponent()
    {}
};

// --------------------------------------------------------------------------
// Function writeComponents()
// --------------------------------------------------------------------------

template<typename TSeq>
bool
writeComponents(CharString & fileName, std::map<ContigId, ContigComponent<TSeq> > & components)
{
    typedef typename Size<TSeq>::Type TSize;
    typedef std::map<ContigId, ContigComponent<TSeq> > TComponents;

    // Open the output file.
    std::fstream outputStream(toCString(fileName), std::ios::out);
    if (!outputStream.good())
    {
        std::cerr << "ERROR: Could not open component output file " << fileName << std::endl;
        return 1;
    }
    
    // Write the components.
    TSize i = 0;
    for (typename TComponents::iterator it = components.begin(); it != components.end(); ++it, ++i)
    {
        if (it->second.alignedPairs.begin()->second.size() == 0) continue;
        for (typename std::map<TSize, std::set<TSize> >::iterator pairIt = it->second.alignedPairs.begin(); pairIt != it->second.alignedPairs.end(); ++pairIt)
        {
            if (pairIt->second.size() == 0) continue;

            typename std::set<TSize>::iterator alignedTo = pairIt->second.begin();
            outputStream << " [" << pairIt->first << "->" << *alignedTo;
            if (alignedTo != pairIt->second.end()) ++alignedTo;
            for (; alignedTo != pairIt->second.end(); ++alignedTo)
                outputStream << "," << *alignedTo;

            outputStream << "]";
        }
        outputStream << "\n";
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function readComponents()
// --------------------------------------------------------------------------

template<typename TSeq, typename TSize>
bool
readComponents(String<ContigComponent<TSeq> > & components, UnionFind<int> & uf, CharString & fileName, TSize len, bool verbose)
{
    typedef String<ContigComponent<TSeq> > TComponents;
    
    unsigned numComps = length(components);
    
    // Open the input file and initialize record reader.
    std::fstream stream(toCString(fileName), std::ios::in);
    
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open components input file " << fileName << std::endl;
        return 1;
    }
    
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    CharString buffer;
    
    // Read the components line by line.
    while (!atEnd(reader))
    {
        if (skipChar(reader, ' ') != 0)
        {
            std::cerr << "ERROR: File format error. Reading whitespace from " << fileName << " failed." << std::endl;
            return 1;
        }

        ContigComponent<TSeq> component;
        while (value(reader) == '[')
        {
            TSize key, val, key_rev, val_rev;
            skipChar(reader, '[');
            
            clear(buffer);
            if (readDigits(buffer, reader) != 0)
            {
                std::cerr << "ERROR: File format error. Reading key from " << fileName << " failed." << std::endl;
                return 1;
            }
            lexicalCast2<TSize>(key, buffer);
            if (key < len) key_rev = key + len;
            else key_rev = key - len;
            
            if (skipChar(reader, '-') != 0 || skipChar(reader, '>') != 0)
            {
                std::cerr << "ERROR: File format error. Reading '->' from " << fileName << " failed." << std::endl;
                return 1;
            }
            
            clear(buffer);
            if (readDigits(buffer, reader) != 0)
            {
                std::cerr << "ERROR: File format error. Reading value from " << fileName << " failed." << std::endl;
                return 1;
            }
            lexicalCast2<TSize>(val, buffer);
            if (val < len) val_rev = val + len;
            else val_rev = val - len;
            
            // Add the aligned pairs.
            component.alignedPairs[key].insert(val);

            // Join sets of key and value.
            joinSets(uf, findSet(uf, key), findSet(uf, val));
            SEQAN_ASSERT_EQ(findSet(uf, key), findSet(uf, val));
            joinSets(uf, findSet(uf, key_rev), findSet(uf, val_rev));
            SEQAN_ASSERT_EQ(findSet(uf, key_rev), findSet(uf, val_rev));
            
            while (value(reader) == ',')
            {
                skipChar(reader, ',');
                clear(buffer);
                if (readDigits(buffer, reader) != 0)
                {
                    std::cerr << "ERROR: File format error. Reading value from " << fileName << " failed." << std::endl;
                    return 1;
                }
                lexicalCast2<TSize>(val, buffer);
                if (val < len) val_rev = val + len;
                else val_rev = val - len;
                
                // Add the aligned pairs.
                component.alignedPairs[key].insert(val);

                // Join sets of key and value.
                joinSets(uf, findSet(uf, key), findSet(uf, val));
                SEQAN_ASSERT_EQ(findSet(uf, key), findSet(uf, val));
                joinSets(uf, findSet(uf, key_rev), findSet(uf, val_rev));
                SEQAN_ASSERT_EQ(findSet(uf, key_rev), findSet(uf, val_rev));
            }
            
            if (skipChar(reader, ']') != 0)
            {
                std::cerr << "ERROR: File format error. Reading ']' from " << fileName << " failed." << std::endl;
                return 1;
            }
            
            if (value(reader) == ' ') skipChar(reader, ' ');
        }
        appendValue(components, component);
        skipLine(reader);
    }
    
    if (verbose) std::cerr << "[" << time(0) << "] " << "Loaded " << fileName << ": " << (length(components) - numComps) << " components." << std::endl;
    
    return 0;
}

#endif  // #ifndef CONTIG_COMPONENT_H_
