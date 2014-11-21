//#include <seqan/bam_io.h>
#include <seqan/align.h>
#include <seqan/index.h>

#ifndef ADAPTER_REMOVAL_H_
#define ADAPTER_REMOVAL_H_


struct NoAdapters_;
typedef Tag<NoAdapters_> NoAdapters;

struct HiSeqAdapters_;
typedef Tag<HiSeqAdapters_> HiSeqAdapters;

struct HiSeqXAdapters_;
typedef Tag<HiSeqXAdapters_> HiSeqXAdapters;


inline StringSet<Dna5String>
complementUniversalOneError()
{
    typedef Dna5String TSeq;
    typedef Size<TSeq>::Type TSize;
    typedef Dna TAlphabet;

    // Append the error-free sequence to sequence set.
    StringSet<TSeq> adaptSeqs;
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC");

    complement(adaptSeqs[0]);

    // Append all sequences with one mismatch to sequence set.
    TSeq errSeq = adaptSeqs[0];
    for (TSize i = 0; i < length(adaptSeqs[0]); ++i)
    {
        for (TSize c = 0; c < ValueSize<TAlphabet>::VALUE; ++c)
        {
            errSeq[i] = c;
            appendValue(adaptSeqs, errSeq);
        }
        errSeq[i] = adaptSeqs[0][i];
    }
    
    return adaptSeqs;
}


inline StringSet<Dna5String>
truSeqs(NoAdapters &)
{
    // Nothing to be done.
    return StringSet<Dna5String>();
}
    

inline StringSet<Dna5String>
truSeqs(HiSeqAdapters)
{
    StringSet<Dna5String> adaptSeqs;
    
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG");
    
    return adaptSeqs;
}

inline StringSet<Dna5String>
truSeqs(HiSeqXAdapters)
{
    StringSet<Dna5String> adaptSeqs;
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACATAGAGGCACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACCCTATCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACGGCTCTGAACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACAGGCGAAGACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACTAATCTTAACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACCAGGACGTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACACGTACTGACACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTCCGGAGAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGCTCATTATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGATTCCATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCAGAAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAATTCGTATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTGAAGCTATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAATGCGCATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGGCTATGATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTCCGCGAAATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTCTCGCGCATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGCGATAGATCTCGTATGCCGTCTTCTGCTTG");
    
    return adaptSeqs;
}

template<typename TAdapterTag>
inline StringSet<Dna5String>
reverseTruSeqsOneError(TAdapterTag tag)
{
    typedef Dna5String TSeq;
    typedef Size<TSeq>::Type TSize;
    typedef Dna TAlphabet;

    StringSet<TSeq> adaptSeqs = truSeqs(tag);
    
    // Append all sequences with one mismatch to sequence set.
    TSize len = length(adaptSeqs);
    for (TSize s = 0; s < len; ++s)
    {
        reverse(adaptSeqs[s]);
        TSeq errSeq = adaptSeqs[s];
        for (TSize i = 0; i < length(adaptSeqs[s]); ++i)
        {
            for (TSize c = 0; c < ValueSize<TAlphabet>::VALUE; ++c)
            {
                errSeq[i] = c;
                appendValue(adaptSeqs, errSeq);
            }
            errSeq[i] = adaptSeqs[s][i];
        }
    }
    
    return adaptSeqs;
}

// Returns the length of the matching prefix.
template<typename TIndexSeq, typename TSpec, typename TSequence>
inline typename Size<TSequence>::Type
prefixMatchLength(Index<TIndexSeq, TSpec> & suffixIndex, TSequence const & seq)
{
    typedef typename Size<TSequence>::Type TSize;
    typedef typename Position<TSequence>::Type TPos;
    typedef Index<TIndexSeq, TSpec> TIndex;


    typename Iterator<TIndex, TopDown<> >::Type it(suffixIndex);

    TSize len = 0;   // length of matching prefix

    while (goDown(it, seq[repLength(it)]))
    {
        // if end of seq reached
        if (repLength(it) >= length(seq))
        {
            // if all of seq matches the beginning of some suffix
            if (infix(representative(it), parentRepLength(it) + 1, length(seq)) ==
                infix(seq, parentRepLength(it) + 1, length(seq)))
            {
                len = length(seq);
            }
            break;
        }

        // if mismatch in remaining characters along this branch
        if (infix(representative(it), parentRepLength(it) + 1, repLength(it)) !=
            infix(seq, parentRepLength(it) + 1, repLength(it)))
        {
            break;
        }

        // if current prefix is a complete suffix
        if (isRightTerminal(it))
        {
            len = repLength(it);
        }
    }

    return len;
}

template<typename TSize>
String<CigarElement<> >
cigarPrefix(String<CigarElement<> > const & cigar, TSize len)
{
    typedef String<CigarElement<> > TCigar;
    typedef typename Iterator<TCigar const>::Type TIter;
    
    TCigar prefixCigar;
    
    TSize pos = 0;
    
    TIter itEnd = end(cigar);
    for (TIter it = begin(cigar); it != itEnd; ++it)
    {
        pos += (*it).count;
        if (pos < len)
        {
            appendValue(prefixCigar, *it);
        }
        else if (pos == len)
        {
            appendValue(prefixCigar, *it);
            break;
        }
        else
        {
            appendValue(prefixCigar, CigarElement<>((*it).operation, (*it).count + len - pos));
            break;
        }
    }
    
    return prefixCigar;
}

template<typename TSize>
String<CigarElement<> >
cigarSuffix(String<CigarElement<> > const & cigar, TSize len)
{
    typedef String<CigarElement<> > TCigar;
    typedef typename Iterator<TCigar const>::Type TIter;
    
    TCigar suffixCigar;
    
    TSize pos = 0;
    
    TIter it = begin(cigar);
    TIter itEnd = end(cigar);
    for (; it != itEnd; ++it)
    {
        pos += (*it).count;
        if (pos >= len) break;
    }
    
    if (it != itEnd && pos > len)
        appendValue(suffixCigar, CigarElement<>((*it).operation, pos - len));
        
    for (; it != itEnd; ++it)
        appendValue(suffixCigar, *it);
    
    return suffixCigar;
}

template<typename TSequence, typename TTag>
int
removeAdapter(BamAlignmentRecord & record,
              Index<StringSet<TSequence> > & indexUniversal,
              Index<StringSet<TSequence> > & indexTruSeqs,
              unsigned minAdapterLength,
              TTag)
{
    typedef typename Size<TSequence>::Type TSize;
    typedef ModifiedString<const TSequence, ModReverse> TRevSequence;
    typedef ModifiedString<const TSequence, ModComplementDna5> TComplSequence;
    typedef ModifiedString<const TComplSequence, ModReverse> TRevComplSequence;
    
    TSize seqLen = length(record.seq);
    
    if (hasFlagRC(record))
    {
        // Search prefex of complemented read in *TruSeq* index.
        TSequence complSeq = record.seq; reverseComplement(complSeq); reverse(complSeq);
        //TComplSequence complSeq(record.seq); // TODO: This modifier seems to be buggy!
        TSize adaptLen = prefixMatchLength(indexTruSeqs, complSeq);
    
        if (adaptLen == seqLen)
        {
            // Read starts with adapter.
            //std::cerr << "Removing full read        " << record.seq << std::endl;
            return 2;
        }
        else if (adaptLen >= minAdapterLength)
        {
            //std::cerr << "Removing first " << adaptLen << " bases from " << record.seq << std::endl;
            replace(record.seq, 0, adaptLen, "");
            replace(record.qual, 0, adaptLen, "");
            record.cigar = cigarSuffix(record.cigar, adaptLen);
            return 1;
        }
        
        // Search prefix of complemented read in *Universal* index.
        adaptLen = prefixMatchLength(indexUniversal, complSeq);
        
        if (adaptLen == seqLen)
        {
            // Read starts with adapter.
            //std::cerr << "Removing full read        " << record.seq << std::endl;
            return 2;
        }
        else if (adaptLen >= minAdapterLength)
        {
            //std::cerr << "Removing first " << adaptLen << " bases from " << record.seq << std::endl;
            replace(record.seq, 0, adaptLen, "");
            replace(record.qual, 0, adaptLen, "");
            record.cigar = cigarSuffix(record.cigar, adaptLen);
            return 1;
        }
    }
    else
    {
        // Search prefix of reversed read in *TruSeq* index.
        TRevSequence revSeq(record.seq);
        TSize adaptLen = prefixMatchLength(indexTruSeqs, revSeq);
    
        if (adaptLen == seqLen)
        {
            // Read starts with adapter.
            //std::cerr << "Removing full read        " << record.seq << std::endl;
            return 2;
        }
        else if (adaptLen >= minAdapterLength)
        {
            //std::cerr << "Removing last " << adaptLen << " bases from  " << record.seq << std::endl;
            replace(record.seq, seqLen-adaptLen, seqLen, "");
            replace(record.qual, seqLen-adaptLen, seqLen, "");
            record.cigar = cigarPrefix(record.cigar, adaptLen);
            return 1;
        }
    
        // Search prefix of read in *Universal* index.
        adaptLen = prefixMatchLength(indexUniversal, revSeq);
            
        if (adaptLen == seqLen)
        {
            // Read starts with adapter.
            //std::cerr << "Removing full read        " << record.seq << std::endl;
            return 2;
        }
        else if (adaptLen >= minAdapterLength)
        {
            //std::cerr << "Removing last " << adaptLen << " bases from  " << record.seq << std::endl;
            replace(record.seq, seqLen-adaptLen, seqLen, "");
            replace(record.qual, seqLen-adaptLen, seqLen, "");
            record.cigar = cigarPrefix(record.cigar, adaptLen);
            return 1;
        }
    }

    return 0;
}

template<typename TSequence>
int
removeAdapter(BamAlignmentRecord &,
              Index<StringSet<TSequence> > &,
              Index<StringSet<TSequence> > &,
              unsigned,
              NoAdapters)
{
    return 0;
}

#endif // #ifndef ADAPTER_REMOVAL_H_
