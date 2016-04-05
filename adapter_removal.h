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


template<typename TTag>
inline Dna5String
getUniversal(TTag)
{
    return "ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
}

inline Dna5String
getUniversal(HiSeqXAdapters)
{
    return "ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGCCTCTATGTGTAGATCTCGGTGGTCGCCGTATCATT";
}
inline Dna5String
getTruSeqPrefix(HiSeqAdapters)
{
    return "ATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
}

inline Dna5String
getTruSeqSuffix(HiSeqAdapters)
{
    return "ATCTCGTATGCCGTCTTCTGCTTG";
}

inline Dna5String
getTruSeqPrefix(HiSeqXAdapters)
{
    return "AATGATACGGCGACCACCGAGATCTACAC";
}

inline Dna5String
getTruSeqSuffix(HiSeqXAdapters)
{
    return "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
}

template<typename TTag>
inline StringSet<Dna5String>
complementUniversalOneError(TTag tag)
{
    typedef Dna5String TSeq;
    typedef Size<TSeq>::Type TSize;
    typedef Dna TAlphabet;

    // Append the error-free sequence to sequence set.
    StringSet<TSeq> adaptSeqs;
    appendValue(adaptSeqs, getUniversal(tag));

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

template<typename TSequence>
bool
startsWithTruSeq(TSequence & seq, HiSeqXAdapters tag)
{
    TSequence truSeqPre = getTruSeqPrefix(tag);
    int score = globalAlignmentScore(truSeqPre, prefix(seq, length(truSeqPre)), Score<int>(1,0,0), -2, 2);
    if (score > (int)length(truSeqPre) - 5)
    {
        TSequence truSeqSuf = getTruSeqSuffix(tag);
        score += globalAlignmentScore(truSeqSuf, infix(seq, length(truSeqPre)+8, length(truSeqPre)+8+length(truSeqSuf)), Score<int>(1,0,0), -2, 2);
        if (score > (int)length(truSeqPre) + (int)length(truSeqSuf) - 10) return 0;
    }
    else
    {
        truSeqPre = getTruSeqPrefix(HiSeqAdapters());
        score = globalAlignmentScore(truSeqPre, prefix(seq, length(truSeqPre)), Score<int>(1,0,0), -2, 2);
        if (score > (int)length(truSeqPre) - 5)
        {
            TSequence truSeqSuf = getTruSeqSuffix(HiSeqAdapters());
            score += globalAlignmentScore(truSeqSuf, infix(seq, length(truSeqPre)+8, length(truSeqPre)+8+length(truSeqSuf)), Score<int>(1,0,0), -2, 2);
            if (score > (int)length(truSeqPre) + (int)length(truSeqSuf) - 10) return 0;
        }
    }
    return 1;
}

template<typename TSequence>
bool
startsWithTruSeq(TSequence & seq, HiSeqAdapters tag)
{
    TSequence truSeqPre = getTruSeqPrefix(tag);
    int score = globalAlignmentScore(truSeqPre, prefix(seq, length(truSeqPre)), Score<int>(1,0,0), -2, 2);
    if (score > (int)length(truSeqPre) - 5)
    {
        TSequence truSeqSuf = getTruSeqSuffix(tag);
        score += globalAlignmentScore(truSeqSuf, infix(seq, length(truSeqPre)+8, length(truSeqPre)+8+length(truSeqSuf)), Score<int>(1,0,0), -2, 2);
        if (score > (int)length(truSeqPre) + (int)length(truSeqSuf) - 10) return 0;
    }
    return 1;
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
    
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATCACG"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGATGT"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TTAGGC"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TGACCA"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ACAGTG"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GCCAAT"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CAGATC"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ACTTGA"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GATCAG"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TAGCTT"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GGCTAC"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CTTGTA"   "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "AGTCAACA" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "AGTTCCGT" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATGTCAGA" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CCGTCCCG" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTCCGCAC" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTGAAACG" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTGGCCTT" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTTTCGGA" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGTACGTA" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GAGTGGAT" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ACTGATAT" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATTCCTTT" "ATCTCGTATGCCGTCTTCTGCTTG");
    
    return adaptSeqs;
}

inline StringSet<Dna5String>
truSeqs(HiSeqXAdapters)
{
    StringSet<Dna5String> adaptSeqs;
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "TATAGCCT" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "ATAGAGGC" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "CCTATCCT" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "GGCTCTGA" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "AGGCGAAG" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "TAATCTTA" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "CAGGACGT" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    appendValue(adaptSeqs, "AATGATACGGCGACCACCGAGATCTACAC" "GTACTGAC" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATTACTCG" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TCCGGAGA" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGCTCATT" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GAGATTCC" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATTCAGAA" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GAATTCGT" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CTGAAGCT" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TAATGCGC" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGGCTATG" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TCCGCGAA" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TCTCGCGC" "ATCTCGTATGCCGTCTTCTGCTTG");
    appendValue(adaptSeqs, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "AGCGATAG" "ATCTCGTATGCCGTCTTCTGCTTG");

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
              TTag tag)
{
    typedef typename Size<TSequence>::Type TSize;
    typedef ModifiedString<const TSequence, ModReverse> TRevSequence;
//  typedef ModifiedString<const TSequence, ModComplementDna5> TComplSequence;
//  typedef ModifiedString<const TComplSequence, ModReverse> TRevComplSequence; // TODO: This modifier seems to be buggy!
    
    TSize seqLen = length(record.seq);
    
    if (hasFlagRC(record))
    {
        TSequence complSeq = record.seq; reverseComplement(complSeq);
        
        // Check for adapter at begin of read.
        if (hasFlagFirst(record))
        {
            // Compute alignment score to TruSeq (excluding barcode)
            if (startsWithTruSeq(complSeq, tag) == 0) return 2;
        }
        else
        {
            // Compute alignment score to reverse complement of Universal
            TSequence universal = getUniversal(tag);
            int score = globalAlignmentScore(universal, prefix(complSeq, length(universal)), Score<int>(1,0,0), -2, 2);
            if (score > (int)length(universal) - 5) return 2;
        }
        
        // Search prefix of complemented read in *TruSeq* index.
        reverse(complSeq);
        //TComplSequence complSeq(record.seq);
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
        // Check for adapter at begin of read.
        if (hasFlagFirst(record))
        {
            // Compute alignemnt score to TruSeq (excluding barcode)
            if (startsWithTruSeq(record.seq, tag) == 0) return 2;
        }
        else
        {
            // Compute alignment score to reverse complement of Universal
            TSequence universal = getUniversal(tag);
            int score = globalAlignmentScore(universal, prefix(record.seq, length(universal)), Score<int>(1,0,0), -2, 2);
            if (score > (int)length(universal) - 5) return 2;
        }

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
