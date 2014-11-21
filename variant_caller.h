#ifndef VARIANT_CALLER_H
#define VARIANT_CALLER_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <map>
#include <string>
#include <vector>

#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/vcf_io.h>

using namespace seqan;

// Sequence, alignment, and alignment row.
typedef String< Dna5> TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;

// Name store, cache, and BAM I/O context.
typedef StringSet< CharString> TNameStore;
typedef NameStoreCache< TNameStore> TNameStoreCache;
typedef BamIOContext< TNameStore> TBamIOContext;

#define COMPONENT_HAS_NO_END -1

enum component_dir{
  both_dir_forward, both_dir_reverse, left_dir_forward, left_dir_reverse, right_dir_forward, right_dir_reverse
};

class bamF
{
public:
  TBamIOContext hN;
  BamIndex< Bai> baiI; 
  Stream< Bgzf> bamS;
  bamF();
};

inline int
qualToInt( char c )
{
    return (int) c - 33;
}

template <typename T> int
trimReadEnds( T& seq, CharString& qual, int bpQclip, bool verbose)
{
  int beg, end;
  for( beg = 0; (beg < (int)length( seq )) && ( qualToInt( qual[beg] ) < bpQclip ); beg++ ){}
  for( end = (int)length( seq )-1; (end > 0) && ( qualToInt( qual[end] ) < bpQclip  )  ; end-- ){
    if( verbose ) std::cout << seq[end] << " " << qual[end] << " " << (int) qual[end] << " " << bpQclip << std::endl;
  }
  if( verbose ) std::cout << "Trimming " << beg << " " << end << " " << length( seq ) << std::endl; 
  if( beg > end ){
    seq = "";
    qual = "";
  }else{
    seq = infix( seq, beg, end+1 );
    qual = infix( qual, beg, end+1 );
  }
  return 1;
}

void parseComponent( CharString& altR, bool verbose, component_dir& cdir, CharString& componentName, int& beginPos, int& endPos )
{
  /* BUG:  assumes that the length of the variant is 1, this is not the case for general variants 
   */
  if( verbose ) std::cout  << "In parseComponent " << altR << " altR " << length( altR ) << std::endl;
  if( altR[0] == ']' ){
    int cNameEnd = 1;
    for( unsigned i = 1; i < length( altR ); i++ ){
      if( altR[i] == ':' ){
	cNameEnd = i;
      }
    }
    componentName = infix( altR, 1, (cNameEnd-1) );
    // Need to make sure the indexing is correct
    if( altR[cNameEnd-1] == 'f' or altR[cNameEnd-1] == 'F' ){
      cdir = right_dir_forward;
    }else if( altR[cNameEnd-1] == 'r' or altR[cNameEnd-1] == 'R' ){
      cdir = right_dir_reverse;
    }else{
      if( verbose ) std::cout  << "Incorrect component name format " << altR[cNameEnd-1] << " " << cNameEnd << std::endl;
    }
    beginPos = COMPONENT_HAS_NO_END;
    CharString endPosChar = infix( altR, cNameEnd+1, length( altR )-2 );
    endPos = atoi( toCString( endPosChar ) );
  }else if( altR[1] == '[' ){
    int cNameEnd = 2;
    for( unsigned i = 2; i < length( altR ); i++ ){
      if( altR[i] == ':' ){
	cNameEnd = i;
      }
    }
    componentName = infix( altR, 2, cNameEnd-1 );
    if( altR[cNameEnd-1] == 'f' or altR[cNameEnd-1] == 'F' ){
      cdir = left_dir_forward;
    }else if( altR[cNameEnd-1] == 'r' or altR[cNameEnd-1] == 'R' ){
      cdir = left_dir_reverse;
    }else{
      if( verbose ) std::cout  << "Incorrect component name format " << altR[cNameEnd-1] << " " << cNameEnd << std::endl;
    }
    endPos = COMPONENT_HAS_NO_END;
    CharString beginPosChar = infix( altR, (cNameEnd+1), length(altR)-1 );
    beginPos = atoi( toCString( beginPosChar ) );
  }else{
    componentName = infix( altR, 1, length( altR )-1 );
    if( altR[length(altR)-1] == 'f' or altR[length(altR)-1] == 'F' )
      cdir = both_dir_forward;
    else if( altR[length(altR)-1] == 'r' or altR[length(altR)-1] == 'R' )
      cdir = both_dir_reverse;
    else
      if( verbose ) std::cout  << "Incorrect component name format " << std::endl;
    beginPos = COMPONENT_HAS_NO_END;
    endPos = COMPONENT_HAS_NO_END;
  }
  if( verbose ) std::cout << "parseComponent " << altR << " " << cdir << " " << componentName << " " << beginPos << " " << endPos << std::endl;
}

/** Input: A sequence and a VCF entry
    Two alternate sequences, a breakpoint describing where one sequence inserted the other sequence
    and the set of reads belonging to those sequences
    Output: The likelihood of the allelotype of the two sequences
 */
// Note!! We probably want to align to the reverse complement as well

template<typename TOptions>
int alignReadToSeq( TSequence& to, BamAlignmentRecord& bar, TOptions & options)
{
  // this is very similar to 
  //  bamRecordToAlignment( align, to, bar );  
  // which reads the cigar info from the string 
  // std::cout << "In alignReadToSeq " << to << " " << bar.seq << std::endl;
  TAlign align;
  //  std::cout << align;  
  // There is a bug in SeqAn that does not allow the two commands above.
  if (options.verbose) std::cout << "Pre trimming " << bar.seq << std::endl; 
  trimReadEnds( bar.seq, bar.qual, options.bpQclip, options.verbose);
  if (options.verbose) std::cout << "After trimming " << bar.seq << std::endl; 
  if (length(bar.seq) < (unsigned)options.minSeqLen){
    return -100;
  }
  resize( rows( align ), 2 );
  assignSource(row(align,0), to);
  assignSource(row(align,1), bar.seq );
  //score = match,mismatch, gap extend, gap open
  //need to fix these parameters and set them as options into our program
  /* We seem to need to do the global alignment twice */
  int score, score2;
  Score<int, Simple> scoringScheme(options.match, options.mismatch, options.gapExtend, options.gapOpen);
  if (options.fullOverlap)
  {
    score = globalAlignment(align, scoringScheme, AlignConfig< true, true, true, true>());
    if (options.verbose) std::cout << "Score1 " << score << " " << to << " " << bar.seq << std::endl;
    reverseComplement( bar.seq );
    assignSource(row(align,1), bar.seq );
    score2 = globalAlignment(align, scoringScheme, AlignConfig< true, true, true, true>());
    if (options.verbose) std::cout << "Score2 " << score2 << " " << to << " " << bar.seq << std::endl;
    reverseComplement( bar.seq );
    if (score2 > score) score = score2;
  }
  else
  {
    score = globalAlignment(align, scoringScheme, AlignConfig< true, false, true, false>());
    if (options.verbose) std::cout << "Score1 " << score  << " " << to << " " << bar.seq << std::endl;
    score2 = globalAlignment(align, scoringScheme, AlignConfig< false, true, false, true>());
    if (score2 > score) score = score2;
    if (options.verbose) std::cout << "Score2 " << score2 << " " << to << " " << bar.seq << std::endl;
    reverseComplement( bar.seq );
    assignSource(row(align,1), bar.seq );
    score2 = globalAlignment(align, scoringScheme, AlignConfig< true, false, true, false>());
    if (score2 > score) score = score2;
    if (options.verbose) std::cout << "Score3 " << score2 << " " << to << " " << bar.seq << std::endl;
    score2 = globalAlignment(align, scoringScheme, AlignConfig< false, true, false, true>());
    if (options.verbose) std::cout << "Score4 " << score2 << " " << to << " " << bar.seq << std::endl;
    if (score2 > score) score = score2;
    reverseComplement(bar.seq);
  }
  if (options.verbose) std::cout << "Out alignReadToSeq " << to << " " << bar.seq << " " << score <<  std::endl;
  return score;
}

void transformLogLtoP( std::vector< double>& L )
{
  assert( L.size() == 3 );
  double max=L[0];
  if( L[1] > max ) max = L[1];
  if( L[2] > max ) max = L[2];
  L[0] -= max;
  L[1] -= max;
  L[2] -= max;
  double div = exp( L[0] ) + exp( L[1] ) + exp( L[2] );
  L[0] = exp( L[0] )/div;
  L[1] = exp( L[1] )/div;
  L[2] = exp( L[2] )/div;
} 

// Should create a class containing the hN, baiI and bamS
// and another class that contains the reads in the given region
int readBamRegion(TBamIOContext& hN, BamIndex< Bai>& baiI, Stream< Bgzf>& bamS,
		          CharString& chrom, int beg, int end, bool addReadGroup, bool verbose,
		          std::map< CharString, BamAlignmentRecord>& bars1, std::map< CharString, BamAlignmentRecord>& bars2 )
{
  //Possible shift by 1 in position
  if( verbose ) std::cout << "reading Bam region " << chrom << " " << beg << " " << end << std::endl;

  int rID = 0;
  if (!getIdByName( nameStore( hN ), chrom, rID, nameStoreCache( hN ) ))
  {
    if( verbose ) std::cout << "ERROR: Reference sequence named " << chrom << " not known.\n";
    return 1;
  }

  // Jump the BGZF stream to this position.
  bool hasAlignments = false;
  if (!jumpToRegion( bamS, hasAlignments, hN, rID, beg, end, baiI))
  {
    std::cerr << "ERROR: Could not jump to " << rID << " " << chrom << ":" << beg << "-" << end << "\n";
    return 1;
  }
  if( verbose ) std::cout << "hasAlignments " << hasAlignments << std::endl;
  if (!hasAlignments)
    return 0;  // No alignments here.
 
  BamAlignmentRecord record;
  while (!atEnd(bamS))
  {
    if (readRecord(record, hN, bamS, Bam()) != 0)
    {
      std::cerr << "ERROR: Could not read record from BAM file.\n";
      return 1;
    }
    //  if( verbose ) std::cout << "Reading " << record.qName << std::endl;
    // If we are on the next reference or at the end already then we stop.
    if (record.rID == -1 || record.rID > rID || record.beginPos >= end )
      break;
    // If we are left of the selected position then we skip this record.
    if (record.beginPos + getAlignmentLengthInRef(record)  < (unsigned)beg) // We would like to read the read even if the end pos is less than the begin of our region
      continue;

    if( (not hasFlagDuplicate( record )) and (not hasFlagQCNoPass( record )) ){
      if( addReadGroup ){
	BamTagsDict tagsDict(record.tags);
	unsigned idx;
	if (findTagKey(idx, tagsDict, "RG")){
	  CharString readGroup = getTagValue(tagsDict, idx);
  
	  readGroup = infix(readGroup, 1, length(readGroup)-1);
	  readGroup += ":";
	  readGroup += record.qName;
	  record.qName = readGroup;
	}
      }
      if( hasFlagFirst( record )){
	bars1[record.qName] = record;
      }else{
	bars2[record.qName] = record;
      }
    }

    /*    if (write2( std::cout, record, hN, Sam()) != 0)
      {
	cerr << "ERROR: Could not write record to stdout.\n";
	return 1;
	} */
  }

  //  std::cerr << "finished reading bam" << std::endl;
  return 0;
}

template<typename TOptions>
int addBARsToVC(std::map< CharString, BamAlignmentRecord>&  bars,  TSequence& refSeq, TSequence& altSeq, TOptions & options, std::vector< double>& vC)
{
  for( auto i = bars.begin(); i != bars.end(); i++ ){
    double asRef = alignReadToSeq( refSeq, (*i).second, options );
    double asAlt = alignReadToSeq( altSeq, (*i).second, options );
    if( asRef <= options.minAlignScore ) asRef = options.minAlignScore;
    if( asAlt <= options.minAlignScore ) asAlt = options.minAlignScore;
    if( asRef != options.minAlignScore || asAlt != options.minAlignScore ){
      double pa,pr;
      if( asRef > asAlt + log((1-options.minReadProb)/options.minReadProb) ){
	pr = log( 1-options.minReadProb );
	pa = log( options.minReadProb );
      }else if( asAlt > asRef + log( (1-options.minReadProb)/options.minReadProb) ){
	pa = log( 1-options.minReadProb );
	pr = log( options.minReadProb );
      }else{
	double div = exp( asRef ) + exp( asAlt );
	pa = log( exp( asAlt )/div );
	pr = log( exp( asRef )/div );
      }
      vC[0] += pr;
      vC[1] += -log( 2.0);
      vC[2] += pa;
    }
  }
  return 0;
}

/** Input: a vector of sequences; the different alleles of a particular variant.  A set of BamAlignmentRecords overlapping the 
    Output: the log likelihood of the reads having been generated by each one of the sequences is updated, vs the alternate that the read
    were generated by one of the other sequences
 */
template<typename TOptions>
int addBARsToSLs(std::map< CharString, BamAlignmentRecord>& bars, std::vector< TSequence>& seqs, TOptions & options,
                 std::vector< std::vector< double > >& seqPairLs)
{
  assert( seqs.size() == seqPairLs.size() );
  for( auto i = bars.begin(); i != bars.end(); i++ ){
    std::vector< double> aS( seqs.size() );
    bool alignmentFound = false;
    double div = 0.0;
    for( unsigned j = 0; j < seqs.size(); j++ ){
      aS[j] = alignReadToSeq( seqs[j], (*i).second, options );
      if( aS[j] > options.minAlignScore )
	alignmentFound = true;
      div += exp( aS[j] );
    }
    if( alignmentFound ){
      for( unsigned j = 0; j < seqs.size(); j++ ){
        aS[j] = exp( aS[j] )/div;
        if( aS[j] > 1-options.minReadProb ) aS[j] = 1-options.minReadProb;
        else if( aS[j] < options.minReadProb ) aS[j] = options.minReadProb;
      }
      for( unsigned j = 0; j < seqs.size(); j++ ){
        for( unsigned k = 0; k < seqs.size(); k++ ){
          seqPairLs[j][k] += log( aS[j] + aS[k] ) - log( 2.0 );
        }
      }
    }
  }  
  return 0;
}

/** 
 *
 */

template<typename TOptions>
int countBARalignments(std::map< CharString, BamAlignmentRecord>& bars, std::vector< TSequence>& seqs, TOptions & options, 
			           std::vector< int>& seqCounts, std::vector< std::vector< int > >& seqPairCounts)
{
  assert( seqs.size() == seqPairCounts.size() );
  for( auto i = bars.begin(); i != bars.end(); i++ ){
    std::vector< double> aS( seqs.size() );
    for( unsigned j = 0; j < seqs.size(); j++ ){
      aS[j] = alignReadToSeq( seqs[j], (*i).second, options );
      if( aS[j] >= options.minAlignScore ) seqCounts[j]++;
    }
    for( unsigned j = 0; j < seqs.size(); j++ ){
      for( unsigned k = 0; k < seqs.size(); k++ ){
        if( aS[j] >= options.minAlignScore || aS[k] >= options.minAlignScore )
          seqPairCounts[j][k]++;
      }
    }
  }  
  return 0;
}

enum barPairType{ refIns, refOverlap 
};

int addBARPairToVC( barPairType barPT, double minReadProb, std::vector< double>& vC )
{
  if( barPT == refIns ){
    vC[0] += log( minReadProb );
    vC[1] += -log( 2.0 );
    vC[2] += log( 1-minReadProb);
  }else if( barPT == refOverlap ){
    vC[0] += log( 1-minReadProb);
    vC[1] += -log( 2.0 );
    vC[2] += log( minReadProb );
  }
  return 0;
}

void
parseInfoField( CharString& infoField, bool verbose, int& refDl, int& refDr ){
  refDl = 0;
  refDr = -1;
  for( unsigned i = 0; i < length( infoField ); i++ ){
    if (infoField[i] == '=' and i >= 4 and infoField[i-1] == 'D' and infoField[i-2] == 'F' and infoField[i-3] == 'E' and infoField[i-4] == 'R'){
      unsigned j = i+1;
      while (infoField[j] != ',') j++;
      CharString leftStr = infix( infoField, i+1, j);
      if (verbose ) std::cout << "parseInfoField left " << leftStr << " " << i << " " << j << std::endl;
      refDl = atoi( toCString( leftStr ) );
      unsigned k = j+1;
      while ( k < length( infoField ) and infoField[k] != ';') k++;
      CharString rightStr = infix( infoField, j+1, k);
      if (verbose) std::cout << "parseInfoField right " << rightStr << " " << j << " " << k <<  std::endl;
      refDr = atoi( toCString( rightStr ) );
    }
  }
}

/** Input:  A fasta index, bamStream and a VCF entry
    Output: vC, with vC[i] = probability of i copies of the alternate given the reads in the region
*/
template<typename TOptions>
int variantCallRegionReadPair(CharString & chrom, CharString & componentName,
                              int beginPos, int devL, int devR, component_dir componentDir,
                              BamIndex< Bai>& baiI, Stream< Bgzf>& bamS, TBamIOContext& hN,
                              BamIndex< Bai>& baiIAlt, Stream< Bgzf>& bamSAlt, TBamIOContext& hNAlt,
                              TOptions & options, std::vector< double>& vC)
{
  std::map< CharString, BamAlignmentRecord> bars1L;
  std::map< CharString, BamAlignmentRecord> bars1R;
  std::map< CharString, BamAlignmentRecord> bars1Ins;
  std::map< CharString, BamAlignmentRecord> bars2L;
  std::map< CharString, BamAlignmentRecord> bars2R;
  std::map< CharString, BamAlignmentRecord> bars2Ins;

  readBamRegion( hN, baiI, bamS, chrom, beginPos-options.maxInsertSize + devL, beginPos + devL, options.addReadGroup, options.verbose, bars1L, bars2L );
  readBamRegion( hN, baiI, bamS, chrom, beginPos+devR, beginPos + devR + options.maxInsertSize, options.addReadGroup, options.verbose, bars1R, bars2R );
  readBamRegion( hNAlt, baiIAlt, bamSAlt, componentName, 0, 1e9, options.addReadGroup, options.verbose, bars1Ins, bars2Ins );
  if( options.verbose ) std::cout << "variantCallRegionReadPair " << bars1L.size() << " " << bars2L.size() << " " << bars1R.size() << " " << bars2R.size() << " " << bars1Ins.size() << " " << bars2Ins.size() << std::endl;
  if( componentDir == left_dir_forward or componentDir == left_dir_reverse ){
    for( auto i = bars1L.begin(); i != bars1L.end(); i++ ){
      if( bars2Ins.count( (*i).first ) != 0 ) addBARPairToVC( refIns, options.minReadProb, vC );
      if( bars2R.count( (*i).first ) != 0 ) addBARPairToVC( refOverlap, options.minReadProb, vC );					      
    }
    for( auto i = bars2L.begin(); i != bars2L.end(); i++ ){
      if( bars1Ins.count( (*i).first ) != 0 ) addBARPairToVC( refIns, options.minReadProb, vC );
      if( bars1R.count( (*i).first ) != 0 ) addBARPairToVC( refOverlap, options.minReadProb, vC );					      
    }
  }else{
    for( auto i = bars1R.begin(); i != bars1R.end(); i++ ){
      if( bars2Ins.count( (*i).first ) != 0 ) addBARPairToVC( refIns, options.minReadProb, vC );
      if( bars2L.count( (*i).first ) != 0 ) addBARPairToVC( refOverlap, options.minReadProb, vC );					      
    }
    for( auto i = bars2R.begin(); i != bars2R.end(); i++ ){
      if( bars1Ins.count( (*i).first ) != 0 ) addBARPairToVC( refIns, options.minReadProb, vC );
      if( bars1L.count( (*i).first ) != 0 ) addBARPairToVC( refOverlap, options.minReadProb, vC );					      
    }
  }
  return 0;
}


/** Input:  A fasta index, bamStream and a VCF entry
    Output: vC, with vC[i] = probability of i copies of the alternate given the 
    reads in the region
*/
template<typename TOptions>
int variantCallRegion(VcfRecord & variant, VcfHeader & vcfH,
                      FaiIndex & faiI, FaiIndex & faiIAlt,
                      BamIndex<Bai> & baiI, Stream<Bgzf> & bamS, TBamIOContext & hN,
                      BamIndex<Bai> & baiIAlt, Stream<Bgzf> & bamSAlt, TBamIOContext & hNAlt,
                      TOptions & options, std::vector< double> & vC)
{
    // Look at all the reads that map to this location, compute the likelihood of the different mappings
    TSequence ref = variant.ref;
    CharString alt = variant.alt;

    //  if( readUntilChar( alt, ',' ) !=  EOF_BEFORE_SUCCESS ){
    //  std::cerr << "cannot handle multiallelic markers ";
    //  return 1;
      // This can probably be handled with TokenizeResult and tokenize.h
      // We would then have to come up with a formula for the multiple alleles
    //  }
    //StringSet< CharString, Owner< ConcatDirect< CharString> > > alts = alt;
    //split( alts, "," );
    __int32 rID = variant.rID;
    __int32 beginPos = variant.beginPos;
    unsigned idx = 0;
    CharString chrom = vcfH.sequenceNames[rID];
    if (!getIdByName(faiI, chrom, idx)){
        if( options.verbose ) std::cout << "rID " << chrom << " " << beginPos << " ERROR: reference FAI index has no entry for rID in Ref mapped.\n";
    }
    TSequence refSeq;
    if (readRegion( refSeq, faiI, idx, beginPos-options.regionWindowSize, beginPos+length(ref)+options.regionWindowSize) != 0){
        if( options.verbose ) std::cout << "ERROR: Could not load reference sequence.\n";
    }    
    if( options.verbose ) std::cout << "refSeq " << refSeq << " " << chrom << " " << beginPos << std::endl;

    CharString componentName;
    component_dir componentDir;
    int beginPosC, endPosC;
    parseComponent( variant.alt, options.verbose, componentDir, componentName, beginPosC, endPosC );

    int devL, devR;
    parseInfoField( variant.info, options.verbose, devL, devR );
    if( devR >= 0 || options.callBoth ){
        if( options.verbose ) std::cout << "variantCallRegionReadPair " << devL << " " << devR << " " << std::endl; 
        variantCallRegionReadPair( chrom, componentName, variant.beginPos, devL, devR, componentDir, baiI, bamS, hN, baiIAlt, bamSAlt, hNAlt, options, vC);
        if( not options.callBoth ){
            transformLogLtoP( vC ); 
            return 0;
        }
    }

    unsigned idxAlt = 0;
    if(!getIdByName( faiIAlt, componentName, idxAlt ) && options.verbose ){
        if( options.verbose ) std::cout << "rID " << componentName << " " << beginPosC << " " << endPosC << std::endl;
        if( options.verbose ) std::cout << "ERROR: reference FAI index has no entry for rID in Alt mapped.\n";
    }

    int seqLenAlt = sequenceLength( faiIAlt, idxAlt );
    int altRegBeg = 0;
    int altRegEnd = seqLenAlt;
    if( componentDir == left_dir_forward ){
        altRegBeg = beginPosC;
        altRegEnd = beginPosC+options.regionWindowSize+length(ref);
        //    readBamRegion( hNAlt, baiIAlt, bamSAlt, componentName, beginPosC, beginPosC+regionWindowSize+length(ref), options.addReadGroup, options.verbose, bars1, bars2 );
    }else if( componentDir == left_dir_reverse ){
        altRegBeg = seqLenAlt-(beginPosC+options.regionWindowSize+length(ref));
        altRegEnd = seqLenAlt-beginPosC;
        //    readBamRegion( hNAlt, baiIAlt, bamSAlt, componentName, seqLenAlt-(beginPosC+regionWindowSize+length(ref)), options.addReadGroup, options.verbose, seqLenAlt-beginPosC, bars1, bars2 );
    }else if( componentDir == right_dir_forward ){
        altRegBeg = endPosC - options.regionWindowSize;
        altRegEnd = endPosC;
        //    readBamRegion( hNAlt, baiIAlt, bamSAlt, componentName, endPosC-regionWindowSize, endPosC, options.addReadGroup, options.verbose, bars1, bars2 );
    }else if( componentDir == right_dir_reverse ){
        altRegBeg = seqLenAlt - endPosC;
        altRegEnd = seqLenAlt - (endPosC-options.regionWindowSize);
        //    readBamRegion( hNAlt, baiIAlt, bamSAlt, componentName, seqLenAlt - endPosC, seqLenAlt - (endPosC-regionWindowSize), options.addReadGroup, options.verbose, bars1, bars2 );
    }else{
        // Nothing to be done.
    }
    int begNs = 0;
    if( altRegBeg < 0 ){
        begNs = -altRegBeg;
        altRegBeg = 0;
    }
    int endNs = 0;
    if( altRegEnd > seqLenAlt ){
        endNs = altRegEnd - seqLenAlt;
        altRegEnd = seqLenAlt;
    }
    // These are no longer optimal, should be replaced by N's
    if( options.verbose ) std::cout << "seqLenAlt " << seqLenAlt << " " << altRegBeg << " " << altRegEnd << " " << componentDir << std::endl;

    TSequence altSeq;
    if( componentDir == left_dir_forward )
    { 
        if (readRegion(altSeq, faiI, idx, beginPos-options.regionWindowSize, beginPos ) != 0 && options.verbose)
        std::cout << "ERROR: Could not load sequence.\n";
        TSequence post;
        if (readRegion(post, faiIAlt, idxAlt, altRegBeg, altRegEnd) != 0 && options.verbose)
            std::cout << "ERROR: Could not load sequence.\n";
        for( int bp = 0; bp < begNs; bp++ ) append( altSeq, 'N' );
        append( altSeq, post );
        if( options.verbose ) std::cout << "AltSeq left_dir_forward " << altSeq << std::endl;
    }
    else if( componentDir == left_dir_reverse )
    { 
        if (readRegion(altSeq, faiI, idx, beginPos-options.regionWindowSize, beginPos ) != 0)
            std::cerr << "ERROR: Could not load sequence.\n";
        TSequence post;
        // We need to read from the end 
        if (readRegion(post, faiIAlt, idxAlt, altRegBeg, altRegEnd ) != 0 && options.verbose)
            std::cout << "ERROR: Could not load sequence.\n";
        reverseComplement( post );
        for( int bp = 0; bp < endNs; bp++ ) append( altSeq, 'N' );
        append( altSeq, post );
        if( options.verbose ) std::cout << "AltSeq left_dir_reverse " << altSeq << std::endl;
    }
    else if( componentDir == right_dir_forward )
    {
        if (readRegion(altSeq, faiIAlt, idxAlt, altRegBeg, altRegEnd ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        for( int bp = 0; bp < endNs; bp++ ) append( altSeq, 'N' );
        TSequence post;
        if (readRegion(post, faiI, idx, beginPos, beginPos+options.regionWindowSize+length(ref) ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        append( altSeq, post );
        if( options.verbose ) std::cout << "AltSeq right_dir_forward " << altSeq << std::endl;
    }
    else if( componentDir == right_dir_reverse )
    {
        //    int seqLenAlt = sequenceLength( faiIAlt, idxAlt ); 
        if (readRegion(altSeq, faiIAlt, idxAlt, altRegBeg, altRegEnd ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        reverseComplement( altSeq );
        for( int bp = 0; bp < begNs; bp++ ) append( altSeq, 'N' );
        TSequence post;
        if (readRegion(post, faiI, idx, beginPos, beginPos+options.regionWindowSize+length(ref) ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        append( altSeq, post );
        if( options.verbose ) std::cout << "AltSeq right_dir_reverse " << altSeq << std::endl;
    }
    else if( componentDir == both_dir_forward )
    {
        if (readRegion(altSeq, faiI, idx, beginPos-options.regionWindowSize, beginPos ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        TSequence post;
        if (readRegion(post, faiIAlt, idxAlt, 0, 1e9 ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        append( altSeq, post );
        if (readRegion(post, faiI, idx, beginPos, beginPos+options.regionWindowSize+length( ref ) ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        append( altSeq, post );
        if( options.verbose ) std::cout << "AltSeq both_dir_forward " << altSeq << std::endl;
    }
    else if( componentDir == both_dir_reverse )
    {
        if (readRegion(altSeq, faiI, idx, beginPos-options.regionWindowSize, beginPos ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        TSequence post;
        if (readRegion(post, faiIAlt, idxAlt, 0, 1e9 ) != 0 && options.verbose)
            std::cout << "ERROR: Could not load sequence.\n";
        reverseComplement( post );
        append( altSeq, post );
        if (readRegion(post, faiI, idx, beginPos, beginPos+options.regionWindowSize+length( ref ) ) != 0 && options.verbose )
            std::cout << "ERROR: Could not load sequence.\n";
        append( altSeq, post );
        if( options.verbose ) std::cout << "AltSeq both_dir_reverse " << altSeq << std::endl;
    }
  
    std::map< CharString, BamAlignmentRecord> bars1;
    std::map< CharString, BamAlignmentRecord> bars2;

    // Need to get chromosome from VCF file
    readBamRegion( hN, baiI, bamS, chrom, beginPos-options.regionWindowSize, beginPos+options.regionWindowSize, options.addReadGroup, options.verbose, bars1, bars2 );

    if( altRegEnd > altRegBeg )
        readBamRegion( hNAlt, baiIAlt, bamSAlt, componentName, altRegBeg, altRegEnd, options.addReadGroup, options.verbose, bars1, bars2 );
    addBARsToVC( bars1, refSeq, altSeq, options, vC );
    addBARsToVC( bars2, refSeq, altSeq, options, vC );

    transformLogLtoP( vC ); 
    return 0;
}


int initializeBam(char* fileName, TBamIOContext & context, BamIndex<Bai> & bamIndex, Stream<Bgzf> & bamStream)
{
    if (!open(bamStream, fileName, "r"))
    {
        std::cerr << "ERROR: Could not open " << fileName << " for reading.\n";
        return 1;
    }
  
    char baifileName[strlen(fileName) + 10];
    strcpy(baifileName, fileName);
    strcat(baifileName, ".bai");
    if (read(bamIndex, baifileName) != 0)
    {
        std::cerr << "ERROR: Could not read BAI index file " << fileName << "\n";
        return 1;
    }

    BamHeader header;
    if (readRecord(header, context, bamStream, Bam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from BAM file " << fileName << "\n";
        return 1;
    }
    return 0;
}

void readBam(char* fileName, std::map<CharString, BamAlignmentRecord> & bars1, std::map<CharString, BamAlignmentRecord> & bars2)
{
    BamStream bamStreamIn(fileName);
    BamAlignmentRecord record;
    while (!atEnd(bamStreamIn))
    {
        readRecord(record, bamStreamIn);
        if((not hasFlagDuplicate(record)) and (not hasFlagQCNoPass(record)))
        {
        	if(hasFlagFirst(record))
	            bars1[record.qName] = record;
        	else
	            bars2[record.qName] = record;
        }
    }
}

#endif // #ifndef VARIANT_CALLER_H
