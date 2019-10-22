#include <iostream>
#include <seqan/file.h>
#include <seqan/bam_io.h>
#include <set>
#include <map>
#include <string>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <sys/stat.h>
#include <ctime>
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace seqan;

namespace std
{
    template <>
    struct hash<seqan::String<Dna5> >
    {
      std::size_t
      operator()(seqan::String<Dna5> const & s) const
      {
          std::size_t seed = 42;

          for (char c : s)
              seed ^= c + 0x9e3779b9 + (c << 6) + (c >> 2);

          return seed;
      }
    };
}

//Structure to store marker information
struct STRinfo {
    CharString chrom;
    int STRstart;
    int STRend;
    Dna5String motif;
    float refRepeatNum; //Number of repeats in reference sequence
    Dna5String refBf;
    Dna5String refAf;
    float refRepPurity;
    unsigned minFlankLeft;
    unsigned minFlankRight;
    std::unordered_set<Dna5String> hash8beforeSet;
    std::unordered_set<Dna5String> hash8afterSet;
    vector<float> motifFrqs;
} ;

//Structure to store marker information
struct STRinfoSmall {
    CharString chrom;
    int STRstart;
    int STRend;
    CharString motif;
    float refRepeatNum; //Number of repeats in reference sequence
} ;

//So I can map from STRinfoSmall in finalMap
bool operator<(const STRinfoSmall & Left, const STRinfoSmall & Right)
{
    return Left.STRstart < Right.STRstart;
}

bool operator ==(const STRinfoSmall & Left, const STRinfoSmall & Right)
{
    return Left.chrom == Right.chrom && Left.STRstart == Right.STRstart && Left.STRend == Right.STRend && Left.motif == Right.motif && Left.refRepeatNum == Right.refRepeatNum;
}

//structure to store read information
struct ReadInfo {
    int STRend;
    CharString motif;
    float refRepeatNum; //Number of repeats in reference sequence
    float numOfRepeats;
    float ratioBf;
    float ratioAf;
    unsigned locationShift;
    float purity;
    float ratioOver20In;
    unsigned mateEditDist;
    CharString repSeq; //Repeat sequence in read
} ;

//structure to store read information
struct ReadPairInfo {
    float numOfRepeats;
    float ratioBf;
    float ratioAf;
    unsigned locationShift;
    float purity;
    float ratioOver20In;
    unsigned mateEditDist;
    CharString repSeq; //Repeat sequence in read
} ;

namespace std
{

template <>
struct hash<seqan::String<char> >
{
  std::size_t
  operator()(seqan::String<char> const & s) const
  {
      std::size_t seed = 42;

      for (char c : s)
          seed ^= c + 0x9e3779b9 + (c << 6) + (c >> 2);

      return seed;
  }
};

template <>
struct hash<seqan::Pair<int> >
{
  std::size_t
  operator()(seqan::Pair<int> const & p) const
  {
      std::size_t seed = 42;

      seed ^= p.i1 + 0x9e3779b9 + (p.i1 << 6) + (p.i1 >> 2);
      seed ^= p.i2 + 0x9e3779b9 + (p.i2 << 6) + (p.i2 >> 2);

      return seed;
  }
};

template <>
struct hash<STRinfoSmall >
{
  std::size_t
  operator()(STRinfoSmall const & p) const
  {
      std::size_t seed = 42;

      for (char c : p.chrom)
          seed ^= c + 0x9e3779b9 + (c << 6) + (c >> 2);
      seed ^= p.STRstart + 0x9e3779b9 + (p.STRstart << 6) + (p.STRstart >> 2);
      seed ^= p.STRend + 0x9e3779b9 + (p.STRend << 6) + (p.STRend >> 2);

      return seed;
  }
};

template <>
struct hash<seqan::Triple<CharString, CharString, int> >
{
  std::size_t
  operator()(seqan::Triple<CharString, CharString, int> const & p) const
  {
      std::size_t seed = 42;

      for (char c : p.i1)
          seed ^= c + 0x9e3779b9 + (c << 6) + (c >> 2);
      for (char c : p.i2)
          seed ^= c + 0x9e3779b9 + (c << 6) + (c >> 2);
      seed ^= p.i3 + 0x9e3779b9 + (p.i3 << 6) + (p.i3 >> 2);

      return seed;
  }
};

} // namespace std

//Vector to check how many repeats I need to find in a read w.r.t. motif length
std::vector<unsigned> repeatNumbers (6);
//Vector to store my repeat purity demands w.r.t motif length.
std::vector<Pair<float> > purityDemands (6);

string chroms[23]  = { "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX" };

//map storing a set of all hamming distance 1 sequences for each motif
unordered_map<Dna5String, std::unordered_set<Dna5String> > motifToPermutations;

//For repeating a motif n times
Dna5String repeat(Dna5String s, int n) {
    Dna5String ret;
    for (int i = 0; i < n; i++) {
        append(ret,s);
    }
    return ret;
}

//Create all one error permutations of motif passed
std::unordered_set<Dna5String> createPermutations(Dna5String motifRepeat)
{
    std::unordered_set<Dna5String> permutations;
    permutations.insert(motifRepeat);
    Dna5String bases = "ATCGN";
    for (unsigned i = 0; i < length(motifRepeat); ++i)
    {
        Dna5String motifCopy = motifRepeat;
        for(unsigned j = 0; j < length(bases); ++j)
        {
            motifCopy[i] = bases[j];
            if(motifCopy[i] != motifRepeat[i])
                permutations.insert(motifCopy);
        }
        erase(motifCopy,i);
        permutations.insert(motifCopy);
    }
    return permutations;
}

//Primitive pattern search, just searches for first and last occurence of #repeats*motif or hammingDistance1(#repeats*motif) without considering what's between them
Pair<Pair<int>, float> findPatternPrim(Dna5String motif, BamAlignmentRecord& record, const unsigned motifLength)
{
    //unsigned motifLength = length(motif);
    Dna5String pattern = repeat(motif,repeatNumbers[length(motif)-1]);
    std::unordered_set<Dna5String> permutations = motifToPermutations[motif];
    unsigned patternLength = length(pattern);
    int startCoordinate = length(record.seq);
    int endCoordinate = 0;
    unsigned index = 0;
    while(index+patternLength <= length(record.seq))
    {
        if (permutations.find(infixWithLength(record.seq, index, patternLength)) != permutations.end())
        {
            if (startCoordinate == length(record.seq))
                startCoordinate = index;
            endCoordinate = index + patternLength;
            index = index + motifLength;
        }
        else
        {
            if (permutations.find(infixWithLength(record.seq, index, patternLength-1)) != permutations.end())
            {
                if (startCoordinate == length(record.seq))
                    startCoordinate = index;
                endCoordinate = index + patternLength - 1;
                if (prefix(pattern, motifLength) == prefix(infixWithLength(record.seq, index, patternLength-1), motifLength))
                    index = index + motifLength;
                else
                    index = index + motifLength - 1;
            }
            else
                index = index + 1;
        }
    }
    float numOfRepeats = (float)(endCoordinate-startCoordinate)/(float)motifLength;
    if (endCoordinate == length(record.seq))
        endCoordinate -= 1;
    return Pair<Pair<int>, int>(Pair<int>(startCoordinate,endCoordinate), numOfRepeats);
}

//Check if read sequence contains motif and return start and end coordinates and number of repeats
Pair<Pair<int>, int> findPattern(Dna5String & pattern, Dna5String & readSequence, int motifLength)
{
    unsigned patternLength = length(pattern);
    std::unordered_set<Dna5String> permutations = createPermutations(pattern);
    int startCoordinate = length(readSequence);
    int endCoordinate = 0;
    unsigned moves = 0;
    int currentStart = length(readSequence);
    int currentEnd = 0;
    unsigned currentMoves = 0;
    unsigned index = 0;
    unsigned errors = 0;
    while(index+patternLength <= length(readSequence))
    {
        Dna5String theSubString = infixWithLength(readSequence, index, patternLength);
        Dna5String theSubStringMini = infixWithLength(readSequence, index, patternLength-1);
        if (permutations.find(theSubString) != permutations.end())
        {
            if (startCoordinate == length(readSequence))
            {
                startCoordinate = index;
                errors = 0;
                moves = 0;
            }
            endCoordinate = index + patternLength;
            index = index + motifLength;
            errors = 0;
            moves += 1;
        }
        else
        {
            if (permutations.find(theSubStringMini) != permutations.end())
            {
                if (startCoordinate == length(readSequence))
                {
                    startCoordinate = index;
                    errors = 0;
                    moves = 0;
                }
                endCoordinate = index + patternLength - 1;
                if (prefix(pattern, motifLength) == prefix(theSubStringMini, motifLength))
                    index = index + motifLength;
                else
                    index = index + motifLength - 1;
                errors = 0;
                moves += 1;
            }
            else
            {
                index = index + 1;
                ++errors;
            }
        }
        if (startCoordinate != length(readSequence) && errors > floor((float)motifLength/(float)2))
        {
            if (endCoordinate - startCoordinate > currentEnd - currentStart)
            {
                currentStart = startCoordinate;
                currentEnd = endCoordinate;
                currentMoves = moves;
            }
            startCoordinate = length(readSequence);
            endCoordinate = 0;
        }
    }
    unsigned numOfRepeats = moves + repeatNumbers[motifLength-1] - 1;
    if (endCoordinate - startCoordinate > currentEnd - currentStart)
        return Pair<Pair<int>, int>(Pair<int>(startCoordinate,endCoordinate), numOfRepeats);
    else
        return Pair<Pair<int>, int>(Pair<int>(currentStart,currentEnd), numOfRepeats);
}

//Compute the actual number of repeats/expected number of repeats based on length
float getPurity(Dna5String & motif, CharString & STRsequence)
{
    //cout << "getPurity( " << motif << "," << STRsequence << ")\n";
    unsigned motifLength = length(motif);
    unsigned expectReps = length(STRsequence)/motifLength;
    unsigned result = 0;
    unsigned index = 0;
    while(index+motifLength <= length(STRsequence))
    {
        //Dna5String theSubString = infixWithLength(STRsequence, index, motifLength);
        //cout << "Checking: " << theSubString << "\n";
        if (infixWithLength(STRsequence, index, motifLength) == motif)
        {
            result++;
            index = index + motifLength;
        }
        else
        {
            index = index + 1;
        }
    }
    return min (1.0f, (float)result/(float)expectReps);
}

//Finds ratio of bases in a sequence with PHRED score higher than 20
float findRatioOver20(CharString sequence)
{
    unsigned numOver20 = 0;
    int numVal;
    for (unsigned i = 0; i<length(sequence); ++i)
    {
        numVal = sequence[i] - 33;
        if(numVal>=20)
            ++numOver20;
    }
    return (float)numOver20/(float)(length(sequence));
}

//Get value of a tag (ATH! does not work if tag-value is not numeric)
int getTagValue(BamAlignmentRecord& record, CharString tagName)
{
    int returnValue;
    unsigned myIdx = 0;
    BamTagsDict tagsDict(record.tags);
    bool keyFound = findTagKey(myIdx, tagsDict, tagName);
    if (!keyFound)
    {
        cerr << "ERROR: Unknown key!" << tagName << " in read: " << record.qName << endl;
        returnValue = 999;
        return returnValue;
    }
    else
    {
        bool ok = extractTagValue(returnValue, tagsDict, myIdx);
        if (!ok)
        {
            cerr << "ERROR: There was an error extracting" << tagName << "from tags!\n";
            returnValue = 999;
            return returnValue;
        }
        else
            return returnValue;
    }
}

unsigned checkForAndRemoveSoftClippingAfter(BamAlignmentRecord& record)
{
    unsigned nRemoved = 0;
    String<CigarElement<> > cigarString = record.cigar;
    CharString cigarOperation = cigarString[length(cigarString)-1].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("S")==0)
    {
        nRemoved = cigarString[length(cigarString)-1].count;
        for (unsigned i = 0; i<cigarString[length(cigarString)-1].count; ++i)
        {
            eraseBack(record.seq);
            eraseBack(record.qual);
        }
        eraseBack(record.cigar);
    }
    return nRemoved;
}

unsigned checkForAndRemoveSoftClippingBefore(BamAlignmentRecord& record)
{
    unsigned nRemoved = 0;
    String<CigarElement<> > cigarString = record.cigar;
    CharString cigarOperation = cigarString[0].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("S")==0)
    {
        nRemoved = cigarString[0].count;
        erase(record.seq, 0, cigarString[0].count);
        erase(record.qual, 0, cigarString[0].count);
        erase(record.cigar, 0);
    }
    return nRemoved;
}

//Computes all sorts of quality indicators for the given read w.r.t. the given microsatellite
unsigned noFrontAlign = 0, noBackAlign = 0, frontAlign = 0, backAlign = 0;

Pair<Triple<CharString, CharString, int>,ReadInfo> computeReadInfo(BamAlignmentRecord& record, STRinfo& markerInfo, Pair<Pair<int>, float> coordinates, int minFlank, int maxRepeatLength)
{
    //Type definition for alignment structure
    typedef Dna5String TSequence;
    typedef Align<TSequence,ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    //cout << "Looking at read: " << record.qName << endl;
    //cout << "coordinates.i1.i1: " << coordinates.i1.i1 << " coordinates.i1.i2: " << coordinates.i1.i2 << "\n";

    //Variables for sequence-parts
    Dna5String before, after, before_8, after_8;
    CharString qualString = record.qual;
    ReadInfo mapValue;
    int oldStartCoord = coordinates.i1.i1;
    unsigned motifLength = length(markerInfo.motif), readLength = length(record.seq);
    //Create key for storing readInfo in map
    Triple<CharString, CharString, int> mapKey = Triple<CharString, CharString, int>(record.qName, markerInfo.chrom, markerInfo.STRstart);

    //Insert values into mapValue in returnPair
    mapValue.STRend = markerInfo.STRend;
    mapValue.motif = markerInfo.motif;
    mapValue.refRepeatNum = markerInfo.refRepeatNum;

    int scoreBf, scoreAf, startCoord, endCoord, leftFlank, rightFlank;
    float rBf, rAf;
    //check if 8-mer in front of and behind repeat in read match reference with edit-distance max 1
    //This is the best case scenario, then I skip flanking alignment
    before_8 = infix(record.seq, max(0,coordinates.i1.i1-8), coordinates.i1.i1);
    //cout << "Made before guy of length: " << length(before_8) << ".\n";
    after_8 = infix(record.seq, min(coordinates.i1.i2+1, (int)(readLength-1)),min(coordinates.i1.i2+9,(int)readLength));
    //cout << "Made after guy of length: " << length(after_8) << ".\n";
    bool doBeforeAlign = true, doAfterAlign = true;
    if (length(before_8) == 8)
    {
        //check this marker's 8mer edit distance 1 hash table for before_8
        if (markerInfo.hash8beforeSet.find(before_8) != markerInfo.hash8beforeSet.end())
        {
            //Skip alignment and set startCoord same as coordinates
            doBeforeAlign = false;
            startCoord = coordinates.i1.i1;
            scoreBf = 8;
            before = before_8;
            leftFlank = coordinates.i1.i1;
            rBf = 1.0f;
            //cout << "Found 8 mer in before map.\n";
        }
    }
    if (length(after_8) == 8)
    {
        //check this marker's 8mer edit distance 1 hash table for after_8
        if (markerInfo.hash8afterSet.find(after_8) != markerInfo.hash8afterSet.end())
        {
            //Skip alignment and set endCoord same as coordinates
            doAfterAlign = false;
            endCoord = coordinates.i1.i2 - coordinates.i1.i1;
            scoreAf = 8;
            after = after_8;
            rightFlank = readLength - coordinates.i1.i2;
            rAf = 1.0f;
            //cout << "Found 8 mer in after map.\n";
        }
    }
    if (markerInfo.minFlankLeft > 8)
        doBeforeAlign = true;
    if (markerInfo.minFlankRight > 8)
        doAfterAlign = true;
    //If I need to realign in front of repeat
    //cout << "doBeforeAlign: " << doBeforeAlign << " do AfterAlign: " << doAfterAlign << "\n";
    if (doBeforeAlign)
    {
        TAlign alignBefore;
        resize(rows(alignBefore), 2);
        before = prefix(record.seq, coordinates.i1.i2+1);
        assignSource(row(alignBefore,0), suffix(markerInfo.refBf,1000-length(before)-1));
        assignSource(row(alignBefore,1),before);
        scoreBf = globalAlignment(alignBefore, Score<int,Simple>(1,-2,-1,-5), AlignConfig<true, false, true, false>());
        int viewPosition = toViewPosition(row(alignBefore,0),length(source(row(alignBefore,0)))-1);
        startCoord = toSourcePosition(row(alignBefore,1),viewPosition)+1;
        if (startCoord == 1 && isGap(row(alignBefore,1),viewPosition))
            startCoord -= 1;
        leftFlank = startCoord;

        //Check wether the alignment score fails our criteria, if it does then we see if we can remove soft clipped sequence and realign.
        /*if ((float)scoreBf/(float)(startCoord+1) < 0.6 && !hasFlagUnmapped(record))
        {
            unsigned nRemoved = checkForAndRemoveSoftClippingBefore(record);
            coordinates.i1.i2 -= nRemoved;
            coordinates.i1.i1 -= nRemoved;
            oldStartCoord = max((int)0, (int)oldStartCoord - (int)nRemoved);
            if (nRemoved>0 && nRemoved <= coordinates.i1.i1+nRemoved)
            {
                before = prefix(record.seq, coordinates.i1.i2+1);
                assignSource(row(alignBefore,1),before);
                scoreBf = globalAlignment(alignBefore, Score<int,Simple>(1,-2,-1,-5), AlignConfig<true, false, true, false>());
                viewPosition = toViewPosition(row(alignBefore,0),length(source(row(alignBefore,0)))-1);
                startCoord = toSourcePosition(row(alignBefore,1),viewPosition)+1;
                if (startCoord == 1 && isGap(row(alignBefore,1),viewPosition))
                    startCoord -= 1;
                leftFlank = startCoord;
            }
        }*/
        //cout << "Finished before alignment\n";
        rBf = (float)scoreBf/(float)(startCoord+1);
        ++frontAlign;

        /*CharString ref_cs = suffix(markerInfo.refBf,1000-length(before)-1);
        const string ref = toCString(ref_cs);
        CharString query_cs = before;
        const string query = toCString(query_cs);
        int32_t maskLen = strlen(query.c_str())/2;
        maskLen = maskLen < 15 ? 15 : maskLen;

        // Declares a default Aligner
        StripedSmithWaterman::Aligner aligner(1,2,5,1);
        // Declares a default filter
        StripedSmithWaterman::Filter filter;
        // Declares an alignment that stores the result
        StripedSmithWaterman::Alignment alignment;
        // Aligns the query to the ref
        aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
        cout << "SeqAn\t" << scoreBf << "\tSW\t" << alignment.sw_score << "\tdelta\t" << alignment.sw_score - scoreBf<< "\n";
        cout << "SeqAnSt\t" << startCoord << "\tSW\t" << alignment.query_end << "\tdelta\t" << alignment.query_end - startCoord<< "\n";
        cout << "findPattern\t" << coordinates.i1.i1 << "-" << coordinates.i1.i2 << "\n";*/
    }
    else
        ++noFrontAlign;

    //If I need to realign behind repeat
    if (doAfterAlign)
    {
        TAlign alignAfter;
        resize(rows(alignAfter), 2);
        after = suffix(record.seq, coordinates.i1.i1);
        assignSource(row(alignAfter,0), prefix(markerInfo.refAf,length(after)+1));
        assignSource(row(alignAfter,1),after);
        scoreAf = globalAlignment(alignAfter, Score<int,Simple>(1,-2,-1,-5), AlignConfig<false, true, false, true>());
        endCoord = toSourcePosition(row(alignAfter,1),toViewPosition(row(alignAfter,0),0))-1;
        rightFlank = length(source(row(alignAfter,1)))-endCoord - 1;

        //Check wether the alignment score fails our criteria, if it does then we see if we can remove soft clipped sequence and realign.
        /*if ((float)scoreAf/(float)(length(after)-endCoord) < 0.6 && !hasFlagUnmapped(record))
        {
            unsigned readLength = length(record.seq);
            unsigned nRemoved = checkForAndRemoveSoftClippingAfter(record);
            if (nRemoved>0 && readLength - nRemoved > coordinates.i1.i2)
            {
                //redefine after region and realign
                after = suffix(record.seq, coordinates.i1.i1);
                assignSource(row(alignAfter,1),after);
                scoreAf = globalAlignment(alignAfter, Score<int,Simple>(1,-2,-1,-5), AlignConfig<false, true, false, true>());
                endCoord = toSourcePosition(row(alignAfter,1),toViewPosition(row(alignAfter,0),0))-1;
                rightFlank = length(source(row(alignAfter,1)))-endCoord - 1;
            }
        }*/
        rAf = (float)scoreAf/(float)(length(after)-endCoord-1);
        ++backAlign;
    }
    else
        ++noBackAlign;

    /*bool debug = true;
    if (debug)
    {
        int flankSum = leftFlank + rightFlank;
        cout << "FlankSum: " << flankSum << " leftFlank: " << leftFlank << " rightFlank: " << rightFlank << endl;
        if (leftFlank > 0)
            cout << "AlignmentScore before: " << (float)scoreBf/(float)startCoord << endl;
        if (rightFlank > 0)
            cout << "AlignmentScore after: " << (float)scoreAf/(float)(length(after)-endCoord-1) << endl;
        cout << "Start coord: " << startCoord << endl;
        cout << "Old start coord: " << oldStartCoord << endl;
        cout << "End coord: " << endCoord << endl;
        //cout << "Infix command is: infix(" << startCoord << "," << oldStartCoord+endCoord+1 << ")" << endl;
        //cout << "Repeat purity: " << getPurity(markerInfo.motif,infix(record.seq, startCoord, oldStartCoord+endCoord+1)) << endl;
    }*/

    float refRepPurity = markerInfo.refRepPurity;
    bool startOk = false;
    bool endOk = false;
    bool purityOk = false;
    int flankSum = leftFlank + rightFlank;
    if ((startCoord >= oldStartCoord + endCoord) || (((oldStartCoord+endCoord+1)-startCoord < maxRepeatLength) && (leftFlank < markerInfo.minFlankLeft || rightFlank < markerInfo.minFlankRight)))
    {
        //cout << "Failing flanking and start/end coord check." << endl;
        coordinates.i1.i1 = startCoord;
        coordinates.i1.i2 = endCoord;
        mapValue.ratioBf = 0;
        mapValue.ratioAf = 0;
        mapValue.numOfRepeats = 666;
        //Debugging code
        //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl;
        mapValue.ratioOver20In = 0;
        mapValue.purity = 0;
        mapValue.repSeq = "";
        mapValue.locationShift = 100;
        return Pair<Triple<CharString, CharString, int>,ReadInfo>(mapKey,mapValue);
    }
    mapValue.repSeq = infix(record.seq, startCoord, oldStartCoord+endCoord+1);
    float readPurity = getPurity(markerInfo.motif, mapValue.repSeq);
    unsigned startCoord1 = startCoord, endCoord1 = oldStartCoord+endCoord+1;
    /*cout << "Read purity: " << readPurity << "\n";
    cout << "Reference purity: " << refRepPurity << "\n";
    cout << "Purity minimum: " << purityDemands[motifLength-1].i1 << "*" << refRepPurity << "=" << purityDemands[motifLength-1].i1*refRepPurity << "\n";*/
    if (readPurity <= purityDemands[motifLength-1].i1*refRepPurity)
    {
        //cout << "Failing read purity check" << endl;
        coordinates.i1.i1 = startCoord;
        coordinates.i1.i2 = endCoord;
        mapValue.ratioBf = 0;
        mapValue.ratioAf = 0;
        mapValue.numOfRepeats = 666;
        //Debugging code
        //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl;
        mapValue.ratioOver20In = 0;
        mapValue.purity = 0;
        mapValue.repSeq = "";
        mapValue.locationShift = 100;
        return Pair<Triple<CharString, CharString, int>,ReadInfo>(mapKey,mapValue);
    }
    else
        purityOk = true;
    //Check if a minimum number of repeats have been found and whether the flanking area is sufficient for both start and end coordinates
    if ((readLength - leftFlank - rightFlank >=length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)-1]) && (leftFlank>=minFlank))
        startOk = true;
    if ((endCoord >= length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)-1]-1)&&(rightFlank >= minFlank))
        endOk = true;
    //I allow only marker specific minimum number of aligning bases on either side if I have more than 2*minFlank aligned bases in total.
    if (flankSum >= 2*minFlank && leftFlank >= markerInfo.minFlankLeft && rightFlank >= markerInfo.minFlankRight)
    {
        startOk = true;
        endOk = true;
    }
    /*if (debug)
    {
        cout << "Start ok: " << startOk << endl;
        cout << "End ok: " << endOk << endl;
        cout << "Purity ok: " << purityOk << endl;
    }*/
    if (((oldStartCoord+endCoord+1)-startCoord >= maxRepeatLength) && (leftFlank < markerInfo.minFlankLeft) && (flankSum>=2*minFlank) && (rAf>0.7) && (readPurity>(purityDemands[motifLength-1].i2*refRepPurity)))
    {
        //cout << "Am making greater than allele on left end." << endl;
        coordinates.i1.i1 = startCoord;
        coordinates.i1.i2 = endCoord;
        mapValue.ratioBf = 0;
        mapValue.ratioAf = rAf;
        //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl;
        mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
        mapValue.ratioOver20In = findRatioOver20(infix(qualString, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1));
        //cout << "getPurity( " << markerInfo.motif << "," << infix(record.seq, startCoord, oldStartCoord+endCoord+1) << ")\n";
        if (startCoord1!=coordinates.i1.i1 || endCoord1 != oldStartCoord+coordinates.i1.i2+1)
        {
            mapValue.repSeq = infix(record.seq, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1);
            mapValue.purity = getPurity(markerInfo.motif,mapValue.repSeq);
        }
        else
            mapValue.purity = readPurity;
    }
    else
    {
        if (((oldStartCoord+endCoord+1)-startCoord >= maxRepeatLength) && (rightFlank < markerInfo.minFlankRight) && (flankSum>=2*minFlank) && (rBf>0.7) && (readPurity>(purityDemands[motifLength-1].i2*refRepPurity)))
        {
            //cout << "Am making greater than allele on right end." << endl;
            coordinates.i1.i1 = startCoord;
            coordinates.i1.i2 = endCoord;
            mapValue.ratioBf = rBf;
            mapValue.ratioAf = 0;
            //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl;
            mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
            mapValue.ratioOver20In = findRatioOver20(infix(qualString, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1));
            if (startCoord1!=coordinates.i1.i1 || endCoord1 != oldStartCoord+coordinates.i1.i2+1)
            {
                mapValue.repSeq = infix(record.seq, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1);
                mapValue.purity = getPurity(markerInfo.motif,mapValue.repSeq);
            }
            else
                mapValue.purity = readPurity;
        }
        else
        {
            if (rightFlank < markerInfo.minFlankRight && leftFlank < markerInfo.minFlankLeft && readPurity>(purityDemands[motifLength-1].i2*refRepPurity) && (oldStartCoord+endCoord+1)-startCoord >= maxRepeatLength && readLength>=maxRepeatLength)
            {
                //cout << "Am making a SUPER allele." << endl;
                coordinates.i1.i1 = 0;
                coordinates.i1.i2 = readLength-1;
                mapValue.ratioBf = 0.0;
                mapValue.ratioAf = 0.0;
                //cout << "Infix command is: infix(" << 0 << "," << length(record.seq)-1 << ")" << endl;
                mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
                mapValue.ratioOver20In = findRatioOver20(record.seq);
                if (length(mapValue.repSeq) != endCoord1-startCoord1)
                {
                    mapValue.repSeq = record.seq;
                    mapValue.purity = getPurity(markerInfo.motif,mapValue.repSeq);
                }
                else
                    mapValue.purity = readPurity;
            }
            else
            {
                //If both coordinates are ok, the distance between them is > motifLength * min#ofMotifs and alignment scores on both end are ok then the read is useful.
                if (startOk && endOk && purityOk && startCoord < oldStartCoord + endCoord && rBf>0.6 && rAf>0.6)
                {
                    coordinates.i1.i1 = startCoord;
                    coordinates.i1.i2 = endCoord;
                    mapValue.ratioBf = rBf;
                    mapValue.ratioAf = rAf;
                    //Debugging code
                    //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl;
                    mapValue.repSeq = infix(record.seq, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1);
                    mapValue.numOfRepeats = (float)length(mapValue.repSeq)/(float)length(markerInfo.motif);
                    if (length(mapValue.repSeq)>=maxRepeatLength)
                        mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
                    //If I find less repeats than the required minimum (depending on the motif length) then I don't use the read
                    if (ceil(mapValue.numOfRepeats) < repeatNumbers[length(markerInfo.motif)-1])
                    {
                        //cout << "Not enough repeats." << endl;
                        mapValue.numOfRepeats = 666;
                    }
                    mapValue.ratioOver20In = findRatioOver20(infix(qualString, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1));
                    if (startCoord1 != coordinates.i1.i1 || endCoord1 != oldStartCoord+coordinates.i1.i2+1)
                        mapValue.purity = getPurity(markerInfo.motif,mapValue.repSeq);
                    else
                        mapValue.purity = readPurity;
                    //cout << "Processed read: " << record.qName << " into map and reported: " << mapValue.numOfRepeats << ".\n";
                }
                //Otherwise I can't use the read so I set numOfRepeats to 666
                else
                {
                    //cout << "Failing alignment score requirements" << endl;
                    //cout << "Setting num of repeats to 666 in the last else statement." << endl;
                    coordinates.i1.i1 = startCoord;
                    coordinates.i1.i2 = endCoord;
                    mapValue.ratioBf = 0;
                    mapValue.ratioAf = 0;
                    mapValue.numOfRepeats = 666;
                    //Debugging code
                    //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl;
                    mapValue.repSeq = "";
                    mapValue.ratioOver20In = 0;
                    mapValue.purity = 0;
                }
            }
        }
    };

    //Debugging code
    /*cout << "Motif: " << mapValue.motif << " Start: " << markerInfo.STRstart << endl;
    cout << "Before: " << prefix(before, coordinates.i1.i1) << endl;
    cout << "Aligned to: " << suffix(source(row(alignBefore,0)),toSourcePosition(row(alignBefore,0),toViewPosition(row(alignBefore,1),0))) << endl;
    cout << "Repeat: " << mapValue.repSeq << endl;
    cout << "After: " << suffix(after, coordinates.i1.i2+1) << endl;
    cout << "Aligned to: " << prefix(source(row(alignAfter,0)),toSourcePosition(row(alignAfter,0),toViewPosition(row(alignAfter,1),length(source(row(alignAfter,1)))))) << endl;
    cout << "Number of repeats: " << mapValue.numOfRepeats << endl;
    if (mapValue.numOfRepeats == 666)
    {
        if (!startOk)
        {
            cout << "Enough repeats according to start? " << ((length(source(row(alignBefore,1)))-(startCoord+1)) >= (length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)-1])) << endl;
            cout << "Enough flanking area in front? " << (startCoord>=minFlank) << endl;
        }
        if (!endOk)
        {
            cout << "Enough repeats according to end? " << (endCoord >= length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)-1]-1) << endl;
            cout << "Enough flanking area behind ?" << (length(source(row(alignAfter,1)))-endCoord >= minFlank) << endl;
        }
        cout << "Start and end not overlapping? " << ((startCoord + length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)-1])< (startCoord + endCoord)) << endl;
    }*/

    //Check location shift
    mapValue.locationShift = abs(oldStartCoord - startCoord);
    /*if (mapValue.numOfRepeats != 666)
        cout << "Used this read and gave: " << mapValue.numOfRepeats << "repeats.\n";
    else
        cout << "Didn't use this read." << endl;*/
    return Pair<Triple<CharString, CharString, int>,ReadInfo>(mapKey,mapValue);
}

String<STRinfo> readMarkerinfo(CharString & markerInfoFile, int minFlank, CharString & attributeDirectory, unordered_map<Pair<int>, Pair<String<long int> ,FILE*> > & startAndEndToStreamAndOffsets, unsigned nPns)
{
    String<STRinfo> markers;
    float freq;
    STRinfo currInfo;
    //Read all markers into memory
    ifstream markerFile(toCString(markerInfoFile));
    while (!markerFile.eof())
    {
        string chromString;
        markerFile >> chromString;
        if (chromString.length() == 0)
            break;
        currInfo.chrom = chromString;
        markerFile >> currInfo.STRstart;
        markerFile >> currInfo.STRend;
        string motifString;
        markerFile >> motifString;
        currInfo.motif = motifString;
        motifToPermutations[currInfo.motif] = createPermutations(repeat(currInfo.motif,repeatNumbers[length(currInfo.motif)-1]));
        markerFile >> currInfo.refRepeatNum;
        string refBfString;
        markerFile >> refBfString;
        currInfo.refBf = refBfString;
        string refAfString;
        markerFile >> refAfString;
        currInfo.refAf = refAfString;
        string refRepSeq;
        markerFile >> refRepSeq;
        currInfo.refRepeatNum = (float)refRepSeq.length()/(float)motifString.length();
        markerFile >> currInfo.minFlankLeft;
        markerFile >> currInfo.minFlankRight;
        markerFile >> currInfo.refRepPurity;
        markerFile >> freq;
        currInfo.motifFrqs.push_back(freq);
        markerFile >> freq;
        currInfo.motifFrqs.push_back(freq);
        markerFile >> freq;
        currInfo.motifFrqs.push_back(freq);
        markerFile >> freq;
        currInfo.motifFrqs.push_back(freq);
        currInfo.hash8beforeSet = createPermutations(suffix(currInfo.refBf,992));
        currInfo.hash8afterSet = createPermutations(prefix(currInfo.refAf,8));

        //Make output file for current marker
        CharString currAttributeDirectory = attributeDirectory;
        append(currAttributeDirectory, to_string(currInfo.STRstart));
        append(currAttributeDirectory, "_");
        append(currAttributeDirectory, currInfo.motif);
        //cout << "Outputfile: " << currAttributeDirectory << "\n";
        Pair<int> mapKey = Pair<int>(currInfo.STRstart, currInfo.STRend);
        FILE *fp;
        fp = fopen(toCString(currAttributeDirectory), "w+");
        startAndEndToStreamAndOffsets[mapKey].i2 = fp;
        if (startAndEndToStreamAndOffsets[mapKey].i2 == NULL)
            cerr << "Couldn't make ofstream for marker @ " << currAttributeDirectory << "\n";
        else
        {
            if (fseek(startAndEndToStreamAndOffsets[mapKey].i2, nPns*(11), SEEK_SET) != 0)
                cerr << "Could not reserve index space at front of " << currAttributeDirectory << "\n";
            else
            {
                resize(startAndEndToStreamAndOffsets[mapKey].i1, nPns, 0);
                fprintf(startAndEndToStreamAndOffsets[mapKey].i2, "\r\n");
            }
        }
        appendValue(markers, currInfo);
    }
    return markers;
}

String<Pair<CharString> > readBamList(CharString & bamListFile)
{
    String<Pair<CharString> > PnsAndBams;
    Pair<CharString> currPnAndBam;
    //Read all Pns and their bamFiles
    ifstream bamList(toCString(bamListFile));
    while (!bamList.eof())
    {
        string PN_ID;
        bamList >> PN_ID;
        if (PN_ID.length() == 0)
            break;
        currPnAndBam.i1 = PN_ID;
        string bamPath;
        bamList >> bamPath;
        currPnAndBam.i2 = bamPath;
        appendValue(PnsAndBams, currPnAndBam);
    }
    return PnsAndBams;
}

vector<string> splitString(string &s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string next;
    while(getline(ss, next, delim)) {
        elems.push_back(next);
    }
    return elems;
}

void readExpansions(CharString & expansionFile, unordered_map<string, String<Triple<vector<float>,unsigned,unsigned> > > & chromToMotifAndIntervals)
{
    Triple<vector<float>, unsigned, unsigned> motifStartEnd;
    string chrom, length, line;
    vector<string> words;
    ifstream expansions(toCString(expansionFile));
    while (getline(expansions, line))
    {
        words = splitString(line, ' ');
        chrom = words[0];
        motifStartEnd.i2 = stoi (words[1]);
        motifStartEnd.i3 = stoi (words[2]);
        motifStartEnd.i1.resize(words.size()-5);
        transform (words.begin() + 5, words.end(), motifStartEnd.i1.begin(), [](std::string const & str){return std::stof(str);});
        appendValue(chromToMotifAndIntervals[chrom], motifStartEnd);
    }
    cout << "Finished reading expansions.\n";
    return;
}

bool compareMotifs(const std::vector<float> & expansionBpFreqs, const std::vector<float> & motifBpFreqs)
{
    unsigned nMotifs = expansionBpFreqs.size() / 4;
    for (unsigned i = 0; i < nMotifs; i++)
    {
        if (std::equal(expansionBpFreqs.cbegin()+(i*4), expansionBpFreqs.cbegin()+(i*4)+4, motifBpFreqs.cbegin()))
            return true;
        if (std::equal(expansionBpFreqs.cbegin()+(i*4), expansionBpFreqs.cbegin()+(i*4)+4, motifBpFreqs.crbegin()))
            return true;
    }
    return false;
}

bool chromIsReal(CharString chrom)
{
    string * p = std::find (chroms, chroms+22, toCString(chrom));
    if (p != chroms+22)
        return true;
    else
        return false;
}

/*unsigned storeReadsAtExpansions(CharString & bamPath, const char* reference, map<string, String<Triple<string,unsigned,unsigned> > > & chromToMotifAndIntervals, map <CharString, Triple<BamAlignmentRecord, string, unsigned> > & readNameToExpansion)
{
    HtsFile hts_file(toCString(bamPath), "r", reference);
    if (!loadIndex(hts_file))
    {
        std::cerr << "ERROR: Could not read index file for " << bamPath << "\n";
        return 1;
    }

    for (auto chrom: chromToMotifAndIntervals)
    {
        for (auto interval: chrom.second)
        {
            Triple<string, unsigned, unsigned> key = Triple<string, unsigned, unsigned>(chrom.first, interval.i2, interval.i3);
            string motif = interval.i1;
            CharString chromosome = chrom.first;
            cout << "Checking possibly expanded reads at " << chromosome << ":" << interval.i2 << "-" << interval.i3 << endl;
            if (!setRegion(hts_file, toCString(chromosome), (int) interval.i2, (int) interval.i3))
            {
                std::cerr << "ERROR: Could not jump to " << chromosome << ":" << interval.i2 << "\n";
                return 1;
            }
            BamAlignmentRecord record;
            while (readRegion(record, hts_file))
            {
                if (hasFlagQCNoPass(record) || hasFlagDuplicate(record) || hasFlagNextUnmapped(record) || hasFlagUnmapped(record) || chromIsReal(hts_file.hdr->target_name[record.rNextId]))
                    continue;
                if (record.rID != record.rNextId || abs(record.tLen) > 100000)
                {
                    Triple<BamAlignmentRecord, string, unsigned> value = Triple<BamAlignmentRecord, string, unsigned>(record, motif, interval.i3 - interval.i2);
                    readNameToExpansion[record.qName] = value;
                }
            }
        }
    }
    return 0;
}*/

BamAlignmentRecord findMate(BamAlignmentRecord & record, CharString & bamPath, const char* reference)
{
    BamAlignmentRecord mate;
    HtsFile hts_file(toCString(bamPath), "r", reference);
    if (!loadIndex(hts_file))
    {
        std::cerr << "ERROR: Could not read index file for " << bamPath << "\n";
        return mate;
    }
    CharString chromosome = hts_file.hdr->target_name[record.rNextId];
    if (!setRegion(hts_file, toCString(chromosome), record.pNext-10, record.pNext+10))
    {
        std::cerr << "ERROR: Could not jump to " << chromosome << ":" << record.pNext << "\n";
        return mate;
    }
    while (readRegion(mate, hts_file))
    {
        if (mate.qName == record.qName)
            break;
    }
    return mate;
}

bool lookForExpansion(bam_hdr_t * header,
                      unsigned chrIdx,
                      unsigned pos,
                      unsigned seqLength,
                      const vector<float> & bpFreqMotif,
                      unsigned markerSize,
                      unordered_map<string, String<Triple<vector<float>,unsigned,unsigned> > > & chromToMotifAndIntervals)
{
    Triple<vector<float>,unsigned,unsigned> pos_triple;
    pos_triple.i3 = pos;
    //std::vector<float>bpFreqRCMotif(bpFreqMotif.rbegin(), bpFreqMotif.rend()) ;
    unsigned endPos = pos + seqLength - 1;
    string chr = toCString(header->target_name[chrIdx]);
    if (chromToMotifAndIntervals.count(chr) == 0)
        return false;

    String<Triple<vector<float>,unsigned,unsigned> > const & intervals = chromToMotifAndIntervals[chr];
    auto it = std::lower_bound(begin(intervals), end(intervals), pos_triple,
                     [](Triple<vector<float>,unsigned,unsigned> const & a,
                        Triple<vector<float>,unsigned,unsigned> const & b)
                     {return a.i3 < b.i3;}
      );
    for (; it != end(intervals); ++it)
    {
        if (std::min((int) endPos, (int) it->i3)-std::max((int) pos, (int) it->i2) >= 50 && it->i3 - it->i2 > markerSize)
        {
            if (compareMotifs(it->i1, bpFreqMotif))// || compareMotifs(it->i1, bpFreqRCMotif))
                return true;
            else
                return false;
        }
        else if (it->i2 > endPos)
        {
            return false;
        }
    }
    return false;
}

int main(int argc, char * argv[])
{
    time_t begin = time(0);
    //Check arguments.
    if (argc != 10)
    {
        cerr << "USAGE: " << argv[0] << " bamFiles outputDirectory markerInfoFile minFlankLength maxRepeatLength chrom reference expansions jumpToExpansions[Y/N]\n";
        return 1;
    }
    //Store maximum repeat length
    int maxRepeatLength = lexicalCast<int>(argv[5]), minFlank = lexicalCast<int>(argv[4]);
    //Store parameters
    CharString bamListFile = argv[1], markerInfoFile = argv[3], attributeDirectory = argv[2], chrom = argv[6], expansionFile = argv[8];
    const char* reference = argv[7];
    bool jumpAround = false;
    string arg9 = argv[9];
    if (arg9 == "Y")
        jumpAround = true;
    //Read pn info
    String<Pair<CharString> > PnsAndBams = readBamList(bamListFile);
    unsigned nPns = length(PnsAndBams);
    if (nPns==0)
    {
        cerr << "The PN and BAM file list supplied is empty, please supply a file containing at least one PN and its BAM.\n";
    }
    cout << "Finished reading pn info, number of PNs: " << nPns << "\n";

    //Map from chrom to motif and interval
    unordered_map<string, String<Triple<vector<float>,unsigned,unsigned> > > chromToMotifAndIntervals;
    //read know expansions in reference
    readExpansions(expansionFile, chromToMotifAndIntervals);
    //Map to store possibly expanded reads if jumpAround = true
    //map <CharString, Triple<BamAlignmentRecord, string, unsigned> > readNameToExpansion;

    struct stat st2;
    if (stat(toCString(attributeDirectory),&st2) != 0)
    {
        cerr << "Output directory does not exist: " << attributeDirectory << "\n";
        return 1;
    }
    append(attributeDirectory, "/attributes/");
    struct stat st;
    if (stat(toCString(attributeDirectory),&st) != 0)
        mkdir(toCString(attributeDirectory),0777);
    append(attributeDirectory, chrom);
    struct stat st3;
    if(stat(toCString(attributeDirectory),&st3) != 0)
        mkdir(toCString(attributeDirectory),0777);
    append(attributeDirectory, "/");
    //Set up a map from marker to String of long int values to store offset of each PN in each markerFile.
    unordered_map<STRinfoSmall, String<long int> > markerToPnOffsets;
    //Map from start and end coordinate to output stream
    unordered_map<Pair<int>, Pair<String<long int>, FILE* > > startAndEndToStreamAndOffsets;
    //Set up how many repeats I require for each motif length
    repeatNumbers[0]=10;
    repeatNumbers[1]=4;
    repeatNumbers[2]=3;
    repeatNumbers[3]=2;
    repeatNumbers[4]=2;
    repeatNumbers[5]=2;
    //Read marker info
    String<STRinfo> markers = readMarkerinfo(markerInfoFile, minFlank, attributeDirectory, startAndEndToStreamAndOffsets, nPns);
    if (length(markers)==0)
    {
        cerr << "The markerInfo file is empty, please supply a file containing at least one marker.\n";
        return 1;
    }
    cout << "Finished reading marker Info, number of markers: " << length(markers) << "\n";
    unsigned finalMarkerIdx = length(markers) - 1;

    //set up the purity I require for each motif length
    purityDemands[0] = Pair<float>(0.9,0.95);
    purityDemands[1] = Pair<float>(0.85,0.9);
    purityDemands[2] = Pair<float>(0.8,0.85);
    purityDemands[3] = Pair<float>(0.8,0.85);
    purityDemands[4] = Pair<float>(0.75,0.85);
    purityDemands[5] = Pair<float>(0.75,0.85);

    //START FOR LOOP OVER PNS HERE
    for (unsigned i=0; i<length(PnsAndBams); ++i)
    {
        //Time process for each PN
        time_t pnStart = time(0);

        //Set up hts file and jump to correct chromosome
        CharString PN_ID = PnsAndBams[i].i1;
        HtsFile hts_file(toCString(PnsAndBams[i].i2), "r", reference);

        int jumpStart = std::max(0,markers[0].STRstart - 1000);
        int jumpEnd = jumpStart + 300000000;
        if (!loadIndex(hts_file))
        {
            std::cerr << "ERROR: Could not read index file for " << PnsAndBams[i].i2 << "\n";
            return 1;
        }

        if (!setRegion(hts_file, toCString(markers[0].chrom), jumpStart, jumpEnd))
        {
            cerr << "ERROR: Could not jump to " << markers[0].chrom << ":" << 0 << "\n";
            return 1;
        }

        //Variables for the start and end coordinates of reads and their mates, length of read, how many repeats to look for and index into string storing marker information
        unsigned bamStart, bamEnd, mateStart, mateEnd, readLength, markerIndex = 0, motifLength = length(markers[currentMarker].motif), numToLook = repeatNumbers[motifLength-1];
        //Map from read name, marker chromosome and marker start to info on read-pair with that read name
        unordered_map<Triple<CharString, CharString, int>, ReadInfo> myMap;
        BamAlignmentRecord record;;
        while (readRegion(record, hts_file))
        {
            //If the read is a duplicate or doesn't pass the quality check I move on
            if (hasFlagQCNoPass(record) || hasFlagDuplicate(record) || !hasFlagMultiple(record))
                continue;
            bamStart = record.beginPos;
            readLength = length(record.seq);
            if (!hasFlagUnmapped(record))
                bamEnd = bamStart + getAlignmentLengthInRef(record);
            else
                bamEnd = bamStart + readLength;
            mateStart = record.pNext;
            mateEnd = mateStart + readLength; //Can't get alignment length in ref for mate, use sequence length
            //Have passed the current marker? -> update markerIndex
            while (bamStart > markers[markerIndex].STRend + 1000)
            {
                ++markerIndex;
                motifLength = length(markers[currentMarker].motif);
                numToLook = repeatNumbers[motifLength-1];
                //If the markerIndex has exceeded the length of the marker string I can't check the while condition (markers[markerIndex].STRend doesn't exist) so I break
                if (markerIndex > finalMarkerIdx)
                    break;
            }
            //If the markerIndex has exceeded the length of the marker string I break the BAM-reading loop
            if(markerIndex > finalMarkerIdx)
                break;
            int currentMarker = markerIndex;

            //Loop for comparing current read to all possible markers, could be useful for more than one marker
            while (true)
            {
                // Have not reached beginning of the interval? -> next read
                if (bamEnd < max(markers[currentMarker].STRstart-1000,0))
                    break;
                //If the read is not long enough for the current marker I need to check the next marker
                if (markers[currentMarker].STRend-markers[currentMarker].STRstart > readLength)
                {
                    ++currentMarker;
                    if (currentMarker > finalMarkerIdx)
                        break;
                    motifLength = length(markers[currentMarker].motif);
                    numToLook = repeatNumbers[motifLength-1];
                    continue;
                }
                // Does the read intersect the microsatellite without being unaligned and the mate is within a reasonable distance? -> Then I look for the current repeat motif
                if (!hasFlagUnmapped(record) && record.rID == record.rNextId && abs(record.tLen) < 100000 && std::min((int) bamEnd, markers[currentMarker].STRend)-std::max((int) bamStart, markers[currentMarker].STRstart) > 0 )
                {
                    Pair<Pair<int>, float> coordinates = findPatternPrim(markers[currentMarker].motif, record, motifLength);
                    int startCoordinate = coordinates.i1.i1;
                    int endCoordinate = coordinates.i1.i2;
                    Pair<Triple<CharString, CharString, int>,ReadInfo> keyValuePair;
                    if (endCoordinate>startCoordinate)
                        keyValuePair.i2.repSeq = infixWithLength(record.seq, startCoordinate, endCoordinate-startCoordinate+1);
                    //Does the read contain all of the microsatellite and flanking regions larger than the set minimum? -> Then I compute read attributes
                    if (((endCoordinate-startCoordinate+1 < readLength-2*minFlank) && (endCoordinate > startCoordinate)) || ((endCoordinate > startCoordinate)&&(getPurity(markers[currentMarker].motif, keyValuePair.i2.repSeq)>0.75*markers[currentMarker].refRepPurity)))
                    {
                        keyValuePair = computeReadInfo(record, markers[currentMarker], coordinates, minFlank, maxRepeatLength);
                        if (myMap.count(keyValuePair.i1) == 0)
                        {
                            if (keyValuePair.i2.numOfRepeats == 666)
                            {
                                myMap[keyValuePair.i1].numOfRepeats = 666;
                                myMap[keyValuePair.i1].mateEditDist = getTagValue(record, "NM");
                            }
                            else
                            {
                                keyValuePair.i2.mateEditDist = 666;
                                myMap[keyValuePair.i1] = keyValuePair.i2;
                            }
                        }
                        else
                        {
                            if (myMap[keyValuePair.i1].numOfRepeats == 666)
                            {
                                keyValuePair.i2.mateEditDist = myMap[keyValuePair.i1].mateEditDist;
                                myMap[keyValuePair.i1] = keyValuePair.i2;
                            }
                            else
                                myMap[keyValuePair.i1].mateEditDist = getTagValue(record, "NM");
                        }
                        ++currentMarker;
                        //If I've reached the end of the marker string I break this loop and take the next read
                        if (currentMarker > finalMarkerIdx)
                            break;
                        motifLength = length(markers[currentMarker].motif);
                        numToLook = repeatNumbers[motifLength-1];
                        continue;
                    }
                }
                // Does the read's mate intersect the microsatellite and the read is not unaligned? -> Then I want the read's edit distance!
                if (!hasFlagUnmapped(record) && record.rID == record.rNextId && abs(record.tLen) < 100000 && std::min((int) mateEnd, markers[currentMarker].STRend)-std::max((int) mateStart, markers[currentMarker].STRstart) > 0)
                {
                    Triple<CharString, CharString, int> mapKey = Triple<CharString, CharString, int>(record.qName, markers[currentMarker].chrom, markers[currentMarker].STRstart);
                    if (myMap.count(mapKey) == 0)
                        myMap[mapKey].numOfRepeats = 666;
                    myMap[mapKey].mateEditDist = getTagValue(record, "NM");
                    ++currentMarker;
                    //If I've reached the end of the marker string I break this loop and take the next read
                    if (currentMarker > finalMarkerIdx)
                        break;
                    motifLength = length(markers[currentMarker].motif);
                    numToLook = repeatNumbers[motifLength-1];
                    continue;
                }
                //Is the read aligned to either 1000 bp before or after the microsatellite-> Then I check it out!
                if ((bamEnd >= max(markers[currentMarker].STRstart-1000,0) && bamEnd <= markers[currentMarker].STRstart) || (bamStart >= markers[currentMarker].STRend && bamStart <= markers[currentMarker].STRend+1000))
                {
                    //Is the reads mate unaligned? -> Then I want the reads edit distance!
                    if (hasFlagNextUnmapped(record))
                    {
                        Triple<CharString, CharString, int> mapKey = Triple<CharString, CharString, int>(record.qName, markers[currentMarker].chrom, markers[currentMarker].STRstart);
                        if (myMap.count(mapKey) == 0)
                            myMap[mapKey].numOfRepeats = 666;
                        myMap[mapKey].mateEditDist = getTagValue(record, "NM");
                        ++currentMarker;
                        //If I've reached the end of the marker string I break this loop and take the next read
                        if (currentMarker > finalMarkerIdx)
                            break;
                        motifLength = length(markers[currentMarker].motif);
                        numToLook = repeatNumbers[motifLength-1];
                        continue;
                    }
                    //Is the read itself unaligned -> Then I look for the current repeat motif
                    if (hasFlagUnmapped(record))
                    {
                        Pair<Pair<int>, float> coordinates = findPatternPrim(markers[currentMarker].motif, record, motifLength);
                        int startCoordinate = coordinates.i1.i1;
                        int endCoordinate = coordinates.i1.i2;
                        Pair<Triple<CharString, CharString, int>,ReadInfo> keyValuePair;
                        if (endCoordinate>startCoordinate)
                            keyValuePair.i2.repSeq = infixWithLength(record.seq, startCoordinate, endCoordinate-startCoordinate+1);
                        //Does the read contain all of the microsatellite and flanking regions larger than the set minimum? -> Then I compute read attributes
                        if (((endCoordinate-startCoordinate+1 < length(record.seq)-2*minFlank) && (endCoordinate > startCoordinate)) || ((endCoordinate > startCoordinate)&&(getPurity(markers[currentMarker].motif, keyValuePair.i2.repSeq)>0.75*markers[currentMarker].refRepPurity)))
                        {
                            keyValuePair = computeReadInfo(record, markers[currentMarker], coordinates, minFlank, maxRepeatLength);
                            if (myMap.count(keyValuePair.i1) == 0)
                            {
                                keyValuePair.i2.mateEditDist = 666;
                                myMap[keyValuePair.i1] = keyValuePair.i2;
                            }
                            else
                            {
                                keyValuePair.i2.mateEditDist = myMap[keyValuePair.i1].mateEditDist;
                                myMap[keyValuePair.i1] = keyValuePair.i2;
                            }
                            ++currentMarker;
                            //If I've reached the end of the marker string I break this loop and take the next read
                            if (currentMarker > finalMarkerIdx)
                                break;
                            motifLength = length(markers[currentMarker].motif);
                            numToLook = repeatNumbers[motifLength-1];
                            continue;
                        }
                    }
                }
                //Read's mate aligned to another chromosome
                if ((record.rID != record.rNextId || abs(record.tLen) > 100000) &&
                    std::min((int) bamEnd, markers[currentMarker].STRend+1000) -
                    std::max((int) bamStart, markers[currentMarker].STRstart-1000) > readLength-1 &&
                    length(markers[currentMarker].motif) > 1)
                {
                    //Find mate and check if it has repeat tract.
                    if (lookForExpansion(hts_file.hdr,
                                         record.rNextId,
                                         record.pNext,
                                         readLength,
                                         markers[currentMarker].motifFrqs,
                                         markers[currentMarker].STRend - markers[currentMarker].STRstart,
                                         chromToMotifAndIntervals))
                    {
                        if (jumpAround)
                        {
                            BamAlignmentRecord mateAtRepeat = findMate(record, PnsAndBams[i].i2, reference);
                            Pair<Pair<int>, float> coordinates = findPatternPrim(markers[currentMarker].motif, mateAtRepeat, motifLength);
                            int startCoordinate = coordinates.i1.i1;
                            int endCoordinate = coordinates.i1.i2;
                            if (endCoordinate-startCoordinate+1 >= readLength-2*minFlank && endCoordinate > startCoordinate)
                            {
                                Pair<Triple<CharString, CharString, int>,ReadInfo> keyValuePair;
                                keyValuePair.i1 = Triple<CharString, CharString, int>(record.qName, markers[currentMarker].chrom, markers[currentMarker].STRstart);
                                keyValuePair.i2.mateEditDist = getTagValue(record, "NM");
                                keyValuePair.i2.STRend = markers[currentMarker].STRend;
                                keyValuePair.i2.motif = markers[currentMarker].motif;
                                keyValuePair.i2.refRepeatNum = markers[currentMarker].refRepeatNum;
                                keyValuePair.i2.numOfRepeats = (float)maxRepeatLength/(float)length(markers[currentMarker].motif);
                                keyValuePair.i2.ratioBf = 0.0;
                                keyValuePair.i2.ratioAf = 0.0;
                                keyValuePair.i2.locationShift = 0;
                                keyValuePair.i2.repSeq = infix(mateAtRepeat.seq, startCoordinate, endCoordinate);
                                keyValuePair.i2.purity = getPurity(markers[currentMarker].motif, keyValuePair.i2.repSeq);
                                keyValuePair.i2.ratioOver20In = findRatioOver20(keyValuePair.i2.repSeq);
                                myMap[keyValuePair.i1] = keyValuePair.i2;
                                ++currentMarker;
                                //If I've reached the end of the marker string I break this loop and take the next read
                                if (currentMarker > finalMarkerIdx)
                                    break;
                                motifLength = length(markers[currentMarker].motif);
                                numToLook = repeatNumbers[motifLength-1];
                                continue;
                            }
                            else
                            {
                                Dna5String revCompMotif = markers[currentMarker].motif;
                                reverseComplement(revCompMotif);
                                Pair<Pair<int>, float> coordinates = findPatternPrim(revCompMotif, mateAtRepeat, motifLength);
                                int startCoordinate = coordinates.i1.i1;
                                int endCoordinate = coordinates.i1.i2;
                                if (endCoordinate-startCoordinate+1 >= readLength-2*minFlank && endCoordinate > startCoordinate)
                                {
                                    Pair<Triple<CharString, CharString, int>,ReadInfo> keyValuePair;
                                    keyValuePair.i1 = Triple<CharString, CharString, int>(record.qName, markers[currentMarker].chrom, markers[currentMarker].STRstart);
                                    keyValuePair.i2.mateEditDist = getTagValue(record, "NM");
                                    keyValuePair.i2.STRend = markers[currentMarker].STRend;
                                    keyValuePair.i2.motif = markers[currentMarker].motif;
                                    keyValuePair.i2.refRepeatNum = markers[currentMarker].refRepeatNum;
                                    keyValuePair.i2.numOfRepeats = (float)maxRepeatLength/(float)length(markers[currentMarker].motif);
                                    keyValuePair.i2.ratioBf = 0.0;
                                    keyValuePair.i2.ratioAf = 0.0;
                                    keyValuePair.i2.locationShift = 0;
                                    keyValuePair.i2.repSeq = infix(mateAtRepeat.seq, startCoordinate, endCoordinate);
                                    reverseComplement(keyValuePair.i2.repSeq);
                                    keyValuePair.i2.purity = getPurity(markers[currentMarker].motif, keyValuePair.i2.repSeq);
                                    keyValuePair.i2.ratioOver20In = findRatioOver20(keyValuePair.i2.repSeq);
                                    myMap[keyValuePair.i1] = keyValuePair.i2;
                                    ++currentMarker;
                                    //If I've reached the end of the marker string I break this loop and take the next read
                                    if (currentMarker > finalMarkerIdx)
                                        break;
                                    motifLength = length(markers[currentMarker].motif);
                                    numToLook = repeatNumbers[motifLength-1];
                                    continue;
                                }
                            }
                        }
                        else
                        {
                            Pair<Triple<CharString, CharString, int>,ReadInfo> keyValuePair;
                            keyValuePair.i1 = Triple<CharString, CharString, int>(record.qName, markers[currentMarker].chrom, markers[currentMarker].STRstart);
                            keyValuePair.i2.mateEditDist = getTagValue(record, "NM");
                            keyValuePair.i2.STRend = markers[currentMarker].STRend;
                            keyValuePair.i2.motif = markers[currentMarker].motif;
                            keyValuePair.i2.refRepeatNum = markers[currentMarker].refRepeatNum;
                            keyValuePair.i2.numOfRepeats = (float)maxRepeatLength/(float)length(markers[currentMarker].motif);
                            keyValuePair.i2.ratioBf = 0.0;
                            keyValuePair.i2.ratioAf = 0.0;
                            keyValuePair.i2.locationShift = 0;
                            keyValuePair.i2.repSeq = std::string(readLength, 'N');
                            keyValuePair.i2.purity = 0.5;
                            keyValuePair.i2.ratioOver20In = findRatioOver20(record.seq);
                            myMap[keyValuePair.i1] = keyValuePair.i2;
                            ++currentMarker;
                            //If I've reached the end of the marker string I break this loop and take the next read
                            if (currentMarker > finalMarkerIdx)
                                break;
                            motifLength = length(markers[currentMarker].motif);
                            numToLook = repeatNumbers[motifLength-1];
                            continue;
                        }
                    }
                }
                ++currentMarker;
                if (currentMarker > finalMarkerIdx)
                    break;
                motifLength = length(markers[currentMarker].motif);
                numToLook = repeatNumbers[motifLength-1];
            }
        }
        //cout << "Skipped front alignment: " << noFrontAlign << " times and performed it " << frontAlign << " times.\n";
        //cout << "Skipped back alignment: " << noBackAlign << " times and performed it " << backAlign << " times.\n";
        noFrontAlign = frontAlign = noBackAlign = backAlign = 0;
        unordered_map<STRinfoSmall, Triple<std::set<float>,vector<float>, String<ReadPairInfo> > > finalMap; //Stores String of ReadPairInfo for each marker
        unordered_map<Triple<CharString, CharString, int>, ReadInfo>::const_iterator ite = myMap.end();
        for(unordered_map<Triple<CharString, CharString, int>, ReadInfo>::const_iterator it = myMap.begin(); it != ite; ++it)
        {
            //If this condition holds then only one member of the read pair has fulfilled the conditions and I can't use the pair.
            if ((it->second.numOfRepeats == 666) || (it->second.mateEditDist == 666))
            {
                //cout << "No: " << it->first.i1 << " "<< it->first.i2 << " " << it->first.i3 << "\n";
                continue;
            }
            //Create key
            //cout << "Yes: " << it->first.i1 << " "<< it->first.i2 << " " << it->first.i3 << "\n";
            STRinfoSmall currentSTR;
            currentSTR.chrom = it->first.i2;
            currentSTR.STRstart = it->first.i3;
            currentSTR.STRend = it->second.STRend;
            currentSTR.motif = it->second.motif;
            currentSTR.refRepeatNum = it->second.refRepeatNum;
            //Create value
            ReadPairInfo currentReadPair;
            currentReadPair.numOfRepeats = it->second.numOfRepeats;
            currentReadPair.ratioBf = it->second.ratioBf;
            currentReadPair.ratioAf = it->second.ratioAf;
            currentReadPair.locationShift = it->second.locationShift;
            currentReadPair.purity = it->second.purity;
            currentReadPair.ratioOver20In = it->second.ratioOver20In;
            currentReadPair.mateEditDist = it->second.mateEditDist;
            currentReadPair.repSeq = it->second.repSeq;
            //Put the lime in the coconut
            appendValue(finalMap[currentSTR].i3, currentReadPair);
            finalMap[currentSTR].i1.insert(currentReadPair.numOfRepeats);
            finalMap[currentSTR].i2.push_back(currentReadPair.numOfRepeats);
        }
        //Check whether any reads have been found
        if (finalMap.empty())
        {
            cout << "No reads found for supplied markers. Finished: " << PN_ID << "\n";
            continue;
        }
        //Set for storing allele-types and vector for storing reported alleles, count occurences in vector for all elements in set to get frequency of each allele
        std::set<float> presentAlleles;
        vector<float> allAlleles;
        int winnerFreq, secondFreq, currentFreq;
        float winner, second;
        //Loop over map of markers and look at all reads for each of them
        unordered_map<STRinfoSmall, Triple<std::set<float>,vector<float>, String<ReadPairInfo> > >::const_iterator ite2 = finalMap.end();
        unordered_map<STRinfoSmall, Triple<std::set<float>,vector<float>, String<ReadPairInfo> > >::iterator it = finalMap.begin();
        for(; it != ite2; ++it)
        {
            Pair<int> mapKey = Pair<int>(it->first.STRstart, it->first.STRend);
            if (startAndEndToStreamAndOffsets.count(mapKey)>0)
            {
                fflush(startAndEndToStreamAndOffsets[mapKey].i2);
                startAndEndToStreamAndOffsets[mapKey].i1[i] = ftell(startAndEndToStreamAndOffsets[mapKey].i2);
                fprintf(startAndEndToStreamAndOffsets[mapKey].i2, "%s\n", toCString(PN_ID));
            }
            else
            {
                cout << "Trying to write to a file that doesn't exist!" << endl;
                return 1;
            }
            String<ReadPairInfo> readPairs = it->second.i3;
            presentAlleles = it->second.i1;
            allAlleles = it->second.i2;
            winnerFreq = 0;
            secondFreq = 0;
            //Loop over set of alleles to consider and count occurences of each to determine initial labelling
            std::set<float>::iterator end = presentAlleles.end();
            for (std::set<float>::iterator allIt = presentAlleles.begin(); allIt!=end; ++allIt)
            {
                if (*allIt < 0)
                    continue;
                currentFreq = count(allAlleles.begin(), allAlleles.end(), *allIt);
                if ( currentFreq > winnerFreq)
                {
                    secondFreq = winnerFreq;
                    second = winner;
                    winnerFreq = currentFreq;
                    winner = *allIt;
                }
                else
                {
                    if (currentFreq > secondFreq)
                    {
                        secondFreq = currentFreq;
                        second = *allIt;
                    }
                    else
                    {
                        if(currentFreq == secondFreq)
                            second = max(second,*allIt);
                    }
                }
            }
            if (secondFreq < 0.10*winnerFreq)
                second = winner;
            //Write attributes and initial labelling to output file
            fprintf(startAndEndToStreamAndOffsets[mapKey].i2,"%s\t%u\t%u\t%s\t%.1f\t%u\t%.1f\t%.1f\n",toCString(it->first.chrom),it->first.STRstart,it->first.STRend,toCString(it->first.motif),it->first.refRepeatNum,length(readPairs),winner,second);
            for (unsigned i=0; i < length(readPairs); ++i)
            {
                ReadPairInfo printMe = readPairs[i];
                //cout << "Printing read number " << i << " which reports " << printMe.numOfRepeats << " repeats." << "\n";
                if (printMe.numOfRepeats < 0)
                    continue;
                //Print attributes to output file.
                fprintf(startAndEndToStreamAndOffsets[mapKey].i2,"%.1f\t%.1f\t%.1f\t%u\t%u\t%.2f\t%.2f\t%s\n",printMe.numOfRepeats,printMe.ratioBf,printMe.ratioAf,printMe.locationShift,printMe.mateEditDist,printMe.purity,printMe.ratioOver20In,toCString(printMe.repSeq));
            }
        }
        //Flush last stream before starting next PN
        //fflush(startAndEndToStreamAndOffsets[Pair<int>(it->first.STRstart, it->first.STRend)].i2);
        time_t pnEnd = time(0);
        cout << "Finished generation of output for " << PN_ID << " in " << pnEnd - pnStart << " seconds.\n";
        myMap.clear();
        finalMap.clear();
    }
    //END FOR LOOP OVER PNS HERE
    //Loop over marker map and write offset vector to the beginning of each file
    for (auto& marker: startAndEndToStreamAndOffsets)
    {
        //Find index of first pn with available reads
        unsigned idx = 0;
        while (idx < length(marker.second.i1) && marker.second.i1[idx] == 0)
            ++idx;
        //Check if any pn has available reads
        if (marker.second.i1[idx] == 0 || idx >= length(marker.second.i1))
        {
            rewind(marker.second.i2);
            fprintf(marker.second.i2, "%i ", -69 );
            continue;
        }
        rewind(marker.second.i2);
        for (unsigned i=0; i<length(marker.second.i1); ++i)
        {
            fprintf(marker.second.i2, "%u ", marker.second.i1[i]);
            fflush(marker.second.i2);
            //Chech if I am writing passed the reserved space at front
            if (ftell(marker.second.i2) > marker.second.i1[idx])
            {
                cerr << "writing pn-offsets over pnData @: " << marker.first.i1 << endl;
                cerr << "Pn idx: " << i << " and offset: " << marker.second.i1[i] << endl;
                return 1;
            }
        }
    }
    time_t end = time(0);
    cout << "Total time: " << end - begin << "\n";
    return 0;
}