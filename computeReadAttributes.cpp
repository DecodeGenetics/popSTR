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

using namespace std;
using namespace seqan;

//Structure to store marker information
struct STRinfo {   
    CharString chrom;
    int STRstart;
    int STRend;
    Dna5String motif;
    float refRepeatNum; //Number of repeats in reference sequence 
    Dna5String refBf;
    Dna5String refAf;
    Dna5String refRepSeq;
    double refRepPurity;
} ;

//Structure to store marker information
struct STRinfoSmall {   
    CharString chrom;
    int STRstart;
    int STRend;
    Dna5String motif;
    float refRepeatNum; //Number of repeats in reference sequence 
    Dna5String refRepSeq;
} ;

//So I can map from STRinfoSmall in finalMap
bool operator<(const STRinfoSmall & Left, const STRinfoSmall & Right)
{
    return Left.STRstart < Right.STRstart;
}

//structure to store read information
struct ReadInfo {
    int STRend;
    Dna5String motif;
    float refRepeatNum; //Number of repeats in reference sequence 
    float numOfRepeats;
    float ratioBf; 
    float ratioAf; 
    unsigned locationShift; 
    float purity; 
    float ratioOver20In;
    float ratioOver20After;
    unsigned sequenceLength;
    unsigned mateEditDist;
    bool wasUnaligned;
    Dna5String repSeq; //Repeat sequence in read
    Dna5String refRepSeq; //Repeat sequence in reference
} ; 

//structure to store read information
struct ReadPairInfo {
    float numOfRepeats;
    float ratioBf; 
    float ratioAf; 
    unsigned locationShift; 
    float purity; 
    float ratioOver20In;
    float ratioOver20After;
    unsigned sequenceLength;
    unsigned mateEditDist;
    bool wasUnaligned;
    Dna5String repSeq; //Repeat sequence in read
} ; 

//Map to check how many repeats I need to find in a read w.r.t. motif length
map<int,int> repeatNumbers;

//For repeating a motif n times
Dna5String repeat(Dna5String s, int n) {
    Dna5String ret;
    for (int i = 0; i < n; i++) {
        append(ret,s);
    }
    return ret;
}

//Create all one error permutations of motif passed 
std::set<Dna5String> createPermutations(Dna5String motif)
{
    std::set<Dna5String> permutations;
    permutations.insert(motif);
    Dna5String bases = "ATCGN";
    for (unsigned i = 0; i < length(motif); ++i)
    {
        Dna5String motifCopy = motif; 
        for(unsigned j = 0; j < length(bases); ++j)
        {
            motifCopy[i] = bases[j];
            if(motifCopy[i] != motif[i])
                permutations.insert(motifCopy);
        }
        erase(motifCopy,i);
        permutations.insert(motifCopy);  
    }
    return permutations;
}

//Primitive pattern search, just searches for first and last occurence of #repeats*motif without considering what's between them
Pair<Pair<int>, float> findPatternPrim(Dna5String pattern, Dna5String readSequence, int motifLength)
{
    unsigned patternLength = length(pattern);
    std::set<Dna5String> permutations = createPermutations(pattern);
    int startCoordinate = length(readSequence);
    int endCoordinate = 0;
    unsigned index = 0;
    while(index+patternLength <= length(readSequence))
    {
        Dna5String theSubString = infixWithLength(readSequence, index, patternLength);
        Dna5String theSubStringMini = infixWithLength(readSequence, index, patternLength-1);
        if (permutations.find(theSubString) != permutations.end())
        {
            if (startCoordinate == length(readSequence))
                startCoordinate = index;
            endCoordinate = index + patternLength;
            index = index + motifLength;      
        }
        else
        {
            if (permutations.find(theSubStringMini) != permutations.end())
            {
                if (startCoordinate == length(readSequence))
                    startCoordinate = index;
                endCoordinate = index + patternLength - 1;
                if (prefix(pattern, motifLength) == prefix(theSubStringMini, motifLength))
                    index = index + motifLength;
                else 
                    index = index + motifLength - 1;
            }
            else
                index = index + 1;
        }
    }
    float numOfRepeats = (float)(endCoordinate-startCoordinate)/(float)motifLength;
    return Pair<Pair<int>, int>(Pair<int>(startCoordinate,endCoordinate), numOfRepeats);
}

//Check if read sequence contains motif and return start and end coordinates and number of repeats
Pair<Pair<int>, int> findPattern(Dna5String pattern, Dna5String readSequence, int motifLength)
{
    unsigned patternLength = length(pattern);
    std::set<Dna5String> permutations = createPermutations(pattern);
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
    unsigned numOfRepeats = moves + repeatNumbers[motifLength] - 1;
    if (endCoordinate - startCoordinate > currentEnd - currentStart)
        return Pair<Pair<int>, int>(Pair<int>(startCoordinate,endCoordinate), numOfRepeats);
    else
        return Pair<Pair<int>, int>(Pair<int>(currentStart,currentEnd), numOfRepeats);
}

//Compute the actual number of repeats/expected number of repeats based on length
float getPurity(Dna5String motif, Dna5String STRsequence)
{
    unsigned motifLength = length(motif);
    float expectReps = (float)length(STRsequence)/(float)motifLength;
    unsigned result = 0;
    unsigned index = 0;
    while(index < length(STRsequence))
    {
        Dna5String theSubString = infixWithLength(STRsequence, index, motifLength);
        if (theSubString == motif)
        {
            result++;
            index = index + motifLength;
        }
        else
            index = index + 1;    
    }
    return (float)result/(float)expectReps;
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

//Computes all sorts of quality indicators for the given read w.r.t. the given microsatellite 
Pair<Pair<Pair<CharString>,int>,ReadInfo> computeReadInfo(BamAlignmentRecord& record, STRinfo& markerInfo, Pair<Pair<int>, float> coordinates, int minFlank, int maxRepeatLength)
{
    //Type definition for alignment structure
    typedef Dna5String TSequence;              
    typedef Align<TSequence,ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    
    //Create alignment structures, before and after
    TAlign alignBefore;
    TAlign alignAfter;
    resize(rows(alignBefore), 2);
    resize(rows(alignAfter), 2);
    assignSource(row(alignBefore,0),markerInfo.refBf);
    assignSource(row(alignAfter,0),markerInfo.refAf); 
    
    //Variables for sequence-parts
    Dna5String before, repeatRegion, after;
    CharString qualString = record.qual;
    ReadInfo mapValue;
    int oldStartCoord = coordinates.i1.i1;
    
    //Create index to map and store sequence length, STRend and motif
    Pair<CharString> readNameAndChrom = Pair<CharString>(record.qName,markerInfo.chrom);
    Pair<Pair<CharString>,int> mapKey = Pair<Pair<CharString>,int>(readNameAndChrom, markerInfo.STRstart); 
    
    //Insert values into mapValue in returnPair
    mapValue.sequenceLength = length(record.seq);
    mapValue.STRend = markerInfo.STRend;
    mapValue.motif = markerInfo.motif;
    mapValue.refRepeatNum = markerInfo.refRepeatNum;
    mapValue.refRepSeq = markerInfo.refRepSeq;
    
    //Split read into 2 parts, both containing the alleged repeat
    before = prefix(record.seq, coordinates.i1.i2+1);
    after = suffix(record.seq, coordinates.i1.i1);
    
    //Align part of read coming before repeat to reference before repeat and rest to reference after repeat
    assignSource(row(alignBefore,1),before);
    assignSource(row(alignAfter,1),after);
    int scoreBf = globalAlignment(alignBefore, Score<int,Simple>(1,-5,-5), AlignConfig<true, false, true, false>());
    int scoreAf = globalAlignment(alignAfter, Score<int,Simple>(1,-5,-5), AlignConfig<false, true, false, true>());
    
    int viewPosition = toViewPosition(row(alignBefore,0),length(source(row(alignBefore,0)))-1);
    int startCoord = toSourcePosition(row(alignBefore,1),viewPosition)+1;
    if (startCoord == 1 && isGap(row(alignBefore,1),viewPosition))
        startCoord -= 1;
    int endCoord = toSourcePosition(row(alignAfter,1),toViewPosition(row(alignAfter,0),0))-1;
    int flankSum = startCoord + length(source(row(alignAfter,1)))-endCoord - 1;
    int leftFlank = startCoord;
    int rightFlank = length(source(row(alignAfter,1)))-endCoord - 1;
    /*if (record.qName == "ST-E00133:149:H3GKMCCXX:5:2218:19077:41497")
    {
        cout << "FlankSum: " << flankSum << " leftFlank: " << leftFlank << " rightFlank: " << rightFlank << endl;
        if (leftFlank > 0)
            cout << "AlignmentScore before: " << (float)scoreBf/(float)startCoord << endl;
        if (rightFlank > 0)
            cout << "AlignmentScore after: " << (float)scoreAf/(float)(length(after)-endCoord-1) << endl;
        cout << "Start coord: " << startCoord << endl;
        cout << "Old start coord: " << oldStartCoord << endl;
        cout << "End coord: " << endCoord << endl;
        cout << "Infix command is: infix(" << startCoord << "," << oldStartCoord+endCoord+1 << ")" << endl;
        cout << "Repeat purity: " << getPurity(markerInfo.motif,infix(record.seq, startCoord, oldStartCoord+endCoord+1)) << endl;
        cout << "Reference repeat purity: " << getPurity(markerInfo.motif, markerInfo.refRepSeq) << endl;
    }*/
    double refRepPurity = markerInfo.refRepPurity;
    bool startOk = false;
    bool endOk = false;
    bool purityOk = false;
    if (startCoord >= oldStartCoord + endCoord || ((oldStartCoord+endCoord+1)-startCoord < maxRepeatLength && (leftFlank < 4 || rightFlank < 4)))
    {
        //cout << "Start coordinate is larger than end coordinate, can't use this read." << endl;
        coordinates.i1.i1 = startCoord;
        coordinates.i1.i2 = endCoord;
        mapValue.ratioBf = 0;
        mapValue.ratioAf = 0;
        mapValue.numOfRepeats = 666;
        //Debugging code
        //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl; 
        repeatRegion = "";
        mapValue.ratioOver20In = 0;
        mapValue.ratioOver20After = 0;
        mapValue.purity = 0;
        mapValue.repSeq = repeatRegion;
        mapValue.locationShift = 100;
        return Pair<Pair<Pair<CharString>,int>,ReadInfo>(mapKey,mapValue);
    }
    double readPurity = getPurity(markerInfo.motif,infix(record.seq, startCoord, oldStartCoord+endCoord+1));
    //Check if a minimum number of repeats have been found and whether the flanking area is sufficient for both start and end coordinates
    if ((length(source(row(alignBefore,1)))-(startCoord) >=length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)]) && (startCoord>=minFlank))
        startOk = true;
    if ((endCoord >= length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)]-1)&&(length(source(row(alignAfter,1)))-endCoord > minFlank))
        endOk = true;
    if (readPurity>(0.75*refRepPurity))
        purityOk = true;
    //I allow only 4 aligning bases on either side if I have more than 2*minFlank aligned bases in total.
    if (flankSum >= 2*minFlank && leftFlank >= 4 && rightFlank >= 4)
    {
        startOk = true;
        endOk = true;
    } 
    /*if (record.qName == "ST-E00133:149:H3GKMCCXX:5:2218:19077:41497")
    {
        cout << "Start ok: " << startOk << endl;
        cout << "End ok: " << endOk << endl;
        cout << "Purity ok: " << purityOk << endl;
    }*/
    if (leftFlank < 4 && flankSum>=2*minFlank && (float)scoreAf/(float)(length(after)-endCoord)>0.7 && readPurity>(0.8*refRepPurity) && (oldStartCoord+endCoord+1)-startCoord >= maxRepeatLength)
    {
        //cout << "Am making greater than allele on left end." << endl;
        coordinates.i1.i1 = startCoord;
        coordinates.i1.i2 = endCoord;
        mapValue.ratioBf = 0;
        mapValue.ratioAf = (float)scoreAf/(float)(length(after)-endCoord-1);
        //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl; 
        repeatRegion = infix(record.seq, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1);
        mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
        mapValue.ratioOver20In = findRatioOver20(infix(qualString, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1));
        mapValue.ratioOver20After = findRatioOver20(suffix(suffix(qualString, oldStartCoord),coordinates.i1.i2+1));
        mapValue.purity = getPurity(markerInfo.motif,repeatRegion);   
    }
    else
    {
        if (rightFlank < 4 && flankSum>=2*minFlank && (float)scoreBf/(float)startCoord>0.7 && readPurity>(0.8*refRepPurity) && (oldStartCoord+endCoord+1)-startCoord >= maxRepeatLength)
        {
            //cout << "Am making greater than allele on right end." << endl;
            coordinates.i1.i1 = startCoord;
            coordinates.i1.i2 = endCoord;
            mapValue.ratioBf = (float)scoreBf/((float)startCoord+1);
            mapValue.ratioAf = 0;
            //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl; 
            repeatRegion = infix(record.seq, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1);
            mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
            mapValue.ratioOver20In = findRatioOver20(infix(qualString, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1));
            mapValue.ratioOver20After = 0;
            mapValue.purity = getPurity(markerInfo.motif,repeatRegion);
        }
        else
        {
            if (rightFlank < 4 && leftFlank < 4 && getPurity(markerInfo.motif,record.seq)>(0.85*refRepPurity) && (oldStartCoord+endCoord+1)-startCoord >= maxRepeatLength)
            {
                //cout << "Am making a SUPER allele." << endl;
                coordinates.i1.i1 = 0;
                coordinates.i1.i2 = length(record.seq)-1;
                mapValue.ratioBf = 0.0;
                mapValue.ratioAf = 0.0;
                //cout << "Infix command is: infix(" << 0 << "," << length(record.seq)-1 << ")" << endl;
                repeatRegion = record.seq;
                mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
                mapValue.ratioOver20In = findRatioOver20(record.seq);
                mapValue.ratioOver20After = 0.0;
                mapValue.purity = getPurity(markerInfo.motif,repeatRegion);
            }
            else
            {            
                //If both coordinates are ok and the distance between them is > motifLength * min#ofMotifs then the read is useful. 
                if (startOk && endOk && purityOk && startCoord < oldStartCoord + endCoord)
                {
                    coordinates.i1.i1 = startCoord;
                    coordinates.i1.i2 = endCoord;
                    mapValue.ratioBf = (float)scoreBf/((float)startCoord+1);
                    mapValue.ratioAf = (float)scoreAf/(float)(length(after)-endCoord-1);
                    //Debugging code
                    //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl; 
                    repeatRegion = infix(record.seq, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1);
                    mapValue.numOfRepeats = (float)length(repeatRegion)/(float)length(markerInfo.motif);
                    if (length(repeatRegion)>=maxRepeatLength)
                        mapValue.numOfRepeats = (float)maxRepeatLength/(float)length(markerInfo.motif);
                    //If I find less repeats than the required minimum (depending on the motif length) then I don't use the read
                    if (ceil(mapValue.numOfRepeats) < repeatNumbers[length(markerInfo.motif)])
                    {
                        //cout << "Not enough repeats." << endl;
                        mapValue.numOfRepeats = 666;
                    }
                    mapValue.ratioOver20In = findRatioOver20(infix(qualString, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1));
                    mapValue.ratioOver20After = findRatioOver20(suffix(suffix(qualString, oldStartCoord),coordinates.i1.i2+1));
                    mapValue.purity = getPurity(markerInfo.motif,repeatRegion);
                }
                //Otherwise I can't use the read so I set numOfRepeats to 666
                else
                {
                    //cout << "Setting num of repeats to 666 in the last else statement." << endl;
                    coordinates.i1.i1 = startCoord;
                    coordinates.i1.i2 = endCoord;
                    mapValue.ratioBf = 0;
                    mapValue.ratioAf = 0;
                    mapValue.numOfRepeats = 666;
                    //Debugging code
                    //cout << "Infix command is: infix(" << coordinates.i1.i1 << "," << oldStartCoord+coordinates.i1.i2+1 << ")" << endl; 
                    repeatRegion = infix(record.seq, coordinates.i1.i1, oldStartCoord+coordinates.i1.i2+1);
                    mapValue.ratioOver20In = 0;
                    mapValue.ratioOver20After = 0;
                    mapValue.purity = 0;
                }
            }
        }
    } 
    mapValue.repSeq = repeatRegion;
    
    //Debugging code
    /*cout << "Motif: " << mapValue.motif << " Start: " << markerInfo.STRstart << endl;
    cout << "Before: " << prefix(before, coordinates.i1.i1) << endl;
    cout << "Aligned to: " << suffix(source(row(alignBefore,0)),toSourcePosition(row(alignBefore,0),toViewPosition(row(alignBefore,1),0))) << endl;
    cout << "Repeat: " << repeatRegion << endl;
    cout << "After: " << suffix(after, coordinates.i1.i2+1) << endl;
    cout << "Aligned to: " << prefix(source(row(alignAfter,0)),toSourcePosition(row(alignAfter,0),toViewPosition(row(alignAfter,1),length(source(row(alignAfter,1)))))) << endl;
    cout << "Number of repeats: " << mapValue.numOfRepeats << endl;
    if (mapValue.numOfRepeats == 666)
    {
        if (!startOk)
        {
            cout << "Enough repeats according to start? " << ((length(source(row(alignBefore,1)))-(startCoord+1)) >= (length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)])) << endl;
            cout << "Enough flanking area in front? " << (startCoord>=minFlank) << endl;
        }
        if (!endOk)
        {
            cout << "Enough repeats according to end? " << (endCoord >= length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)]-1) << endl; 
            cout << "Enough flanking area behind ?" << (length(source(row(alignAfter,1)))-endCoord >= minFlank) << endl;
        }
        cout << "Start and end not overlapping? " << ((startCoord + length(markerInfo.motif)*repeatNumbers[length(markerInfo.motif)])< (startCoord + endCoord)) << endl; 
    }*/
    
    TRow &row2B = row(alignBefore,1);
    
    //Check location shift
    mapValue.locationShift = abs(record.beginPos - (markerInfo.STRstart - 1000 + toViewPosition(row2B, 0)));
    return Pair<Pair<Pair<CharString>,int>,ReadInfo>(mapKey,mapValue);
}

bool qualityClipBegin(BamAlignmentRecord& record, int windowSize)
{        
    int qualSum = 0;
    for (unsigned i = 0; i<windowSize; ++i)
        qualSum += (record.qual[i]-33);
    int averageQual = round((double)qualSum/(double)windowSize);
    if (averageQual>=25)
        return true;
    int index = 0;
    while (averageQual<25 && index+windowSize < length(record.qual))
    {
        qualSum -= (record.qual[index]-33);
        qualSum += (record.qual[index+windowSize]-33);
        averageQual = round((double)qualSum/(double)windowSize);
        ++index;
    }
    erase(record.seq, 0, index-1);
    erase(record.qual, 0, index-1);
    record.beginPos += (index - 1);
    if (length(record.seq)>=30)
        return true;
    else 
        return false;
}

bool qualityClipEnd(BamAlignmentRecord& record, int windowSize)
{        
    int qualSum = 0;
    for (unsigned i = length(record.qual)-1; i>=length(record.qual)-windowSize; --i)
        qualSum += (record.qual[i]-33);
    int averageQual = round((double)qualSum/(double)windowSize);
    if (averageQual>=25)
        return true;
    int index = 0;
    while (averageQual<25 && index < length(record.qual))
    {
        ++index;
        qualSum -= (record.qual[length(record.qual)-index]-33);
        qualSum += (record.qual[length(record.qual)-windowSize-index]-33);
        averageQual = round((double)qualSum/(double)windowSize);        
    }
    erase(record.seq, length(record.seq)-index, length(record.seq)-1);
    erase(record.qual,length(record.qual)-index, length(record.qual)-1);
    if (length(record.seq)>=30)
        return true;
    else 
        return false;
}

int main(int argc, char const ** argv)
{          
    //Check arguments.
    if (argc != 8)
    {
        cerr << "USAGE: " << argv[0] << " IN.bam IN.bam.bai outputDirectory markerInfoFile minFlankLength maxRepeatLength PN-id\n";
        return 1;
    }
    
    //Save PN-id
    int maxRepeatLength = lexicalCast<int>(argv[6]);
    CharString PN_ID = argv[7];
    int windowSize = 5;
    
    //min-flanking area
    int minFlank = lexicalCast<int>(argv[5]);
    
    //To store marker info
    String<STRinfo> markers;
    STRinfo currInfo;
    //Read all markers into memory
    ifstream markerFile(argv[4]);
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
        markerFile >> currInfo.refRepeatNum; 
        string refBfString;
        markerFile >> refBfString;
        currInfo.refBf = refBfString;
        string refAfString;
        markerFile >> refAfString;
        currInfo.refAf = refAfString;
        string refRepSeq;
        markerFile >> refRepSeq;
        currInfo.refRepSeq = refRepSeq;
        currInfo.refRepeatNum = (float)refRepSeq.length()/(float)motifString.length();
        currInfo.refRepPurity = getPurity(currInfo.motif, currInfo.refRepSeq);
        if (currInfo.STRend - currInfo.STRstart < 151 - 2*minFlank)
            appendValue(markers, currInfo);
    }
    cout << "Finished reading marker Info, number of markers: " << length(markers) << endl; 
    
    //Create output stream
    CharString attributeDirectory = argv[3];
    append(attributeDirectory, "/attributes/");
    struct stat st;
    if(stat(toCString(attributeDirectory),&st) != 0)
        mkdir(toCString(attributeDirectory),0777);
    if (length(currInfo.chrom) > 2)    
        append(attributeDirectory, currInfo.chrom);
    else
    {
        append(attributeDirectory, "chr");
        append(attributeDirectory, currInfo.chrom);
    }
    struct stat st3;
    if(stat(toCString(attributeDirectory),&st3) != 0)
        mkdir(toCString(attributeDirectory),0777);
    append(attributeDirectory, "/");
	append(attributeDirectory, PN_ID);
    ofstream outputFile(toCString(attributeDirectory));
    
    //Set up how many repeats I require for each motif length
    repeatNumbers[2]=4;
    repeatNumbers[3]=3;
    repeatNumbers[4]=2;
    repeatNumbers[5]=2;
    repeatNumbers[6]=2;
    
    // Setup name store, cache, and BAM I/O context.
    typedef StringSet<CharString> TNameStore;
    typedef NameStoreCache<TNameStore>   TNameStoreCache;
    typedef BamIOContext<TNameStore>     TBamIOContext;
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);
    
    // Open BAM Stream for reading.
    Stream<Bgzf> inStream;
    if (!open(inStream, argv[1], "r"))
    {
        cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }
    
    // Read header.
    BamHeader header;
    if (readRecord(header, context, inStream, Bam()) != 0)
    {
        cerr << "ERROR: Could not read header from BAM file " << argv[1] << "\n";
        return 1;
    }

    // Read BAI index.
    BamIndex<Bai> baiIndex;
    if (read(baiIndex, argv[2]) != 0)
    {
        cerr << "ERROR: Could not read BAI index file " << argv[2] << "\n";
        return 1;
    }
    
    //Variables for the start and end coordinates of reads and their mates
    int bamStart, bamEnd, mateStart, mateEnd;
    
    //Get rID of chromosome
    int rID = 0;
    if (!getIdByName(nameStore, markers[0].chrom, rID, nameStoreCache))
    {
        std::cerr << "ERROR: Reference sequence named " << markers[0].chrom << " not known.\n";
        return 1;
    }
    
    //Jump to beginning of chromosome, put 100000000 as end to make sure I find alignments.
    bool hasAlignments = false;
    int jumpStart = std::max(0,markers[0].STRstart - 2000);
    int jumpEnd = jumpStart + 100000000;
    if (!jumpToRegion(inStream, hasAlignments, context, rID, jumpStart, jumpEnd, baiIndex))
    {
        cerr << "ERROR: Could not jump to " << currInfo.chrom << ":" << 0 << "\n";
        return 1;
    }
    if (!hasAlignments)
    {
        //no alignments in interval
        cout << "No alignments found!" << endl;
        return 0;
    } 
    
    //Map from read name, marker chromosome and marker start to info on read-pair with that read name
    map<Pair<Pair<CharString>, int>, ReadInfo> myMap;
    //Index into string storing marker information
    unsigned markerIndex = 0;
    BamAlignmentRecord record;
    unsigned numToLook;
    while (!atEnd(inStream))
    {
        if (readRecord(record, context, inStream, Bam()) != 0)
        {
            cerr << "ERROR: Could not read record from BAM file.\n";
            return 1;
        }
        //If the read is a duplicate or doesn't pass the quality check I move on
        if (hasFlagQCNoPass(record) || hasFlagDuplicate(record))
            continue;
        //Clip away low quality ends and return false if the read becomes too short after this
        //if (!qualityClipBegin(record, windowSize) || !qualityClipEnd(record, windowSize))
        //    continue;
        bamStart = record.beginPos;
        bamEnd = bamStart + length(record.seq);
        mateStart = record.pNext;
        mateEnd = mateStart + length(record.seq);
        //Have passed the current marker? -> update markerIndex   
        while (bamStart > markers[markerIndex].STRend + 1000)
        {
            ++markerIndex;
            //If the markerIndex has exceeded the length of the marker string I can't check the while condition (markers[markerIndex].STRend doesn't exist) so I break
            if (markerIndex > length(markers)-1)
                break;
        }
        //If the markerIndex has exceeded the length of the marker string I break the BAM-reading loop
        if(markerIndex > length(markers)-1)
            break;
        int currentMarker = markerIndex;
        numToLook = repeatNumbers[length(markers[currentMarker].motif)];
        // At end or at next chromosome? -> done
        if (record.rID == -1 || record.rID > rID)
            break;
        //Loop for comparing current read to all possible markers, could be useful for more than one marker 
        while (true)
        {
            // Have not reached beginning of the interval? -> next read
            if (bamEnd < max(markers[currentMarker].STRstart-1000,0))
                break;
            //If the read is not long enough for the current marker I need to check the next marker
            if (markers[currentMarker].STRend-markers[currentMarker].STRstart > length(record.seq)-2*minFlank)
            {
                ++currentMarker;
                if (currentMarker > length(markers)-1)
                    break;
                numToLook = repeatNumbers[length(markers[currentMarker].motif)];
                continue;
            }
            // Does the read intersect the microsatellite without being unaligned? -> Then I look for the current repeat motif   
            if (((bamStart >= markers[currentMarker].STRstart && bamStart <= markers[currentMarker].STRend)||(bamStart < markers[currentMarker].STRstart && bamEnd >= markers[currentMarker].STRstart)) && !hasFlagUnmapped(record))
            {
                Pair<Pair<int>, float> coordinates = findPatternPrim(repeat(markers[currentMarker].motif,numToLook), record.seq, length(markers[currentMarker].motif));
                int startCoordinate = coordinates.i1.i1;
                int endCoordinate = coordinates.i1.i2;
                //Does the read contain all of the microsatellite and flanking regions larger than the set minimum? -> Then I compute read attributes
                if ((startCoordinate != length(record.seq)) && (endCoordinate > startCoordinate))
                {                
                    Pair<Pair<Pair<CharString>,int>,ReadInfo> keyValuePair = computeReadInfo(record, markers[currentMarker], coordinates, minFlank, maxRepeatLength);
                    keyValuePair.i2.wasUnaligned = false; 
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
                }   
                ++currentMarker;
                //If I've reached the end of the marker string I break this loop and take the next read
                if (currentMarker > length(markers)-1)
                    break; 
                numToLook = repeatNumbers[length(markers[currentMarker].motif)];
                continue;
            }
            // Does the reads mate intersect the microsatellite and the read is not unaligned? -> Then I want the reads edit distance!
            if (((mateStart >= markers[currentMarker].STRstart && mateStart <= markers[currentMarker].STRend)||(mateStart < markers[currentMarker].STRstart && mateEnd >= markers[currentMarker].STRstart)) && !hasFlagUnmapped(record))
            {
                Pair<CharString> readNameAndChrom = Pair<CharString>(record.qName,markers[currentMarker].chrom);
                Pair<Pair<CharString>,int> mapKey = Pair<Pair<CharString>,int>(readNameAndChrom, markers[currentMarker].STRstart); 
                if (myMap.count(mapKey) == 0)
                    myMap[mapKey].numOfRepeats = 666;
                myMap[mapKey].mateEditDist = getTagValue(record, "NM");
                ++currentMarker;
                //If I've reached the end of the marker string I break this loop and take the next read
                if (currentMarker > length(markers)-1)
                    break;
                numToLook = repeatNumbers[length(markers[currentMarker].motif)];
                continue;
            }
            //Is the read aligned to either 1000 bp before or after the microsatellite-> Then I check it out!
            if ((bamEnd >= max(markers[currentMarker].STRstart-1000,0) && bamEnd <= markers[currentMarker].STRstart) ||(bamStart >= markers[currentMarker].STRend && bamStart <= markers[currentMarker].STRend+1000))
            {
                //Is the reads mate unaligned? -> Then I want the reads edit distance!
                if (hasFlagNextUnmapped(record))
                {
                    Pair<CharString> readNameAndChrom = Pair<CharString>(record.qName,markers[currentMarker].chrom);
                    Pair<Pair<CharString>,int> mapKey = Pair<Pair<CharString>,int>(readNameAndChrom, markers[currentMarker].STRstart);
                    if (myMap.count(mapKey) == 0)
                        myMap[mapKey].numOfRepeats = 666; 
                    myMap[mapKey].mateEditDist = getTagValue(record, "NM");
                    ++currentMarker;
                    //If I've reached the end of the marker string I break this loop and take the next read
                    if (currentMarker > length(markers)-1)
                        break;
                    numToLook = repeatNumbers[length(markers[currentMarker].motif)];
                    continue;       
                }
                //Or, is the read itself unaligned -> Then I look for the current repeat motif
                if (hasFlagUnmapped(record))
                {
                    Pair<Pair<int>, float> coordinates = findPatternPrim(repeat(markers[currentMarker].motif,numToLook), record.seq, length(markers[currentMarker].motif));
                    int startCoordinate = coordinates.i1.i1;
                    int endCoordinate = coordinates.i1.i2;
                    //Does the read contain all of the microsatellite and flanking regions larger than the set minimum? -> Then I compute read attributes
                    if ((startCoordinate != length(record.seq)) && (endCoordinate > startCoordinate))
                    {
                        Pair<Pair<Pair<CharString>,int>,ReadInfo> keyValuePair = computeReadInfo(record, markers[currentMarker], coordinates, minFlank, maxRepeatLength);
                        keyValuePair.i2.wasUnaligned = true;
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
                    }
                    ++currentMarker;
                    //If I've reached the end of the marker string I break this loop and take the next read
                    if (currentMarker > length(markers)-1)
                        break;
                    numToLook = repeatNumbers[length(markers[currentMarker].motif)];
                    continue;
                }
            }
            ++currentMarker;  
            if (currentMarker > length(markers)-1)
                break;
            numToLook = repeatNumbers[length(markers[currentMarker].motif)];     
        }
    }
    map<STRinfoSmall, Pair<Pair<std::set<float>,vector<float> >,String<ReadPairInfo> > > finalMap; //Stores String of ReadPairInfo for each marker
    map<Pair<Pair<CharString>, int>, ReadInfo>::const_iterator ite = myMap.end();
    for(map<Pair<Pair<CharString>, int>, ReadInfo>::const_iterator it = myMap.begin(); it != ite; ++it)
    {
        //If this condition holds then only one member of the read pair has fulfilled the conditions and I can't use the pair.
        if ((it->second.numOfRepeats == 666) || (it->second.mateEditDist == 666))
            continue;
        //Create key     
        STRinfoSmall currentSTR; 
        currentSTR.chrom = it->first.i1.i2;
        currentSTR.STRstart = it->first.i2;
        currentSTR.STRend = it->second.STRend;
        currentSTR.motif = it->second.motif;
        currentSTR.refRepeatNum = it->second.refRepeatNum;
        currentSTR.refRepSeq = it->second.refRepSeq;
        //Create value 
        ReadPairInfo currentReadPair;
        currentReadPair.numOfRepeats = it->second.numOfRepeats;
        currentReadPair.ratioBf = it->second.ratioBf;
        currentReadPair.ratioAf = it->second.ratioAf;
        currentReadPair.locationShift = it->second.locationShift;
        currentReadPair.purity = it->second.purity;
        currentReadPair.ratioOver20In = it->second.ratioOver20In;
        currentReadPair.ratioOver20After = it->second.ratioOver20After;
        currentReadPair.sequenceLength = it->second.sequenceLength;
        currentReadPair.mateEditDist = it->second.mateEditDist;
        currentReadPair.repSeq = it->second.repSeq;
        currentReadPair.wasUnaligned = it->second.wasUnaligned;
        //Put the lime in the coconut 
        appendValue(finalMap[currentSTR].i2, currentReadPair);
        finalMap[currentSTR].i1.i1.insert(currentReadPair.numOfRepeats);
        finalMap[currentSTR].i1.i2.push_back(currentReadPair.numOfRepeats);
    }
    
    //I write PN-id and check whether any reads have been found, if not I exit.
    outputFile << PN_ID << endl;
    if (finalMap.empty())
    {
        cout << "Finished: " << PN_ID << endl;
        return 0;
    }

    //Set for storing allele-types and vector for storing reported alleles, count occurences in vector for all elements in set to get frequency of each allele
    std::set<float> presentAlleles; 
    vector<float> allAlleles;
    vector<float> biggerAlleles;
    int winnerFreq, secondFreq, currentFreq;
    float winner, second;
    //Loop over map of markers and look at all reads for each of them
    map<STRinfoSmall, Pair<Pair<std::set<float>,vector<float> >,String<ReadPairInfo> > >::const_iterator ite2 = finalMap.end();          
    for(map<STRinfoSmall, Pair<Pair<std::set<float>,vector<float> >,String<ReadPairInfo> > >::iterator it = finalMap.begin(); it != ite2; ++it)
    {        
        String<ReadPairInfo> readPairs = it->second.i2;        
        presentAlleles = it->second.i1.i1;                
        allAlleles = it->second.i1.i2;        
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
        if (secondFreq < 0.15*winnerFreq)
            second = winner;
        //Write attributes and initial labelling to output file
        outputFile << it->first.chrom << "\t" << it->first.STRstart << "\t" << it->first.STRend << "\t" << it->first.motif << "\t" << setprecision(1) << fixed << it->first.refRepeatNum << "\t" << length(readPairs) << "\t" << it->first.refRepSeq << "\t" << setprecision(1) << fixed << winner << "\t" << setprecision(1) << fixed << second << endl;
        for (unsigned i=0; i < length(readPairs); ++i)
        {
            ReadPairInfo printMe = readPairs[i];
            //cout << "Printing read number " << i << " which reports " << printMe.numOfRepeats << " repeats." << endl;
            if (printMe.numOfRepeats < 0)
                continue;
            //Print attributes to output file.
            outputFile << setprecision(1) << fixed << printMe.numOfRepeats << flush << "\t" << setprecision(2) << fixed << printMe.ratioBf << "\t" << setprecision(2) << fixed << printMe.ratioAf << "\t" << printMe.locationShift << "\t" << printMe.mateEditDist << "\t" << setprecision(2) << fixed << std::min((float)1.00, printMe.purity) << "\t" << setprecision(2) << fixed << printMe.ratioOver20In << "\t" << setprecision(2) << fixed << printMe.ratioOver20After << "\t" << printMe.sequenceLength << "\t" << printMe.wasUnaligned << "\t" << printMe.repSeq << endl; 
        }        
    }
    cout << "Finished: " << PN_ID << endl;
    return 0;
}
