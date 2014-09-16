#include <iostream>
#include <set>
#include <map>
#include <string>
#include <ctime>
#include <math.h> 
#include <sstream>
#include <sys/stat.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/vcf_io.h>
#include <seqan/sequence.h>
#include "/home/snaedisk/work/liblinear/linear.h"
#include "/home/snaedisk/work/liblinear/linear.cpp"
#include "/home/snaedisk/work/liblinear/tron.h"
#include "/home/snaedisk/work/liblinear/tron.cpp"
#include "/home/snaedisk/work/liblinear/blas/blas.h"
#include "/home/snaedisk/work/liblinear/blas/blasp.h"
#include "/home/snaedisk/work/liblinear/blas/daxpy.c"
#include "/home/snaedisk/work/liblinear/blas/ddot.c"
#include "/home/snaedisk/work/liblinear/blas/dnrm2.c"
#include "/home/snaedisk/work/liblinear/blas/dscal.c"
#include <boost/math/distributions/poisson.hpp>
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;
using namespace seqan;

//structure to store read information
struct AttributeLine {
    string PnId;
    float numOfRepeats;
    float ratioBf; 
    float ratioAf; 
    unsigned locationShift; 
    unsigned mateEditDist;
    float purity; 
    float ratioOver20In;
    float ratioOver20After;
    unsigned sequenceLength;
    bool wasUnaligned;
    int label;
    double pValue;
    string repSeq;
} ; 

//For storing marker information
struct Marker {
    string chrom;
    int start;
    int end;
    string motif;
    float refRepeatNum;
    string refRepSeq;
} ;

//For storing number of members in each class and their pValue-sum
struct LabelProps {
    Pair<int, double> p1;
    Pair<int, double> p2;
    Pair<int, double> p3;
} ;

//For storing all possible genotypes, their pValues, the chosen genotype, its pValue, all alleles with their frequencies and number of reads available for the decision 
struct GenotypeInfo {
    String<Pair<float> > genotypes;
    String<long double> pValues;
    Pair<float> genotype;
    double pValue;
    int numOfReads;
    map<float, int> alleleToFreq; //This maps from reported alleles to their frequencies
} ;

//So I can map from Markers in mapPerMarker
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//Sums the pValues of reads and counts the number of reads for each type of label at every marker
map<Marker, LabelProps> markerToSize;
//Parameter, problem and model structs to use in training of logistic regression model and computing pValues
parameter param;
problem prob; 
problem probBig; //problem including reads with label = 2 to use for predict function
model* model_;
double bias = -1;

//For repeating a motif n times and n can be a float
CharString repeat(CharString s, float n) {
    CharString ret;
    int m = floor(n);
    for (int i = 0; i < m; i++) {
        append(ret,s);
    }
    float ratio = 0;
    float fraction = n - m;
    unsigned index = 0;
    while (ratio < fraction)
    {
        appendValue(ret,s[index]);
        ++index;
        ratio = (float)index/(float)length(s);
    }
    return ret;
}

//Fills in the x-part of a problem structure from an AttributeLine structure
void fillProblemX(int idx, AttributeLine currentLine, problem& myProb)
{
    myProb.x[idx][0].index = 1;                              
    myProb.x[idx][0].value = currentLine.ratioBf;
    myProb.x[idx][1].index = 2;
    myProb.x[idx][1].value = currentLine.ratioAf;
    myProb.x[idx][2].index = 3;
    myProb.x[idx][2].value = currentLine.locationShift;
    myProb.x[idx][3].index = 4;
    myProb.x[idx][3].value = currentLine.mateEditDist;
    myProb.x[idx][4].index = 5;
    myProb.x[idx][4].value = currentLine.purity;
    myProb.x[idx][5].index = 6;
    myProb.x[idx][5].value = currentLine.ratioOver20In;
    myProb.x[idx][6].index = 7;
    myProb.x[idx][6].value = currentLine.ratioOver20After;
    myProb.x[idx][7].index = 8;
    myProb.x[idx][7].value = currentLine.sequenceLength;
    myProb.x[idx][8].index = 9;
    myProb.x[idx][8].value = currentLine.wasUnaligned;
    myProb.x[idx][9].index = -1; // This is to indicate that there aren't any more attributes.
    myProb.x[idx][9].value = 0;    
}

//Parses one line from output file of computeReadAttributes by filling up and returning an AttributeLine
AttributeLine parseNextLine(float winner, float second, ifstream& attributeFile, Marker& marker)
{
    AttributeLine currentLine;
    attributeFile >> currentLine.numOfRepeats;
    attributeFile >> currentLine.ratioBf;
    attributeFile >> currentLine.ratioAf;
    attributeFile >> currentLine.locationShift;
    attributeFile >> currentLine.mateEditDist;
    attributeFile >> currentLine.purity;
    attributeFile >> currentLine.ratioOver20In;
    attributeFile >> currentLine.ratioOver20After;
    attributeFile >> currentLine.sequenceLength;
    attributeFile >> currentLine.wasUnaligned;
    attributeFile >> currentLine.repSeq;
    //Initialize pValue for all reads as 0.95
    currentLine.pValue = 0.95;
    //Determining the initial label of the read
    if (currentLine.numOfRepeats == winner || currentLine.numOfRepeats == second)
    {
        currentLine.label = 1;
        ++markerToSize[marker].p1.i1;
        markerToSize[marker].p1.i2 += currentLine.pValue;
    }
    else 
    {
        if ((currentLine.numOfRepeats == winner - 1) || (currentLine.numOfRepeats == second - 1))
        {
            currentLine.label = 2;
            ++markerToSize[marker].p2.i1;
            markerToSize[marker].p2.i2 += currentLine.pValue; 
        }
        else
        { 
            currentLine.label = -1;
            ++markerToSize[marker].p3.i1;
            markerToSize[marker].p3.i2 += currentLine.pValue;
        }
    }
    return currentLine;
}

void combinationUtil(String<float> values, float data[], int start, int end, int index, int r, string PnId, map<string, Pair<String<float>, String<Pair<float> > > >& PnToAlleles)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        appendValue(PnToAlleles[PnId].i2, Pair<int>(data[0],data[1]));
        return;
    }
 
    // replace index with all possible elements. The condition
    // "end-i+1 >= r-index" makes sure that including one element
    // at index will make a combination with remaining elements
    // at remaining positions
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = values[i];
        combinationUtil(values, data, i+1, end, index+1, r, PnId, PnToAlleles);
    }
}

void createHetero(string PnId, map<string, Pair<String<float>, String<Pair<float> > > >& PnToAlleles)
{
    float data[2];
    String<float> values = PnToAlleles[PnId].i1;
    combinationUtil(values, data, 0, length(values)-1, 0, 2, PnId, PnToAlleles);
}

void createHomo(string PnId, map<string, Pair<String<float>, String<Pair<float> > > >& PnToAlleles)
{
    String<float> values = PnToAlleles[PnId].i1;
    for(unsigned i=0; i<length(values); ++i)
        appendValue(PnToAlleles[PnId].i2,Pair<float>(values[i],values[i]));
}

String<Pair<float> > makeGenotypes(std::set<float> alleles)
{
    String<Pair<float> > genotypes;
    String<float> alleleString;
    std::set<float>::iterator allelesEnd = alleles.end();
    for (std::set<float>::iterator alleleIt = alleles.begin(); alleleIt!=allelesEnd; ++alleleIt)
        appendValue(alleleString, *alleleIt);  
    for (unsigned i=0; i<length(alleleString); ++i)
    {
        appendValue(genotypes,Pair<float>(alleleString[i],alleleString[i]));
        if (i == (length(alleleString)-1))
            break;
        for (unsigned j=i+1; j<length(alleleString); ++j)
            appendValue(genotypes,Pair<float>(alleleString[i],alleleString[j])); 
    }
    return genotypes;
}

bool findElement(String<float> alleles, float allele)
{
    if (length(alleles) == 0) 
        return false;
    for (unsigned i=0; i<length(alleles); ++i)
    {
        if (alleles[i] == allele)
            return true;
    }
    return false;
}

int findMaxIndex(String<long double> probs)
{
    int maxIndex = 0;
    long double maxValue = probs[0];
    for (unsigned i = 1; i<length(probs); ++i)
    {
        if (probs[i]>maxValue)
        {
            maxIndex=i;
            maxValue=probs[i];
        }
    }
    return maxIndex;
}

Pair<GenotypeInfo, bool> determineGenotype(String<AttributeLine> reads, float markerSlippage, String<Pair<float> > genotypes, int numberOfAlleles)
{
    GenotypeInfo returnValue;
    returnValue.genotypes = genotypes;
    returnValue.numOfReads = length(reads);
    boost::math::poisson_distribution<> myPoiss(markerSlippage);
    Pair<float> genotypeToCheck;
    AttributeLine readToCheck;
    std::set<float> currentGenotype;
    std::set<float> newGenotypeSet; 
    String<long double> probs;
    long double probSum = 0;
    resize(probs, length(genotypes));
    bool isHomo;
    float posNegSlipp = 1;
    float posNegSlipp2 = 1;
    float diff;
    float diff2;
    int indexOfWinner;
    for (unsigned i=0; i<length(genotypes); ++i)
    {
        probs[i] = 1;
        genotypeToCheck = genotypes[i];
        isHomo = genotypeToCheck.i1 == genotypeToCheck.i2;
        for (unsigned j=0; j<length(reads); ++j)
        {
            posNegSlipp = 1;
            posNegSlipp2 = 1;
            readToCheck = reads[j];
            if (i == 0)
            {
                ++returnValue.alleleToFreq[readToCheck.numOfRepeats];
                //Debugging code
                //cout << "P-value of read " << j << " with " << readToCheck.numOfRepeats << " repeats:  " << readToCheck.pValue << endl;
            }
            if (readToCheck.label == 1)
                currentGenotype.insert(readToCheck.numOfRepeats);
            if (isHomo)
            {
                if (readToCheck.numOfRepeats < genotypeToCheck.i1)
                    posNegSlipp = 0.85;
                if (readToCheck.numOfRepeats > genotypeToCheck.i1)
                    posNegSlipp = 0.15;
                diff = readToCheck.numOfRepeats - genotypeToCheck.i1;
                //Debugging code
                //cout << "Diff for homozygous genotype " << genotypeToCheck.i1 << " from read with " << readToCheck.numOfRepeats << " repeats is: " << diff << endl;
                //cout << "Update of genotype probability: " << probs[i] << " *= (" << readToCheck.pValue << "*" << pdf(myPoiss, abs(diff)) << "*" << posNegSlipp << "+ (" << (double)(1.0-readToCheck.pValue)/(double)numberOfAlleles << "))";
                probs[i] *= (readToCheck.pValue * pdf(myPoiss, fabs(diff)) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                //cout << " = " << probs[i] << endl;
            }
            else
            {
                if (readToCheck.numOfRepeats < genotypeToCheck.i1)
                    posNegSlipp = 0.85;
                if (readToCheck.numOfRepeats > genotypeToCheck.i1)
                    posNegSlipp = 0.15;
                if (readToCheck.numOfRepeats < genotypeToCheck.i2)
                    posNegSlipp2 = 0.85;
                if (readToCheck.numOfRepeats > genotypeToCheck.i2)
                    posNegSlipp2 = 0.15;
                diff = readToCheck.numOfRepeats - genotypeToCheck.i1;
                diff2 = readToCheck.numOfRepeats - genotypeToCheck.i2;
                //Debugging code
                //cout << "Diffs for heterozygous genotype " << genotypeToCheck.i1 << "/" << genotypeToCheck.i2 << " from read with " << readToCheck.numOfRepeats << " repeats are: " << diff << " and " << diff2 << endl;
                //cout << "Update of genotype probability: " << probs[i] << " *= (" << readToCheck.pValue << " * (0.5 * " << pdf(myPoiss, abs(diff)) << "*" << posNegSlipp << "+ 0.5 * "<< pdf(myPoiss, abs(diff2)) << "*" << posNegSlipp2 << ") + (" << (double)(1.0-readToCheck.pValue)/(double)numberOfAlleles << "))"; 
                probs[i] *= (readToCheck.pValue * (0.5 * pdf(myPoiss, fabs(diff)) * posNegSlipp + 0.5 * pdf(myPoiss, fabs(diff2)) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(float)numberOfAlleles));
                //cout << " = " << probs[i] << endl;
            }
        }
        //Debugging code
        //cout << "Genotype: " << genotypeToCheck.i1 << "/" << genotypeToCheck.i2 << " with probability: " << probs[i] << endl;
    }
    for (unsigned i=0; i<length(probs); ++i)
        probSum = probSum+probs[i];
    for (unsigned i=0; i<length(probs); ++i)
        probs[i] = probs[i]/probSum;
    returnValue.pValues = probs;
    indexOfWinner = findMaxIndex(probs);
    returnValue.genotype = genotypes[indexOfWinner];
    returnValue.pValue = probs[indexOfWinner];
    newGenotypeSet.insert(returnValue.genotype.i1);
    newGenotypeSet.insert(returnValue.genotype.i2);    
    if (newGenotypeSet == currentGenotype)
        return Pair<GenotypeInfo, bool>(returnValue, false);
    else
        return Pair<GenotypeInfo, bool>(returnValue, true);
}

void relabelReads(String<AttributeLine>& readsToRelabel, int start, int end, Pair<float> newGenotype, map<Marker, LabelProps>& markerToSize, Marker marker)
{
    for (unsigned i=start; i<end; ++i)
    {
        if ((fabs(readsToRelabel[i].numOfRepeats - newGenotype.i1)<=0.1) || (fabs(readsToRelabel[i].numOfRepeats - newGenotype.i2)<=0.1))
        {
            readsToRelabel[i].label = 1;
            ++markerToSize[marker].p1.i1;
            markerToSize[marker].p1.i2 += readsToRelabel[i].pValue;
        }
        else 
        {
            if ((fabs(readsToRelabel[i].numOfRepeats - (newGenotype.i1 - 1))<=0.1) || (fabs(readsToRelabel[i].numOfRepeats - (newGenotype.i2 - 1))<=0.1))
            {
                readsToRelabel[i].label = 2;
                ++markerToSize[marker].p2.i1;
                markerToSize[marker].p2.i2 += readsToRelabel[i].pValue; 
            }
            else
            { 
                readsToRelabel[i].label = -1;
                ++markerToSize[marker].p3.i1;
                markerToSize[marker].p3.i2 += readsToRelabel[i].pValue;
            }
        }
    }
}
 
//Write all sorts of info to the header of the vfc file I pass to the function
void makeVcfHeader(VcfStream& out, String<string> PnIds, string chrom)
{
    appendValue(out.header.sequenceNames, chrom);
    //Add IDs of all PNs to the header
    for (unsigned i = 0; i<length(PnIds); ++i)
        appendValue(out.header.sampleNames, PnIds[i]);
    
    //Complicated way of getting todays date
    time_t rawtime;
    tm* timeinfo;
    char buffer [80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,10,"%Y%m%d",timeinfo);
    string date = buffer;
    appendValue(out.header.headerRecords, VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("fileDate", date));
    appendValue(out.header.headerRecords, VcfHeaderRecord("source", "genotyperV1.0"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("reference", "/nfs/gpfs/data/reference/Homo_sapiens-NCBI-build37/Sequence/WholeGenomeFasta/genome.fa"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=MOTIF,Number=1,Type=String,Description=\"Repeat motif\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=REF,Number=1,Type=Float,Description=\"Copy number in reference\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=RL,Number=1,Type=Integer,Description=\"Length of STR in reference in bp\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=VT,Number=String,Type=Flag,Description=\"Variant type\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FILTER", "<ID=N,Description=\"N = number of alleles in population\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">"));
    //appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
    //appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"));
    //appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=AF,Number=1,Type=String,Description=\"All reported alleles along with their frequency: a1|f1;a2|f2;...\">"));
    //appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"));
}

//Clear a stringstream and a string, I use this a lot in fillRecord so a function seemed appropriate
void stringClear(stringstream& ss, CharString& str)
{
    ss.str("");  
    ss.clear(); 
    clear(str);
}

VcfRecord fillRecordMarker(Marker marker, std::set<float> allelesAtThisMarker)
{
    int refLength = marker.end - (marker.start-1);
    VcfRecord record;
    record.rID = 0;
    record.beginPos = marker.start-1;
    record.id = ".";
    stringstream ss;
    std::set<float>::iterator smallestIt = allelesAtThisMarker.begin();
    float smallestAllele = *smallestIt;
    ss << smallestAllele;
    CharString str = ss.str();
    record.ref = str;
    stringClear(ss,str);
    record.qual = 0; //Which qual goes here, one for each call??
    record.info = "END=";
    ss << marker.end;
    str = ss.str();
    append(record.info,str);
    stringClear(ss,str);
    append(record.info,";");
    append(record.info,"MOTIF=");
    append(record.info,marker.motif);
    append(record.info,";");
    append(record.info,"REF=");
    ss << marker.refRepeatNum;
    str = ss.str();
    append(record.info,str);
    stringClear(ss,str);
    append(record.info,";");
    append(record.info,"RL=");
    ss << refLength;
    str = ss.str();
    append(record.info,str);
    stringClear(ss,str);
    append(record.info,";");
    append(record.info,"VT=STR");
    record.format = "GT:PL";
    //Loop over all alternative alleles and add them to ALT field in record
    allelesAtThisMarker.erase(smallestAllele);
    float currentAllele;
    std::set<float>::iterator allEnd = allelesAtThisMarker.end();    
    for (std::set<float>::iterator allIt = allelesAtThisMarker.begin(); allIt!=allEnd; ++allIt)
    {            
        currentAllele = *allIt;
        ss << currentAllele;
        str = ss.str();
        append(record.alt, str);
        append(record.alt, ",");
        stringClear(ss,str);
    }
    if (allelesAtThisMarker.size() > 0)
        eraseBack(record.alt);
    else 
        record.alt = ".";
    return record;
}

//Fill up the PN-specific fields in the VCF-record.
void fillRecordPn(GenotypeInfo genotype, VcfRecord& record, String<Pair<float> > genotypesAtThisMarker)
{
    stringstream ss;
    CharString str = ss.str();
    CharString gtInfo; //First I make the string containing the genotype info 
    ss << genotype.genotype.i1;
    str = ss.str();
    append(gtInfo,str);
    stringClear(ss,str);
    append(gtInfo,"/");
    ss << genotype.genotype.i2;
    str = ss.str();
    append(gtInfo,str);
    stringClear(ss,str);
    append(gtInfo,":");
    Pair<float> genotypeToLookFor;
    Pair<float> genotypeToCompare;
    int index;
    long double numerator;
    long double denominator;
    int pl;
    for (unsigned i=0; i<length(genotypesAtThisMarker); ++i)
    {
        index = -1;
        genotypeToLookFor = genotypesAtThisMarker[i];
        for (unsigned j=0; j<length(genotype.genotypes); ++j)
        {
            genotypeToCompare=genotype.genotypes[j];
            if (genotypeToCompare == genotypeToLookFor)
            {
                index = j;
                break;
            }
        }
        if (index == -1)
        {
            ss << 255;
            str = ss.str();
            append(gtInfo,str);
        }
        else
        {
            numerator = genotype.pValues[index];
            denominator = genotype.pValue;
            pl = round(-10*log10((long double)numerator/(long double)denominator));
            pl = std::min(255,pl);
            if (pl == -2147483648)
                cout << "Numerator: " << numerator << " Denominator: " << denominator << endl; 
            ss << pl;
            str = ss.str();
            append(gtInfo,str);
        }
        stringClear(ss,str);
        append(gtInfo,",");
    }
    eraseBack(gtInfo);
    /*ss << numberOfReads;
    str = ss.str();
    append(gtInfo,str);
    append(gtInfo,":");
    stringClear(ss,str);
    map<float, int>::const_iterator ite = genotype.alleleToFreq.end();
    for(map<float, int>::const_iterator it = genotype.alleleToFreq.begin(); it != ite; ++it)
    {
        ss << it->first;
        str = ss.str();
        append(gtInfo,str);
        append(gtInfo,"|");
        stringClear(ss,str);
        ss << it->second;
        str = ss.str();
        append(gtInfo,str);
        append(gtInfo,";");
        stringClear(ss,str);
    }
    eraseBack(gtInfo);
    append(gtInfo,":");
    ss << q;
    str = ss.str();
    append(gtInfo,str);
    append(gtInfo,":");
    stringClear(ss,str);
    ss << gq;
    str = ss.str();
    append(gtInfo,str);*/
    //When the genotype info string is ready I append it to the stringset of charstrings
    appendValue(record.genotypeInfos, gtInfo);
}

double computeAlleleDist(String<Pair<float> > genotypes, map<float,int> allelesToFreq, int pnsAtMarker, int numOfAlleles)
{
    float genotype1, genotype2;
    double freq1, freq2;
    double currVal;
    double totalSum = 0;
    int indicator;
    for (unsigned i=0; i<length(genotypes); ++i)
    {
        indicator = 0;
        genotype1 = genotypes[i].i1;
        genotype2 = genotypes[i].i2;
        float n1 = fabs(genotype1-genotype2);
        float n2 = ceil(fabs(genotype1-genotype2));
        float n3 = n1-n2;
        if (n3 == 0)
            indicator = 1;
        freq1 = (double)allelesToFreq[genotype1]/(double)(2*pnsAtMarker);
        freq2 = (double)allelesToFreq[genotype2]/(double)(2*pnsAtMarker);
        currVal = (double)indicator*freq1*freq2;
        if (genotype1 != genotype2)
            currVal = 2*currVal;
        totalSum += currVal;
    }
    return totalSum;
}

int main(int argc, char const ** argv)
{   
    //Check arguments.
    if (argc != 4)
    {
        cerr << "USAGE: " << argv[0] << " attributeFile initialLabellingFile vcfOutputDirectory\n";
        return 1;
    }
    
    //Make input streams for attribute file and initial labelling file 
    ifstream attributeFile(argv[1]);
    ifstream initialLabels(argv[2]);
    
    string PnId, chrom, motif, nextWord, refRepSeq;
    String<string> PnIds;
    std::set<Marker> markers;
    int start, end, numberOfReads;
    float refRepeatNum;
    float winner, second;
    
    //Map from marker to all reads covering it 
    map<Marker, String<AttributeLine> > mapPerMarker;
    //Count the total number of alleles in the population for each marker -- by checking the size of the set
    map<Marker, std::set<float> > markerToAlleles;
    
    Marker marker;
    AttributeLine currentLine;
    attributeFile >> PnId;
    //Debugging code
    cout << "PnId: " << PnId << endl;
    appendValue(PnIds,PnId);
    attributeFile >> chrom;
    
    //Read training data into maps
    while (!attributeFile.eof() && !initialLabels.eof())
    {
        attributeFile >> start;        
        attributeFile >> end;
        attributeFile >> motif;
        attributeFile >> refRepeatNum;
        attributeFile >> numberOfReads;
        attributeFile >> refRepSeq;
        marker.chrom = chrom;
        marker.start = start;
        marker.end = end;
        marker.motif = motif;
        marker.refRepeatNum = refRepeatNum;
        marker.refRepSeq = refRepSeq;
        markers.insert(marker);
        initialLabels >> winner;
        initialLabels >> second; 
        for (unsigned i = 0; i < numberOfReads; ++i)
        {
            currentLine = parseNextLine(winner, second, attributeFile, marker);
            currentLine.PnId = PnId;
            if (currentLine.label == 1)
                markerToAlleles[marker].insert(currentLine.numOfRepeats);
            appendValue(mapPerMarker[marker],currentLine);
        }
        attributeFile >> nextWord;
        if (nextWord != chrom)
        {
            PnId = nextWord;
            //Debugging code
            cout << "PnId: " << PnId << endl;
            appendValue(PnIds, PnId);
            attributeFile >> chrom;
        }
    }
    //Debugging code 
    cout << "Number of markers: " << mapPerMarker.size() << endl;
    cout << "Number of pns: " << length(PnIds) << endl;
    
    double *prob_estimates = NULL;
    double predict_label;
    float changedRatio;
    int PnsAtMarker, updatedPns;
    unsigned numOfAlleles;
    //Initialize parameter object for logistic regression
    param.solver_type = 0; 
    param.C = 1;
    param.eps = 0.1;
    param.p = 0.1;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    //Initialize problem objects for logistic regression training and predicting
    prob.n = 9;
    probBig.n = 9;
    prob.bias = bias;
    probBig.bias = bias;
    //Map to store alleles present in each individual and the possible genotypes. Is cleared for each marker.
    map<string, Pair<String<float>, String<Pair<float> > > > PnToAlleles;
    //Map to store a map from alleles to their frequencies in the population for each marker, used for estimating the probability that the distance between alleles at the marker is an integer
    map<Marker, Pair<map<float,int>,int> > markerToAlleleFreqs;
    //Map to store current genotype(and lots of other things) of a person for each marker maps from pnId and Marker-struct to GenotypeInfo struct
    map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype;
    //Loop over map from Marker to string<AttributeLine> and train model for each marker and use it to determine genotype
    map<Marker, String<AttributeLine> >::iterator itEnd = mapPerMarker.end();
    int z = 1;
    float markerSlippage;
    for (map<Marker, String<AttributeLine> >::iterator it = mapPerMarker.begin(); it != itEnd; ++it)
    {        
        cout << "Starting marker number: " << z << endl;
        changedRatio=1;
        PnsAtMarker = 1;
        int loops = 0;
        String<AttributeLine>& currentMarker = it->second;
        probBig.l = length(currentMarker);
        //convergence condition = less than 0.5% of Pns get a new genotype 
        while (changedRatio > 0.005)
        {
            //Backup if convergence has not been reached within 10 iterations 
            if (loops > 10)
                break;
            ++loops;
            numOfAlleles = markerToAlleles[it->first].size();
            markerToAlleles[it->first].clear();
            markerToAlleleFreqs[it->first].i1.clear();
            //Set the marker slippage
            markerSlippage = std::max(0.01,(double)(markerToSize[it->first].p2.i2+markerToSize[it->first].p3.i2)/(2*(double)(markerToSize[it->first].p1.i2+markerToSize[it->first].p2.i2+markerToSize[it->first].p3.i2)));
            
            //Debugging code
            /*cout << "Marker slippage: " << markerSlippage << endl;
            boost::math::poisson_distribution<> myPoiss(markerSlippage);
            cout << "Testing poisson distribution: " << endl;
            cout << "pdf(0): " << pdf(myPoiss, 0) << " pdf(1): " << pdf(myPoiss, 1) << " pdf(2): " << pdf(myPoiss, 2)<< endl; */
            
            //Reads with label 2 are not included in training
            prob.l = length(currentMarker) - markerToSize[it->first].p2.i1;
            //Now I "nullSet" the markerToSize map for the marker I am looking at so I can update it when I relabel the reads 
            markerToSize[it->first].p1 = Pair<int,float>(0,0);
            markerToSize[it->first].p2 = Pair<int,float>(0,0);
            markerToSize[it->first].p3 = Pair<int,float>(0,0);
            int idx = 0;
            prob.y = Malloc(double,prob.l);
            prob.x = (feature_node **) malloc(prob.l * sizeof(feature_node *));
            probBig.x = (feature_node **) malloc(probBig.l * sizeof(feature_node *));
            PnId = currentMarker[0].PnId;            
            updatedPns = 0;
            for (unsigned i = 0; i<length(currentMarker); ++i)
            { 
                currentLine = currentMarker[i];
                //I just do this in the first iteration because this map will remain the same despite changes in the labels
                if (loops == 1)
                {
                    if (currentLine.PnId != PnId)
                    {
                        createHomo(PnId, PnToAlleles);
                        createHetero(PnId, PnToAlleles);
                        PnId = currentLine.PnId;
                        ++PnsAtMarker;
                    } 
                    if (!findElement(PnToAlleles[currentLine.PnId].i1, currentLine.numOfRepeats))
                        appendValue(PnToAlleles[currentLine.PnId].i1, currentLine.numOfRepeats);
                }
                probBig.x[i] = (feature_node *) malloc(10 * sizeof(feature_node));
                fillProblemX(i,currentLine,probBig);
                //Just include reads with label = 1 or label = -1 in logistic regression training
                if (currentLine.label == 1 || currentLine.label == -1) 
                {
                    prob.x[idx] = (feature_node *) malloc(10 * sizeof(feature_node));
                    fillProblemX(idx,currentLine,prob);
                    prob.y[idx] = currentLine.label;
                    ++idx;
                }
                
            }
            //Create PnToAlleles entry for last Pn but only in first iteration
            if (loops == 1)
            {
                createHomo(PnId, PnToAlleles); 
                createHetero(PnId, PnToAlleles);
            }
            const char *error_msg;
            error_msg = check_parameter(&prob,&param);
            if (error_msg != NULL)
                cout << "Error message: " << error_msg << endl;
            //Train logistic regression model
            model_ = train(&prob, &param);
            prob_estimates = (double *) malloc(2*sizeof(double));
            Pair<GenotypeInfo, bool> changed;
            PnId = currentMarker[0].PnId;
            String<AttributeLine> reads;
            for (unsigned i = 0; i < length(currentMarker); ++i) 
            {
                currentLine = currentMarker[i];
                if (currentLine.PnId != PnId)
                {
                    //make decision about genotype for PnId at the current marker.
                    //Debugging code
                    //cout << "Decision for: " << PnId << endl;
                    changed = determineGenotype(reads, markerSlippage, PnToAlleles[PnId].i2, numOfAlleles);
                    if (changed.i2)
                        ++updatedPns;
                    relabelReads(currentMarker, i-length(reads), i, changed.i1.genotype, markerToSize, it->first);
                    PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;
                    markerToAlleles[it->first].insert(changed.i1.genotype.i1);
                    markerToAlleles[it->first].insert(changed.i1.genotype.i2);
                    ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
                    ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];
                    PnId = currentLine.PnId;
                    clear(reads);
                }
                //Use logistic regression model to get pValue for all reads at marker
                predict_label = predict_probability(model_,probBig.x[i],prob_estimates);
                free(probBig.x[i]);
                if (i<idx)
                    free(prob.x[i]);
                mapPerMarker[it->first][i].pValue = prob_estimates[0];
                appendValue(reads, mapPerMarker[it->first][i]);
            }
            free(prob.y);
            free(prob.x);
            free(probBig.x);
            free(prob_estimates);
            //cout << "Decision for: " << PnId << endl;
            changed = determineGenotype(reads, markerSlippage, PnToAlleles[PnId].i2, numOfAlleles);
            if (changed.i2)
                ++updatedPns;
            relabelReads(currentMarker, length(currentMarker)-length(reads), length(currentMarker), changed.i1.genotype, markerToSize, it->first);
            PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;  
            markerToAlleles[it->first].insert(changed.i1.genotype.i1);
            markerToAlleles[it->first].insert(changed.i1.genotype.i2);   
            ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
            ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];       
            changedRatio = (float)updatedPns/(float)PnsAtMarker;
            cout << "Changed: " << changedRatio << endl;
        }
        markerToAlleleFreqs[it->first].i2 = PnsAtMarker;
        PnToAlleles.clear();
        mapPerMarker.erase(it->first);
        cout << "Finished marker number: " << z << endl;
        ++z;                 
    }
    cout << "Finished determining genotypes" << endl;
    
    //Here I need code to initialize things common to the vcf files
    GenotypeInfo genotype;
    string thisPn;
    Marker thisMarker;
    VcfStream out;
    VcfRecord record;
    CharString outputDirectory = argv[3];
    append(outputDirectory, chrom);
    append(outputDirectory, ".vcf");
    bool outOk = open(out,toCString(outputDirectory), VcfStream::WRITE);
    makeVcfHeader(out, PnIds, chrom);    
    stringstream ss;
    CharString str;
    std::set<float> allelesAtThisMarker;    
    String<Pair<float> > genotypesAtThisMarker;
    double alleleDistance;
    int pnsAtMarker;
    
    //Loop over markers in outer loop
    std::set<Marker>::iterator markersEnd = markers.end();
    cout << "Starting to write VCF file " << endl;
    for (std::set<Marker>::iterator markerIt = markers.begin(); markerIt!=markersEnd; ++markerIt)
    {
        thisMarker = *markerIt;
        //Determine which alleles are present in the population for this marker
        allelesAtThisMarker = markerToAlleles[thisMarker];
        //Check how many PNs have been typed for this marker
        pnsAtMarker = markerToAlleleFreqs[thisMarker].i2;
        //Check the number of alleles in the population, shouldn't be above ~20.
        if ((allelesAtThisMarker.size() > 20) || (allelesAtThisMarker.size() < 2))
            continue;      
        //Make a String<Pair<float> > which contains a list of genotypes
        genotypesAtThisMarker = makeGenotypes(allelesAtThisMarker);
        //Compute abs(allele1-allele2)*allele1Freq*allele2Freq for all genotypes and return average of those, estimate of distance between alleles.
        alleleDistance = computeAlleleDist(genotypesAtThisMarker, markerToAlleleFreqs[thisMarker].i1, pnsAtMarker, allelesAtThisMarker.size()); 
        cout << "Average allele distance for marker: " << alleleDistance << endl;
        cout << "PNs at marker: " << pnsAtMarker << endl;
        cout << "Alleles at marker: " << allelesAtThisMarker.size() << endl;     
        //First fill marker specific fields of vcfRecord
        record = fillRecordMarker(thisMarker, allelesAtThisMarker);       
        //Loop over Pns in inner loop and fill in PN specific fields of vcfRecord for each PN
        for (unsigned i = 0; i<length(PnIds); ++i)
        {
            thisPn = PnIds[i];
            //If a decision has been made for thisPn at thisMarker I add it to the genotypeInfos stringSet
            if (PnAndMarkerToGenotype.count(Pair<string,Marker>(thisPn, thisMarker)) != 0)
            {
                genotype = PnAndMarkerToGenotype[Pair<string,Marker>(thisPn, thisMarker)];
                fillRecordPn(genotype, record, genotypesAtThisMarker);
                PnAndMarkerToGenotype.erase(Pair<string,Marker>(thisPn, thisMarker));
            }
            //If a decision has not been made I add a CharString with no decision(0:0,0,0,0....etc) to the set to maintain order of Pns vs genotypeInfos in output
            else 
            {         
                CharString gtInfo = "0:";       
                for (unsigned i=0; i<length(genotypesAtThisMarker); ++i)
                    append(gtInfo,"0,");
                eraseBack(gtInfo);
                appendValue(record.genotypeInfos, gtInfo);
            }
        }                
        ss << allelesAtThisMarker.size();
        str = ss.str();
        record.filter = str;
        stringClear(ss,str);
        //After adding info for all PNs at thisMarker I write the record to the vcf output file
        if (writeRecord(out, record) != 0)
            cerr << "ERROR: Problem writing VCF-output file.";
        clear(record);
        allelesAtThisMarker.clear();
    }
    return 0;
}
