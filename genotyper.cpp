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
    double numOfRepeats;
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
    double refRepeatNum;
    string refRepSeq;
} ;

//For storing number of members in each class and their pValue-sum
struct LabelProps {
    Pair<int, double> p1;
    Pair<int, double> p2;
    Pair<int, double> p3;
} ;

//For storing all possible genotypes, their pValues, the chosen genotype, its pValue and number of reads available for the decision 
struct GenotypeInfo {
    String<Pair<double> > genotypes;
    String<double> pValues;
    Pair<double> genotype;
    double pValue;
    int numOfReads;
} ;

//So I can map from Markers in mapPerMarker
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//Sums the pValues of reads and counts the number of reads for each type of label at every marker
map<Marker, LabelProps> markerToSize;
//Count the total number of alleles in the population for each marker -- by checking the size of the set
map<Marker, std::set<double> > markerToAlleles;
//Parameter, problem and model structs to use in training of logistic regression model and computing pValues
parameter param;
problem prob; 
problem probBig; //problem including reads with label = 2 to use for predict function
model* model_;
double bias = -1;

//For repeating a motif n times and n can be a double
CharString repeat(CharString s, double n) {
    CharString ret;
    int m = floor(n);
    for (int i = 0; i < m; i++) {
        append(ret,s);
    }
    double ratio = 0;
    double fraction = n - m;
    unsigned index = 0;
    while (ratio < fraction)
    {
        appendValue(ret,s[index]);
        ++index;
        ratio = (double)index/(double)length(s);
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
AttributeLine parseNextLine(double winner, double second, ifstream& attributeFile, Marker& marker)
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
    markerToAlleles[marker].insert(currentLine.numOfRepeats);
    return currentLine;
}

void combinationUtil(String<double> values, double data[], int start, int end, int index, int r, string PnId, map<string, Pair<String<double>, String<Pair<double> > > >& PnToAlleles)
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

void createHetero(string PnId, map<string, Pair<String<double>, String<Pair<double> > > >& PnToAlleles)
{
    double data[2];
    String<double> values = PnToAlleles[PnId].i1;
    combinationUtil(values, data, 0, length(values)-1, 0, 2, PnId, PnToAlleles);
}

void createHomo(string PnId, map<string, Pair<String<double>, String<Pair<double> > > >& PnToAlleles)
{
    String<double> values = PnToAlleles[PnId].i1;
    for(unsigned i=0; i<length(values); ++i)
        appendValue(PnToAlleles[PnId].i2,Pair<double>(values[i],values[i]));
}

bool findElement(String<double> alleles, double allele)
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

int findMaxIndex(String<double> probs)
{
    int maxIndex = 0;
    double maxValue = probs[0];
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

//Pair<Pair<double, bool>, Pair<double> > determineGenotype(String<AttributeLine> reads, double markerSlippage, String<Pair<double> > genotypes, int numberOfAlleles)
Pair<GenotypeInfo, bool> determineGenotype(String<AttributeLine> reads, double markerSlippage, String<Pair<double> > genotypes, int numberOfAlleles)
{
    GenotypeInfo returnValue;
    returnValue.genotypes = genotypes;
    returnValue.numOfReads = length(reads);
    boost::math::poisson_distribution<> myPoiss(markerSlippage);
    Pair<double> genotypeToCheck;
    AttributeLine readToCheck;
    std::set<double> currentGenotype;
    std::set<double> newGenotypeSet; 
    String<double> probs; 
    double probSum = 0; 
    resize(probs, length(genotypes));
    bool isHomo;
    double posNegSlipp = 1;
    double posNegSlipp2 = 1;
    int diff;
    int diff2;
    int indexOfWinner;
    for (unsigned i=0; i<length(genotypes); ++i)
    {
        probs[i] = 1;
        genotypeToCheck = genotypes[i];
        isHomo = genotypeToCheck.i1 == genotypeToCheck.i2;
        for (unsigned j=0; j<length(reads); ++j)
        {
            readToCheck = reads[j];
            if (readToCheck.label == 1)
                currentGenotype.insert(readToCheck.numOfRepeats);
            if (isHomo)
            {
                if (readToCheck.numOfRepeats < genotypeToCheck.i1)
                    posNegSlipp = 0.85;
                if (readToCheck.numOfRepeats > genotypeToCheck.i1)
                    posNegSlipp = 0.15;
                diff = round(readToCheck.numOfRepeats - genotypeToCheck.i1);
                probs[i] *= (readToCheck.pValue * pdf(myPoiss, abs(diff)) * posNegSlipp + (1-(float)readToCheck.pValue/(float)numberOfAlleles));
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
                diff = round(readToCheck.numOfRepeats - genotypeToCheck.i1);
                diff2 = round(readToCheck.numOfRepeats - genotypeToCheck.i2);
                probs[i] *= (readToCheck.pValue * (0.5 * pdf(myPoiss, abs(diff)) * posNegSlipp + 0.5 * pdf(myPoiss, abs(diff2)) * posNegSlipp2) + (1-(float)readToCheck.pValue/(float)numberOfAlleles));
            }
        }
        probSum += probs[i];
    }
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

void relabelReads(String<AttributeLine>& readsToRelabel, int start, int end, Pair<double> newGenotype, map<Marker, LabelProps>& markerToSize, Marker marker)
{
    for (unsigned i=start; i<end; ++i)
    {
        if ((abs(readsToRelabel[i].numOfRepeats - newGenotype.i1)<=0.1) || (abs(readsToRelabel[i].numOfRepeats - newGenotype.i2)<=0.1))
        {
            readsToRelabel[i].label = 1;
            ++markerToSize[marker].p1.i1;
            markerToSize[marker].p1.i2 += readsToRelabel[i].pValue;
        }
        else 
        {
            if ((abs(readsToRelabel[i].numOfRepeats - (newGenotype.i1 - 1))<=0.1) || (abs(readsToRelabel[i].numOfRepeats - (newGenotype.i2 - 1))<=0.1))
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

//Function to check if a vcf file already exists for the PN in the given directory, then I just append genotyping results to it instead of making a new one.
bool exists(const string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

//Write all sorts of info to the header of the vfc file I pass to the function
void makeVcfHeader(VcfStream& out, string PnId, string chrom)
{
    clear(out.header.sequenceNames);
    clear(out.header.sampleNames);
    clear(out.header.headerRecords);
    appendValue(out.header.sequenceNames, chrom);
    appendValue(out.header.sampleNames, PnId);
    
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
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=RPA,Number=A,Type=Float,Description=\"Repeats per allele\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=MOTIF,Number=1,Type=String,Description=\"Repeat motif\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=REF,Number=1,Type=Float,Description=\"Copy number in reference\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=RL,Number=1,Type=Integer,Description=\"Length of STR in reference in bp\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=VT,Number=String,Type=Flag,Description=\"Variant type\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FILTER", "<ID=U5,Description=\"Number of reads used for Genotyping<5\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FILTER", "<ID=U10,Description=\"Number of reads used for Genotyping<10\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"));
}

//Function to update the vcf header if a file aldready exists for the current pn containing genotyping data for another chromosome.
void updateVcfHeader(VcfStream& out, string chrom)
{
    appendValue(out.header.sequenceNames, chrom);
}

//Fill up the VCF-record before I write it to the outputStream
VcfRecord fillRecord(GenotypeInfo genotype, Marker marker)
{
    bool isHomo = genotype.genotype.i1 == genotype.genotype.i2;
    int numberOfReads = genotype.numOfReads;
    double pValue = genotype.pValue;
    double gq = -10*log10(std::max(numeric_limits< double >::min(),1-pValue));
    int refLength = marker.end - marker.start;
    VcfRecord record;
    record.rID = 0;
    record.beginPos = marker.start;
    record.id = ".";
    record.ref = marker.refRepSeq;
    if (isHomo)
        if (genotype.genotype.i1 == marker.refRepeatNum)
            record.alt = ".";
        else
            record.alt = repeat(marker.motif,genotype.genotype.i1);
    else
    {
        if (genotype.genotype.i1 == marker.refRepeatNum)
            record.alt = repeat(marker.motif,genotype.genotype.i2);
        else
        {
            if (genotype.genotype.i2 == marker.refRepeatNum)
                record.alt = repeat(marker.motif,genotype.genotype.i1);
            else
            {
                record.alt = repeat(marker.motif,genotype.genotype.i1);
                append(record.alt,",");
                CharString temp = repeat(marker.motif,genotype.genotype.i2);
                append(record.alt,temp);
            }
        }
    } 
    record.qual = 29; //Which qual goes here??
    if (numberOfReads < 5) 
        record.filter = "U5"; 
    else
    {
        if (numberOfReads < 10)
            record.filter = "U10";
        else 
            record.filter = "PASS";
    }
    record.info = "END=";
    stringstream ss;
    ss << marker.end;
    CharString str = ss.str();
    append(record.info,str);
    ss.str("");  
    ss.clear(); 
    clear(str);
    append(record.info,";");
    append(record.info,"MOTIF=");
    append(record.info,marker.motif);
    append(record.info,";");
    append(record.info,"REF=");
    ss << marker.refRepeatNum;
    str = ss.str();
    append(record.info,str);
    ss.str("");  
    ss.clear(); 
    clear(str);
    append(record.info,";");
    append(record.info,"RPA=");
    if (isHomo)
    {
        ss << genotype.genotype.i1;
        str = ss.str();
        append(record.info,str);
    }
    else
    {
        ss << genotype.genotype.i1;
        str = ss.str();
        append(record.info,str);
        ss.str("");  
        ss.clear(); 
        clear(str);
        append(record.info,",");
        ss << genotype.genotype.i2;
        str = ss.str();
        append(record.info,str);
    }
    ss.str("");  
    ss.clear(); 
    clear(str);
    append(record.info,";");
    append(record.info,"RL=");
    ss << refLength;
    str = ss.str();
    append(record.info,str);
    ss.str("");  
    ss.clear(); 
    clear(str);
    append(record.info,";");
    append(record.info,"VT=STR"); 
    record.format = "GT:DP:GL:GQ";
    CharString gtInfo; //First I make the string containing the genotype info 
    if (isHomo)
    {
        ss << genotype.genotype.i1;
        str = ss.str();
        append(gtInfo,str);
    }
    else
    {
        ss << genotype.genotype.i1;
        str = ss.str();
        append(gtInfo,str);
        ss.str("");  
        ss.clear(); 
        clear(str);
        append(gtInfo,"/");
        ss << genotype.genotype.i2;
        str = ss.str();
        append(gtInfo,str);
    }
    ss.str("");  
    ss.clear(); 
    clear(str);
    append(gtInfo,":");
    ss << numberOfReads;
    str = ss.str();
    append(gtInfo,str);
    append(gtInfo,":");
    ss.str("");  
    ss.clear(); 
    clear(str);
    append(gtInfo,".:");
    ss << gq;
    str = ss.str();
    append(gtInfo,str);
    appendValue(record.genotypeInfos, gtInfo); //When the genotype info string is ready I append it to the stringset of charstrings
    return record;
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
    double refRepeatNum;
    double winner, second;
    
    //Map from marker to all reads covering it 
    map<Marker, String<AttributeLine> > mapPerMarker;
    
    Marker marker;
    AttributeLine currentLine;
    attributeFile >> PnId;
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
            appendValue(mapPerMarker[marker],currentLine);
        }
        attributeFile >> nextWord;
        if (nextWord != chrom)
        {
            PnId = nextWord;
            appendValue(PnIds, PnId);
            attributeFile >> chrom;
        }
    }
    cout << "Number of markers: " << mapPerMarker.size() << endl;
    
    double *prob_estimates = NULL;
    double predict_label;
    double changedPns;
    int PnsAtMarker, updatedPns;
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
    map<string, Pair<String<double>, String<Pair<double> > > > PnToAlleles;
    //Map to store current genotype of person for each marker 
    map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype;
    //Loop over map from Marker to string<AttributeLine> and train model for each marker
    map<Marker, String<AttributeLine> >::iterator itEnd = mapPerMarker.end();
    for (map<Marker, String<AttributeLine> >::iterator it = mapPerMarker.begin(); it != itEnd; ++it)
    {
        changedPns=1;
        int loops = 0;
        String<AttributeLine>& currentMarker = it->second;
        while (changedPns > 0.005)
        {
            if (loops > 10)
                break;
            ++loops;
            //Set the marker slippage
            double markerSlippage = std::max(0.01,(double)(markerToSize[it->first].p2.i2+markerToSize[it->first].p3.i2)/(2*(double)(markerToSize[it->first].p1.i2+markerToSize[it->first].p2.i2+markerToSize[it->first].p3.i2)));
            //Reads with label 2 are not included in training
            prob.l = length(currentMarker) - markerToSize[it->first].p2.i1;
            //Now I "nullSet" the markerToSize map for the marker I am looking at so I can update it when I relabel the reads 
            markerToSize[it->first].p1 = Pair<int,double>(0,0);
            markerToSize[it->first].p2 = Pair<int,double>(0,0);
            markerToSize[it->first].p3 = Pair<int,double>(0,0);
            probBig.l = length(currentMarker);
            int idx = 0;
            prob.y = Malloc(double,prob.l);
            prob.x = (feature_node **) malloc(prob.l * sizeof(feature_node *));
            probBig.x = (feature_node **) malloc(probBig.l * sizeof(feature_node *));
            PnId = currentMarker[0].PnId;
            PnsAtMarker = 1;
            updatedPns = 0;
            for (unsigned i = 0; i<length(currentMarker); ++i)
            { 
                currentLine = currentMarker[i];
                if (currentLine.PnId != PnId)
                {
                    createHomo(PnId, PnToAlleles);
                    createHetero(PnId, PnToAlleles);
                    PnId = currentLine.PnId;
                    ++PnsAtMarker;
                } 
                if (!findElement(PnToAlleles[currentLine.PnId].i1, currentLine.numOfRepeats))
                    appendValue(PnToAlleles[currentLine.PnId].i1, currentLine.numOfRepeats);
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
            createHomo(PnId, PnToAlleles); 
            createHetero(PnId, PnToAlleles);
            const char *error_msg;
            error_msg = check_parameter(&prob,&param);
            if (error_msg != NULL)
                cout << "Error message: " << error_msg << endl;
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
                    changed = determineGenotype(reads, markerSlippage, PnToAlleles[PnId].i2, markerToAlleles[it->first].size());
                    if (changed.i2)
                        ++updatedPns;
                    relabelReads(currentMarker, i-length(reads), i, changed.i1.genotype, markerToSize, it->first);
                    PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;
                    PnId = currentLine.PnId;
                    clear(reads);
                }
                predict_label = predict_probability(model_,probBig.x[i],prob_estimates);
                mapPerMarker[it->first][i].pValue = prob_estimates[0];
                appendValue(reads, mapPerMarker[it->first][i]);
            }
            changed = determineGenotype(reads, markerSlippage, PnToAlleles[PnId].i2, markerToAlleles[it->first].size());
            if (changed.i2)
                ++updatedPns;
            relabelReads(currentMarker, length(currentMarker)-length(reads), length(currentMarker), changed.i1.genotype, markerToSize, it->first);
            PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;
            PnToAlleles.clear();
            changedPns = (double)updatedPns/(double)PnsAtMarker;
            cout << "Changed: " << changedPns << endl;
        }                 
    }
    
    //Here I need code to initialize things common to the vcf files
    GenotypeInfo genotype;
    string thisPn;
    Marker thisMarker;
    VcfStream out;
    VcfRecord record;
    
    //Loop over Pns in outer loop
    for (unsigned i = 0; i<length(PnIds); ++i)
    {
        CharString outputDirectory = argv[3];
        thisPn = PnIds[i];
        append(outputDirectory,thisPn);
        append(outputDirectory,".vcf"); 
        //If the file doesn't exist I make a new output vcf stream and a header for it. 
        if (!exists(toCString(outputDirectory)))
        {
            bool outOk = open(out,toCString(outputDirectory), VcfStream::WRITE);
            makeVcfHeader(out, thisPn, chrom);
            //Loop over markers in inner loop
            std::set<Marker>::iterator end = markers.end();
            for (std::set<Marker>::iterator markerIt = markers.begin(); markerIt!=end; ++markerIt)
            {
                thisMarker = *markerIt;
                //Check if a decision has been made for this person at this marker.
                if (PnAndMarkerToGenotype.count(Pair<string,Marker>(thisPn, thisMarker)) != 0)
                {
                    genotype = PnAndMarkerToGenotype[Pair<string,Marker>(thisPn, thisMarker)];
                    record = fillRecord(genotype, thisMarker);
                    if (writeRecord(out, record) != 0)
                        cerr << "ERROR: Problem writing VCF-output file.";
                    clear(record);
                }
            }
        }
        //If the file exists I read its header to a temp header, update it and then make a new output vcfStream and copy the temp header to it.
        else 
        {
            //Make a string of vcfRecords to store the ones in the existing file and then write them back out. 
            String<VcfRecord> oldRecords;
            //open the existing vcf file to read it 
            VcfStream in(toCString(outputDirectory));
            //copy the header of the existing file to a temp header
            VcfHeader tempHeader = in.header;
            //Structure to store the old record
            VcfRecord oldRecord;
            //go through the existing file and store all records to a string
            while (!atEnd(in))
            {
                if (readRecord(oldRecord, in) != 0)
                {
                    cerr << "ERROR: Problem reading from existing vcf file\n";
                    return 1;
                }
                appendValue(oldRecords,oldRecord);
            }
            //Close the stream when I finish reading from it
            int closeIn = close(in);
            //Open it again to write to it
            bool outOk = open(out, toCString(outputDirectory), VcfStream::WRITE);
            //Assign the temp header to the new VCF stream and add the chromosome I'm working on to sequenceNames
            out.header = tempHeader; 
            updateVcfHeader(out, chrom);
            //First write old records to the stream
            for (unsigned i = 0; i<length(oldRecords); ++i)
            {
                if (writeRecord(out, oldRecords[i]) != 0)
                {
                    cerr << "ERROR: Problem writing to existing vcf file.\n";
                    return 1;
                }
            }
            //Loop over markers in inner loop (write new records to the stream)
            std::set<Marker>::iterator end = markers.end();
            for (std::set<Marker>::iterator markerIt = markers.begin(); markerIt!=end; ++markerIt)
            {
                thisMarker = *markerIt;
                //Check if a decision has been made for this person at this marker.
                if (PnAndMarkerToGenotype.count(Pair<string,Marker>(thisPn, thisMarker)) != 0)
                {
                    genotype = PnAndMarkerToGenotype[Pair<string,Marker>(thisPn, thisMarker)];
                    record = fillRecord(genotype, thisMarker);
                    if (writeRecord(out, record) != 0)
                        cerr << "ERROR: Problem writing VCF-output file.";
                    clear(record);
                }
            }
        }
        //Close the output stream and clear output directory string        
        close(out);
        clear(outputDirectory);
    }
    return 0;
}
