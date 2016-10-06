#include <iostream>
#include <set>
#include <map>
#include <string>
#include <ctime>
#include <math.h> 
#include <limits.h>
#include <vector>
#include <numeric>
#include <sstream>
#include <sys/stat.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/vcf_io.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <liblinear-2.01/linear.h>
#include <liblinear-2.01/linear.cpp>
#include <liblinear-2.01/tron.h>
#include <liblinear-2.01/tron.cpp>
#include <liblinear-2.01/blas/blas.h>
#include <liblinear-2.01/blas/blasp.h>
#include <liblinear-2.01/blas/daxpy.c>
#include <liblinear-2.01/blas/ddot.c>
#include <liblinear-2.01/blas/dnrm2.c>
#include <liblinear-2.01/blas/dscal.c>
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
    double pValueSum;
    double fullMotifSlippageSum;
} ;

//For storing command line arguments 
struct MsGenotyperOptions
{
    CharString attDirChromNum, pnSlippageFile, markerSlippageFile, vcfOutputDirectory, vcfFileName, intervalIndex;    
    unsigned int startCoordinate;
    unsigned long int endCoordinate;
    
    MsGenotyperOptions() :
        startCoordinate(0), endCoordinate(ULONG_MAX)
    {}
} ;

struct MakeGenotypesRet {
	std::set<Pair<float> > genotypesSet;
	String<Pair<float> > genotypes;
};

//So I can map from Markers
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//Sums the pValues of reads and counts the number of reads for each type of label at every marker, also stores the current slippage rate value for each marker
map<Marker, Pair<LabelProps,double> > markerToLabelsAndSlipp;
//Stores the slippage rate value for each PN - read from input file.
map<string, double> pnToSize;
//Maps from pnId to reads and labels for one marker, is cleared for each marker 
map<string, String<Pair<float> > > pnToLabels;
//Maps from pnIds to the sequencing methods
map<string, string> pnToSeqMethod;
//Map from marker to (average stepsize)%period
map<Marker, double> markerToStepSum;
//for debugging purposes
bool useGeom = true;

/*Map to store sum of pValues for each alleleLength at each marker
map<Marker, map<int,double> > MarkToAllLenToPvalSum;
//Map to store sum of pValue*alleleLength for all reads at each marker
map<Marker, double> MarkToPvalAllLenSum;*/

//Parameter, problem and model structs to use in training of logistic regression model and computing pValues
parameter param;
problem prob; 
problem probBig; //problem structure including reads with label = 2 to use for predict function
model* model_;
double bias = -1;

ArgumentParser::ParseResult parseCommandLine(MsGenotyperOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("msGenotyperDefault");
    setShortDescription(parser, "Microsatellite genotyper");
    setVersion(parser, "1.3");
    setDate(parser, "February 2016");
    addUsageLine(parser, "\\fI-ADCN\\fP attributesDirectory/chromNum \\fI-PNS\\fP pnSlippageFile \\fI-MS\\fP markerSlippageFile \\fI-VD\\fP vcfOutputDirectory \\fI-VN\\fP vcfFileName \\fI-R\\fP start end \\fI-I\\fP intervalIndex ");
    addDescription(parser, "This program performs genptyping on a per chromosome basis for all PNs in the pnSlippageFile given that it can find an attribute file for the PN. The genotypes are written to a file specified by the user.");
    
    addOption(parser, ArgParseOption("ADCN", "attributesDirectory/chromNum", "Path to attributes files for the chromosome being genotyped.", ArgParseArgument::INPUTFILE, "IN-DIR"));
    setRequired(parser, "attributesDirectory/chromNum");
    
    addOption(parser, ArgParseOption("PNS", "pnSlippageFile", "A file containing slippage rates for the pns to be genotyped.", ArgParseArgument::INPUTFILE, "IN-FILE"));
    setRequired(parser, "pnSlippageFile");
    
    addOption(parser, ArgParseOption("MS", "markerSlippageFile", "A file containing slippage rates for the microsatellites.", ArgParseArgument::OUTPUTFILE, "OUT-FILE"));
    setRequired(parser, "markerSlippageFile");
    
    addOption(parser, ArgParseOption("VD", "vcfOutputDirectory", "A directory to write the vcf file to.", ArgParseArgument::OUTPUTFILE, "OUT-DIR"));
    setRequired(parser, "vcfOutputDirectory");
    
    addOption(parser, ArgParseOption("VN", "vcfFileName", "Name of vcf output file.", ArgParseArgument::OUTPUTFILE, "OUT-FILE"));
    setRequired(parser, "vcfFileName");
    
    addOption(parser, ArgParseOption("R", "range", "The range to genotype.", ArgParseArgument::INTEGER, "BEGIN END", false, 2));
    setRequired(parser, "range");
    
    addOption(parser, ArgParseOption("I", "intervalIndex", "Index of the interval being processed", ArgParseArgument::STRING, "INDEX"));
    setRequired(parser, "intervalIndex");    
	setDefaultValue(parser, "intervalIndex", "1");
	
	ArgumentParser::ParseResult res = parse(parser, argc, argv);
	
	if (res != ArgumentParser::PARSE_OK)
	    return res;
	    
	getOptionValue(options.attDirChromNum, parser, "attributesDirectory/chromNum");
	getOptionValue(options.pnSlippageFile, parser, "pnSlippageFile");
	getOptionValue(options.markerSlippageFile, parser, "markerSlippageFile");
	getOptionValue(options.vcfOutputDirectory, parser, "vcfOutputDirectory");
	getOptionValue(options.vcfFileName, parser, "vcfFileName");
	getOptionValue(options.startCoordinate, parser, "range", 0);
	getOptionValue(options.endCoordinate, parser, "range", 1);
	getOptionValue(options.intervalIndex, parser, "intervalIndex");
	
	return ArgumentParser::PARSE_OK;
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
    myProb.x[idx][9].index = -1; // This is to indicate that there aren't any more attributes to read in.
    myProb.x[idx][9].value = 0;    
}

//Parses one line from attribute file by filling up and returning an AttributeLine, also initializes markerToLabelsAndSlipp map using the labels 
AttributeLine parseNextLine(float winner, float second, ifstream& attributeFile, Marker& marker, string PnId, map<Pair<string,Marker>, GenotypeInfo>& PnAndMarkerToGenotype, String<string> firstLine, bool useFirstLine, bool enoughReads)
{
    PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].genotype = Pair<float>(winner,second);
    AttributeLine currentLine;
    currentLine.PnId = PnId;
    if (useFirstLine)
    {
        currentLine.numOfRepeats = lexicalCast<float>(firstLine[0]);
        currentLine.ratioBf = lexicalCast<float>(firstLine[1]);
        currentLine.ratioAf = lexicalCast<float>(firstLine[2]);
        currentLine.locationShift = lexicalCast<unsigned int>(firstLine[3]);
        currentLine.mateEditDist = lexicalCast<unsigned int>(firstLine[4]);
        currentLine.purity = lexicalCast<float>(firstLine[5]);
        currentLine.ratioOver20In = lexicalCast<float>(firstLine[6]);
        currentLine.ratioOver20After = lexicalCast<float>(firstLine[7]);
        currentLine.sequenceLength = lexicalCast<unsigned int>(firstLine[8]);
        currentLine.wasUnaligned = lexicalCast<bool>(firstLine[9]);        
        currentLine.repSeq = firstLine[10];        
    }
    else
    {
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
    }
    currentLine.pValue = 0.95;    
    if (enoughReads)    
        PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].pValueSum += currentLine.pValue;
    //Check if the read is a result of a full motif slippage
    float diff1 = fabs(winner-currentLine.numOfRepeats), diff2 = fabs(second-currentLine.numOfRepeats);
    if (std::min(diff1,diff2)>=0.9)
        PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].fullMotifSlippageSum += currentLine.pValue;
    //Determining the initial label of the read
    if ((fabs(currentLine.numOfRepeats - winner) <= 0.05) || (fabs(currentLine.numOfRepeats - second) <= 0.05))
    {
        currentLine.label = 1;
        ++markerToLabelsAndSlipp[marker].i1.p1.i1;        
        if (enoughReads)           
            markerToLabelsAndSlipp[marker].i1.p1.i2 += currentLine.pValue;
    }
    else 
    {
        if ((fabs(currentLine.numOfRepeats - (winner - 1)) <= 0.05) || (fabs(currentLine.numOfRepeats - (second - 1)) <= 0.05))
        {
            currentLine.label = 2;            
            ++markerToLabelsAndSlipp[marker].i1.p2.i1;
            if (enoughReads)
                markerToLabelsAndSlipp[marker].i1.p2.i2 += currentLine.pValue;
            markerToStepSum[marker] += (float)marker.motif.size();
        }
        else
        {             
            currentLine.label = -1;
            ++markerToLabelsAndSlipp[marker].i1.p3.i1;
            if (enoughReads)
                markerToLabelsAndSlipp[marker].i1.p3.i2 += currentLine.pValue;
            float diff1, diff2;
            diff1 = fabs(currentLine.numOfRepeats - winner);
            diff2 = fabs(currentLine.numOfRepeats - second);
            if (fmod(diff1,1.0)<0.05)
            {
                if (fmod(diff2,1.0)<0.05)
                    markerToStepSum[marker] += round(std::min(diff1,diff2) * (float)marker.motif.size());
                else 
                    markerToStepSum[marker] += round(diff1 * (float)marker.motif.size());
            }
            else 
            {
                if (fmod(diff2,1.0)<0.05)
                    markerToStepSum[marker] += round(diff2 * (float)marker.motif.size());
                else
                    markerToStepSum[marker] += round(std::min(diff1,diff2) * (float)marker.motif.size());
            }
        }
    }    
    return currentLine;
}

//TODO:Fix this so it creates the genotypes in the right order for the vcf file. Need to pass the reference allele too.
MakeGenotypesRet makeGenotypes(std::set<float> alleles)
{
	String<Pair<float> > genotypes;
    String<float> alleleString;
    std::set<Pair<float> > genotypeSet;
    std::set<float>::reverse_iterator allelesBegin = alleles.rend();
    for (std::set<float>::reverse_iterator alleleIt = alleles.rbegin(); alleleIt!=allelesBegin; ++alleleIt)
        appendValue(alleleString, *alleleIt);  
    for (unsigned i=0; i<length(alleleString); ++i)
    {
        appendValue(genotypes,Pair<float>(alleleString[i],alleleString[i]));
        genotypeSet.insert(Pair<float>(alleleString[i],alleleString[i]));
        if (i == (length(alleleString)-1))
            break;
        for (unsigned j=i+1; j<length(alleleString); ++j)
        {
            appendValue(genotypes,Pair<float>(alleleString[j],alleleString[i]));
            genotypeSet.insert(Pair<float>(alleleString[j],alleleString[i]));
        }
    }
    reverse(genotypes);
    MakeGenotypesRet returnValue;
    returnValue.genotypes = genotypes;
    returnValue.genotypesSet = genotypeSet;
    return returnValue;
}

int findMinIndex(String<long double> probs)
{
    int minIndex = 0;
    long double minValue = probs[0];
    for (unsigned i = 1; i<length(probs); ++i)
    {
        if (probs[i]<minValue)
        {
            minIndex=i;
            minValue=probs[i];
        }
    }
    return minIndex;
}

float dgeom(int diff, double psucc) 
{
    if (diff < 0) 
        return 0;
    double p = psucc;
    for (int i = 0; i < diff; i++) 
    {
        p = p*(1-psucc);
    }
    return p;
}

float dpois(int step, float mean) {
  if (step < 0) 
    return 0;
  float p = exp(-1*mean);
  for (int i = 0; i < step; i++) {
    p = p*mean;
    p = p/(i+1);
  }
  return p;  
}

Pair<GenotypeInfo, Pair<bool> > determineGenotype(String<AttributeLine> reads, double markerSlippage, String<Pair<float> > genotypes, int numberOfAlleles, int motifLength, double psucc)
{
    GenotypeInfo returnValue;
    returnValue.pValueSum = 0;
    returnValue.fullMotifSlippageSum = 0;
    returnValue.genotypes = genotypes;
    returnValue.numOfReads = length(reads);
    Pair<float> genotypeToCheck;
    AttributeLine readToCheck;
    std::set<float> currentGenotype, newGenotypeSet; 
    String<long double> probs;
    double errorProbSum = 0, lengthSlippage;
    resize(probs, length(genotypes));
    bool isHomo, enoughDistance = true;
    float posNegSlipp = 1, posNegSlipp2 = 1, lambda = std::max((double)0.001,markerSlippage), diff, diff2;
    int indexOfWinner, indexOfSecond;    
    for (unsigned i=0; i<length(genotypes); ++i)
    {
        probs[i] = 0;
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
                if (length(reads) >= 10)
                    returnValue.pValueSum += readToCheck.pValue;
            }
            if (readToCheck.label == 1)
                currentGenotype.insert(readToCheck.numOfRepeats);
            if (isHomo)
            {
                if (readToCheck.numOfRepeats < genotypeToCheck.i1)
                    posNegSlipp = 0.8;
                if (readToCheck.numOfRepeats > genotypeToCheck.i1)
                    posNegSlipp = 0.2;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                if (useGeom)
                    probs[i] += -(double)10*log10(readToCheck.pValue * dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * dpois(floor(diff), lambda) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                else                
                    probs[i] += -(double)10*log10(readToCheck.pValue * dpois(ceil(diff), lambda) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
            }
            else
            {
                if (readToCheck.numOfRepeats < genotypeToCheck.i1)
                    posNegSlipp = 0.8;
                if (readToCheck.numOfRepeats > genotypeToCheck.i1)
                    posNegSlipp = 0.2;
                if (readToCheck.numOfRepeats < genotypeToCheck.i2)
                    posNegSlipp2 = 0.8;
                if (readToCheck.numOfRepeats > genotypeToCheck.i2)
                    posNegSlipp2 = 0.2;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                diff2 = fabs(readToCheck.numOfRepeats - genotypeToCheck.i2);
                if (useGeom)
                    probs[i] += -(double)10*log10(readToCheck.pValue * 0.5 * (dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * dpois(floor(diff), lambda) * posNegSlipp + dgeom(static_cast<int>((diff2-(float)floor(diff2))*motifLength), psucc) * dpois(floor(diff2), lambda) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));    
                else
                    probs[i] += -(double)10*log10(readToCheck.pValue * 0.5 * (dpois(ceil(diff), lambda) * posNegSlipp + dpois(ceil(diff2), lambda) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));                
            }
        }
    } 
    returnValue.pValues = probs;    
    indexOfWinner = findMinIndex(probs);
    if (length(probs)>1)
    {
        String<long double> probsCopy = probs;
        erase(probsCopy,indexOfWinner);
        indexOfSecond = findMinIndex(probsCopy);
        enoughDistance = round(probsCopy[indexOfSecond]-probs[indexOfWinner]) > 10;
    }
    returnValue.genotype = genotypes[indexOfWinner];        
    returnValue.pValue = probs[indexOfWinner];    
    newGenotypeSet.insert(returnValue.genotype.i1);
    newGenotypeSet.insert(returnValue.genotype.i2);
    for (unsigned j=0; j<length(reads); ++j)
    {
        readToCheck = reads[j];
        if ((readToCheck.numOfRepeats != returnValue.genotype.i1) && (readToCheck.numOfRepeats != returnValue.genotype.i2))
            errorProbSum += readToCheck.pValue;
        float diff1 = fabs(returnValue.genotype.i1-readToCheck.numOfRepeats), diff2 = fabs(returnValue.genotype.i2-readToCheck.numOfRepeats);
        if (std::min(diff1,diff2)>=0.9)
            returnValue.fullMotifSlippageSum += readToCheck.pValue;
    }
    if (newGenotypeSet == currentGenotype)
        return Pair<GenotypeInfo, Pair<bool> >(returnValue, Pair<bool>(false,enoughDistance));
    else
        return Pair<GenotypeInfo, Pair<bool> >(returnValue, Pair<bool>(true,enoughDistance));    
}

void relabelReads(String<AttributeLine>& readsToRelabel, int start, int end, Pair<float> newGenotype, Marker marker)
{    
    int numOfReads = end - start; 
    for (unsigned i=start; i<end; ++i)
    {
        string pnId = readsToRelabel[i].PnId;        
        if ((fabs(readsToRelabel[i].numOfRepeats - newGenotype.i1)<=0.05) || (fabs(readsToRelabel[i].numOfRepeats - newGenotype.i2)<=0.05))
        {
            readsToRelabel[i].label = 1;
            ++markerToLabelsAndSlipp[marker].i1.p1.i1;
            if (numOfReads >= 10)
                markerToLabelsAndSlipp[marker].i1.p1.i2 += readsToRelabel[i].pValue;                
        }
        else 
        {
            if ((fabs(readsToRelabel[i].numOfRepeats - (newGenotype.i1 - 1))<=0.05) || (fabs(readsToRelabel[i].numOfRepeats - (newGenotype.i2 - 1))<=0.05))
            {
                readsToRelabel[i].label = 2;
                ++markerToLabelsAndSlipp[marker].i1.p2.i1;
                markerToStepSum[marker] += (float)marker.motif.size();
                if (numOfReads >= 10)
                    markerToLabelsAndSlipp[marker].i1.p2.i2 += readsToRelabel[i].pValue;                    
            }
            else
            { 
                readsToRelabel[i].label = -1;
                ++markerToLabelsAndSlipp[marker].i1.p3.i1;
                if (numOfReads >= 10)             
                    markerToLabelsAndSlipp[marker].i1.p3.i2 += readsToRelabel[i].pValue;
                float diff1, diff2;
                diff1 = fabs(readsToRelabel[i].numOfRepeats - newGenotype.i1);
                diff2 = fabs(readsToRelabel[i].numOfRepeats - newGenotype.i2);
                if (fmod(diff1,1.0)<0.05)
                {
                    if (fmod(diff2,1.0)<0.05)
                        markerToStepSum[marker] += round(std::min(diff1,diff2) * (float)marker.motif.size());
                    else 
                        markerToStepSum[marker] += round(diff1 * (float)marker.motif.size());
                }
                else 
                {
                    if (fmod(diff2,1.0)<0.05)
                        markerToStepSum[marker] += round(diff2 * (float)marker.motif.size());
                    else
                        markerToStepSum[marker] += round(std::min(diff1,diff2) * (float)marker.motif.size());
                }
            }
        }
    }
}

unsigned getChrLength(string chrom)
{
	if (chrom.compare("chr1") == 0)
		return 248956422;
	if (chrom.compare("chr2") == 0)
		return 242193529;
	if (chrom.compare("chr3") == 0)
		return 198295559;
	if (chrom.compare("chr4") == 0)
		return 190214555;
	if (chrom.compare("chr5") == 0)
		return 181538259;
	if (chrom.compare("chr6") == 0)
		return 170805979;
	if (chrom.compare("chr7") == 0)
		return 159345973;
	if (chrom.compare("chr8") == 0)
		return 145138636;
	if (chrom.compare("chr9") == 0)
		return 138394717;
	if (chrom.compare("chr10") == 0)
		return 133797422;
	if (chrom.compare("chr11") == 0)
		return 135086622;
	if (chrom.compare("chr12") == 0)
		return 133275309;
	if (chrom.compare("chr13") == 0)
		return 114364328;
	if (chrom.compare("chr14") == 0)
		return 107043718;
	if (chrom.compare("chr15") == 0)
		return 101991189;
	if (chrom.compare("chr16") == 0)
		return 90338345;
	if (chrom.compare("chr17") == 0)
		return 83257441;
	if (chrom.compare("chr18") == 0)
		return 80373285;
	if (chrom.compare("chr19") == 0)
		return 58617616;
	if (chrom.compare("chr20") == 0)
		return 64444167;
	if (chrom.compare("chr21") == 0)
		return 46709983;
	if (chrom.compare("chr22") == 0)
		return 50818468;
	else
		return 0;
}
 
//Write all sorts of info to the header of the vfc file I pass to the function
void makeVcfHeader(VcfStream& out, String<string> PnIds, string chrom)
{
    appendValue(out.header.sequenceNames, chrom);
    //Add IDs of all PNs to the header
    for (unsigned i = 0; i<length(PnIds); ++i)
        appendValue(out.header.sampleNames, PnIds[i]);
    unsigned chromLength = getChrLength(chrom);
    string contigString = "<ID=" + chrom + ", length=" + to_string(chromLength) + ">";
    //Complicated way of getting todays date
    time_t rawtime;
    tm* timeinfo;
    char buffer [80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,10,"%Y%m%d",timeinfo);
    string date = buffer;
    appendValue(out.header.headerRecords, VcfHeaderRecord("fileformat", "VCFv4.2"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("fileDate", date));
    appendValue(out.header.headerRecords, VcfHeaderRecord("source", "PopSTR"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("source_bin", "/odinn/tmp/bjarnih/Genotyping/160205/bin/computeReadAttributes"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("source_bin", "/odinn/tmp/bjarnih/Genotyping/160205/bin/computePnSlippage"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("source_bin", "/odinn/tmp/bjarnih/Genotyping/160205/bin/msGenotyper"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("reference", "/odinn/data/reference/Homo_sapiens-deCODE-hg38/Sequence/WholeGenomeFasta/genome.fa"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("contig", contigString));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=RefLen,Number=A,Type=Integer,Description=\"Length of the reference allele\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=Motif,Number=1,Type=String,Description=\"Microsatellite repeat motif\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">"));
}

//Clear a stringstream and a string, I use this a lot in fillRecordMarker/Pn so a function seemed appropriate
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
    stringstream ss;
    ss << record.beginPos+1;
    CharString str = ss.str();
    record.id = marker.chrom + ":" + toCString(str) + ":M";
    stringClear(ss,str);
    //Set smallest allele as reference allele
    std::set<float>::iterator smallestIt = allelesAtThisMarker.begin();
    float smallestAllele = *smallestIt;
    ss << smallestAllele;
    str = ss.str();
    record.ref = str;
    stringClear(ss,str);
    record.qual = record.MISSING_QUAL();
    record.info = "RefLen=";
    ss << refLength;
	str = ss.str();
	append(record.info,str);
	stringClear(ss,str);
    append(record.info,";");
    append(record.info,"Motif=");
    append(record.info,marker.motif);
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
void fillRecordPn(GenotypeInfo genotype, VcfRecord& record, MakeGenotypesRet genotypesAtThisMarker)
{
    unsigned indexOfCurrGt;
    while (genotypesAtThisMarker.genotypesSet.find(genotype.genotype)== genotypesAtThisMarker.genotypesSet.end() && length(genotype.genotypes)>1)
    {
        indexOfCurrGt = findMinIndex(genotype.pValues);
        erase(genotype.pValues,indexOfCurrGt);
        erase(genotype.genotypes,indexOfCurrGt);
        indexOfCurrGt = findMinIndex(genotype.pValues);
        genotype.genotype = genotype.genotypes[indexOfCurrGt];
        genotype.pValue = genotype.pValues[indexOfCurrGt];
    }
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
    if (genotype.genotype.i1 == 0)
    {
        for (unsigned i=0; i<length(genotypesAtThisMarker.genotypes); ++i)
                    append(gtInfo,"0,");
        eraseBack(gtInfo);    
    }
    else
    {
        Pair<float> genotypeToLookFor;
        Pair<float> genotypeToCompare;
        int index;
        long double numerator;
        long double denominator;
        int pl;
        for (unsigned i=0; i<length(genotypesAtThisMarker.genotypes); ++i)
        {
            index = -1;
            genotypeToLookFor = genotypesAtThisMarker.genotypes[i];
            if (genotypeToLookFor == genotype.genotype)
            {
                ss << 0;
                str = ss.str();
                append(gtInfo,str);
            }
            else
            {
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
                    pl = round(numerator-denominator);
                    pl = std::min(255,pl);
                    if (pl<0)
                        pl = 255;
                    ss << pl;
                    str = ss.str();
                    append(gtInfo,str);
                }
            }
            stringClear(ss,str);
            append(gtInfo,",");
        }
        eraseBack(gtInfo);
    }
    //When the genotype info string is ready I append it to the stringset of charstrings
    appendValue(record.genotypeInfos, gtInfo);
}

double computeAlleleDist(String<Pair<float> > genotypes, map<float,int> allelesToFreq, int pnsAtMarker)
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

//Count number of words in a sentence, use to parse input from attribute file
Pair<int, String<string> > countNumberOfWords(string sentence)
{
    int numberOfWords = 0;
    String<string> words;
    resize(words, 11);
    int currentWordLength;
    
    if (!isspace(sentence[0])) 
    {
        numberOfWords++;
        words[0] = sentence[0];
    }

    for (unsigned i = 1; i < sentence.length(); i++) 
    {
        if ((!isspace(sentence[i])) && (isspace(sentence[i-1]))) 
        {
            numberOfWords++;
            words[numberOfWords-1] = sentence[i];
        }
        else
        {
            if (!isspace(sentence[i]))            
                words[numberOfWords-1].push_back(sentence[i]);            
        }
    }
    
    resize(words, numberOfWords);
    return Pair<int, String<string> >(numberOfWords, words);
} 

void readPnSlippage(ifstream& pnSlippageFile)
{
    string PnId;
    int nMarkers;
    double currPnSlipp;
    while (!pnSlippageFile.eof())
    {
        pnSlippageFile >> PnId;
        pnSlippageFile >> currPnSlipp;
        pnToSize[PnId]= currPnSlipp;
    }
    pnSlippageFile.close();
}

inline bool exists(const std::string& name) 
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

Pair<double, int> estimateSlippage(String<string> PnIds, map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype, Marker marker, int currItNum)
{
    int nMissing = 0, nAvailable;
    vector<double> weights;
    vector<double> slippFragments;
    double currPnSlipp, currPvalSum = 0, weightSum = 0, fullMotifSlippageSum = 0, result, currMarkSlipp;
    if (currItNum == 1)
    {
        for (unsigned i = 0; i<length(PnIds); ++i)
        {
            if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
                continue;
            fullMotifSlippageSum += PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].fullMotifSlippageSum;
            currPvalSum += PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum;
        }
        currMarkSlipp = 0.5*(fullMotifSlippageSum/currPvalSum);
    }
    else
        currMarkSlipp = markerToLabelsAndSlipp[marker].i2;
    for (unsigned i = 0; i<length(PnIds); ++i)
    {                
        if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
        {
            weights.push_back(0);
            ++nMissing;
            continue;
        }
        if (pnToSize[PnIds[i]] == 0)
            currPnSlipp = 0.001;
        else 
            currPnSlipp = pnToSize[PnIds[i]];
        currPvalSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum;
        weights.push_back(currPvalSum/((currPnSlipp+currMarkSlipp)*(1-(currPnSlipp+currMarkSlipp))));
    }
    weightSum = accumulate(weights.begin(),weights.end(),0.0);
    nAvailable = length(PnIds) - nMissing;
    for (unsigned i = 0; i<length(PnIds); ++i)
    {                        
        if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
        {
            slippFragments.push_back(0);            
            continue;
        }
        if (pnToSize[PnIds[i]] == 0)
            currPnSlipp = 0.001;
        else 
            currPnSlipp = pnToSize[PnIds[i]];
        fullMotifSlippageSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].fullMotifSlippageSum;
        currPvalSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum;
        slippFragments.push_back((weights[i]/weightSum)*((fullMotifSlippageSum/currPvalSum)-currPnSlipp));
    }
    result = std::max(0.0,accumulate(slippFragments.begin(),slippFragments.end(),0.0));
    cout << "Slippage rate: " << result << endl;
    return Pair<double, int>(result, nAvailable);
}

int main(int argc, char const ** argv)
{   
    //Check arguments.
    MsGenotyperOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
	    return res == seqan::ArgumentParser::PARSE_ERROR;
    
    //Assign parameters     
    int startCoord = options.startCoordinate, endCoord = options.endCoordinate;
    CharString attributePath = options.attDirChromNum, intervalIndex = options.intervalIndex, markerSlippageFile = options.markerSlippageFile;
    ifstream pnSlippageFile(toCString(options.pnSlippageFile));
    if(pnSlippageFile.fail())        
        cout << "Unable to locate pnSlippageFile @ " << options.pnSlippageFile << endl;
    append(markerSlippageFile, "_");
    append(markerSlippageFile, intervalIndex);
    ofstream markerSlippageOut(toCString(markerSlippageFile));
    string PnId, chrom, motif, nextWord, refRepSeq;
    String<string> PnIds;
    std::set<Marker> markers;    
    int start, end, numberOfReads;    
    float refRepeatNum, winner, second;     
    bool enoughReads = true;
    //Read the slippage rate for all PNs into the pnToSize map.
    readPnSlippage(pnSlippageFile);
    cout << "Finished reading pnSlipp." << endl;
    //Map from marker to all reads covering it 
    map<Marker, String<AttributeLine> > mapPerMarker;
    //Count the total number of alleles in the population for each marker -- by checking the size of the set
    map<Marker, std::set<float> > markerToAlleles;
    //Map to store current genotype(and lots of other things) of a person for each marker maps from pnId and Marker-struct to GenotypeInfo struct
    map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype;    
    
    //Open vcf stream and make header if the estimateMarkerSlippage switch is off
    VcfStream out;
    CharString outputDirectory = options.vcfOutputDirectory;
    CharString outputFileName = options.vcfFileName;
    append(outputDirectory, "/");
    append(outputDirectory, outputFileName);
    append(outputDirectory, "_");
    append(outputDirectory, intervalIndex);
    append(outputDirectory, ".vcf");
    bool outOk = open(out,toCString(outputDirectory), VcfStream::WRITE);
    if (!outOk)
    {
        cerr << "Cannot create vcf file at: " << outputDirectory << endl;
        return 1;  
    }
    
    Marker marker;
    string nextLine;
    AttributeLine currentLine;    
    Pair<int, String<string> > numberOfWordsAndWords;
    
    //Iterate over all Pns I have slippage for and read from attribute files in the given interval
    map<string, double>::const_iterator pnEnd = pnToSize.end();
    int nProcessedPns = 0;
    for (map<string, double>::iterator pnStart = pnToSize.begin(); pnStart != pnEnd; ++pnStart)
    {
        PnId = pnStart->first;
        append(attributePath, "/");
        append(attributePath, PnId);
        ifstream attributeFile(toCString(attributePath));
        if (attributeFile.fail())
        {
            cout << "Unable to locate attribute file for " << PnId << " at " << attributePath << endl;
            attributePath = options.attDirChromNum;
            continue;
        }
        ++nProcessedPns;
        if (nProcessedPns % 1000==0)
            cout << "Working on pn number: " << nProcessedPns << endl;
        while (!attributeFile.eof())
        {
            getline (attributeFile,nextLine);
            if (nextLine.length() == 0) 
                continue;
            numberOfWordsAndWords = countNumberOfWords(nextLine);
            if (numberOfWordsAndWords.i1 == 1)
            {
                PnId = nextLine;
                appendValue(PnIds,PnId);         
            }
            if (numberOfWordsAndWords.i1 == 9) 
            {                      
                marker.chrom = numberOfWordsAndWords.i2[0];
                marker.start = lexicalCast<int>(numberOfWordsAndWords.i2[1]);                          
                marker.end = lexicalCast<int>(numberOfWordsAndWords.i2[2]);
                marker.motif = numberOfWordsAndWords.i2[3];
                marker.refRepeatNum = lexicalCast<float>(numberOfWordsAndWords.i2[4]);
                numberOfReads = lexicalCast<int>(numberOfWordsAndWords.i2[5]);
                marker.refRepSeq = numberOfWordsAndWords.i2[6];
                //If I am in front of the interval, I move to the next marker.
                if (marker.start < startCoord)
                {
                    for (unsigned i = 0; i < numberOfReads; ++i)
                    {
                        getline (attributeFile,nextLine);
                    }
                    continue;    
                }
                //If I have passed the interval, I move on to the next pn.
                if (marker.start > endCoord)
                    break;
                //Otherwise the marker is in the interval and I process it.
                markers.insert(marker);
                winner = lexicalCast<float>(numberOfWordsAndWords.i2[7]);
                second = lexicalCast<float>(numberOfWordsAndWords.i2[8]);
                markerToAlleles[marker].insert(winner);
                markerToAlleles[marker].insert(second);
                
            }
            if (numberOfWordsAndWords.i1 == 11) 
            {
                if (numberOfReads < 10)
                {
                    PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].pValueSum = 0;
                    enoughReads = false;
                }
                for (unsigned i = 0; i < numberOfReads; ++i)
                {
                    if (i == 0)
                        currentLine = parseNextLine(winner, second, attributeFile, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, true, enoughReads);
                    else 
                        currentLine = parseNextLine(winner, second, attributeFile, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, false, enoughReads);
                    appendValue(mapPerMarker[marker],currentLine);
                }
                enoughReads = true;                                    
            }
            if (numberOfWordsAndWords.i1 != 1 && numberOfWordsAndWords.i1 != 9 && numberOfWordsAndWords.i1 != 11) 
                cerr << "Format error in attribute file!" << endl;
        }            
        attributePath = options.attDirChromNum;
        attributeFile.close();
    }
    chrom = marker.chrom;
    cout << "Reading data from input complete." << endl;    
    
    //Make header for vcf file
    makeVcfHeader(out, PnIds, chrom);
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
    //Map to store alleles present in each individual at a given marker. Is cleared for each marker.
    map<string, std::set<float> > PnToAlleles;
    //Store set of all alleles in population for marker being considered and use it to generate all possible genotypes and compute their p-values
    std::set<float> allelesAtMarker;
    std::set<float> allelesToConsider;
    String<Pair<float> > genotypesToConsider;
    float currAllele;
    //Map to store a map from alleles to their frequencies in the population for each marker, used for estimating the probability that the distance between alleles at the marker is an integer and number of PNs at each marker.
    map<Marker, Pair<map<float,int>,int> > markerToAlleleFreqs;
    int z = 1, nAvailable;
    
    //Stuff for vcf file
    GenotypeInfo genotype;
    string thisPn;
    Marker thisMarker;
    VcfRecord record;
    stringstream ss;
    CharString str;
    MakeGenotypesRet genotypesAtThisMarker;
    double alleleDistance;
    Pair<double, int> slippAndNavail;
    
    //Loop over map from Marker to string<AttributeLine> and train model for each marker and use it to determine genotype
    map<Marker, String<AttributeLine> >::iterator itEnd = mapPerMarker.end();
    for (map<Marker, String<AttributeLine> >::iterator it = mapPerMarker.begin(); it != itEnd; ++it)
    {        
        thisMarker = it->first;        
        changedRatio=1;
        PnsAtMarker = 1;
        int loops = 0;
        String<AttributeLine>& currentMarker = it->second;
        probBig.l = length(currentMarker);
        //convergence condition = less than 0.5% of Pns get a new genotype 
        cout << "Starting marker number: " << z << " with start coordinate: " << thisMarker.start << endl;
        while (changedRatio > 0.005)
        {
            //Backup if convergence has not been reached within 10 iterations 
            if (loops > 10)
                break;
            ++loops;
            /*cout << "Iteration number: " << loops << endl;            
            cout << "Average step size modulo the period: fmod(" << markerToStepSum[thisMarker] << "/" << length(currentMarker) << ",1.0)" << endl;
            cout << "P for geometric distribution is: " << "1.0/(" << fmod(markerToStepSum[thisMarker]/(float)length(currentMarker),1.0) << "+1.0) = " << 1.0/(fmod(markerToStepSum[thisMarker]/(float)length(currentMarker),1.0)+1.0) << endl;*/
            double geomP = 1/(fmod(markerToStepSum[thisMarker]/(float)length(currentMarker),1.0)+1);
            allelesAtMarker = markerToAlleles[thisMarker];                        
            numOfAlleles = allelesAtMarker.size();
            //cout << "Number of alleles at the marker: " << numOfAlleles << endl;
            markerToAlleles[thisMarker].clear();
            markerToAlleleFreqs[thisMarker].i1.clear();
            //estimate the marker slippage - have to subtract PN (and allele length slippage)
            //Estimate marker slippage and update in markerToLabelsAndSlipp map            
            slippAndNavail = estimateSlippage(PnIds, PnAndMarkerToGenotype, it->first, loops);
            markerToLabelsAndSlipp[it->first].i2 = std::max((double)0.0,slippAndNavail.i1);
            nAvailable = slippAndNavail.i2;                            
            //Reads with label 2 are not included in training
            prob.l = length(currentMarker) - markerToLabelsAndSlipp[it->first].i1.p2.i1;
            //Now I "nullSet" the markerToLabelsAndSlipp map for the marker I am looking at so I can update it when I relabel the reads 
            markerToLabelsAndSlipp[it->first].i1.p1 = Pair<int,double>(0,0);
            markerToLabelsAndSlipp[it->first].i1.p2 = Pair<int,double>(0,0);
            markerToLabelsAndSlipp[it->first].i1.p3 = Pair<int,double>(0,0);
            //Nullset the averageStepSize map 
            markerToStepSum[it->first] = 0.0;          
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
                        PnId = currentLine.PnId;
                        ++PnsAtMarker;
                    } 
                    PnToAlleles[currentLine.PnId].insert(currentLine.numOfRepeats);
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
            const char *error_msg;
            error_msg = check_parameter(&prob,&param);
            if (error_msg != NULL)
                cout << "Error message: " << error_msg << endl;
            //Train logistic regression model 
            model_ = train(&prob, &param);
            prob_estimates = (double *) malloc(2*sizeof(double));
            Pair<GenotypeInfo, Pair<bool> > changed;            
            PnId = currentMarker[0].PnId;
            String<AttributeLine> reads;
            for (unsigned i = 0; i < length(currentMarker); ++i) 
            {                
                currentLine = currentMarker[i];
                if (currentLine.PnId != PnId)
                {
                    //Need to reset this set for each PN so I can add alleles from their reads
                    allelesToConsider = allelesAtMarker;
                    //Add alleles in reads from this PN to the set of alleles present at the marker.
                    std::set<float>::iterator pnAllsEnd = PnToAlleles[PnId].end();
                    for (std::set<float>::iterator pnAlls = PnToAlleles[PnId].begin(); pnAlls!=pnAllsEnd; ++pnAlls)
                    {
                        currAllele = *pnAlls;
                        allelesToConsider.insert(currAllele);
                    }
                    genotypesToConsider = makeGenotypes(allelesToConsider).genotypes;
                    //make decision about genotype for PnId at the current marker.        
                    changed = determineGenotype(reads, markerToLabelsAndSlipp[it->first].i2+pnToSize[PnId], genotypesToConsider, numOfAlleles, it->first.motif.size(), geomP);                   
                    if (changed.i2.i1)
                        ++updatedPns;
                    relabelReads(currentMarker, i-length(reads), i, changed.i1.genotype, it->first);
                    PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1; 
                    if (changed.i2.i2 || !changed.i2.i2)
                    {
                        markerToAlleles[it->first].insert(changed.i1.genotype.i1);
                        markerToAlleles[it->first].insert(changed.i1.genotype.i2);
                    }
                    ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
                    ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];
                    //If I am estimating the marker slippage then I should update map from Pn to labels. (before I update PnId to currentLine.PnId)
                    if (loops>1)
                        eraseBack(pnToLabels[PnId]);                       
                    append(pnToLabels[PnId],Pair<float>(changed.i1.genotype.i1, changed.i1.genotype.i2));
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
            //Make decision about genotype for last PnId at the current marker.
            //Need to reset this set for last PN so I can add alleles from his reads
            allelesToConsider = allelesAtMarker;
            //Add alleles in reads from last PN to the set of alleles present at the marker.
            std::set<float>::iterator pnAllsEnd = PnToAlleles[PnId].end();
            for (std::set<float>::iterator pnAlls = PnToAlleles[PnId].begin(); pnAlls!=pnAllsEnd; ++pnAlls)
            {
                currAllele = *pnAlls;
                allelesToConsider.insert(currAllele);
            }
            genotypesToConsider = makeGenotypes(allelesToConsider).genotypes;
            changed = determineGenotype(reads, markerToLabelsAndSlipp[it->first].i2+pnToSize[PnId], genotypesToConsider, numOfAlleles, it->first.motif.size(), geomP);
            if (changed.i2.i1)
                ++updatedPns;
            relabelReads(currentMarker, length(currentMarker)-length(reads), length(currentMarker), changed.i1.genotype, it->first);
            PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;
            if (changed.i2.i2 || !changed.i2.i2)
            {
                markerToAlleles[it->first].insert(changed.i1.genotype.i1);
                markerToAlleles[it->first].insert(changed.i1.genotype.i2);
            }
            ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
            ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];
            if (loops>1)
                eraseBack(pnToLabels[PnId]);
            append(pnToLabels[PnId],Pair<float>(changed.i1.genotype.i1, changed.i1.genotype.i2));
            changedRatio = (float)updatedPns/(float)PnsAtMarker;
            cout << "Changed: " << changedRatio << endl;
        }
        //Write vcf record for the marker I just finished.    
        //Check the number of alleles in the population, should be 2 or higher to be considered polymorphic.
        if (markerToAlleles[thisMarker].size() < 2)
        {
            cout << "Not enough alleles at marker " << z << endl;
            ++z;
            continue;      
        }
        //Make a Pair<std::set<Pair<float> > String<Pair<float> > > which contains a list and set of genotypes
        genotypesAtThisMarker = makeGenotypes(markerToAlleles[thisMarker]);
        //Compute abs(allele1-allele2)*allele1Freq*allele2Freq for all genotypes and return average of those, estimate of distance between alleles.
        alleleDistance = computeAlleleDist(genotypesAtThisMarker.genotypes, markerToAlleleFreqs[it->first].i1, PnsAtMarker);
        //First fill marker specific fields of vcfRecord
        record = fillRecordMarker(thisMarker, markerToAlleles[thisMarker]);  
        //Loop over Pns and fill in PN specific fields of vcfRecord for each PN
        for (unsigned i = 0; i<length(PnIds); ++i)
        {
            thisPn = PnIds[i];
            //If a decision has been made for thisPn at thisMarker I add it to the genotypeInfos stringSet
            if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(thisPn, thisMarker)) != 0))
            {
                genotype = PnAndMarkerToGenotype[Pair<string,Marker>(thisPn, thisMarker)];                   
                fillRecordPn(genotype, record, genotypesAtThisMarker);
                PnAndMarkerToGenotype.erase(Pair<string,Marker>(thisPn, thisMarker));
            }
            //If a decision has not been made I add a CharString with no decision(0:0,0,0,0....etc) to the set to maintain order of Pns vs genotypeInfos in output
            else 
            {         
                CharString gtInfo = "0/0:";       
                for (unsigned i=0; i<length(genotypesAtThisMarker.genotypes); ++i)
                    append(gtInfo,"0,");
                eraseBack(gtInfo);
                appendValue(record.genotypeInfos, gtInfo);
            }
        }               
        ss << markerToAlleles[thisMarker].size();
        str = ss.str();
        record.filter = ".";
        stringClear(ss,str);
        //After adding info for all PNs at thisMarker I write the record to the vcf output file
        if (writeRecord(out, record) != 0)
            cerr << "ERROR: Problem writing VCF-output file." << endl;
        clear(record);
        markerToAlleleFreqs[it->first].i2 = PnsAtMarker;
        PnToAlleles.clear();
        markerSlippageOut << thisMarker.chrom << "\t" << thisMarker.start << "\t" << thisMarker.end << "\t" << thisMarker.motif << "\t" << thisMarker.refRepeatNum << "\t" << thisMarker.refRepSeq << "\t" << setprecision(4) << fixed << markerToLabelsAndSlipp[it->first].i2 << "\t" << nAvailable<< endl;
        cout << thisMarker.start << " totalSlipp: " << setprecision(4) << fixed << markerToLabelsAndSlipp[it->first].i2 << endl;
        cout << "Finished marker number: " << z << endl;
        ++z;                 
    }
    cout << "Finished determining genotypes" << endl;
    return 0;
}
