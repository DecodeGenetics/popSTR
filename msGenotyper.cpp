#include <iostream>
#include <set>
#include <map>
#include <string>
#include <ctime>
#include <math.h>
#include <vector>
#include <numeric>
#include <sstream>
#include <fstream>
#include <ios>
#include <sys/stat.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/vcf_io.h>
#include <seqan/sequence.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
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

namespace io = boost::iostreams;
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
map<Marker, Pair<Pair<LabelProps,double> ,model* > > markerToSizeAndModel;
//Store stutter rate and nAlleles for each marker
map<Marker, Pair<int, double> > markerToNallelesAndStutter;
//Stores the slippage rate value for each PN - read from input file.
map<string, Pair<double, int> > pnToSize;
//Maps from pnId to reads and labels for one marker, is cleared for each marker
map<string, String<Pair<float> > > pnToLabels;
//Map from marker to (average stepsize)%period
map<Marker, double> markerToStepSum;
//Map from marker to alleles that I have found more than 3 reads of with freq > 25% in a single PN
map<Marker, std::set<float> > markerToTrueAlleles;
//for debugging purposes
bool useGeom = true;

//Parameter, problem and model structs to use in training of logistic regression model and computing pValues
parameter param;
problem prob;
problem probBig; //problem structure including reads with label = 2 to use for predict function
model* model_;
double bias = -1;

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

double getPval(Marker marker, AttributeLine currentLine)
{
    double predict_label;
    model* model_ = markerToSizeAndModel[marker].i2;
    double *prob_estimates = (double *) malloc(2*sizeof(double));
    problem prob_;
    prob_.bias = -1;
    prob_.l = 1;
    prob_.n = 9;
    prob_.x = (feature_node **) malloc(prob_.l * sizeof(feature_node *));
    prob_.x[0] = (feature_node *) malloc(10 * sizeof(feature_node));
    fillProblemX(0, currentLine, prob_);
    predict_label = predict_probability(model_,prob_.x[0],prob_estimates);
    return prob_estimates[0];
}

//Parses one line from attribute file by filling up and returning an AttributeLine, also initializes markerToSizeAndModel map using the labels
AttributeLine parseNextLine(float winner, float second, io::filtering_istream& attributeFile, Marker& marker, string PnId, map<Pair<string,Marker>, GenotypeInfo>& PnAndMarkerToGenotype, String<string> firstLine, bool useFirstLine, bool useModelAndLabels, bool enoughReads)
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
    if (!useModelAndLabels)
        currentLine.pValue = 0.95;
    else
        currentLine.pValue = getPval(marker, currentLine);
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
        ++markerToSizeAndModel[marker].i1.i1.p1.i1;
        if (enoughReads)
            markerToSizeAndModel[marker].i1.i1.p1.i2 += currentLine.pValue;
    }
    else
    {
        if ((fabs(currentLine.numOfRepeats - (winner - 1)) <= 0.05) || (fabs(currentLine.numOfRepeats - (second - 1)) <= 0.05))
        {
            currentLine.label = 2;
            ++markerToSizeAndModel[marker].i1.i1.p2.i1;
            if (enoughReads)
                markerToSizeAndModel[marker].i1.i1.p2.i2 += currentLine.pValue;
            markerToStepSum[marker] += (float)marker.motif.size();
        }
        else
        {
            currentLine.label = -1;
            ++markerToSizeAndModel[marker].i1.i1.p3.i1;
            if (enoughReads)
                markerToSizeAndModel[marker].i1.i1.p3.i2 += currentLine.pValue;
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

MakeGenotypesRet makeGenotypes(std::set<float> alleles, float refAllele)
{
    String<Pair<float> > genotypes;
    String<float> alleleString;
    std::set<Pair<float> > genotypeSet;
    std::set<float>::reverse_iterator allelesBegin = alleles.rend();
    float closestToRef = *alleles.rbegin();
    float minDistToRef = fabs(closestToRef-refAllele);
    unsigned refIdx = 0;
    for (std::set<float>::reverse_iterator alleleIt = alleles.rbegin(); alleleIt!=allelesBegin; ++alleleIt)
    {
        appendValue(alleleString, *alleleIt);
        if (fabs(*alleleIt-refAllele)<minDistToRef)
        {
            closestToRef = *alleleIt;
            minDistToRef = fabs(*alleleIt-refAllele);
            refIdx = length(alleleString)-1;
        }
    }
    //Put reference allele at end of alleleString
    erase(alleleString, refIdx);
    appendValue(alleleString, closestToRef);
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

float dpois(int step, float lambda) {
  if (step < 0)
      return 0;
  float p = exp(-1*lambda);
  for (int i = 0; i < step; i++) {
    p = p*lambda;
    p = p/(i+1);
  }
  return p;
}

Pair<GenotypeInfo, bool> determineGenotype(String<AttributeLine> reads, double markerSlippage, String<Pair<float> > genotypes, int numberOfAlleles, int motifLength, double psucc)
{
    GenotypeInfo returnValue;
    returnValue.pValueSum = 0;
    returnValue.genotypes = genotypes;
    returnValue.numOfReads = length(reads);
    Pair<float> genotypeToCheck;
    AttributeLine readToCheck;
    String<long double> probs;
    double errorProbSum = 0;
    resize(probs, length(genotypes));
    bool isHomo, enoughDistance = true;
    float diff, diff2, lambda = std::max((double)0.001,markerSlippage), posNegSlipp = 1, posNegSlipp2 = 1;
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
            if (isHomo)
            {
                if (readToCheck.numOfRepeats < genotypeToCheck.i1)
                    posNegSlipp = 0.8;
                if (readToCheck.numOfRepeats > genotypeToCheck.i1)
                    posNegSlipp = 0.2;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                if (useGeom)
                    probs[i] += -(double)10*log10(readToCheck.pValue * dgeom(static_cast<int>(roundf((diff-(float)floor(diff))*motifLength)), psucc) * dpois(floor(diff), lambda) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
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
                    probs[i] += -(double)10*log10(readToCheck.pValue * 0.5 * (dgeom(static_cast<int>(roundf((diff-(float)floor(diff))*motifLength)), psucc) * dpois(floor(diff), lambda) * posNegSlipp + dgeom(static_cast<int>(roundf((diff2-(float)floor(diff2))*motifLength)), psucc) * dpois(floor(diff2), lambda) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
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
    for (unsigned j=0; j<length(reads); ++j)
    {
       readToCheck = reads[j];
       if ((readToCheck.numOfRepeats != returnValue.genotype.i1) && (readToCheck.numOfRepeats != returnValue.genotype.i2))
        errorProbSum += readToCheck.pValue;
    }
    return Pair<GenotypeInfo, bool>(returnValue, enoughDistance);
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

//Write all sorts of info to the header of the vfc file.
void makeVcfHeader(VcfStream& out, String<string> PnIds, string chrom)
{
    appendValue(out.header.sequenceNames, chrom);
    //Add IDs of all PNs to the header
    for (unsigned i = 0; i<length(PnIds); ++i)
        appendValue(out.header.sampleNames, PnIds[i]);
    unsigned chromLength = getChrLength(chrom);
    string contigString = "<ID=" + chrom + ",length=" + to_string((long long unsigned int)chromLength) + ">";
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
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">"));
}

//Clear a stringstream and a string, I use this a lot in fillRecordMarker/Pn, so a function seemed appropriate
void stringClear(stringstream& ss, CharString& str)
{
    ss.str("");
    ss.clear();
    clear(str);
}

Pair<int> findGenotypeIndices(GenotypeInfo& genotype, std::set<float>& allelesAtMarker)
{
    Pair<int> gtIdxs = Pair<int>(0,0);
    unsigned idx = 1;
    float currentAllele;
    std::set<float>::iterator allEnd = allelesAtMarker.end();
    for (std::set<float>::iterator allIt = allelesAtMarker.begin(); allIt!=allEnd; ++allIt)
    {
        currentAllele = *allIt;
        if (currentAllele == genotype.genotype.i1)
            gtIdxs.i1 = idx;
        if (currentAllele == genotype.genotype.i2)
            gtIdxs.i2 = idx;
        idx++;
    }
    return gtIdxs;
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
    //Get ref seq using crop_fasta and set as record.ref
    string refSeq;
    FILE * stream;
    char buffer[256];
    stringstream cmdStream;
    cmdStream << "/odinn/users/snaedisk/bin/crop_fasta /odinn/data/reference/Homo_sapiens-deCODE-hg38/Sequence/WholeGenomeFasta/genome.fa " << marker.chrom << ":" << toCString(str) << "-" << marker.end;
    string cmd = cmdStream.str();
    cmd.append(" 2>&1");
    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
        if (fgets(buffer, 256, stream) != NULL)
        {
            buffer[strcspn(buffer, "\n")] = 0;
            refSeq.append(buffer);
        }
        pclose(stream);
    }
    //Set reference allele and delete it from allele set
    stringClear(ss,str);
    record.ref = refSeq;
    allelesAtThisMarker.erase(marker.refRepeatNum);
    record.qual = record.MISSING_QUAL();
    record.info = "RefLen=";
    ss << refLength;
	str = ss.str();
	append(record.info,str);
	stringClear(ss,str);
    append(record.info,";");
    append(record.info,"Motif=");
    append(record.info,marker.motif);
    record.format = "GT:AD:DP:PL";
    //Loop over all alternative alleles and add them to ALT field in record
    float currentAllele;
    std::set<float>::iterator allEnd = allelesAtThisMarker.end();
    for (std::set<float>::iterator allIt = allelesAtThisMarker.begin(); allIt!=allEnd; ++allIt)
    {
        append(record.alt, "<");
        currentAllele = *allIt;
        ss << currentAllele;
        str = ss.str();
        append(record.alt, str);
        append(record.alt, ">,");
        stringClear(ss,str);
    }
    if (allelesAtThisMarker.size() > 0)
        eraseBack(record.alt);
    else
        record.alt = ".";
    return record;
}

//Fill up the PN-specific fields in the VCF-record.
void fillRecordPn(GenotypeInfo genotype, VcfRecord& record, MakeGenotypesRet genotypesAtThisMarker, std::set<float>& allelesAtMarker, Marker& marker)
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
    float refAllele = marker.refRepeatNum;
    //Remove ref allele from set.
    allelesAtMarker.erase(refAllele);
    stringstream ss;
    CharString str = ss.str();
    CharString gtInfo; //First I make the string containing the genotype info
    Pair<int> gtIdxs = findGenotypeIndices(genotype, allelesAtMarker);
    ss << gtIdxs.i1;
    str = ss.str();
    append(gtInfo,str);
    stringClear(ss,str);
    append(gtInfo,"/");
    ss << gtIdxs.i2;
    str = ss.str();
    append(gtInfo,str);
    stringClear(ss,str);
    append(gtInfo,":");
    //Add alleleDepth count for the ref allele
    if (genotype.alleleToFreq.count(refAllele)==0)
        append(gtInfo,"0,"); //No reads supporting ref allele
    else
    {
        ss << genotype.alleleToFreq[refAllele];
        str = ss.str();
        append(gtInfo,str);
        stringClear(ss,str);
        append(gtInfo,",");
    }
    //Remove ref allele from set before looping over
    allelesAtMarker.erase(refAllele);
    //Add allele depth counts for other alleles
    for (auto const & allele : allelesAtMarker)
    {
        if (genotype.alleleToFreq.count(allele) == 0)
            append(gtInfo,"0,");
        else
        {
            ss << genotype.alleleToFreq[allele];
            str = ss.str();
            append(gtInfo,str);
            stringClear(ss,str);
            append(gtInfo,",");
        }
    }
    //put ref allele back in set after looping over
    allelesAtMarker.insert(refAllele);
    unsigned possibleGenotypes = 0;
    for (unsigned i=0; i<allelesAtMarker.size(); ++i)
    {
        possibleGenotypes += i+1;
    }
    if (length(genotypesAtThisMarker.genotypes) != possibleGenotypes)
        cerr << "ATTENTION: number of genotypes considered is wrong. It should be: " << possibleGenotypes << " but it is:" << length(genotypesAtThisMarker.genotypes) << endl;
    eraseBack(gtInfo);
    append(gtInfo,":");
    //Add readDepth to gtInfo
    ss << genotype.numOfReads;
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

//Read in data from markerSlippageFile
void readMarkerSlippage(CharString markerSlippageFile, map<Marker, Pair<Pair<LabelProps,double>, model* > >& markerToSizeAndModel, int startCoord, int endCoord)
{
    ifstream markerSlippageIn(toCString(markerSlippageFile));
    double currMarkSlipp, currMarkStutt;
    int nPns, nAlleles;
    Marker currMarker;
    while (true)
    {
        markerSlippageIn >> currMarker.chrom;
        markerSlippageIn >> currMarker.start;
        markerSlippageIn >> currMarker.end;
        markerSlippageIn >> currMarker.motif;
        markerSlippageIn >> currMarker.refRepeatNum;
        markerSlippageIn >> currMarker.refRepSeq;
        markerSlippageIn >> currMarkSlipp;
        markerSlippageIn >> nPns;
        markerSlippageIn >> nAlleles;
        markerSlippageIn >> currMarkStutt;
        if (markerSlippageIn.eof())
            break;
        if (currMarker.start < startCoord)
            continue;
        if (currMarker.start > endCoord)
            break;
        markerToNallelesAndStutter[currMarker].i1 = nAlleles;
        markerToNallelesAndStutter[currMarker].i2 = currMarkStutt;
        markerToSizeAndModel[currMarker].i1.i2 = currMarkSlipp;
        markerToSizeAndModel[currMarker].i2 = NULL;
    }
    cout << "Finished reading marker slippage." << endl;
    markerSlippageIn.close();
}

void readPnSlippage(ifstream& pnSlippageFile)
{
    string PnId;
    int nMarkers;
    double currPnSlipp;
    while (true)
    {
        pnSlippageFile >> PnId;
        pnSlippageFile >> currPnSlipp;
        pnSlippageFile >> nMarkers;
        if (pnSlippageFile.eof() || PnId.length() == 0)
            break;
        pnToSize[PnId].i1= currPnSlipp;
        pnToSize[PnId].i2 = nMarkers;
    }
    cout << "Finished reading pn Slippage, number of pns:" << pnToSize.size() << endl;
    pnSlippageFile.close();
}

inline bool exists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

Pair<double, int> estimateSlippage(String<string> PnIds, map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype, Marker marker, map<Marker, Pair<Pair<LabelProps,double> ,model* > > currMarkerToSizeAndModel, int currItNum)
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
        currMarkSlipp = currMarkerToSizeAndModel[marker].i1.i2;
    for (unsigned i = 0; i<length(PnIds); ++i)
    {
        if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
        {
            weights.push_back(0);
            ++nMissing;
            continue;
        }
        if (pnToSize[PnIds[i]].i1 == 0)
            currPnSlipp = 0.001;
        else
            currPnSlipp = pnToSize[PnIds[i]].i1;
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
        if (pnToSize[PnIds[i]].i1 == 0)
            currPnSlipp = 0.001;
        else
            currPnSlipp = pnToSize[PnIds[i]].i1;
        fullMotifSlippageSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].fullMotifSlippageSum;
        currPvalSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum;
        slippFragments.push_back((weights[i]/weightSum)*((fullMotifSlippageSum/currPvalSum)-currPnSlipp));
    }
    result = std::max(0.0,accumulate(slippFragments.begin(),slippFragments.end(),0.0));
    cout << "Slippage rate: " << result << endl;
    return Pair<double, int>(result, nAvailable);
}

void appendChromName(CharString& dir, CharString chromName)
{
    append(dir, "/");
    append(dir, chromName);
    append(dir, "/");
}

Pair<bool> genotypeIsConfident(GenotypeInfo& genotypeInfo)
{
    bool A1ok = false, A2ok = false;
    if (genotypeInfo.alleleToFreq[genotypeInfo.genotype.i1]>=3 && (double)genotypeInfo.alleleToFreq[genotypeInfo.genotype.i1]/(double)genotypeInfo.numOfReads >=0.25)
        A1ok = true;
    if (genotypeInfo.alleleToFreq[genotypeInfo.genotype.i2]>=3 && (double)genotypeInfo.alleleToFreq[genotypeInfo.genotype.i2]/(double)genotypeInfo.numOfReads >=0.25)
        A2ok = true;
    return Pair<bool>(A1ok, A2ok);
}

int main(int argc, char const ** argv)
{
    //Check arguments.
    if (argc != 10 && argc != 12)
    {
        cerr << "USAGE: " << argv[0] << " attDir PN-slippageFile startCoordinate endCoordinate intervalIndex markerSlippageDir modelAndLabelDir iterationNumber chromosomeName vcfOutputDirectory vcfFileName \n";
        return 1;
    }

    //Parse parameters
    int startCoord = lexicalCast<int>(argv[3]), endCoord = lexicalCast<int>(argv[4]), currItNum = lexicalCast<int>(argv[8]), prevItNum = lexicalCast<int>(argv[8]) - 1;
    CharString attributePath = argv[1], pnSlippagePath = argv[2], intervalIndex = argv[5], markerSlippageDir = argv[6], modelAndLabelDir = argv[7], currItNumStr = argv[8], prevItNumStr = to_string((long long int)prevItNum), chromName = argv[9];
    append(pnSlippagePath, prevItNumStr);
    append(attributePath, "/attributes");
    appendChromName(attributePath, chromName);
    appendChromName(modelAndLabelDir, chromName);
    appendChromName(markerSlippageDir, chromName);
    CharString modelAndLabelDirBase = modelAndLabelDir;
    CharString markerSlippageDirBase = markerSlippageDir;
    CharString attributePathBase = attributePath;
    ifstream pnSlippageFile(toCString(pnSlippagePath));
    ofstream markerSlippageOut; //Will use if I am estimating marker slippage.

    string PnId, chrom, motif, nextWord, refRepSeq;
    String<string> PnIds;
    std::set<Marker> markers;
    bool writeVcf = false, loadModAndLab = true, enoughReads = true;
    int start, end, numberOfReads, finalItNum = 5;
    if (currItNum == finalItNum)
        writeVcf = true;
    float refRepeatNum, winner, second;
    //If this is not the first iteration, I read the marker slippage values from the previous iteration into the markerToSizeAndModel map.
    if (currItNum > 1)
    {
        append(markerSlippageDir,"markerSlippage");
        append(markerSlippageDir, prevItNumStr);
        cout << "Path to marker slippage file: " << markerSlippageDir << endl;
        readMarkerSlippage(markerSlippageDir, markerToSizeAndModel, startCoord, endCoord);
        markerSlippageDir = markerSlippageDirBase;
    }
    //Else, this is the first iteration and I don't have any regression models or labels to load.
    else
    {
        loadModAndLab = false;
        //Check if markerSlipps folder exist and chr folder, otherwise create them
        struct stat st;
        if(stat(toCString(markerSlippageDirBase),&st) != 0)
            mkdir(toCString(markerSlippageDirBase),0777);
        if(stat(toCString(modelAndLabelDirBase),&st) != 0)
            mkdir(toCString(modelAndLabelDirBase),0777);
    }
    append(markerSlippageDir,"markerSlippage");
    append(markerSlippageDir, currItNumStr);
    append(markerSlippageDir, "_");
    append(markerSlippageDir, intervalIndex);
    markerSlippageOut.open(toCString(markerSlippageDir));
    //Read the slippage rate for all PNs into the pnToSize map.
    readPnSlippage(pnSlippageFile);
    //Map from marker to all reads covering it
    map<Marker, String<AttributeLine> > mapPerMarker;
    //Count the total number of alleles in the population for each marker -- by checking the size of the set
    map<Marker, std::set<float> > markerToAlleles;
    //Map to store current genotype(and lots of other things) of a person for each marker maps from pnId and Marker-struct to GenotypeInfo struct
    map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype;

    Marker marker;
    string nextLine;
    AttributeLine currentLine;
    Pair<int, String<string> > numberOfWordsAndWords;
    ifstream labelFile;

    //Iterate over all Pns I have slippage for and read from attribute files in the given interval
    map<string, Pair<double, int> >::const_iterator pnEnd = pnToSize.end();
    int nProcessedPns = 0;
    for (map<string, Pair<double, int> >::iterator pnStart = pnToSize.begin(); pnStart != pnEnd; ++pnStart)
    {
        PnId = pnStart->first;
        append(attributePath, PnId);
        append(attributePath, ".gz");
        std::ifstream attributeFileGz(toCString(attributePath), std::ios_base::in | std::ios_base::binary);
        io::filtering_istream attributeFile;
        attributeFile.push(io::gzip_decompressor());
        attributeFile.push(attributeFileGz);
        if (attributeFile.fail())
        {
            cout << "Unable to locate attribute file for " << PnId << " at " << attributePath << endl;
            attributePath = attributePathBase;
            continue;
        }
        ++nProcessedPns;
        if (nProcessedPns % 1000==0)
            cout << "Working on pn number: " << nProcessedPns << endl;
        while (getline (attributeFile,nextLine))
        {
            if (nextLine.length() == 0)
                continue;
            numberOfWordsAndWords = countNumberOfWords(nextLine);
            if (numberOfWordsAndWords.i1 == 1)
            {
                PnId = nextLine;
                appendValue(PnIds,PnId);
                if (loadModAndLab)
                {
                    if (labelFile.is_open())
                        labelFile.close();
                    append(modelAndLabelDir,PnId);
                    append(modelAndLabelDir,"labels");
                    append(modelAndLabelDir, prevItNumStr);
                    labelFile.open(toCString(modelAndLabelDir));
                    if(labelFile.fail())
                        cout << "Could not open label file: " << modelAndLabelDir << endl;
                    modelAndLabelDir = modelAndLabelDirBase;
                }
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
                    if (loadModAndLab)
                        getline(labelFile, nextLine);
                    continue;
                }
                //If I have passed the interval, I move on to the next pn.
                if (marker.start > endCoord)
                    break;
                //Otherwise the marker is in the interval and I process it.
                markers.insert(marker);
                if (loadModAndLab)
                {
                    labelFile >> winner;
                    labelFile >> second;
                    if (markerToSizeAndModel[marker].i2 == NULL || markerToSizeAndModel.count(marker) == 0)
                    {
                        append(modelAndLabelDir, "/");
                        append(modelAndLabelDir, marker.chrom);
                        append(modelAndLabelDir, "_");
                        stringstream startStr;
                        startStr << marker.start;
                        append(modelAndLabelDir, startStr.str());
                        startStr.clear();
                        startStr.str("");
                        const char *model_in_file = toCString(modelAndLabelDir);
                        markerToSizeAndModel[marker].i2 = load_model(model_in_file);
                        modelAndLabelDir = modelAndLabelDirBase;
                    }
                }
                else
                {
                    winner = lexicalCast<float>(numberOfWordsAndWords.i2[7]);
                    second = lexicalCast<float>(numberOfWordsAndWords.i2[8]);
                }
                markerToAlleles[marker].insert(winner);
                markerToAlleles[marker].insert(second);
            }
            if (numberOfWordsAndWords.i1 == 11)
            {
                if (numberOfReads < 10)
                {
                    PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].pValueSum = 0;
                    PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].fullMotifSlippageSum = 0;
                    enoughReads = false;
                }
                for (unsigned i = 0; i < numberOfReads; ++i)
                {
                    if (i == 0)
                        currentLine = parseNextLine(winner, second, attributeFile, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, true, loadModAndLab, enoughReads);
                    else
                        currentLine = parseNextLine(winner, second, attributeFile, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, false, loadModAndLab, enoughReads);
                    appendValue(mapPerMarker[marker],currentLine);
                }
                enoughReads = true;
            }
            if (numberOfWordsAndWords.i1 != 1 && numberOfWordsAndWords.i1 != 9 && numberOfWordsAndWords.i1 != 11)
                cerr << "Format error in attribute file!" << endl;
        }
        attributePath = attributePathBase;
        labelFile.close();
    }
    chrom = marker.chrom;
    cout << "Reading data from input complete." << endl;

    //Open vcf stream and make header if the writeVcf switch is on
    VcfStream out;
    if (writeVcf)
    {
        CharString outputDirectory = argv[10];
        CharString outputFileName = argv[11];
        append(outputDirectory, "/");
        append(outputDirectory, outputFileName);
        append(outputDirectory, "_");
        append(outputDirectory, intervalIndex);
        append(outputDirectory, ".vcf");
        bool outOk = open(out,toCString(outputDirectory), VcfStream::WRITE);
        makeVcfHeader(out, PnIds, chrom);
    }

    cout << "Number of markers: " << mapPerMarker.size() << endl;
    cout << "Number of pns: " << length(PnIds) << endl;

    double *prob_estimates = NULL;
    double predict_label;
    int PnsAtMarker;
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
    //map<Marker, String<AttributeLine> >::iterator itEnd = mapPerMarker.end();
    //for (map<Marker, String<AttributeLine> >::iterator it = mapPerMarker.begin(); it != itEnd; ++it)
    map<Marker, String<AttributeLine> >::iterator it = mapPerMarker.begin();
    while (mapPerMarker.size()>0)
    {
        thisMarker = it->first;
        PnsAtMarker = 1;
        String<AttributeLine>& currentMarker = it->second;
        probBig.l = length(currentMarker);
        cout << "Starting marker number: " << z << " with start coordinate: " << thisMarker.start << endl;
        double geomP = 1/(fmod(markerToStepSum[thisMarker]/(float)length(currentMarker),1.0)+1);
        allelesAtMarker = markerToAlleles[thisMarker];
        numOfAlleles = allelesAtMarker.size();
        markerToNallelesAndStutter[thisMarker] = Pair<int,double>(numOfAlleles, geomP);
        markerToAlleles[thisMarker].clear();
        markerToAlleleFreqs[thisMarker].i1.clear();
        //Estimate marker slippage and update in markerToSizeAndModel map
        slippAndNavail = estimateSlippage(PnIds, PnAndMarkerToGenotype, it->first, markerToSizeAndModel, currItNum);
        markerToSizeAndModel[it->first].i1.i2 = std::max((double)0.0,slippAndNavail.i1);
        nAvailable = slippAndNavail.i2;
        //Reads with label 2 are not included in training
        prob.l = length(currentMarker) - markerToSizeAndModel[it->first].i1.i1.p2.i1;
        //Now I "nullSet" the markerToSizeAndModel map for the marker I am looking at so I can update it when I relabel the reads
        markerToSizeAndModel[it->first].i1.i1.p1 = Pair<int,double>(0,0);
        markerToSizeAndModel[it->first].i1.i1.p2 = Pair<int,double>(0,0);
        markerToSizeAndModel[it->first].i1.i1.p3 = Pair<int,double>(0,0);
        //Nullset the averageStepSize map
        markerToStepSum[it->first] = 0.0;
        int idx = 0;
        prob.y = Malloc(double,prob.l);
        prob.x = (feature_node **) malloc(prob.l * sizeof(feature_node *));
        probBig.x = (feature_node **) malloc(probBig.l * sizeof(feature_node *));
        PnId = currentMarker[0].PnId;
        for (unsigned i = 0; i<length(currentMarker); ++i)
        {
            currentLine = currentMarker[i];
            if (currentLine.PnId != PnId)
            {
                PnId = currentLine.PnId;
                ++PnsAtMarker;
            }
            PnToAlleles[currentLine.PnId].insert(currentLine.numOfRepeats);
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
        Pair<GenotypeInfo, bool> changed;
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
                genotypesToConsider = makeGenotypes(allelesToConsider, thisMarker.refRepeatNum).genotypes;
                //make decision about genotype for PnId at the current marker.
                changed = determineGenotype(reads, markerToSizeAndModel[it->first].i1.i2+pnToSize[PnId].i1, genotypesToConsider, numOfAlleles, it->first.motif.size(), geomP);
                PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;
                Pair<bool> alleleConfidence = genotypeIsConfident(changed.i1);
                if (alleleConfidence.i1)
                    markerToTrueAlleles[it->first].insert(changed.i1.genotype.i1);
                if (alleleConfidence.i2)
                    markerToTrueAlleles[it->first].insert(changed.i1.genotype.i2);
                if (changed.i2)
                {
                    markerToAlleles[it->first].insert(changed.i1.genotype.i1);
                    markerToAlleles[it->first].insert(changed.i1.genotype.i2);
                }
                ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
                ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];
                //If I am estimating the marker slippage then I should update map from Pn to labels. (before I update PnId to currentLine.PnId)
                if (!writeVcf)
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
        genotypesToConsider = makeGenotypes(allelesToConsider, thisMarker.refRepeatNum).genotypes;
        changed = determineGenotype(reads, markerToSizeAndModel[it->first].i1.i2+pnToSize[PnId].i1, genotypesToConsider, numOfAlleles, it->first.motif.size(), geomP);
        PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;
        Pair<bool> alleleConfidence = genotypeIsConfident(changed.i1);
        if (alleleConfidence.i1)
            markerToTrueAlleles[it->first].insert(changed.i1.genotype.i1);
        if (alleleConfidence.i2)
            markerToTrueAlleles[it->first].insert(changed.i1.genotype.i2);
        if (changed.i2)
        {
            markerToAlleles[it->first].insert(changed.i1.genotype.i1);
            markerToAlleles[it->first].insert(changed.i1.genotype.i2);
        }
        ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
        ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];
        //If I am estimating the marker slippage then here is where I update the pn to labels for the last PN.
        if (!writeVcf)
            append(pnToLabels[PnId],Pair<float>(changed.i1.genotype.i1, changed.i1.genotype.i2));
        //Save logistic regression model to output file so I can use it in pn-slippage estimation
        stringstream zstr;
        zstr << thisMarker.chrom;
        append(modelAndLabelDir, zstr.str());
        append(modelAndLabelDir, "_");
        zstr.str("");
        zstr.clear();
        zstr << thisMarker.start;
        append(modelAndLabelDir, zstr.str());
        const char *model_out_file = toCString(modelAndLabelDir);
        if (save_model(model_out_file,model_) != 0)
                cout << "Unable to save model for marker number " << z << endl;
        modelAndLabelDir = modelAndLabelDirBase;
        //Write vcf record for the marker I just finished.
        if (writeVcf)
        {
            //Check the number of alleles in the population, should be 2 or higher to be considered polymorphic.
            if (markerToTrueAlleles[thisMarker].size() < 2)
            {
                cout << "Not enough alleles at marker " << z << endl;
                mapPerMarker.erase(thisMarker);
                it = mapPerMarker.begin();
                ++z;
                continue;
            }
            //Make a String<Pair<float> > which contains a list of genotypes
            cout << "Number of verified alleles: " << markerToTrueAlleles[thisMarker].size() << endl;
            //Have to add ref allele to true allele set in case no one has it.
            markerToTrueAlleles[thisMarker].insert(thisMarker.refRepeatNum);
            genotypesAtThisMarker = makeGenotypes(markerToTrueAlleles[thisMarker], thisMarker.refRepeatNum);
            //Compute abs(allele1-allele2)*allele1Freq*allele2Freq for all genotypes and return average of those, estimate of distance between alleles.
            alleleDistance = computeAlleleDist(genotypesAtThisMarker.genotypes, markerToAlleleFreqs[it->first].i1, PnsAtMarker);
            //First fill marker specific fields of vcfRecord
            record = fillRecordMarker(thisMarker, markerToTrueAlleles[thisMarker]);
            //Loop over Pns and fill in PN specific fields of vcfRecord for each PN
            for (unsigned i = 0; i<length(PnIds); ++i)
            {
                thisPn = PnIds[i];
                //If a decision has been made for thisPn at thisMarker I add it to the genotypeInfos stringSet
                if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(thisPn, thisMarker)) != 0))
                {
                    genotype = PnAndMarkerToGenotype[Pair<string,Marker>(thisPn, thisMarker)];
                    fillRecordPn(genotype, record, genotypesAtThisMarker,markerToTrueAlleles[thisMarker], thisMarker);
                    PnAndMarkerToGenotype.erase(Pair<string,Marker>(thisPn, thisMarker));
                }
                //If a decision has not been made I add a CharString with no decision(0:0,0,0,0....etc) to the set to maintain order of Pns vs genotypeInfos in output
                else
                {
                    CharString gtInfo = "./.:";
                    //Adding alleleDepth zeros
                    for (unsigned i=0; i<markerToTrueAlleles[thisMarker].size(); ++i)
                        append(gtInfo,"0,");
                    eraseBack(gtInfo);
                    //Adding read depth zero
                    append(gtInfo,":0:");
                    //Adding phred likelihood zeros
                    for (unsigned i=0; i<length(genotypesAtThisMarker.genotypes); ++i)
                        append(gtInfo,"0,");
                    eraseBack(gtInfo);
                    appendValue(record.genotypeInfos, gtInfo);
                }
            }
            ss << markerToTrueAlleles[thisMarker].size();
            str = ss.str();
            record.filter = ".";
            stringClear(ss,str);
            //After adding info for all PNs at thisMarker I write the record to the vcf output file
            if (writeRecord(out, record) != 0)
                cerr << "ERROR: Problem writing VCF-output file." << endl;
            clear(record);
        }
        markerToAlleleFreqs[it->first].i2 = PnsAtMarker;
        PnToAlleles.clear();
        markerSlippageOut << thisMarker.chrom << "\t" << thisMarker.start << "\t" << thisMarker.end << "\t" << thisMarker.motif << "\t" << thisMarker.refRepeatNum << "\t" << thisMarker.refRepSeq << "\t" << setprecision(4) << fixed << markerToSizeAndModel[thisMarker].i1.i2 << "\t" << nAvailable << "\t" << markerToNallelesAndStutter[thisMarker].i1 << "\t" << markerToNallelesAndStutter[thisMarker].i2 << endl;
        cout << thisMarker.start << " totalSlipp: " << setprecision(4) << fixed << markerToSizeAndModel[thisMarker].i1.i2 << endl;
        cout << "Finished marker number: " << z << endl;
        mapPerMarker.erase(thisMarker);
        it = mapPerMarker.begin();
        ++z;
    }
    String<Pair<float> > labels;
    if (!writeVcf)
    {
        ofstream labelsOut;
        for (unsigned i = 0; i<length(PnIds); ++i)
        {
            labels = pnToLabels[PnIds[i]];
            if (length(labels) == 0)
                continue;
            append(modelAndLabelDir, PnIds[i]);
            append(modelAndLabelDir, "labels");
            append(modelAndLabelDir, currItNumStr);
            append(modelAndLabelDir, "_");
            append(modelAndLabelDir, intervalIndex);
            labelsOut.open(toCString(modelAndLabelDir));
            for (unsigned j = 0; j<length(labels); ++j)
                labelsOut << labels[j].i1 << "\t" << labels[j].i2 << endl;
            /*if (currItNum > 1)
            {
                modelAndLabelDir = modelAndLabelDirBase;
                append(modelAndLabelDir, PnIds[i]);
                append(modelAndLabelDir, "labels");
                append(modelAndLabelDir, prevItNumStr);
                append(modelAndLabelDir, "_");
                append(modelAndLabelDir, intervalIndex);
                if (remove(toCString(modelAndLabelDir)) !=0)
                    cout << "Remove operation of old labelFile failed with code " << errno << endl;
            }*/
            modelAndLabelDir = modelAndLabelDirBase;
            labelsOut.close();
        }
    }
    cout << "Finished determining genotypes" << endl;
    return 0;
}
