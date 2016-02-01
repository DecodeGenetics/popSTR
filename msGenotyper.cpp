#include <iostream>
#include <set>
#include <map>
#include <string>
#include <ctime>
#include <math.h> 
#include <vector>
#include <numeric>
#include <sstream>
#include <sys/stat.h>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/vcf_io.h>
#include <seqan/sequence.h>
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
} ;

//So I can map from Markers
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//Sums the pValues of reads and counts the number of reads for each type of label at every marker, also stores the current slippage rate value for each marker
map<Marker, Pair<Pair<LabelProps,double> ,model* > > markerToSizeAndModel;
map<Marker, Pair<Pair<LabelProps,double> ,model* > > markerToSizeAndModelHiSeq;
map<Marker, Pair<Pair<LabelProps,double> ,model* > > markerToSizeAndModelHiSeqX;
//Stores the slippage rate value for each PN - read from input file.
map<string, Pair<double, int> > pnToSize;
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
AttributeLine parseNextLine(float winner, float second, ifstream& attributeFile, Marker& marker, string PnId, map<Pair<string,Marker>, GenotypeInfo>& PnAndMarkerToGenotype, String<string> firstLine, bool useFirstLine, bool useModelAndLabels, bool enoughReads)
{
    PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].genotype = Pair<float>(winner,second);
    AttributeLine currentLine;
    currentLine.PnId = PnId;
    string seqMethod;
    if (pnToSeqMethod.count(PnId) != 0)
        seqMethod = pnToSeqMethod[PnId];
    else
        seqMethod = "mergedBAMfile"; 
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
        //Determining the initial label of the read
    if ((fabs(currentLine.numOfRepeats - winner) <= 0.05) || (fabs(currentLine.numOfRepeats - second) <= 0.05))
    {
        currentLine.label = 1;
        ++markerToSizeAndModel[marker].i1.i1.p1.i1;        
        if (enoughReads)
        {            
            markerToSizeAndModel[marker].i1.i1.p1.i2 += currentLine.pValue;
            if (seqMethod.compare("HiSeq")==0)
                markerToSizeAndModelHiSeq[marker].i1.i1.p1.i2 += currentLine.pValue;
            if (seqMethod.compare("HiSeqX")==0)
                markerToSizeAndModelHiSeqX[marker].i1.i1.p1.i2 += currentLine.pValue;
        }
    }
    else
    {
        if ((fabs(currentLine.numOfRepeats - (winner - 1)) <= 0.05) || (fabs(currentLine.numOfRepeats - (second - 1)) <= 0.05))
        {
            currentLine.label = 2;            
            ++markerToSizeAndModel[marker].i1.i1.p2.i1;
            if (enoughReads)
            {
                markerToSizeAndModel[marker].i1.i1.p2.i2 += currentLine.pValue;
                if (seqMethod.compare("HiSeq")==0)
                    markerToSizeAndModelHiSeq[marker].i1.i1.p2.i2 += currentLine.pValue;
                if (seqMethod.compare("HiSeqX")==0)
                    markerToSizeAndModelHiSeqX[marker].i1.i1.p2.i2 += currentLine.pValue;
            }
            markerToStepSum[marker] += (float)marker.motif.size();
        }
        else
        {             
            currentLine.label = -1;
            ++markerToSizeAndModel[marker].i1.i1.p3.i1;
            if (enoughReads)
            {
                markerToSizeAndModel[marker].i1.i1.p3.i2 += currentLine.pValue;
                if (seqMethod.compare("HiSeq")==0)
                    markerToSizeAndModelHiSeq[marker].i1.i1.p3.i2 += currentLine.pValue;
                if (seqMethod.compare("HiSeqX")==0)
                    markerToSizeAndModelHiSeqX[marker].i1.i1.p3.i2 += currentLine.pValue;
            }
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

String<Pair<float> > makeGenotypes(std::set<float> alleles)
{
    String<Pair<float> > genotypes;
    String<float> alleleString;
    std::set<float>::reverse_iterator allelesBegin = alleles.rend();
    for (std::set<float>::reverse_iterator alleleIt = alleles.rbegin(); alleleIt!=allelesBegin; ++alleleIt)
        appendValue(alleleString, *alleleIt);  
    for (unsigned i=0; i<length(alleleString); ++i)
    {
        appendValue(genotypes,Pair<float>(alleleString[i],alleleString[i]));
        if (i == (length(alleleString)-1))
            break;
        for (unsigned j=i+1; j<length(alleleString); ++j)
            appendValue(genotypes,Pair<float>(alleleString[j],alleleString[i])); 
    }
    reverse(genotypes);
    return genotypes;
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

/* same as R dgeom */
static float dgeom(int diff, double psucc) 
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

Pair<GenotypeInfo, Pair<bool> > determineGenotype(String<AttributeLine> reads, double markerSlippage, String<Pair<float> > genotypes, int numberOfAlleles, int motifLength, double psucc)
{
    GenotypeInfo returnValue;
    returnValue.pValueSum = 0;
    returnValue.genotypes = genotypes;
    returnValue.numOfReads = length(reads);
    Pair<float> genotypeToCheck;
    AttributeLine readToCheck;
    std::set<float> currentGenotype;
    std::set<float> newGenotypeSet; 
    String<long double> probs;
    double errorProbSum = 0;
    double lengthSlippage;
    resize(probs, length(genotypes));
    bool isHomo, enoughDistance = true;
    float posNegSlipp = 1;
    float posNegSlipp2 = 1;
    float diff;
    float diff2;
    int indexOfWinner, indexOfSecond;    
    for (unsigned i=0; i<length(genotypes); ++i)
    {
        probs[i] = 1;
        genotypeToCheck = genotypes[i];
        isHomo = genotypeToCheck.i1 == genotypeToCheck.i2;
        for (unsigned j=0; j<length(reads); ++j)
        {
            posNegSlipp = 1;
            posNegSlipp2 = 1;
            boost::math::poisson_distribution<> myPoiss(std::max((double)0.01,markerSlippage));
            readToCheck = reads[j];
            if (i == 0)
            {
                ++returnValue.alleleToFreq[readToCheck.numOfRepeats];
                if (length(reads) >= 10)
                    returnValue.pValueSum += readToCheck.pValue;
                //Debugging code
                //cout << "P-value of read " << j << " with " << readToCheck.numOfRepeats << " repeats:  " << readToCheck.pValue << endl;
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
                //Debugging code
                //cout << "Diff for homozygous genotype " << genotypeToCheck.i1 << " from read with " << readToCheck.numOfRepeats << " repeats is: " << diff << endl;
                //cout << "Update of genotype probability: " << probs[i] << " *= (" << readToCheck.pValue << "*" << pdf(myPoiss, abs(diff)) << "*" << posNegSlipp << "+ (" << (double)(1.0-readToCheck.pValue)/(double)numberOfAlleles << "))";
                if (fmod(diff,1.0)>0.0 && useGeom)
                {
                    if (diff < 1.0)
                        probs[i] *= (readToCheck.pValue * dgeom(static_cast<int>(diff*motifLength), psucc) * pdf(myPoiss, 0) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    else                    
                        probs[i] *= (readToCheck.pValue * dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * pdf(myPoiss, floor(diff)) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));                    
                }
                else
                {
                    if (!(fmod(diff,1.0)>0.0) && useGeom)
                        probs[i] *= (readToCheck.pValue * pdf(myPoiss, diff) * dgeom(0, psucc) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    else
                        probs[i] *= (readToCheck.pValue * pdf(myPoiss, diff) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                }
                //cout << " = " << probs[i] << endl;
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
                //Debugging code
                //cout << "Diffs for heterozygous genotype " << genotypeToCheck.i1 << "/" << genotypeToCheck.i2 << " from read with " << readToCheck.numOfRepeats << " repeats are: " << diff << " and " << diff2 << endl;
                //cout << "Update of genotype probability: " << probs[i] << " *= (" << readToCheck.pValue << " * (0.5 * " << pdf(myPoiss, abs(diff)) << "*" << posNegSlipp << "+ 0.5 * "<< pdf(myPoiss, abs(diff2)) << "*" << posNegSlipp2 << ") + (" << (double)(1.0-readToCheck.pValue)/(double)numberOfAlleles << "))";
                if (fmod(diff,1.0)>0.0 && fmod(diff2,1.0)>0.0 && useGeom)
                {
                    if (diff < 1.0 && diff2 < 1.0)
                        probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>(diff*motifLength), psucc) * pdf(myPoiss, 0) * posNegSlipp + dgeom(static_cast<int>(diff2*motifLength), psucc) * pdf(myPoiss, 0) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    if (diff < 1.0 && !(diff2 < 1.0))
                        probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>(diff*motifLength), psucc) * pdf(myPoiss, 0) * posNegSlipp + dgeom(static_cast<int>((diff2-(float)floor(diff2))*motifLength), psucc) * pdf(myPoiss, floor(diff2)) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    if (!(diff < 1.0) && diff2 < 1.0)
                        probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * pdf(myPoiss, floor(diff)) * posNegSlipp + dgeom(static_cast<int>(diff2*motifLength), psucc) * pdf(myPoiss, 0) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    else
                        probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * pdf(myPoiss, floor(diff)) * posNegSlipp + dgeom(static_cast<int>((diff2-(float)floor(diff2))*motifLength), psucc) * pdf(myPoiss, floor(diff2)) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                }
                if (fmod(diff,1.0)>0.0 && !fmod(diff2,1.0)>0.0 && useGeom)
                {
                    if (!(diff < 1.0))    
                        probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * pdf(myPoiss, floor(diff)) * posNegSlipp + pdf(myPoiss, diff2) * dgeom(0, psucc) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    else
                        probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>(diff*motifLength), psucc) * posNegSlipp + pdf(myPoiss, diff2) * dgeom(0, psucc) *posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                }
                if (!fmod(diff,1.0)>0.0 && fmod(diff2,1.0)>0.0 && useGeom)
                {
                    if (!(diff2 < 1.0))
                        probs[i] *= (readToCheck.pValue * 0.5 * (pdf(myPoiss, diff) * dgeom(0, psucc) * posNegSlipp + dgeom(static_cast<int>((diff2-(float)floor(diff2))*motifLength), psucc) * pdf(myPoiss, floor(diff2)) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    else
                        probs[i] *= (readToCheck.pValue * 0.5 * (pdf(myPoiss, diff) * dgeom(0, psucc) * posNegSlipp + dgeom(static_cast<int>(diff2*motifLength), psucc) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                }
                else
                {
                    if (useGeom)
                        probs[i] *= (readToCheck.pValue * 0.5 * (pdf(myPoiss, diff) * dgeom(0, psucc) * posNegSlipp + pdf(myPoiss, diff2) * dgeom(0, psucc) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                    else
                        probs[i] *= (readToCheck.pValue * 0.5 * (pdf(myPoiss, diff) * posNegSlipp + pdf(myPoiss, diff2) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                }
            }
        }
        //Debugging code
        //cout << "Genotype: " << genotypeToCheck.i1 << "/" << genotypeToCheck.i2 << " with probability: " << probs[i] << endl;
    } 
    returnValue.pValues = probs;    
    indexOfWinner = findMaxIndex(probs);
    if (length(probs)>1)
    {
        String<long double> probsCopy = probs;
        erase(probsCopy,indexOfWinner);
        indexOfSecond = findMaxIndex(probsCopy);
        enoughDistance = round(-10*log10(probsCopy[indexOfSecond]/probs[indexOfWinner])) > 10;
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
        string pnId = readsToRelabel[i].PnId, seqMethod;
        if (pnToSeqMethod.count(pnId) != 0)
            seqMethod = pnToSeqMethod[pnId];
        else
            seqMethod = "mergedBAMfile";
        if ((fabs(readsToRelabel[i].numOfRepeats - newGenotype.i1)<=0.05) || (fabs(readsToRelabel[i].numOfRepeats - newGenotype.i2)<=0.05))
        {
            readsToRelabel[i].label = 1;
            ++markerToSizeAndModel[marker].i1.i1.p1.i1;
            if (numOfReads >= 10)
            {
                markerToSizeAndModel[marker].i1.i1.p1.i2 += readsToRelabel[i].pValue;
                if (seqMethod.compare("HiSeq")==0)
                    markerToSizeAndModelHiSeq[marker].i1.i1.p1.i2 += readsToRelabel[i].pValue;
                if (seqMethod.compare("HiSeqX")==0)
                    markerToSizeAndModelHiSeqX[marker].i1.i1.p1.i2 += readsToRelabel[i].pValue;
            }
        }
        else 
        {
            if ((fabs(readsToRelabel[i].numOfRepeats - (newGenotype.i1 - 1))<=0.05) || (fabs(readsToRelabel[i].numOfRepeats - (newGenotype.i2 - 1))<=0.05))
            {
                readsToRelabel[i].label = 2;
                ++markerToSizeAndModel[marker].i1.i1.p2.i1;
                markerToStepSum[marker] += (float)marker.motif.size();
                if (numOfReads >= 10)
                {
                    markerToSizeAndModel[marker].i1.i1.p2.i2 += readsToRelabel[i].pValue;
                    if (seqMethod.compare("HiSeq")==0)
                        markerToSizeAndModelHiSeq[marker].i1.i1.p2.i2 += readsToRelabel[i].pValue;
                    if (seqMethod.compare("HiSeqX")==0)
                        markerToSizeAndModelHiSeqX[marker].i1.i1.p2.i2 += readsToRelabel[i].pValue;
                }
            }
            else
            { 
                readsToRelabel[i].label = -1;
                ++markerToSizeAndModel[marker].i1.i1.p3.i1;
                if (numOfReads >= 10)
                {                
                    markerToSizeAndModel[marker].i1.i1.p3.i2 += readsToRelabel[i].pValue;
                    if (seqMethod.compare("HiSeq")==0)
                        markerToSizeAndModelHiSeq[marker].i1.i1.p3.i2 += readsToRelabel[i].pValue;
                    if (seqMethod.compare("HiSeqX")==0)
                        markerToSizeAndModelHiSeqX[marker].i1.i1.p3.i2 += readsToRelabel[i].pValue;
                }
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
    appendValue(out.header.headerRecords, VcfHeaderRecord("reference", "/odinn/data/reference/Homo_sapiens-deCODE-hg38/Sequence/WholeGenomeFasta/genome.fa"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=MOTIF,Number=1,Type=String,Description=\"Repeat motif\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=REF,Number=1,Type=Float,Description=\"Copy number in reference\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=RL,Number=1,Type=Integer,Description=\"Length of STR in reference in bp\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("INFO", "<ID=VT,Number=String,Type=Flag,Description=\"Variant type\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FILTER", "<ID=N,Description=\"N = number of alleles in population\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(out.header.headerRecords, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">"));
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
    ss << marker.start-1;
    CharString str = ss.str();
    record.id = str;
    stringClear(ss,str);
    //Set smallest allele as reference allele
    std::set<float>::iterator smallestIt = allelesAtThisMarker.begin();
    float smallestAllele = *smallestIt;
    ss << smallestAllele;
    str = ss.str();
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
    if (genotype.genotype.i1 == 0)
    {
        for (unsigned i=0; i<length(genotypesAtThisMarker); ++i)
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
                ss << pl;
                str = ss.str();
                append(gtInfo,str);
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

//Read in data from markerSlippageFile and set estimateMarkerSlippage switch to false
void readMarkerSlippage(CharString markerSlippageFile, map<Marker, Pair<Pair<LabelProps,double>, model* > >& markerToSizeAndModel, int startCoord, int endCoord)
{    
    ifstream markerSlippageIn(toCString(markerSlippageFile));
    double currMarkSlipp;
    int nPns;
    Marker currMarker; 
    unsigned counter = 1; 
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
        if (markerSlippageIn.eof())
            break;
        if (currMarker.start < startCoord) 
            continue;
        if (currMarker.start > endCoord)
            break;
        markerToSizeAndModel[currMarker].i1.i2 = currMarkSlipp;
        markerToSizeAndModel[currMarker].i2 = NULL;
        cout << "Finished processing marker number: " << counter << endl;
        ++counter;
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
        if (pnSlippageFile.eof())
            break;
        pnToSize[PnId].i1= currPnSlipp/2.0;        
        pnToSize[PnId].i2 = nMarkers; 
    }
    cout << "Finished reading pn Slippage." << endl;
    pnSlippageFile.close();
}

void readPnSeqMethod(ifstream& pnToSeqMethodFile, map<string, String<string> >& seqMethodToPns)
{
    string pnId, seqMethod;
    while(true)
    {
        pnToSeqMethodFile >> pnId;
        pnToSeqMethodFile >> seqMethod;
        if (pnToSeqMethodFile.eof() )
            break;
        appendValue(seqMethodToPns[seqMethod], pnId);
        pnToSeqMethod[pnId] = seqMethod;
    }
    cout << "Finished reading pn to seqMethod file." << endl;
}

inline bool exists(const std::string& name) 
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

Pair<double, int> estimateSlippage(String<string> PnIds, map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype, Marker marker, map<Marker, Pair<Pair<LabelProps,double> ,model* > > currMarkerToSizeAndModel)
{
    int nMissing = 0, nAvailable;
    double firstPart, secondPart, subSum, finalSub=0, t, result;
    vector<double> subtractions;
    for (unsigned i = 0; i<length(PnIds); ++i)
    {                
        if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
        {
            subtractions.push_back(0);
            ++nMissing;
            continue;
        }
        if (pnToSize[PnIds[i]].i1 == 0)
            t = 0.001;
        else 
            t = pnToSize[PnIds[i]].i1;
        firstPart = 1.0/(t*(1.0-t));                
        secondPart = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum/(currMarkerToSizeAndModel[marker].i1.i1.p1.i2+currMarkerToSizeAndModel[marker].i1.i1.p2.i2+currMarkerToSizeAndModel[marker].i1.i1.p3.i2);                                  
        subtractions.push_back(firstPart*secondPart);
    }
    subSum = accumulate(subtractions.begin(),subtractions.end(),0.0);
    cout << "subSum: " << subSum << endl;
    nAvailable = length(PnIds) - nMissing;
    cout << "Number of pns available for estimating marker slippage: " << nAvailable << endl;
    if (subSum == 0.0)
        subSum = 0.001;
    for (unsigned i = 0; i<length(PnIds); ++i)
    {                
        finalSub += (pnToSize[PnIds[i]].i1)*(subtractions[i]/subSum);
    }
    cout << "(" << currMarkerToSizeAndModel[marker].i1.i1.p2.i2 << "+" << currMarkerToSizeAndModel[marker].i1.i1.p3.i2 << ")/(" << currMarkerToSizeAndModel[marker].i1.i1.p1.i2<< "+" << currMarkerToSizeAndModel[marker].i1.i1.p2.i2 << "+" << currMarkerToSizeAndModel[marker].i1.i1.p3.i2 <<  ")-" << finalSub << endl;
    result = (double)(currMarkerToSizeAndModel[marker].i1.i1.p2.i2+currMarkerToSizeAndModel[marker].i1.i1.p3.i2)/(double)(currMarkerToSizeAndModel[marker].i1.i1.p1.i2+currMarkerToSizeAndModel[marker].i1.i1.p2.i2+currMarkerToSizeAndModel[marker].i1.i1.p3.i2) - (double)finalSub;
    cout << "Slippage rate: " << result << endl;
    return Pair<double, int>(result, nAvailable);    
}

int main(int argc, char const ** argv)
{   
    //Check arguments.
    if (argc != 10 && argc != 12)
    {
        cerr << "USAGE: " << argv[0] << " attDir/chromNum/ PN-slippageFile startCoordinate endCoordinate intervalNumber markerSlippageFile modelAndLabelDir/chromNum/ iterationNumber pnToSeqMethodFile vcfOutputDirectory/ vcfFileName \n";
        return 1;
    }
    
    //Parse parameters     
    int startCoord = lexicalCast<int>(argv[3]), endCoord = lexicalCast<int>(argv[4]), currItNum = lexicalCast<int>(argv[8]), prevItNum = lexicalCast<int>(argv[8]) - 1;
    CharString attributePath = argv[1], pnSlippagePath = argv[2], intervalIndex = argv[5], markerSlippageFile = argv[6], modelAndLabelDir = argv[7], currItNumStr = argv[8], prevItNumStr = to_string(prevItNum);
    append(pnSlippagePath, prevItNumStr);    
    ifstream pnSlippageFile(toCString(pnSlippagePath)), pnToSeqMethodFile(argv[9]);
    ofstream markerSlippageOut; //Will use if I am estimating marker slippage, otherwise I take the file as input
    
    string PnId, chrom, motif, nextWord, refRepSeq;
    String<string> PnIds;
    std::set<Marker> markers;    
    int start, end, numberOfReads, finalItNum = 5;    
    float refRepeatNum, winner, second;     
    bool estimateMarkerSlippage = true, loadModAndLab = true, enoughReads = true;    
    //If this is the final iteration, I read the marker slippage values from the previous iteration into the markerToSizeAndModel map.
    if (currItNum == finalItNum)
    {                
        estimateMarkerSlippage = false;
        append(markerSlippageFile, prevItNumStr);
        cout << "Path to marker slippage file: " << markerSlippageFile << endl;
        readMarkerSlippage(markerSlippageFile, markerToSizeAndModel, startCoord, endCoord);
    }
    //Otherwise I'm going to estimate the marker slippage values and write them to an output file.
    else
    { 
        append(markerSlippageFile, currItNumStr);
        append(markerSlippageFile, "_");
        append(markerSlippageFile, intervalIndex);
        markerSlippageOut.open(toCString(markerSlippageFile)); 
    }
    //If this is the first iteration I don't have any regression models or labels to load.
    if (currItNum == 1)
        loadModAndLab = false;
    //Read the slippage rate for all PNs into the pnToSize map and divide by two.
    readPnSlippage(pnSlippageFile);
    //Read sequencing method for PNs from input into map (will not have method for every PN, some are merged from more than one method)
    map<string, String<string> > seqMethodToPns;
    readPnSeqMethod(pnToSeqMethodFile, seqMethodToPns);
    String<string> hiSeqPns = seqMethodToPns["HiSeq"];
    String<string> hiSeqXPns = seqMethodToPns["HiSeqX"];
    cout << "Number of hiSeqPns: " << length(hiSeqPns) << endl;
    cout << "Number of hiSeqXPns: " << length(hiSeqXPns) << endl;
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
    for (map<string, Pair<double, int> >::iterator pnStart = pnToSize.begin(); pnStart != pnEnd; ++pnStart)
    {
        PnId = pnStart->first;
        append(attributePath, PnId);
        ifstream attributeFile(toCString(attributePath));
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
                    modelAndLabelDir = argv[7];
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
                        append(modelAndLabelDir, marker.chrom);
                        append(modelAndLabelDir, "_");
                        stringstream startStr;
                        startStr << marker.start;
                        append(modelAndLabelDir, startStr.str());
                        startStr.clear();
                        startStr.str("");
                        const char *model_in_file = toCString(modelAndLabelDir);
                        markerToSizeAndModel[marker].i2 = load_model(model_in_file);
                        modelAndLabelDir = argv[7];
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
                    enoughReads = false;                    
                }
                for (unsigned i = 0; i < numberOfReads; ++i)
                {
                    if (i == 0)
                        currentLine = parseNextLine(winner, second, attributeFile, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, true, loadModAndLab, enoughReads);
                    else 
                        currentLine = parseNextLine(winner, second, attributeFile, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, false, loadModAndLab, enoughReads);;
                    if (currentLine.label == 1)
                        markerToAlleles[marker].insert(currentLine.numOfRepeats);
                    appendValue(mapPerMarker[marker],currentLine);
                }
                enoughReads = true;                                    
            }
            if (numberOfWordsAndWords.i1 != 1 && numberOfWordsAndWords.i1 != 9 && numberOfWordsAndWords.i1 != 11) 
                cerr << "Format error in attribute file!" << endl;
        }            
        attributePath = argv[1];
        attributeFile.close();
        labelFile.close();
    }
    chrom = marker.chrom;
    cout << "Reading data from input complete." << endl;    
    
    //Open vcf stream and make header if the estimateMarkerSlippage switch is off
    VcfStream out;
    if (!estimateMarkerSlippage)
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
    
    //Here I estimate initialization of alpha and beta using linear regression on how slippage changes as a function of allele length.
    /*double xSum = 0, ySum = 0, xySum = 0;
    double alpha, beta;
    double allLengthBar, lengthSlippBar;
    allLengthBar = (double)allLengthSum/(double)alleleLengths.size();
    lengthSlippBar = lengthSlippSum/(double)lengthSlipps.size();
    for (unsigned i = 0; i< lengthSlipps.size(); ++i)
    {
       xSum += pow(alleleLengths[i]-allLengthBar,2);
       ySum += pow(lengthSlipps[i]-lengthSlippBar,2);
       xySum += (alleleLengths[i]-allLengthBar) * (lengthSlipps[i]-lengthSlippBar);
    }
    beta = xySum/xSum;
    alpha = lengthSlippBar - beta*allLengthBar;
    alleleLengths.clear();
    lengthSlipps.clear();*/
    
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
    int z = 1, nAvailable, hiSeqAvailable, hiSeqXAvailable;
    
    //Stuff for vcf file 
    GenotypeInfo genotype;
    string thisPn;
    Marker thisMarker;
    VcfRecord record;        
    stringstream ss;
    CharString str;    
    String<Pair<float> > genotypesAtThisMarker;
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
            //Estimate marker slippage and update in markerToSizeAndModel map if it is not provided in an input file
            if (estimateMarkerSlippage)
            {               
                slippAndNavail = estimateSlippage(PnIds, PnAndMarkerToGenotype, it->first, markerToSizeAndModel);
                markerToSizeAndModel[it->first].i1.i2 = std::max((double)0.0,slippAndNavail.i1);
                nAvailable = slippAndNavail.i2;
                slippAndNavail = estimateSlippage(hiSeqPns, PnAndMarkerToGenotype, it->first, markerToSizeAndModelHiSeq);
                markerToSizeAndModelHiSeq[it->first].i1.i2 = std::max((double)0.0,slippAndNavail.i1);
                hiSeqAvailable = slippAndNavail.i2;
                slippAndNavail = estimateSlippage(hiSeqXPns, PnAndMarkerToGenotype, it->first, markerToSizeAndModelHiSeqX);
                markerToSizeAndModelHiSeqX[it->first].i1.i2 = std::max((double)0.0,slippAndNavail.i1);
                hiSeqXAvailable = slippAndNavail.i2;
            }
            //Reads with label 2 are not included in training
            prob.l = length(currentMarker) - markerToSizeAndModel[it->first].i1.i1.p2.i1;
            //Now I "nullSet" the markerToSizeAndModel map for the marker I am looking at so I can update it when I relabel the reads 
            markerToSizeAndModel[it->first].i1.i1.p1 = Pair<int,double>(0,0);
            markerToSizeAndModel[it->first].i1.i1.p2 = Pair<int,double>(0,0);
            markerToSizeAndModel[it->first].i1.i1.p3 = Pair<int,double>(0,0);
            //Do it also for extra maps counting only hiseq and hiseqx pns
            markerToSizeAndModelHiSeq[it->first].i1.i1.p1 = Pair<int,double>(0,0);
            markerToSizeAndModelHiSeq[it->first].i1.i1.p2 = Pair<int,double>(0,0);
            markerToSizeAndModelHiSeq[it->first].i1.i1.p3 = Pair<int,double>(0,0);
            markerToSizeAndModelHiSeqX[it->first].i1.i1.p1 = Pair<int,double>(0,0);
            markerToSizeAndModelHiSeqX[it->first].i1.i1.p2 = Pair<int,double>(0,0);
            markerToSizeAndModelHiSeqX[it->first].i1.i1.p3 = Pair<int,double>(0,0);
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
                    genotypesToConsider = makeGenotypes(allelesToConsider);  
                    //make decision about genotype for PnId at the current marker.        
                    changed = determineGenotype(reads, markerToSizeAndModel[it->first].i1.i2/2.0+pnToSize[PnId].i1, genotypesToConsider, numOfAlleles, it->first.motif.size(), geomP);                   
                    if (changed.i2.i1)
                        ++updatedPns;
                    relabelReads(currentMarker, i-length(reads), i, changed.i1.genotype, it->first);
                    PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1; 
                    if (changed.i2.i2)
                    {
                        markerToAlleles[it->first].insert(changed.i1.genotype.i1);
                        markerToAlleles[it->first].insert(changed.i1.genotype.i2);
                    }
                    ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
                    ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];
                    //If I am estimating the marker slippage then I should update map from Pn to labels. (before I update PnId to currentLine.PnId)
                    if (estimateMarkerSlippage)
                    {
                        if (loops>1)
                            eraseBack(pnToLabels[PnId]);                       
                        append(pnToLabels[PnId],Pair<float>(changed.i1.genotype.i1, changed.i1.genotype.i2));
                    }
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
            genotypesToConsider = makeGenotypes(allelesToConsider);
            changed = determineGenotype(reads, markerToSizeAndModel[it->first].i1.i2+pnToSize[PnId].i1, genotypesToConsider, numOfAlleles, it->first.motif.size(), geomP);
            if (changed.i2.i1)
                ++updatedPns;
            relabelReads(currentMarker, length(currentMarker)-length(reads), length(currentMarker), changed.i1.genotype, it->first);
            PnAndMarkerToGenotype[Pair<string,Marker>(PnId,it->first)] = changed.i1;
            if (changed.i2.i2)
            {
                markerToAlleles[it->first].insert(changed.i1.genotype.i1);
                markerToAlleles[it->first].insert(changed.i1.genotype.i2);
            }
            ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i1];
            ++markerToAlleleFreqs[it->first].i1[changed.i1.genotype.i2];
            //If I am estimating the marker slippage then here is where I update the pn to labels for the last PN.
            if (estimateMarkerSlippage)
            {
                if (loops>1)
                    eraseBack(pnToLabels[PnId]);
                append(pnToLabels[PnId],Pair<float>(changed.i1.genotype.i1, changed.i1.genotype.i2));
            }
            changedRatio = (float)updatedPns/(float)PnsAtMarker;
            cout << "Changed: " << changedRatio << endl;
        }
        //Save logistic regression model to output file so I can use it in pn-slippage estimation, don't do this if I have markerSlippage given.
        if (estimateMarkerSlippage)
        {            
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
            modelAndLabelDir = argv[7];        
        }
        //Write vcf record for the marker I just finished.
        if (!estimateMarkerSlippage) 
        {     
            //Check the number of alleles in the population, shouldn't be above ~20 but should be 2 or higher to be considered polymorphic.
            /*if ((markerToAlleles[thisMarker].size() > 20) || (markerToAlleles[thisMarker].size() < 2))
            {
                cout << "Not enough or too many alleles at marker " << z << endl;
                ++z;
                continue;      
            }*/
            //Make a String<Pair<float> > which contains a list of genotypes
            genotypesAtThisMarker = makeGenotypes(markerToAlleles[thisMarker]);
            //Compute abs(allele1-allele2)*allele1Freq*allele2Freq for all genotypes and return average of those, estimate of distance between alleles.
            alleleDistance = computeAlleleDist(genotypesAtThisMarker, markerToAlleleFreqs[it->first].i1, PnsAtMarker);
            //To make processing of hdf5 file easier I add 150 as a "dummy" allele when population only has 2 alleles (for Agnar)
	        //if (allelesAtMarker.size() == 2)
            //allelesAtMarker.insert(150);
            //genotypesAtThisMarker = makeGenotypes(allelesAtMarker); 
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
                    for (unsigned i=0; i<length(genotypesAtThisMarker); ++i)
                        append(gtInfo,"0,");
                    eraseBack(gtInfo);
                    appendValue(record.genotypeInfos, gtInfo);
                }
            }                
            ss << markerToAlleles[thisMarker].size();
            str = ss.str();
            record.filter = str;
            stringClear(ss,str);
            //After adding info for all PNs at thisMarker I write the record to the vcf output file
            if (writeRecord(out, record) != 0)
                cerr << "ERROR: Problem writing VCF-output file." << endl;
            clear(record);
        }
        markerToAlleleFreqs[it->first].i2 = PnsAtMarker;
        PnToAlleles.clear();
        if (estimateMarkerSlippage)
        {
            markerSlippageOut << thisMarker.chrom << "\t" << thisMarker.start << "\t" << thisMarker.end << "\t" << thisMarker.motif << "\t" << thisMarker.refRepeatNum << "\t" << thisMarker.refRepSeq << "\t" << setprecision(4) << fixed << markerToSizeAndModel[it->first].i1.i2 << "\t" << nAvailable<< endl;
            cout << thisMarker.start << " totalSlipp: " << setprecision(4) << fixed << markerToSizeAndModel[it->first].i1.i2 << " hiSeqSlipp: " << setprecision(4) << fixed << markerToSizeAndModelHiSeq[it->first].i1.i2 << " hiSeqXSlipp: " << setprecision(4) << fixed << markerToSizeAndModelHiSeqX[it->first].i1.i2 << endl;
         }
        cout << "Finished marker number: " << z << endl;
        ++z;                 
    }
    String<Pair<float> > labels;
    if (estimateMarkerSlippage)
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
                modelAndLabelDir = argv[7];
                append(modelAndLabelDir, PnIds[i]);
                append(modelAndLabelDir, "labels");
                append(modelAndLabelDir, prevItNumStr);
                append(modelAndLabelDir, "_");
                append(modelAndLabelDir, intervalIndex);               
                if (remove(toCString(modelAndLabelDir)) !=0)
                    cout << "Remove operation of old labelFile failed with code " << errno << endl;
            }*/
            modelAndLabelDir = argv[7];
            labelsOut.close();
        }
    }
    cout << "Finished determining genotypes" << endl;
    return 0;
}
