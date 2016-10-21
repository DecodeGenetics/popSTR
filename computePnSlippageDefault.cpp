#include <iostream>
#include <set>
#include <map>
#include <string>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <seqan/stream.h>
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

using namespace std;
using namespace seqan;

//For storing number of members in each class and their pValue-sum
struct LabelProps {
    double p1;
    double p2;
    double p3;
} ;

//For storing marker information
struct Marker {
    string chrom;
    int start;
    int end;
    string motif;
} ;

struct AttributeLine {
    float ratioBf; 
    float ratioAf; 
    float numOfRepeats;
    unsigned locationShift; 
    unsigned mateEditDist;
    float purity; 
    float ratioOver20In;
    float ratioOver20After;
    unsigned sequenceLength;
    bool wasUnaligned;
    int label;
    double pValue;
} ; 

struct GenotypeInfo {
    String<Pair<float> > genotypes;
    String<long double> pValues;
    Pair<float> genotype;
    double pValue;
    int numOfReads;
    map<float, int> alleleToFreq; //This maps from reported alleles to their frequencies
    double pValueSum;
} ;

//For storing command line arguments 
struct ComputePnSlippageOptions
{
    CharString attDirChromNumPN, outputFile, markerSlippageFile, modelDirectory;
} ;

//So I can map from Markers
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//Sums the pValues of reads at every marker, also stores the current slippage rate value for each marker
map<Marker, Pair<int, String<double> > > markerToNallelesPSumSlippAndStutt;

//Stores the logistic regression model definition for each marker, is cleared for each chromosome
map<Marker, model*> markerToModel;

//Stores number of pns used in estimating marker slippage for each marker
map<Marker, int> markerToNpns;

//Store AttributeLine for all reads 
map<Marker, String<AttributeLine> > markerToReads;

//Store alleles at each marker
map<Marker, Pair<std::set<float>, String<Pair<float> >> > markerToAllelesAndGenotypes;

ArgumentParser::ParseResult parseCommandLine(ComputePnSlippageOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("computePnSlippageDefault");
    setShortDescription(parser, "Compute individual specific slippage rate for a single individual.");
    setVersion(parser, "1.3");
    setDate(parser, "October 2016");
    addUsageLine(parser, "\\fI-AF\\fP attributesDirectory/chromNum/PN-ID \\fI-OF\\fP outputFile \\fI-MS\\fP markerSlippageFile \\fI-MD\\fP modelDirectory ");
    addDescription(parser, "This program will estimate an individual specific slipppage rate for the individual specified based on the marker slippage rates and models specified.");
    
    addOption(parser, ArgParseOption("AF", "attributesDirectory/chromNum/PN-ID", "Path to attributes file for the individual to estimate a slippage rate for.", ArgParseArgument::INPUTFILE, "IN-FILE"));
    setRequired(parser, "attributesDirectory/chromNum/PN-ID");
    
    addOption(parser, ArgParseOption("OF", "outputFile", "The slippage rate estimated will be appended to this file.", ArgParseArgument::INPUTFILE, "OUT-FILE"));
    setRequired(parser, "outputFile");
    
    addOption(parser, ArgParseOption("MS", "markerSlippageFile", "A file containing slippage rates for the microsatellites.", ArgParseArgument::OUTPUTFILE, "OUT-FILE"));
    setRequired(parser, "markerSlippageFile");
    
    addOption(parser, ArgParseOption("MD", "modelDirectory", "A directory where logistic regression models for all markers in the markerSlippageFile are stored.", ArgParseArgument::OUTPUTFILE, "IN-DIR"));
    setRequired(parser, "modelDirectory");
	
	ArgumentParser::ParseResult res = parse(parser, argc, argv);
	
	if (res != ArgumentParser::PARSE_OK)
	    return res;
	    
	getOptionValue(options.attDirChromNumPN, parser, "attributesDirectory/chromNum/PN-ID");
	getOptionValue(options.outputFile, parser, "outputFile");
	getOptionValue(options.markerSlippageFile, parser, "markerSlippageFile");
	getOptionValue(options.modelDirectory, parser, "modelDirectory");

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

void readMarkerSlippage(ifstream& markerSlippageFile)
{
    Marker currMarker;
    string tempVal;        
    while (!markerSlippageFile.eof())
    {
        markerSlippageFile >> currMarker.chrom;
        markerSlippageFile >> currMarker.start;
        markerSlippageFile >> currMarker.end;
        markerSlippageFile >> currMarker.motif;
        markerSlippageFile >> tempVal;
        markerSlippageFile >> tempVal;
        resize(markerToNallelesPSumSlippAndStutt[currMarker].i2,3);
        markerToNallelesPSumSlippAndStutt[currMarker].i2[0] = -1.0;
        markerSlippageFile >> markerToNallelesPSumSlippAndStutt[currMarker].i2[1];
        markerSlippageFile >> markerToNpns[currMarker];
        markerSlippageFile >> markerToNallelesPSumSlippAndStutt[currMarker].i1;
        markerSlippageFile >> markerToNallelesPSumSlippAndStutt[currMarker].i2[2];
    }
    cout << "Finished reading marker slippage." << endl;
}

//Count number of words in a sentence, use to parse input from attribute file
Pair<int, String<string> > countNumberOfWords(string sentence){
    int numberOfWords = 0;
    String<string> words;
    resize(words, 11);
    
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

double getPval(Marker marker, AttributeLine currentLine)
{
    double predict_label;
    model* model_ = markerToModel[marker];
    double *prob_estimates = (double *) malloc(2*sizeof(double));
    problem prob;
    prob.bias = -1;
    prob.l = 1;
    prob.n = 9;
    prob.x = (feature_node **) malloc(prob.l * sizeof(feature_node *));
    prob.x[0] = (feature_node *) malloc(10 * sizeof(feature_node));
    fillProblemX(0, currentLine, prob);
    predict_label = predict_probability(model_,prob.x[0],prob_estimates);
    return prob_estimates[0];
}

//Parses one line from attribute file by filling up and returning an AttributeLine, also initializes markerToSizeAndModel map using the labels 
AttributeLine parseNextLine(float winner, float second, ifstream& attributeFile, Marker& marker, String<string> firstLine, bool useFirstLine, LabelProps& slippCount)
{
    AttributeLine currentLine;
    string temp;    
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
        attributeFile >> temp;
    }
    currentLine.pValue = getPval(marker, currentLine);
    markerToNallelesPSumSlippAndStutt[marker].i2[0] += currentLine.pValue;
    if (currentLine.numOfRepeats == winner || currentLine.numOfRepeats == second)
    {
        slippCount.p1 += currentLine.pValue;
        currentLine.label = 1;
    }
    else 
    {
        float diff1 = fabs(currentLine.numOfRepeats - winner), diff2 = fabs(currentLine.numOfRepeats - second);
        if (std::min(diff1,diff2)>=0.9)
        {
            slippCount.p2 += currentLine.pValue;
            currentLine.label = 2;
        }
        else
        {
            slippCount.p3 += currentLine.pValue;
            currentLine.label = 3;
        }
    }
    return currentLine;
}

double estimateSlippage(double current_sp)
{
    vector<double> weights;
    vector<double> slippFragments;
    double currMarkSlipp, currPvalSum, weightSum = 0, fullMotifSlippageSum = 0;
    map<Marker, Pair<int, String<double> > >::const_iterator markerEnd =  markerToNallelesPSumSlippAndStutt.end(); 
    for (map<Marker, Pair<int, String<double> > >::iterator markerStart =  markerToNallelesPSumSlippAndStutt.begin(); markerStart != markerEnd; ++markerStart)
    {            
        if ( markerStart->second.i2[0] == -1.0)        
            continue;
        if (markerStart->second.i2[1] == 0)
            currMarkSlipp = 0.001;
        else
            currMarkSlipp = markerStart->second.i2[1];
        currPvalSum = markerStart->second.i2[0];
        weights.push_back(currPvalSum/((current_sp+currMarkSlipp)*(1-(current_sp+currMarkSlipp))));
    }
    weightSum = accumulate(weights.begin(),weights.end(),0.0);
    unsigned index = 0;
    for (map<Marker, Pair<int, String<double> > >::iterator markerStart =  markerToNallelesPSumSlippAndStutt.begin(); markerStart != markerEnd; ++markerStart)
    {
        if ( markerStart->second.i2[0] == -1.0)
            continue;
        String<AttributeLine> readsAtI = markerToReads[markerStart->first];
        for (unsigned i=0; i<length(readsAtI); ++i)
        {
            if (readsAtI[i].label == 2)
                fullMotifSlippageSum += readsAtI[i].pValue;
        }
        slippFragments.push_back((weights[index]/weightSum)*((fullMotifSlippageSum/markerStart->second.i2[0]) - markerStart->second.i2[1]));
        fullMotifSlippageSum = 0;
        ++index;
    }
    double slippage = std::max(0.0,accumulate(slippFragments.begin(),slippFragments.end(),0.0));
    return slippage;
}

String<Pair<float> > makeGenotypes(std::set<float>& alleles)
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

int findMaxIndex(String<long double>& probs)
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

void relabelReads(std::set<float>& newGenotype, String<AttributeLine>& reads)
{
    for (unsigned i = 0; i < length(reads); ++i)
    {
        if (newGenotype.find(reads[i].numOfRepeats) != newGenotype.end())
            reads[i].label = 1;
        else
        {
            std::set<float>::iterator it=newGenotype.begin();
            float allele1 = *it;
            ++it;
            float allele2 = *it;
            float diff1 = fabs(allele1-reads[i].numOfRepeats), diff2 = fabs(allele2-reads[i].numOfRepeats);
            if (std::min(diff1,diff2)>=0.9)
                reads[i].label = 2;
            else
                reads[i].label = 3;
        }
    }
}

bool determineGenotype(String<AttributeLine>& reads, double s_ij, String<Pair<float> > genotypes, int numberOfAlleles, int motifLength, double psucc)
{
    Pair<float> genotypeToCheck;
    AttributeLine readToCheck;
    String<long double> probs;
    std::set<float> currentGenotype, newGenotype;
    resize(probs, length(genotypes));
    bool isHomo;
    float posNegSlipp = 1, posNegSlipp2 = 1, lambda = std::max((double)0.001,s_ij), diff, diff2;
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
            if (readToCheck.label == 1)
                currentGenotype.insert(readToCheck.numOfRepeats);
            if (isHomo)
            {
                if (readToCheck.numOfRepeats < genotypeToCheck.i1)
                    posNegSlipp = 0.85;
                if (readToCheck.numOfRepeats > genotypeToCheck.i1)
                    posNegSlipp = 0.15;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                probs[i] *= (readToCheck.pValue * dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * dpois(floor(diff), lambda) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));                
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
                probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>((diff-(float)floor(diff))*motifLength), psucc) * dpois(floor(diff), lambda) * posNegSlipp + dgeom(static_cast<int>((diff2-(float)floor(diff2))*motifLength), psucc) * dpois(floor(diff2), lambda) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
            }
        }
    } 
    indexOfWinner = findMaxIndex(probs);
    newGenotype.insert(genotypes[indexOfWinner].i1);
    newGenotype.insert(genotypes[indexOfWinner].i2);
    if (newGenotype == currentGenotype)
        return false;
    else
    {
        relabelReads(newGenotype,reads);
        return true;
    }
}

int updateGenotypes(double current_sp)
{
    int nChanged = 0;
    bool changed;
    String<Pair<float> > genotypes;
    map<Marker, Pair<int, String<double> > >::const_iterator markerEnd =  markerToNallelesPSumSlippAndStutt.end(); 
    for (map<Marker, Pair<int, String<double> > >::iterator markerStart =  markerToNallelesPSumSlippAndStutt.begin(); markerStart != markerEnd; ++markerStart)
    {
        if ( markerStart->second.i2[0] == -1.0)
            continue;
        genotypes = markerToAllelesAndGenotypes[markerStart->first].i2;
        changed = determineGenotype(markerToReads[markerStart->first], current_sp+markerStart->second.i2[1], genotypes, markerStart->second.i1, markerStart->first.motif.size(), markerStart->second.i2[2]);
        if (changed)
            ++nChanged;
    }
    return nChanged;
}

int main(int argc, char const ** argv)
{   
    ComputePnSlippageOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
	    return res == seqan::ArgumentParser::PARSE_ERROR;
    
    CharString modelDir = options.modelDirectory;
    ifstream attributeFile(toCString(options.attDirChromNumPN)), slippageFile(toCString(options.markerSlippageFile));
    if(attributeFile.fail())
    {
        cout << "Unable to locate attributes file." << options.attDirChromNumPN << endl;
        return 1;            
    }
    if(slippageFile.fail())
    {
        cout << "Unable to locate slippageFile: " << options.markerSlippageFile << endl;
        return 1;
    }
    ofstream outputFile;
    outputFile.open(toCString(options.outputFile), ios_base::app);
    if(outputFile.fail())
    {
        cout << "Unable to create output file." << endl;
        return 1;
    }
    string pnId;
    LabelProps slippCount;
    slippCount.p1 = 0;
    slippCount.p2 = 0;
    slippCount.p3 = 0; 
    string nextLine;
    Marker marker;
    int numberOfReads, nMarkers = 0, nChanged = 0;
    float winner, second;
    Pair<int, String<string> > numberOfWordsAndWords;            
    attributeFile >> pnId;
    readMarkerSlippage(slippageFile);
    AttributeLine currentLine;
    while (!attributeFile.eof())
    {        
        getline (attributeFile,nextLine);
        if (nextLine.length() == 0) 
            continue;
        numberOfWordsAndWords = countNumberOfWords(nextLine);
        if (numberOfWordsAndWords.i1 == 9) 
        {
            marker.chrom = numberOfWordsAndWords.i2[0];
            marker.start = lexicalCast<int>(numberOfWordsAndWords.i2[1]);
            marker.end = lexicalCast<int>(numberOfWordsAndWords.i2[2]);
            numberOfReads = lexicalCast<int>(numberOfWordsAndWords.i2[5]);
            winner = lexicalCast<float>(numberOfWordsAndWords.i2[7]);
            second = lexicalCast<float>(numberOfWordsAndWords.i2[8]);
            if (numberOfReads >= 10)
            {
                append(modelDir, "/");
                append(modelDir, marker.chrom);
                append(modelDir, "_");
                stringstream startStr;
                startStr << marker.start;
                append(modelDir, startStr.str());
                startStr.clear();
                startStr.str("");
                const char *model_in_file = toCString(modelDir);
                markerToModel[marker] = load_model(model_in_file);
                modelDir = options.modelDirectory;
            }
            else
                markerToNallelesPSumSlippAndStutt[marker].i2[0] = -1.0;
        }
        if (numberOfWordsAndWords.i1 == 11)
        {
            if (numberOfReads >= 10)
            {
                ++nMarkers;
                for (unsigned i = 0; i < numberOfReads; ++i)
                {
                    if (i == 0)
                        currentLine = parseNextLine(winner, second, attributeFile, marker, numberOfWordsAndWords.i2, true, slippCount);
                    else 
                        currentLine = parseNextLine(winner, second, attributeFile, marker, numberOfWordsAndWords.i2, false, slippCount);
                    appendValue(markerToReads[marker],currentLine);
                    markerToAllelesAndGenotypes[marker].i1.insert(currentLine.numOfRepeats);
                }
                markerToAllelesAndGenotypes[marker].i2 = makeGenotypes(markerToAllelesAndGenotypes[marker].i1);                
            }
            else 
            {
                for (unsigned i = 0; i < numberOfReads-1; ++i)
                    getline (attributeFile,nextLine);
            }
        }                
        if (numberOfWordsAndWords.i1 != 9 && numberOfWordsAndWords.i1 != 11)
            cerr << "Format error in attribute file!" << endl;
    }
    attributeFile.close();
    cout << "Finished reading attributes for: " << pnId << endl;
    double current_sp = (0.5*slippCount.p2)/(slippCount.p1 + slippCount.p2 + slippCount.p3);
    if (nMarkers < 1)
    {
        current_sp = 0.0;
        cout << "No markers with more than minimum number of reads for: " << pnId << endl;
        outputFile << pnId << "\t" << current_sp << endl;
        return 0;
    }
    double changed = 1;    
    while (changed > 0.005)
    {
        current_sp = estimateSlippage(current_sp);
        cout << "Estimated slippage." << endl;
        nChanged = updateGenotypes(current_sp);
        cout << "Updated genotypes." << endl;
        changed = (float)nChanged/(float)nMarkers;
        cout << nChanged << " " << nMarkers << endl;
    }      
    cout << "Number of markers available for estimating pnSlippage for " << pnId << " is: " << nMarkers << endl;
    outputFile << pnId << "\t" << current_sp << endl;    
    return 0;    
}
