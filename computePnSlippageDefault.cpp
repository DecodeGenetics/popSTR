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
    int nMarkers;
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

//For storing marker values
struct MarkerStats
{
    model* regressionModel;
    std::map<string, double> pnToPsum;
    double slippage;
    unsigned nAlleles;
    double stutter;
    double fullMotifSlippageSum;
    unsigned nPns;
    double posSlippProb;
    double negSlippProb;
};

//So I can map from Markers
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//stores various marker specific values
map<Marker, MarkerStats> markerToStats;

//Store AttributeLine for all reads
map<Pair<string, Marker>, String<AttributeLine> > markerAndPnToReads;

//Store alleles at each marker
map<Marker, Pair<std::set<float>, String<Pair<float> >> > markerToAllelesAndGenotypes;

//For storing command line arguments
struct ComputePnSlippageOptions
{
    CharString pnList, attributesDirectory, outputFile, markerSlippageFile, modelDirectory;
    unsigned firstPnIdx;
} ;

ArgumentParser::ParseResult parseCommandLine(ComputePnSlippageOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("computePnSlippageDefault");
    setShortDescription(parser, "Compute individual specific slippage rate for the provided individuals.");
    setVersion(parser, "1.4");
    setDate(parser, "April 2019");
    addUsageLine(parser, "\\fI-PL\\fP pnList \\fI-AD\\fP attributesDirectory \\fI-OF\\fP outputFile \\fI-FP\\fP firstPnIdx \\fI-MS\\fP markerSlippageFile \\fI-MD\\fP regressionModelDirectory");
    addDescription(parser, "This program will estimate an individual specific slipppage rate for the individuals specified based on the marker slippage rates and models provided.");

    addOption(parser, ArgParseOption("PL", "pnList", "A list of PNs whose slippage will be estimated.", ArgParseArgument::INPUT_FILE, "IN-FILE"));
    setRequired(parser, "pnList");

    addOption(parser, ArgParseOption("AD", "attributesDirectory", "Path to attributes directory.", ArgParseArgument::INPUT_FILE, "IN-FILE"));
    setRequired(parser, "attributesDirectory");

    addOption(parser, ArgParseOption("OF", "outputFile", "The slippage rate estimated will be appended to this file.", ArgParseArgument::INPUT_FILE, "OUT-FILE"));
    setRequired(parser, "outputFile");

    addOption(parser, ArgParseOption("FP", "firstPnIdx", "Index of first Pn in pnList within the attributeFile.", ArgParseArgument::INTEGER, "INTEGER"));
    setRequired(parser, "firstPnIdx");

    addOption(parser, ArgParseOption("MS", "markerSlippageFile", "A file containing slippage rates for the microsatellites.", ArgParseArgument::OUTPUT_FILE, "OUT-FILE"));
    setRequired(parser, "markerSlippageFile");

    addOption(parser, ArgParseOption("MD", "modelDirectory", "A directory where logistic regression models for all markers in the markerSlippageFile are stored.", ArgParseArgument::OUTPUT_FILE, "IN-DIR"));
    setRequired(parser, "modelDirectory");
	
	ArgumentParser::ParseResult res = parse(parser, argc, argv);
	
	if (res != ArgumentParser::PARSE_OK)
	    return res;
	
    getOptionValue(options.pnList, parser, "pnList");
    getOptionValue(options.attributesDirectory, parser, "attributesDirectory");
	getOptionValue(options.outputFile, parser, "outputFile");
    getOptionValue(options.firstPnIdx, parser, "firstPnIdx"); 
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
    myProb.x[idx][6].index = -1; // This is to indicate that there aren't any more attributes to read in.
    myProb.x[idx][6].value = 0;
}

void readMarkerSlippage(ifstream& markerSlippageFile, CharString regressionModelDirectory)
{
    Marker currMarker;
    string tempVal;
    CharString currMarkerModelDir = regressionModelDirectory;
    //unsigned i = 1;
    while (!markerSlippageFile.eof())
    {
        markerSlippageFile >> currMarker.chrom;
        markerSlippageFile >> currMarker.start;
        markerSlippageFile >> currMarker.end;
        markerSlippageFile >> currMarker.motif;
        markerSlippageFile >> tempVal;
        markerSlippageFile >> markerToStats[currMarker].slippage; //marker slippage rate
        markerSlippageFile >> markerToStats[currMarker].nPns; //how many pns available to estimate the marker slippage
        markerSlippageFile >> markerToStats[currMarker].nAlleles; //read number of alleles
        markerSlippageFile >> markerToStats[currMarker].stutter; //marker stutter rate
        markerSlippageFile >> markerToStats[currMarker].posSlippProb; //probability of adding repeats in slippage
        markerSlippageFile >> markerToStats[currMarker].negSlippProb; //probability of losing repeats in slippage
        append(currMarkerModelDir, "/model_");
        append(currMarkerModelDir, to_string(currMarker.start));
        append(currMarkerModelDir, "_");
        append(currMarkerModelDir, currMarker.motif);
        const char *model_in_file = toCString(currMarkerModelDir);
        markerToStats[currMarker].regressionModel = load_model(model_in_file);
        currMarkerModelDir = regressionModelDirectory;
        //if (i % 1000==0 && i>0)
            //cout << "Working on marker number: " << i << endl;
        //i++;
    }
    cout << "Finished reading marker slippage." << endl;
}

//Count number of words in a sentence, use to parse input from attribute file
Pair<int, String<string> > countNumberOfWords(string sentence)
{
    int numberOfWords = 0;
    String<string> words;
    resize(words, 9);

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
    model* model_ = markerToStats[marker].regressionModel;
    double *prob_estimates = (double *) malloc(2*sizeof(double));
    problem prob;
    prob.bias = -1;
    prob.l = 1;
    prob.n = 6;
    prob.x = (feature_node **) malloc(prob.l * sizeof(feature_node *));
    prob.x[0] = (feature_node *) malloc(10 * sizeof(feature_node));
    fillProblemX(0, currentLine, prob);
    predict_label = predict_probability(model_,prob.x[0],prob_estimates);
    return prob_estimates[0];
}

//Parses one line from attribute file by filling up and returning an AttributeLine, also initializes markerToSizeAndModel map using the labels
AttributeLine parseNextLine(float winner, float second, ifstream& attributeFile, Marker& marker, map<string, LabelProps>& pnToLabelProps, string pnId)
{
    AttributeLine currentLine;
    string temp;
    attributeFile >> currentLine.numOfRepeats;
    attributeFile >> currentLine.ratioBf;
    attributeFile >> currentLine.ratioAf;
    attributeFile >> currentLine.locationShift;
    attributeFile >> currentLine.mateEditDist;
    attributeFile >> currentLine.purity;
    attributeFile >> currentLine.ratioOver20In;
    attributeFile >> temp;
    currentLine.pValue = getPval(marker, currentLine);
    markerToStats[marker].pnToPsum[pnId] += currentLine.pValue;
    if (currentLine.numOfRepeats == winner || currentLine.numOfRepeats == second)
    {
        pnToLabelProps[pnId].p1 += currentLine.pValue;
        currentLine.label = 1;
    }
    else
    {
        float diff1 = fabs(currentLine.numOfRepeats - winner), diff2 = fabs(currentLine.numOfRepeats - second);
        if (std::min(diff1,diff2)>=0.9) //full motif slippage
        {
            pnToLabelProps[pnId].p2 += currentLine.pValue;
            currentLine.label = 2;
        }
        else //stutter
        {
            pnToLabelProps[pnId].p3 += currentLine.pValue;
            currentLine.label = 3;
        }
    }
    return currentLine;
}

double estimateSlippage(double current_sp, string pnId)
{
    vector<double> weights;
    vector<double> slippFragments;
    double currMarkSlipp, currPvalSum, weightSum = 0, fullMotifSlippageSum = 0;
    for (auto& marker: markerToStats)
    {
        if ( marker.second.pnToPsum[pnId] == 0.0)
            continue;
        if (marker.second.slippage == 0)
            currMarkSlipp = 0.001;
        else
            currMarkSlipp = marker.second.slippage;
        currPvalSum = marker.second.pnToPsum[pnId];
        weights.push_back(currPvalSum/((current_sp+currMarkSlipp)*(1-(current_sp+currMarkSlipp))));
    }
    weightSum = accumulate(weights.begin(),weights.end(),0.0);
    unsigned index = 0;
    for (auto& marker: markerToStats)
    {
        if ( marker.second.pnToPsum[pnId] == 0.0)
            continue;
        Pair<string, Marker> mapKey = Pair<string, Marker>(pnId, marker.first); 
        String<AttributeLine> readsAtI = markerAndPnToReads[mapKey];
        for (unsigned i=0; i<length(readsAtI); ++i)
        {
            if (readsAtI[i].label == 2)
                fullMotifSlippageSum += readsAtI[i].pValue;
        }
        slippFragments.push_back((weights[index]/weightSum)*((fullMotifSlippageSum/marker.second.pnToPsum[pnId]) - marker.second.slippage));
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

bool determineGenotype(String<AttributeLine>& reads, double s_ij, String<Pair<float> > genotypes, int numberOfAlleles, int motifLength, double psucc, double posSlippProb, double negSlippProb)
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
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 < -0.9)
                    posNegSlipp = negSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 > 0.9)
                    posNegSlipp = posSlippProb;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                probs[i] *= (readToCheck.pValue * dgeom(static_cast<int>(roundf((diff-(float)floor(diff))*motifLength)), psucc) * dpois(floor(diff), lambda) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
            }
            else
            {
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 < -0.9)
                    posNegSlipp = negSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 > 0.9)
                    posNegSlipp = posSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i2 < -0.9)
                    posNegSlipp2 = negSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i2 > 0.9)
                    posNegSlipp2 = posSlippProb;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                diff2 = fabs(readToCheck.numOfRepeats - genotypeToCheck.i2);
                probs[i] *= (readToCheck.pValue * 0.5 * (dgeom(static_cast<int>(roundf((diff-(float)floor(diff))*motifLength)), psucc) * dpois(floor(diff), lambda) * posNegSlipp + dgeom(static_cast<int>(roundf((diff2-(float)floor(diff2))*motifLength)), psucc) * dpois(floor(diff2), lambda) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
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

int updateGenotypes(double current_sp, string pnId)
{
    int nChanged = 0;
    bool changed;
    String<Pair<float> > genotypes;
    for (auto& marker: markerToStats)
    {
        if ( marker.second.pnToPsum[pnId] == 0.0)
            continue;
        genotypes = markerToAllelesAndGenotypes[marker.first].i2;
        Pair<string, Marker> mapKey = Pair<string, Marker>(pnId, marker.first);
        changed = determineGenotype(markerAndPnToReads[mapKey], current_sp+marker.second.slippage, genotypes, marker.second.nAlleles, marker.first.motif.size(), marker.second.stutter, marker.second.posSlippProb, marker.second.negSlippProb);
        if (changed)
            ++nChanged;
    }
    return nChanged;
}

map<string, LabelProps> readPnList(CharString & pnInfoFile)
{
    map<string, LabelProps> pnToLabelProps;
    ifstream pnList(toCString(pnInfoFile));
    while (!pnList.eof())
    {
        string PN_ID;
        pnList >> PN_ID;
        if (PN_ID.length() == 0 || pnList.eof())
            break;
        LabelProps slippCount;
        slippCount.p1 = 0;
        slippCount.p2 = 0;
        slippCount.p3 = 0;
        pnToLabelProps[PN_ID] = slippCount;
    }
    cout << "Finished reading PnList.\n";
    return pnToLabelProps;
}

long int readOffSets(ifstream & attsFile, unsigned firstPnIdx, unsigned nPns)
{
    long int offset;
    for (unsigned i = 1; i<=firstPnIdx; ++i)
    {
        attsFile >> offset;
        if (offset == -69)
            return 0;
    }
    while (offset == 0 && !attsFile.eof())
    {
        ++firstPnIdx;
        attsFile >> offset;
    }
    //cout << "firstPnIdx: " << firstPnIdx << "\n";
    //cout << "offset: "<< offset << "\n";
    if (firstPnIdx > nPns)
        return 0;
    else
        return offset;
}

void readMarkerData(CharString attributesDirectory, Marker marker, map<string, LabelProps>& pnToLabelProps, unsigned firstPnIdx)
{
    //variables
    int numberOfReads, pnsFound = 0;
    float winner, second, numOfRepeats;
    string nextLine, temp;
    Pair<int, String<string> > numberOfWordsAndWords;
    AttributeLine currentLine;
    //make input stream
    append(attributesDirectory, "/");
    append(attributesDirectory, to_string(marker.start));
    append(attributesDirectory, "_");
    append(attributesDirectory, marker.motif);
    //cout << "Reading data from " << attributesDirectory << endl;
    ifstream attsFile(toCString(attributesDirectory));
    long int offset = readOffSets(attsFile, firstPnIdx, firstPnIdx + pnToLabelProps.size() - 1);
    if (offset != 0)
        attsFile.seekg(offset);
    else 
        return;
    while (!attsFile.eof() && pnsFound < pnToLabelProps.size())
    {
        getline (attsFile,nextLine);
        if (nextLine.length() == 0)
            continue;
        numberOfWordsAndWords = countNumberOfWords(nextLine);
        if (numberOfWordsAndWords.i1 == 1)
        {
            //first check if we passed the last pn in our map
            if (nextLine > pnToLabelProps.rbegin()->first)
                break;
            if (pnToLabelProps.count(nextLine) != 0)
            {
                ++pnsFound;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> numberOfReads;
                attsFile >> winner;
                attsFile >> second;
                //Just use markers where I have more than 10 reads
                if (numberOfReads >= 10)
                {
                    ++pnToLabelProps[nextLine].nMarkers;
                    for (unsigned i = 0; i < numberOfReads; ++i)
                    {
                        currentLine = parseNextLine(winner, second, attsFile, marker, pnToLabelProps, nextLine);
                        Pair<string,Marker> mapKey = Pair<string,Marker>(nextLine, marker);
                        appendValue(markerAndPnToReads[mapKey],currentLine);
                        markerToAllelesAndGenotypes[marker].i1.insert(currentLine.numOfRepeats);
                    }
                }
                //Not enough reads
                else
                {
                    markerToStats[marker].pnToPsum[nextLine] = 0.0;
                    for (unsigned i = 0; i <= numberOfReads; ++i)
                        getline (attsFile,nextLine);
                }
            }
           // Don't want this pn, walk on by
            else
            {
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> numberOfReads;
                attsFile >> temp;
                attsFile >> temp;
                for (unsigned i = 0; i <= numberOfReads; ++i)
                    getline (attsFile,nextLine);
            }
        }
        else
            cerr << "Something went sideways while reading attributes @: " << attributesDirectory << "\n";
    }
    markerToAllelesAndGenotypes[marker].i2 = makeGenotypes(markerToAllelesAndGenotypes[marker].i1);
}

int main_3(int argc, char const ** argv)
{
    ComputePnSlippageOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
	    return res == seqan::ArgumentParser::PARSE_ERROR;

    CharString modelDir = options.modelDirectory;
    ifstream slippageFile(toCString(options.markerSlippageFile));
    if(slippageFile.fail())
    {
        cout << "Unable to locate slippageFile: " << options.markerSlippageFile << endl;
        return 1;
    }
    else
        readMarkerSlippage(slippageFile, options.modelDirectory);
    //Read pn list
    map<string, LabelProps> pnToLabelProps = readPnList(options.pnList);
    //make output file
    ofstream outputFile;
    outputFile.open(toCString(options.outputFile), ios_base::app);
    if(outputFile.fail())
    {
        cout << "Unable to create output file." << endl;
        return 1;
    }
    for (auto &marker: markerToStats)
    {
        readMarkerData(options.attributesDirectory, marker.first, pnToLabelProps, options.firstPnIdx);
    }
    cout << "Finished reading attributes data." << endl;
    for (auto& pn: pnToLabelProps)
    {
        if (pn.second.nMarkers < 1)
        {
            cout << "No markers with more than minimum number of reads for: " << pn.first << endl;
            outputFile << pn.first << "\t" << 0.0 << endl;
            continue;
        }
        double current_sp = (0.5*pn.second.p2)/(pn.second.p1 + pn.second.p2 + pn.second.p3);
        double changed = 1, nChanged = 0;
        while (changed > 0.005)
        {
            current_sp = estimateSlippage(current_sp, pn.first);
            cout << "Estimated slippage." << endl;
            nChanged = updateGenotypes(current_sp, pn.first);
            cout << "Updated genotypes." << endl;
            changed = (float)nChanged/(float)pn.second.nMarkers;
            cout << nChanged << " " << pn.second.nMarkers << endl;
        }
        cout << "Number of markers available for estimating pnSlippage for " << pn.first << " is: " << pn.second.nMarkers << endl;
        outputFile << pn.first << "\t" << current_sp << endl;
    }
    return 0;
}
