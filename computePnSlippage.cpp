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
#include <climits>
#include <liblinear.hpp>


namespace computePnSlippage
{

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

//For storing marker values
struct MarkerStats
{
    model* regressionModel;
    std::map<string, Pair<double> > pnToPandFullMotifSlippSums;
    double slippage;
    unsigned nPns;
};

struct AttributeLine {
    float ratioBf;
    float ratioAf;
    float numOfRepeats;
    unsigned locationShift;
    unsigned mateEditDist;
    float purity;
    float ratioOver20In;
    double pValue;
} ;

//So I can map from Markers
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//stores various marker specific values
map<Marker, MarkerStats> markerToStats;

//For storing command line arguments
struct ComputePnSlippageOptions
{
    CharString markerList, pnList, attributesDirectory, outputFile, iterationNumber, markerSlippageFile, regressionModelDirectory, previousSlippageRates;
    unsigned minPnsPerMarker, firstPnIdx;
} ;
ArgumentParser::ParseResult parseCommandLine(ComputePnSlippageOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("computePnSlippage");
    setShortDescription(parser, "Compute slippage rate for a list of individuals.");
    setVersion(parser, "1.4");
    setDate(parser, "June 2019");
    addUsageLine(parser, "\\fI-ML\\fP MarkerList \\fI-PL\\fP pnList \\fI-AD\\fP attributesDirectory \\fI-OF\\fP outputFile \\fI-IN\\fP iterationNumber \\fI-FP\\fP firstPnIdx \\fI-NPN\\fP minPnsPerMarker \\fI-MS\\fP markerSlippageFile \\fI-MD\\fP regressionModelDirectory");
    addDescription(parser, "This program will estimate an individual specific slipppage rate for the individuals specified based on the marker slippage rates and models specified.");

    addOption(parser, ArgParseOption("ML", "markerList", "List of marker to estimate slippage over.", ArgParseArgument::INPUT_FILE, "IN-FILE"));
    setRequired(parser, "markerList");

    addOption(parser, ArgParseOption("PL", "pnList", "A list of PNs whose slippage will be estimated. Required only when iteration number = 0", ArgParseArgument::INPUT_FILE, "IN-FILE"));

    addOption(parser, ArgParseOption("AD", "attributesDirectory", "Path to attributes files for the markers being used for slippage estimation.", ArgParseArgument::INPUTPREFIX, "IN-DIR"));
    setRequired(parser, "attributesDirectory");

    addOption(parser, ArgParseOption("OF", "outputFile", "The slippage rates estimated will be appended to this file.", ArgParseArgument::OUTPUT_FILE, "OUT-FILE"));
    setRequired(parser, "outputFile");

    addOption(parser, ArgParseOption("IN", "iterationNumber", "Index of the iteration being performed, 0-based.", ArgParseArgument::STRING, "INDEX"));
    setRequired(parser, "iterationNumber");

    addOption(parser, ArgParseOption("FP", "firstPnIdx", "Index of first Pn in pnList within the attributeFile.", ArgParseArgument::INTEGER, "INTEGER"));
    setRequired(parser, "firstPnIdx");

    //Optional parametets (i.e. for when iterationNumber > 0)
    addOption(parser, ArgParseOption("NPN", "minPnsPerMarker", "Minimum number of pns available at marker so it can be included in estimation.", ArgParseArgument::INTEGER, "INTEGER"));

    addOption(parser, ArgParseOption("MS", "markerSlippageFile", "A file containing slippage rates for the microsatellites, supplied when iterationNumber>0.", ArgParseArgument::INPUT_FILE, "IN-FILE"));

    addOption(parser, ArgParseOption("MD", "regressionModelDirectory", "A directory where logistic regression models for all markers in the markerList are stored, supplied when iterationNumber>0.", ArgParseArgument::INPUTPREFIX, "IN-DIR"));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.markerList, parser, "markerList");
    getOptionValue(options.attributesDirectory, parser, "attributesDirectory");
    getOptionValue(options.outputFile, parser, "outputFile");
    getOptionValue(options.iterationNumber, parser, "iterationNumber");
    getOptionValue(options.firstPnIdx, parser, "firstPnIdx");

    if (isSet(parser,"markerSlippageFile"))
    {
        getOptionValue(options.markerList, parser, "minPnsPerMarker");
        getOptionValue(options.markerSlippageFile, parser, "markerSlippageFile");
        getOptionValue(options.regressionModelDirectory, parser, "regressionModelDirectory");
    }
    else
        getOptionValue(options.pnList, parser, "pnList");

    return ArgumentParser::PARSE_OK;
}

//Count number of words in a sentence, use to parse input from attribute file
Pair<int, String<string> > countNumberOfWords(string& sentence){
    int numberOfWords = 0;
    String<string> words;
    resize(words, 9);
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

String<Marker> readMarkerList(CharString & markerInfoFile)
{
    String<Marker> markers;
    Marker currMarker;
    //Read all markers into memory
    ifstream markerFile(toCString(markerInfoFile));
    while (!markerFile.eof())
    {
        markerFile >> currMarker.chrom;
        if (currMarker.chrom.length() == 0 || markerFile.eof())
            break;
        markerFile >> currMarker.start;
        markerFile >> currMarker.end;
        markerFile >> currMarker.motif;

        appendValue(markers, currMarker);
    }
    return markers;
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
    return pnToLabelProps;
}

//Fills in the x-part of a problem structure from an AttributeLine structure
void fillProblemX(int idx, AttributeLine& currentLine, problem& myProb)
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

double getPval(Marker& marker, AttributeLine& currentLine)
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

//Parses one line from attribute file by filling up and returning an AttributeLine, also initializes markerToStats map using the labels
void parseNextLine(float winner, float second, ifstream& attributeFile, Marker& marker, string pn)
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
    markerToStats[marker].pnToPandFullMotifSlippSums[pn].i1 += currentLine.pValue;
    if (currentLine.numOfRepeats != winner && currentLine.numOfRepeats != second)
    {
        float diff1 = fabs(currentLine.numOfRepeats - winner), diff2 = fabs(currentLine.numOfRepeats - second);
        if (std::min(diff1,diff2)>=0.9)
            markerToStats[marker].pnToPandFullMotifSlippSums[pn].i2 += currentLine.pValue;
    }
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
    while (offset == 0)
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

void readMarkerData_level2(CharString attributesDirectory, Marker& marker, map <string, Pair<float> >& pnToLabels, map<string, Pair<double, unsigned> >& pnToPrevSlipp, unsigned minNpns, unsigned firstPnIdx)
{
    //variables
    int numberOfReads, pnsFound = 0;
    float winner, second, numOfRepeats;
    string nextLine, temp;
    Pair<int, String<string> > numberOfWordsAndWords;
    //make input stream
    append(attributesDirectory, "/");
    append(attributesDirectory, to_string(marker.start));
    append(attributesDirectory, "_");
    append(attributesDirectory, marker.motif);
    ifstream attsFile(toCString(attributesDirectory));
    long int offset = readOffSets(attsFile, firstPnIdx, firstPnIdx + pnToLabels.size() - 1);
    if (offset != 0)
        attsFile.seekg(offset);
    else
        return;
    while (!attsFile.eof() && pnsFound < pnToLabels.size())
    {
        getline (attsFile,nextLine);
        if (nextLine.length() == 0)
            continue;
        numberOfWordsAndWords = countNumberOfWords(nextLine);
        if (numberOfWordsAndWords.i1 == 1)
        {
            //first check if we passed the last pn in our map
            if (nextLine > pnToLabels.rbegin()->first)
                break;
            if (pnToLabels.count(nextLine) != 0)
            {
                ++pnsFound;
                winner = pnToLabels[nextLine].i1;
                second = pnToLabels[nextLine].i2;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> temp;
                attsFile >> numberOfReads;
                attsFile >> temp;
                attsFile >> temp;
                //Just use markers where I have more than 10 reads
                if (numberOfReads >= 10 && markerToStats[marker].nPns >= minNpns)
                {
                    for (unsigned i = 0; i < numberOfReads; ++i)
                        parseNextLine(winner, second, attsFile, marker, nextLine);
                }
                //Not enough reads
                else
                {
                    markerToStats[marker].pnToPandFullMotifSlippSums[nextLine].i1 = 0.0;
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
}

void readMarkerData(CharString attributesDirectory, Marker & marker, map<string, LabelProps> & pnToLabelProps, unsigned firstPnIdx)
{
    //make input stream
    append(attributesDirectory, "/");
    append(attributesDirectory, to_string(marker.start));
    append(attributesDirectory, "_");
    append(attributesDirectory, marker.motif);
    ifstream attsFile(toCString(attributesDirectory));
    //cout << "Reading data from: " << attributesDirectory << "\n";
    long int offset = readOffSets(attsFile, firstPnIdx, firstPnIdx + pnToLabelProps.size()-1);
    if (offset < 1)
        return;
    else
        attsFile.seekg(offset);
    //cout <<  "Finished seek command.\n";
    //variable declaration
    int numberOfReads, pnsFound = 0;
    float winner, second, numOfRepeats;
    string nextLine, temp;
    Pair<int, String<string> > numberOfWordsAndWords;
    //Go through file
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
            //check if PN is in pnToLabelProps and process if it is
            if (pnToLabelProps.count(nextLine) != 0)
            {
                ++pnsFound;
                attsFile >> temp;;
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
                        attsFile >> numOfRepeats;
                        attsFile >> temp;
                        attsFile >> temp;
                        attsFile >> temp;
                        attsFile >> temp;
                        attsFile >> temp;
                        attsFile >> temp;
                        attsFile >> temp;
                        if (numOfRepeats == winner || numOfRepeats == second)
                            pnToLabelProps[nextLine].p1 += 0.95;
                        else
                        {
                            float diff1 = fabs(numOfRepeats - winner), diff2 = fabs(numOfRepeats - second);
                            if (std::min(diff1,diff2)>=0.9)
                                pnToLabelProps[nextLine].p2 += 0.95;
                            else
                                pnToLabelProps[nextLine].p3 += 0.95;
                        }
                    }
                }
                //Not enough reads
                else
                {
                    for (unsigned i = 0; i <= numberOfReads; ++i)
                        getline (attsFile,nextLine);
                }
            }
            //Don't want this PN, walk on by
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
}

map<string, Pair<double, unsigned> > readPrevSlipp(CharString & previousSlippageRate, CharString & iterationNumber)
{
    CharString filePath = previousSlippageRate;
    string prevIterIdx;
    int itNum;
    lexicalCast(itNum, iterationNumber);
    --itNum;
    prevIterIdx = to_string(itNum);
    append(filePath, "_");
    append(filePath, prevIterIdx);
    ifstream prevSlipp(toCString(filePath));
    if(prevSlipp.fail())
        cout << "Unable to locate pnSlippageFile @ " << filePath << endl;

    map<string, Pair<double, unsigned> > pnToPrevSlipp;
    string PnId;
    int nMarkers;
    double currPnSlipp;
    while (!prevSlipp.eof())
    {
        prevSlipp >> PnId;
        prevSlipp >> currPnSlipp;
        prevSlipp >> nMarkers;
        pnToPrevSlipp[PnId].i1 = currPnSlipp;
        pnToPrevSlipp[PnId].i2 = nMarkers;
    }
    prevSlipp.close();
    return pnToPrevSlipp;
}

void readMarkerSlippage(CharString & markerSlippFile, CharString & itNumStr, CharString regressionModelDirectory)
{
    Marker currMarker;
    string tempVal;
    append(markerSlippFile, itNumStr);
    ifstream markerSlippageFile(toCString(markerSlippFile));
    CharString currMarkerModelDir = regressionModelDirectory;
    while (!markerSlippageFile.eof())
    {
        markerSlippageFile >> currMarker.chrom;
        markerSlippageFile >> currMarker.start;
        markerSlippageFile >> currMarker.end;
        markerSlippageFile >> currMarker.motif;
        markerSlippageFile >> tempVal;
        markerSlippageFile >> markerToStats[currMarker].slippage; //marker slippage rate
        markerSlippageFile >> markerToStats[currMarker].nPns; //how many pns available to estimate the marker slippage
        markerSlippageFile >> tempVal;
        markerSlippageFile >> tempVal;
        markerSlippageFile >> tempVal;
        markerSlippageFile >> tempVal;
        append(currMarkerModelDir, "/model_");
        append(currMarkerModelDir, to_string(currMarker.start));
        append(currMarkerModelDir, "_");
        append(currMarkerModelDir, currMarker.motif);
        const char *model_in_file = toCString(currMarkerModelDir);
        markerToStats[currMarker].regressionModel = load_model(model_in_file);
        currMarkerModelDir = regressionModelDirectory;
    }
}

double estimatePnSlippage(double current_sp, string pn)
{
    vector<double> weights;
    vector<double> slippFragments;
    double currMarkSlipp, currPvalSum, weightSum = 0, fullMotifSlippageSum = 0;
    for (auto& marker: markerToStats)
    {
        if ( marker.second.pnToPandFullMotifSlippSums[pn].i1 == 0.0)
            continue;
        if (marker.second.slippage == 0)
            currMarkSlipp = 0.001;
        else
            currMarkSlipp = marker.second.slippage;
        currPvalSum = marker.second.pnToPandFullMotifSlippSums[pn].i1;
        weights.push_back(currPvalSum/((current_sp+currMarkSlipp)*(1-(current_sp+currMarkSlipp))));
    }
    weightSum = accumulate(weights.begin(),weights.end(),0.0);
    unsigned index = 0;
    for (auto& marker: markerToStats)
    {
        if ( marker.second.pnToPandFullMotifSlippSums[pn].i1 == 0.0)
            continue;
        fullMotifSlippageSum = marker.second.pnToPandFullMotifSlippSums[pn].i2;
        slippFragments.push_back((weights[index]/weightSum)*((fullMotifSlippageSum/marker.second.pnToPandFullMotifSlippSums[pn].i1) - marker.second.slippage));
        fullMotifSlippageSum = 0;
        ++index;
    }
    double slippage = std::max(0.0,accumulate(slippFragments.begin(),slippFragments.end(),0.0));
    return slippage;
}

void readPnLabels(CharString modelAndLabelDir, CharString iterationNumber, map <string, Pair<float> >& pnToLabels, map<string, Pair<double, unsigned> >& pnToPrevSlipp, Marker& marker)
{
    append(modelAndLabelDir, "/");
    append(modelAndLabelDir, to_string(marker.start));
    append(modelAndLabelDir, "_");
    append(modelAndLabelDir, marker.motif);
    append(modelAndLabelDir, iterationNumber);
    ifstream labelFile(toCString(modelAndLabelDir));
    while (!labelFile.eof())
    {
        string PN_ID;
        float A1, A2;
        labelFile >> PN_ID;
        if (PN_ID.length() == 0 || PN_ID > pnToPrevSlipp.rbegin()->first)
            break;
        labelFile >> A1 >> A2;
        if (pnToPrevSlipp.count(PN_ID)!= 0)
        {
            Pair<float> labels = Pair<float>(A1, A2);
            pnToLabels[PN_ID] = labels;
        }
    }
    labelFile.close();
}

int main(int argc, char const ** argv)
{
    ComputePnSlippageOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    //read markers to use for slippage estimation
    String<Marker> markers = readMarkerList(options.markerList);
    cout << "Finished reading " << length(markers) << " markers.\n";
    //make output stream
    ofstream outputFile;
    CharString outputFilePath = options.outputFile;
    append(outputFilePath, "_");
    append(outputFilePath, options.iterationNumber);
    outputFile.open(toCString(outputFilePath));

    //if iterationNumber == 0 just read data and "initialize" estimation
    if (options.iterationNumber == "0")
    {
        //read pns to estimate slippage for
        map<string, LabelProps> pnToLabelProps = readPnList(options.pnList);
        cout << "Finished reading " << pnToLabelProps.size() << " pns.\n";

        //Loop over supplied markers and look for entries of pns in the pnList
        for (unsigned i = 0; i < length(markers); ++i)
        {
            readMarkerData(options.attributesDirectory, markers[i], pnToLabelProps, options.firstPnIdx);
            if (i % 1000 == 0 && i>0)
                cout << "Finished marker number " << i << endl;
        }

        //Loop over pns and print initial slippRates estimates to output file
        for (auto& pn: pnToLabelProps)
        {
            double slippage = (0.5*pn.second.p2)/(pn.second.p1 + pn.second.p2 + pn.second.p3);
            outputFile << pn.first << "\t" << slippage << "\t" << pn.second.nMarkers << "\n";
        }
    }
    else
    {
        map<string, Pair<float> > pnToLabels; //map to store labels at each marker, must clear for each marker
        //read previous slippage
        map<string, Pair<double, unsigned> > pnToPrevSlipp = readPrevSlipp(options.outputFile, options.iterationNumber);
        cout << "Finished reading " << pnToPrevSlipp.size() << " pns.\n";
        //read marker slippage and models
        readMarkerSlippage(options.markerSlippageFile, options.iterationNumber, options.regressionModelDirectory);
        cout << "Finished reading markerSlippage and regression models." << endl;
        //Loop over markers and accumulate data
        for (unsigned i = 0; i < length(markers); ++i)
        {
            readPnLabels(options.regressionModelDirectory, options.iterationNumber, pnToLabels, pnToPrevSlipp, markers[i]);
            readMarkerData_level2(options.attributesDirectory, markers[i], pnToLabels, pnToPrevSlipp, options.minPnsPerMarker, options.firstPnIdx);
            if (i % 1000 == 0)
                cout << "Finished marker number " << i << endl;
            pnToLabels.clear();
        }
        //Loop over pns estimate slippage and print to output file
        for (auto& pn: pnToPrevSlipp)
        {
            double prevSlipp = pn.second.i1;
            double slippage = estimatePnSlippage(prevSlipp, pn.first);
            outputFile << pn.first << "\t" << slippage << "\t" << pn.second.i2<< "\n";
        }
    }

    return 0;
}

} // namespace computePnSlippage
