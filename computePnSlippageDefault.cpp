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
#include <liblinear/linear.h>
#include <liblinear/linear.cpp>
#include <liblinear/tron.h>
#include <liblinear/tron.cpp>
#include <liblinear/blas/blas.h>
#include <liblinear/blas/blasp.h>
#include <liblinear/blas/daxpy.c>
#include <liblinear/blas/ddot.c>
#include <liblinear/blas/dnrm2.c>
#include <liblinear/blas/dscal.c>

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
    double pValue;
} ; 

//So I can map from Markers
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//Sums the pValues of reads at every marker, also stores the current slippage rate value for each marker
map<Marker, Pair<double> > markerToPSumSlipp;

//Stores the logistic regression model definition for each marker, is cleared for each chromosome
map<Marker, model*> markerToModel;

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
        markerSlippageFile >> tempVal;
        markerSlippageFile >> tempVal;
        markerSlippageFile >> tempVal;
        markerSlippageFile >> markerToPSumSlipp[currMarker].i2;
    }
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
void parseNextLine(float winner, float second, ifstream& attributeFile, Marker& marker, String<string> firstLine, bool useFirstLine, LabelProps& slippCount)
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
    markerToPSumSlipp[marker].i1 += currentLine.pValue;
    if (currentLine.numOfRepeats == winner || currentLine.numOfRepeats == second)    
        slippCount.p1 += currentLine.pValue;    
    else 
    {
        if ((currentLine.numOfRepeats == winner - 1) || (currentLine.numOfRepeats == second - 1))        
            slippCount.p2 += currentLine.pValue;        
        else
        {
            slippCount.p3 += currentLine.pValue;
        }
    }
}

int main(int argc, char const ** argv)
{   
    //Check arguments.
    if (argc != 6)
    {
        cerr << "USAGE: " << argv[0] << " pathToAttributeFile/ PN-id outputFile markerSlippageFile modelDirectory/";
        return 1;
    }
    CharString attDir = argv[1], modelDir = argv[5];
    string pnId = argv[2];
    append(attDir,pnId);
    ifstream attributeFile(toCString(attDir)), slippageFile(argv[4]);
    ofstream outputFile(argv[3]);  
    LabelProps slippCount;
    slippCount.p1 = 0;
    slippCount.p2 = 0;
    slippCount.p3 = 0;
    double slippage; 
    string nextLine;
    Marker marker;
    int numberOfReads;
    float winner, second;
    Pair<int, String<string> > numberOfWordsAndWords;        
    if(attributeFile.fail())
    {
        cout << "Unable to locate attributes file." << endl;
        return 1;            
    }
    attributeFile >> pnId;
    cout << "Computing pnSlippage for: " << pnId << endl;
    readMarkerSlippage(slippageFile);
    cout << "Finished reading the markerSlippage." << endl;
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
                append(modelDir, marker.chrom);
                stringstream startStr;
                startStr << marker.start;
                append(modelDir, startStr.str());
                startStr.clear();
                startStr.str("");
                const char *model_in_file = toCString(modelDir);                    
                markerToModel[marker] = load_model(model_in_file);
                modelDir = argv[5];
            }
            else
                markerToPSumSlipp[marker].i1 = -1.0;
        }
        if (numberOfWordsAndWords.i1 == 11)
        {
            if (numberOfReads >= 10)
            {
                for (unsigned i = 0; i < numberOfReads; ++i)
                {
                    if (i == 0)
                        parseNextLine(winner, second, attributeFile, marker, numberOfWordsAndWords.i2, true, slippCount);
                    else 
                        parseNextLine(winner, second, attributeFile, marker, numberOfWordsAndWords.i2, false, slippCount);
                }
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
    cout << "Finished reading attribute file." << endl;                              
    attributeFile.close();
    markerToModel.clear();
    
    int nMissing = 0;
    vector<double> numerators;
    double currMarkSlipp, currPvalSum, currNumerator, denominator;
    double finalSub = 0;
    map<Marker, Pair<double> >::const_iterator markerEnd = markerToPSumSlipp.end(); 
    for (map<Marker, Pair<double> >::iterator markerStart = markerToPSumSlipp.begin(); markerStart != markerEnd; ++markerStart)
    {            
        if (markerToPSumSlipp[markerStart->first].i1 == -1.0)
        {
            ++nMissing;
            continue;
        }
        if (markerStart->second.i2 == 0)
            currMarkSlipp = 0.001;
        else
            currMarkSlipp = markerStart->second.i2;    
        currPvalSum = markerStart->second.i1;
        currNumerator = (1.0/(currMarkSlipp*(1.0-currMarkSlipp))) * (currPvalSum/(slippCount.p1 + slippCount.p2 + slippCount.p3));
        numerators.push_back(currNumerator);
    }
    denominator = accumulate(numerators.begin(),numerators.end(),0.0);
    unsigned index = 0;
    for (map<Marker, Pair<double> >::iterator markerStart = markerToPSumSlipp.begin(); markerStart != markerEnd; ++markerStart)
    {
        if (markerToPSumSlipp[markerStart->first].i1 == -1.0)
            continue;
        finalSub += markerStart->second.i2 * numerators[index]/denominator;
        ++index;    
    }
    slippage = (slippCount.p2 + slippCount.p3)/(slippCount.p1 + slippCount.p2 + slippCount.p3) - finalSub;
    cout << "Number of markers available for estimating pnSlippage for " << pnId << " is: " << markerToPSumSlipp.size() - nMissing << endl;
    outputFile << pnId << "\t" << slippage << endl;    
    return 0;    
}
