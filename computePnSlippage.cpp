#include <iostream>
#include <set>
#include <map>
#include <string>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;

//For storing number of members in each class and their pValue-sum
struct LabelProps {
    double p1;
    double p2;
    double p3;
} ;

void appendChrAndPnId(CharString& dir, string chromNum, string pnId)
{
    append(dir,"chr");
    append(dir,chromNum);
    append(dir,"/");
    append(dir,pnId);
}

int main(int argc, char const ** argv)
{   
    //Check arguments.
    if (argc != 4)
    {
        cerr << "USAGE: " << argv[0] << " directoryToAttributesAndLabels PN-id outputFile";
        return 1;
    }
        
    ofstream outputFile;
    outputFile.open(argv[3], ios_base::app);
    CharString attDir = argv[1];
    append(attDir, "/attributes/");
    CharString labDir = argv[1];
    append(labDir, "/initialLabels/");
    CharString currAttDir, currLabDir;
    string pnId = argv[2];
    string temp, chrId;
    int numberOfReads;
    float winner, second, numOfRepeats;
    LabelProps slippCount;
    slippCount.p1 = 0;
    slippCount.p2 = 0;
    slippCount.p3 = 0;
    double slippage;
    
    for (unsigned i=1; i<25; ++i)
    {
        if (i<23)
        {
            stringstream ss;
            ss << i;
            chrId = ss.str();
        }
        if (i == 23)
            chrId = "X";
        if (i == 24)
            chrId = "Y";
        currAttDir = attDir;
        currLabDir = labDir;
        appendChrAndPnId(currAttDir, chrId, pnId);
        appendChrAndPnId(currLabDir, chrId, pnId);
        append(currAttDir, "attributes");
        append(currLabDir, "initialLabels");        
        ifstream attributeFile(toCString(currAttDir));
        ifstream initialLabels(toCString(currLabDir));
        if(attributeFile.fail() || initialLabels.fail())
        {
            continue;
        }
        cout << "Starting chromosome: " << chrId << endl;
        attributeFile >> pnId;
        while (!attributeFile.eof() && !initialLabels.eof())
        {            
            attributeFile >> temp;
            attributeFile >> temp;        
            attributeFile >> temp;
            attributeFile >> temp;
            attributeFile >> temp;
            attributeFile >> numberOfReads;
            attributeFile >> temp;
            initialLabels >> winner;
            initialLabels >> second;
            for (unsigned i = 0; i < numberOfReads; ++i)
            {
                attributeFile >> numOfRepeats;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                attributeFile >> temp;
                if (numOfRepeats == winner || numOfRepeats == second)
                    slippCount.p1 += 0.95;
                else 
                {
                    if ((numOfRepeats == winner - 1) || (numOfRepeats == second - 1))
                        slippCount.p2 += 0.95;
                    else
                        slippCount.p3 += 0.95;
                }
            }
        }
        attributeFile.close();
        initialLabels.close();
        cout << "Finished chromosome: " << chrId << endl;
    }
    slippage = (slippCount.p2 + slippCount.p3)/(2.0*(slippCount.p1 + slippCount.p2 + slippCount.p3));
    outputFile << pnId << "\t" << slippage << endl;
    return 0;    
}
