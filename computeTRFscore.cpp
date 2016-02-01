#include <iostream>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>

using namespace std;
using namespace seqan;

//For repeating a motif n times
Dna5String repeat(Dna5String s, float n) {
    Dna5String ret;
    float diff = n - floor(n);
    for (int i = 0; i < floor(n)+1; i++) {
        append(ret,s);
    }
    if (diff > 0)
    {
        int t = round(diff*length(s));
        append(ret,prefix(s,t));
    }
    return ret;
}

int getPurity(Dna5String pureRepeat, Dna5String realRepeat)
{    
    //if (length(pureRepeat)!=length(realRepeat))
        //cout << "Pure and real repeats have unequal lengths! Real: " << length(realRepeat) << " and pure: " << length(pureRepeat) << endl;
    //Type definition for alignment structure
    typedef Dna5String TSequence;              
    typedef Align<TSequence,ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    
    //Create alignment structures, before and after
    TAlign alignStruct;
    resize(rows(alignStruct), 2);
    assignSource(row(alignStruct,0),pureRepeat);
    assignSource(row(alignStruct,1),realRepeat); 
    int alignmentScore = localAlignment(alignStruct, Score<int,Simple>(2,-7,-7));
    return alignmentScore;
}

int main(int argc, char const ** argv)
{ 
    //Check arguments.
    if (argc != 3)
    {
        cerr << "USAGE: " << argv[0] << " (infoFile with: motif refRepeatNum repRepeatSeq) output ";
        return 1;
    }
    ifstream inputFile(argv[1]);
    ofstream outputFile(argv[2]);
    string motifString, repeatString;
    float refRepeatNum;
    int trfScore;
    while (inputFile >> motifString >> refRepeatNum >> repeatString)
    {
        Dna5String pureRepeat = repeat(motifString, refRepeatNum);
        trfScore = getPurity(pureRepeat, repeatString);
        outputFile << trfScore << endl; 
    }
    return 0;
}
