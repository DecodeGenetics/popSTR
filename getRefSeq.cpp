#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

using namespace std;
using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 4)
    {
        std::cerr << "USAGE: refGenome.fa markerFile outputDirectory\n";
        return 1;
    }
    
    // Try to load index and create on the fly if necessary.
    FaiIndex faiIndex;
    if (read(faiIndex, argv[1]) != 0)
    {
        if (build(faiIndex, argv[1]) != 0)
        {
            cerr << "ERROR: Index could not be loaded or built.\n";
            return 1;
        }
        if (write(faiIndex) != 0)  // Name is stored from when reading.
        {
            cerr << "ERROR: Index could not be written do disk.\n";
            return 1;
        }
    }
    
    //Open file containing marker locations
    ifstream markerFile(argv[2]);
    
    //Check which chromosome I'm working on
    string chromString;
    markerFile >> chromString;
    
    // Translate sequence name to index.
    unsigned idx = 0;
    if (!getIdByName(faiIndex, chromString, idx))
    {
        cerr << "ERROR: Index does not know about sequence " << chromString << "\n";
        return 1;
    }
    
    //Create output stream
    string outputDirectory = argv[3];
    outputDirectory.append(chromString);
    outputDirectory.append("markerInfo");
    const char * outFile = outputDirectory.c_str();
    ofstream outputFile(outFile);
    
    while (!markerFile.eof())
    {
        string motifString;
        Dna5String refBf, refAf, refRepSeq;
        int beginPos = 0, endPos = 0;
        double refRepeatNum;                                   
        int LengthOfSeq = sequenceLength(faiIndex, idx);
        // Read marker coordinates, marker repeatMotif and number of repeats in reference.        
        markerFile >> beginPos;
        markerFile >> endPos;        
        markerFile >> motifString;        
        markerFile >> refRepeatNum;
        // Make sure begin and end pos are on the sequence and begin <= end.
        if (beginPos > sequenceLength(faiIndex, idx))
            beginPos = sequenceLength(faiIndex, idx);
        if (endPos > sequenceLength(faiIndex, idx))
            endPos = sequenceLength(faiIndex, idx);
        if (beginPos > endPos)
            endPos = beginPos;
        //Get 1000 bases infront of marker from reference
        if (readRegion(refBf, faiIndex, idx, max(0,beginPos-1001), beginPos-1) != 0)
        {
            cerr << "ERROR: Could not load reference before.\n";
            return 1;
        }
        //Get repeat sequence from reference
        if (readRegion(refRepSeq, faiIndex, idx, beginPos-1, endPos) != 0)
        {
            cerr << "ERROR: Could not load repeat sequence.\n";
            return 1;
        }
        //Get 1000 bases behind marker from reference
        if (readRegion(refAf, faiIndex, idx, endPos, min(LengthOfSeq,endPos+1000)) != 0)
        {
            cerr << "ERROR: Could not load reference after.\n";
            return 1;
        }
        outputFile << chromString << "\t" << beginPos << "\t" << endPos << "\t" << motifString << "\t" << refRepeatNum << "\t" << refBf << "\t" << refAf << "\t" << refRepSeq << endl;
        markerFile >> chromString;
        if (chromString.length() == 0)
            break;  
    }
    return 0;
}
