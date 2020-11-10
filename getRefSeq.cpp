#include <iostream>
#include <SeqAnHTS/include/seqan/sequence.h>
#include <SeqAnHTS/include/seqan/seq_io.h>
#include <SeqAnHTS/include/seqan/stream.h>

using namespace std;
using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 4)
    {
        std::cerr << "USAGE: " << argv[0] << " refGenome.fa markerFile outputDirectory\n";
        return 1;
    }
    
    // Try to load index and create on the fly if necessary.
    FaiIndex faiIndex;
    if (!open(faiIndex, argv[1]))
    {
        if (build(faiIndex, argv[1]) != 0)
        {
            cerr << "ERROR: Index could not be loaded or built.\n";
            return 1;
        }
        /*if (write(faiIndex) != 0)  // Name is stored from when reading.
        {
            cerr << "ERROR: Index could not be written do disk.\n";
            return 1;
        }*/
    }
    
    //Open file containing marker locations
    ifstream markerFile(argv[2]);
    
    //Check which chromosome I'm working on
    string chromString;
    markerFile >> chromString;
    
    // Translate sequence name to index.
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, chromString))
    {
        cerr << "ERROR: Index does not know about sequence " << chromString << "\n";
        return 1;
    }
    //Get length of chromosome
    int LengthOfSeq = sequenceLength(faiIndex, idx);
    
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
        if (beginPos==0 && endPos==0)
        {
            cout << "both positions are 0 " << beginPos << "-" << endPos << "\n";
            return 0;
        }
        if (abs(round(refRepeatNum*motifString.length()) - (endPos-beginPos+1)) > 1)
        {
            cout << "new region is not equal in length to old for " << chromString << ":" << beginPos << "-" << endPos << "\n";
            markerFile >> chromString;
            continue;
        }
        //Get 1000 bases infront of marker from reference
        readRegion(refBf, faiIndex, idx, max(0,beginPos-1001), beginPos-1);
        if (length(refBf) != 1000)
        {
            if (beginPos > 1001)
            {
                cerr << "ERROR: Could not load reference before for " << chromString << ":" << beginPos << "-" << endPos << "\n";
                markerFile >> chromString;
                continue;
            }
            else
            {
                if (length(refBf) < beginPos -1)
                {
                    cerr << "ERROR: Could not load reference before for " << chromString << ":" << beginPos << "-" << endPos << "\n";
                    markerFile >> chromString;
                    continue;
                }
            }
        }
        //Get repeat sequence from reference
        readRegion(refRepSeq, faiIndex, idx, beginPos-1, endPos);
        if (!length(refRepSeq) > 0)
        {
            cerr << "ERROR: Could not load repeat sequence for " << chromString << ":" << beginPos << "-" << endPos << "\n";
            markerFile >> chromString;
            continue;
        }
        //Get 1000 bases behind marker from reference
        readRegion(refAf, faiIndex, idx, endPos, min(LengthOfSeq,endPos+1000));
        if (length(refAf) != 1000)
        {
            if (endPos+1000 < LengthOfSeq)
            {
                cerr << "ERROR: Could not load reference after for " << chromString << ":" << beginPos << "-" << endPos << "\n";
                markerFile >> chromString;
                continue;
            }
            else
            {
                if (length(refAf) < LengthOfSeq - endPos)
                {
                    cerr << "ERROR: Could not load reference before for " << chromString << ":" << beginPos << "-" << endPos << "\n";
                    markerFile >> chromString;
                    continue;
                }
            }
        }
        outputFile << chromString << "\t" << beginPos << "\t" << endPos << "\t" << motifString << "\t" << refRepeatNum << "\t" << refBf << "\t" << refAf << "\t" << refRepSeq << endl;
        markerFile >> chromString;
        if (chromString.length() == 0)
            break;  
    }
    return 0;
}
