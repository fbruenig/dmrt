
#include "dmrtalg.h"
#include "dmrtreader.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc,char* argv[])
{
    if(argc < 6)
    {
        cout << "Not enough arguments given, please enter as follows:\n\ndmrt inputfile.txt outputfile.txt lowCutoff[nm] deltaR[nm] highCutoff[nm]\n" << endl;
        return 0;
    }

    cout << "Starting..." << endl;

    ifstream myfile;
    myfile.open(argv[1]);
    if(!myfile)
    {
        cout << "ERROR cannot open input file!" << endl;
        return 0;
    }

    ofstream outfile;
    outfile.open(argv[2]);
    if(!outfile)
    {
        cout << "ERROR cannot open output file!" << endl;
        return 0;
    }


    dmrtReader reader = dmrtReader(&myfile);
    vector< vector<float> >* vec = reader.read2Dvector();

    if (!vec || vec->size()== 0)
    {
        cout << "ERROR no data was read!" << endl;
        return 0;
    }


    dmrtAlg eval = dmrtAlg();
    vector<vector<float> > * dmrt = eval.getDMRTfrom2DVector(vec,atof(argv[5]),atof(argv[3]),atof(argv[4]));

    if (!dmrt || dmrt->size()== 0)
    {
        cout << "ERROR algorithm did not produce results!" << endl;
        return 0;
    }
    else
    {
        int c = 0;
        for(size_t i = 0; i< dmrt->size(); i++)
        {
            if((*dmrt)[i][1]!=(*dmrt)[i][1]) c++;
        }
        if (c == dmrt->size())
        {
            cout << "ERROR algorithm did not produce results!" << endl;
            return 0;
        }
    }
    reader.print2DVectorToXVG(dmrt,&outfile);

    //vector<vector<float>> * fpts = eval.getFPTfrom2DVector(vec,10.0);

    outfile.close();
    myfile.close();
    return 0;
}

