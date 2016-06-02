#include "dmrtmain.h"
#include "dmrtalg.h"
#include "dmrtalg2.h"
#include "dmrtreader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>

using namespace std;



dmrtMain::dmrtMain(const char *mode, bool verb)
{
    this->mVerb=verb;
    this->mMode=mode;
}

void dmrtMain::execute2(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, const char *input, const char *output, const double start, const double interval, const double end, const int dataColumn)
{

    if(this->mVerb){ cout << "Starting in" << input << endl;}

    ifstream myfile;
    myfile.open(input);
    if(!myfile)
    {
        cout << "ERROR cannot open input file!" << endl;
        return;
    }
    else
    {
        if(this->mVerb){ cout << "Opened input:" << input << endl;}
    }
    ofstream outfile;
    outfile.open(output);
    if(!outfile)
    {
        cout << "ERROR cannot open output file!" << endl;
        return;
    }
    else
    {
        if(this->mVerb){ cout << "Opened ouput:" << output << endl;}
    }

    dmrtReader reader = dmrtReader(&myfile,this->mVerb);
    dmrtalg2 eval = dmrtalg2(this->mMode,this->mVerb,end,start,interval, dataColumn);
    int vecLength = eval.getVecLength();

    // The final vectors get an extra row to save the the radii in, in BINS mode the last column will be empty
    // according to common histogramn convention, size(bins)=size(hist)+1

    //(*finalDmrts) = vector< vector<double> >(vecLength+2,vector<double>(vecLength+1,0.0));
    //(*finalCounts) = vector< vector<int> >(vecLength+1,vector<int>(vecLength+1,0));
    (*finalDmrts) = vector< vector<double> >(vecLength+1,vector<double>(vecLength,0.0));
    (*finalCounts) = vector< vector<int> >(vecLength,vector<int>(vecLength,0));

    vector<double> radii = eval.getRadii();

    for (int i=0;i<vecLength;i++)
    {
        (*finalDmrts)[vecLength+1][i]=radii[i];
    }

    bool success = true;
    int part = 0;
    while(success==true)
    {
        vector< vector<double> >* vec= new vector<vector<double> > ;
        vec = reader.read2DvectorSpace4gb(start-interval,end+interval);
        part++;
        if (!vec || vec->size()!= MAXDOUBLEVEC)
        {
            cout << "Vecsize: " << vec->size()<< endl;
            cout << "Reached end of file in this run" << endl;
            success = false;
        }
        if (mVerb){cout << "Read part "<< part << " of file. Evaluating..." << endl;}
        if (strncmp(this->mMode,"rt",2)==0)
        {
            if (strncmp(this->mMode+2,"bins",4)==0)
            {
                eval.getRTTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
            }
            else if (strncmp(this->mMode+2,"cross",5)==0)
            {
                eval.getRTTfrom2DVectorCross((*finalDmrts),(*finalCounts),vec);
            }
        }
        else if (strncmp(this->mMode,"mftp",4)==0)
        {
            if (strncmp(this->mMode+4,"bins",4)==0)
            {
                eval.getMFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
            }
            else if (strncmp(this->mMode+4,"cross",5)==0)
            {
                eval.getMFPTfrom2DVectorCross((*finalDmrts),(*finalCounts),vec);
            }
        }
        else if (strncmp(this->mMode,"cftp",4)==0)
        {
            if (strncmp(this->mMode+4,"bins",4)==0)
            {
                //eval.getFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
            }
            else if (strncmp(this->mMode+4,"cross",5)==0)
            {
                cout << "Initial config: "<< vecLength << endl;
                eval.getFPTfrom2DVectorCross((*finalCounts),vec);
            }
        }
        delete vec;
    }
    cout << "run complete!"<< endl;


    /*
    if ((strncmp(this->mMode+4,"bins",4)==0) || (strncmp(this->mMode+2,"bins",4)==0))
    {
        for (int i=0;i<vecLength;i++)
        {
            (*finalDmrts)[vecLength+1][i]=radii[i];
        }
    }
    else
    {
        for (int i=0;i<vecLength;i++)
        {
            (*finalDmrts)[vecLength+1][i]=radii[i+1];
        }
    }
    (*finalDmrts)[vecLength+1][vecLength]=(*finalDmrts)[vecLength+1][vecLength-1]+((*finalDmrts)[vecLength+1][2]-(*finalDmrts)[vecLength+1][1]);
    */

    if(this->mVerb){cout << "Finished calculation!" << endl;}
    if ((*finalDmrts).size()== 0)
    {
        cout << "ERROR algorithm did not produce results!" << endl;
        return;
    }
    else
    {
        size_t c = 0;
        for(size_t i = 0; i< (*finalDmrts).size(); i++)
        {
            if((*finalDmrts)[i][1]!=(*finalDmrts)[i][1]) c++;
        }
        if (c == (*finalDmrts).size())
        {
            cout << "ERROR algorithm did not produce results!" << endl;
            return;
        }
    }
    reader.print2DVectorToXVG(finalDmrts,&outfile);
    //vector<vector<double>> * fpts = eval.getFPTfrom2DVector(vec,10.0);

    outfile.close();
    myfile.close();
}

void dmrtMain::executeFly(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, const vector< vector<double> >* vec, const double start, const double interval, const double end, const int dataColumn)
{
    dmrtalg2 eval = dmrtalg2(this->mMode,this->mVerb,end,start,interval, dataColumn);
    int vecLength = eval.getVecLength();

    // The final vectors get an extra row to save the the radii in, in BINS mode the last column will be empty
    // according to common histogramn convention, size(bins)=size(hist)+1

    //(*finalDmrts) = vector< vector<double> >(vecLength+2,vector<double>(vecLength+1,0.0));
    //(*finalCounts) = vector< vector<int> >(vecLength+1,vector<int>(vecLength+1,0));
    (*finalDmrts) = vector< vector<double> >(vecLength+1,vector<double>(vecLength,0.0));
    (*finalCounts) = vector< vector<int> >(vecLength,vector<int>(vecLength,0));

    vector<double> radii = eval.getRadii();

    for (int i=0;i<vecLength;i++)
    {
        (*finalDmrts)[vecLength][i]=radii[i];
    }

    bool success = true;
    int part = 0;
    while(success==true)
    {
        part++;
        if (!vec || vec->size()!= MAXDOUBLEVEC)
        {
            if(mVerb)
            {
                cout << "Vecsize: " << vec->size()<< endl;
                cout << "Reached end of file in this run" << endl;
            }
            success = false;
        }
        if (mVerb){cout << "Read part "<< part << " of file. Evaluating..." << endl;}
        if (strncmp(this->mMode,"rt",2)==0)
        {
            if (strncmp(this->mMode+2,"bins",4)==0)
            {
                eval.getRTTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
            }
            else if (strncmp(this->mMode+2,"cross",5)==0)
            {
                eval.getRTTfrom2DVectorCross((*finalDmrts),(*finalCounts),vec);
            }
        }
        else if (strncmp(this->mMode,"mftp",4)==0)
        {
            if (strncmp(this->mMode+4,"bins",4)==0)
            {
                eval.getMFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
            }
            else if (strncmp(this->mMode+4,"cross",5)==0)
            {
                eval.getMFPTfrom2DVectorCross((*finalDmrts),(*finalCounts),vec);
            }
        }
        else if (strncmp(this->mMode,"cftp",4)==0)
        {
            if (strncmp(this->mMode+4,"bins",4)==0)
            {
                cout << "Warning! Mode not implemented in current version!" << endl;
                //eval.getFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
            }
            else if (strncmp(this->mMode+4,"cross",5)==0)
            {
                cout << "Initial config: "<< vecLength << endl;
                eval.getFPTfrom2DVectorCross((*finalCounts),vec);
            }
        }
        else if (strncmp(this->mMode,"tftp",4)==0)
        {
            if (strncmp(this->mMode+4,"bins",4)==0)
            {
                cout << "Initial config: "<< vecLength << endl;
                //eval.getFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
                eval.getTFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
                cout << (*finalDmrts)[0][0] << " " <<  (*finalDmrts)[10][0] << " " << (*finalDmrts)[20][0] << endl;
            }
            else if (strncmp(this->mMode+4,"cross",5)==0)
            {
                cout << "Initial config: "<< vecLength << endl;
                cout << "TFTPCROSS not implemented (in fact incorrect!): Calling TFTPBINS"<< endl;
                //eval.getTFPTfrom2DVectorCross((*finalDmrts)[0][0],(*finalDmrts)[0][1],(*finalCounts),vec);
                eval.getTFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
                cout << (*finalDmrts)[0][0] << " " <<  (*finalDmrts)[10][0] << " " << (*finalDmrts)[20][0] << endl;
            }
        }
    }

    if(this->mVerb){cout << "Finished calculation!" << endl;}
    if ((*finalDmrts).size()== 0)
    {
        cout << "ERROR algorithm did not produce results!" << endl;
        return;
    }
    else
    {
        size_t c = 0;
        for(size_t i = 0; i< (*finalDmrts).size(); i++)
        {
            if((*finalDmrts)[i][1]!=(*finalDmrts)[i][1]) c++;
        }
        if (c == (*finalDmrts).size())
        {
            cout << "ERROR algorithm did not produce results!" << endl;
            return;
        }
    }
}



int dmrtMain::execute2(int argc, const char *argv[])
{
    if(argc < 6)
    {
        cout << "Not enough arguments given, please enter as follows:\n\ndmrt inputfile.txt outputfile.txt lowCutoff[nm] deltaR[nm] highCutoff[nm]\n" << endl;
        return 0;
    }
    vector<vector<double> >* tes = new vector< vector<double> >;
    vector<vector<int> >* count = new vector< vector<int> >;
    execute2(tes,count,argv[1],argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]));
    dmrtReader reader = dmrtReader();
    return 0;
}
