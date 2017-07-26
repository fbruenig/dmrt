#include "dmrtmain.h"
#include "dmrtalg2.h"
#include "dmrtreader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>

using namespace std;

dmrtMain::dmrtMain(const char *mode, bool verb)
{
    bRt=false;
    bMfpt=false;
    bCftp=false;
    bTftp=false;
    bRate=false;
    bBins=false;
    bCross=false;
    bFull=false;
    bVerb=verb;
    mMode=mode;
    decodeMode();
    cout<< "initialized " << endl;
}

void dmrtMain::decodeMode()
{
    if (strncmp(mMode,"rt",2)==0){bRt=true;}
    else if (strncmp(mMode,"mfpt",4)==0){bMfpt=true;}
    else if (strncmp(mMode,"cftp",4)==0){bCftp=true;}
    else if (strncmp(mMode,"tftp",4)==0){bTftp=true;}
    else if (strncmp(mMode,"rate",4)==0){bRate=true;}

    if (strncmp(mMode+4,"bins",4)==0 || strncmp(mMode+2,"bins",4)==0){this->bBins=true;}
    else if (strncmp(mMode+4,"cross",5)==0 || strncmp(mMode+2,"cross",5)==0){this->bCross=true;}
    else if (strncmp(mMode+4,"full",4)==0){bFull=true;}
}

// Warning this function will not initialize dmrtalg2.mfptDistribution
void dmrtMain::initLocalVectors(const double start, const double interval, const double end, const int dataColumn)
{
    this->eval = dmrtalg2(this->mMode,this->bVerb,end,start,interval, dataColumn);
    eval.initializeLocalVectors();
}

void dmrtMain::initLocalVectors(const double start, const double interval, const double end, const int dataColumn,vector< vector<vector<double> > >* finalDist)
{

  this->eval = dmrtalg2(this->mMode,this->bVerb,end,start,interval, dataColumn);
  (*finalDist) = vector<vector<vector<double> > > (eval.getVecLength(),vector<vector<double>>(eval.getVecLength(),vector<double>(0)));
  eval.mfptDistribution = finalDist;
  eval.initializeLocalVectors();
}

void dmrtMain::initLocalVectors(const vector <double> &radii, const int dataColumn,vector< vector<vector<double> > >* finalDist)
{
  this->eval = dmrtalg2(this->mMode,this->bVerb,radii, dataColumn);
  (*finalDist) = vector<vector<vector<double> > > (eval.getVecLength(),vector<vector<double>>(eval.getVecLength(),vector<double>(0)));
  eval.mfptDistribution = finalDist;
  eval.initializeLocalVectors();
}

void dmrtMain::execute2(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, vector< vector<int> >* finalUpts, const char *input, const char *output, const double start, const double interval, const double end, const int dataColumn)
{

    if(this->bVerb){ cout << "Starting in" << input << endl;}

    ifstream myfile;
    myfile.open(input);
    if(!myfile)
    {
        cout << "ERROR cannot open input file!" << endl;
        return;
    }
    else
    {
        if(this->bVerb){ cout << "Opened input:" << input << endl;}
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
        if(this->bVerb){ cout << "Opened ouput:" << output << endl;}
    }

    dmrtReader reader = dmrtReader(&myfile,this->bVerb);
    dmrtalg2 eval = dmrtalg2(this->mMode,this->bVerb,end,start,interval, dataColumn);
    int vecLength = eval.getVecLength();

    // The final vectors get an extra row to save the radii in, in BINS mode the last column will be empty
    // according to common histogramn convention, size(bins)=size(hist)+1

    (*finalDmrts) = vector< vector<double> >(vecLength+1,vector<double>(vecLength,0.0));
    (*finalCounts) = vector< vector<int> >(vecLength,vector<int>(vecLength,0));
    (*finalUpts) = vector< vector<int> >(vecLength,vector<int>(vecLength,0));

    vector<double> radii = eval.getRadii();

    for (int i=0;i<vecLength;i++)
    {
        (*finalDmrts)[vecLength][i]=radii[i];
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
        if (bVerb){cout << "Read part "<< part << " of file. Evaluating..." << endl;}
        if (strncmp(this->mMode,"rt",2)==0)
        {
            if (strncmp(this->mMode+2,"bins",4)==0)
            {
                eval.getRTTfrom2DVectorBins((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
            else if (strncmp(this->mMode+2,"cross",5)==0)
            {
                eval.getRTTfrom2DVectorCross((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
        }
        else if (strncmp(this->mMode,"mftp",4)==0 || strncmp(this->mMode,"rate",4)==0)
        {
            if (strncmp(this->mMode+4,"bins",4)==0)
            {
                eval.getMFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
            else if (strncmp(this->mMode+4,"cross",5)==0)
            {
                eval.getMFPTfrom2DVectorCross((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
            else if (strncmp(this->mMode+4,"full",5)==0)
            {
                eval.getMFPTfrom2DVectorCross((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
        }
        else if (strncmp(this->mMode,"cftp",4)==0)
        {
            if (strncmp(this->mMode+4,"bins",4)==0)
            {
                cout << "Warning! Mode not implemented in current version!" << endl;
            }
            else if (strncmp(this->mMode+4,"cross",5)==0)
            {
                eval.getFPTfrom2DVectorCross((*finalCounts),vec);
            }
        }
        delete vec;
    }
    cout << "run complete!"<< endl;

    if(this->bVerb){cout << "Finished calculation!" << endl;}
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

void dmrtMain::executeFly(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, vector< vector<int> >* finalUpts,vector< vector<vector<double> > >* finalDist, const vector< vector<double> >* vec, const double start, const double interval, const double end, const int dataColumn)
{
        initLocalVectors(start, interval, end, dataColumn,finalDist);
        executeFly_continue(finalDmrts,finalCounts,finalUpts,vec);
}

void dmrtMain::executeFly(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, vector< vector<int> >* finalUpts,vector< vector<vector<double> > >* finalDist, const vector< vector<double> >* vec, const vector<double> radii, const int dataColumn)
{
        initLocalVectors(radii,dataColumn,finalDist);
        executeFly_continue(finalDmrts,finalCounts,finalUpts,vec);
}

void dmrtMain::executeFly_continue(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, vector< vector<int> >* finalUpts,const vector< vector<double> >* vec)
{

    int vecLength = eval.getVecLength();

    // The final vectors get an extra row to save the the radii in, in BINS mode the last column will be empty
    // according to common histogramn convention, size(bins)=size(hist)+1

    (*finalDmrts) = vector< vector<double> >(vecLength+1,vector<double>(vecLength,0.0));
    (*finalCounts) = vector< vector<int> >(vecLength,vector<int>(vecLength,0));
    (*finalUpts) = vector< vector<int> >(vecLength,vector<int>(vecLength,0));

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
            if(bVerb)
            {
                cout << "Vecsize: " << vec->size()<< endl;
                cout << "Reached end of file in this run" << endl;
            }
            success = false;
        }
        if (bVerb){cout << "Read part "<< part << " of file. Evaluating..." << endl;}
        if (bRt==true)
        {
            if (bBins==true)
            {
                eval.getRTTfrom2DVectorBins((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
            else if (bCross==true)
            {
                eval.getRTTfrom2DVectorCross((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
        }
        else if (bMfpt==true)
        {
            if (bBins==true)
            {
                eval.getMFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
            else if (bCross==true)
            {
                eval.getMFPTfrom2DVectorCross((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
        }
        else if (bCftp==true)
        {
            if (bBins==true)
            {
                cout << "Warning! Mode not implemented in current version!" << endl;
            }
            else if (bCross==true)
            {
                cout << "Initial config: "<< vecLength << endl;
                eval.getFPTfrom2DVectorCross((*finalCounts),vec);
            }
        }
        else if (bRate==true)
        {
            if (bBins==true)
            {
                cout << "Warning! Mode not implemented in current version!" << endl;
            }
            else if (bCross==true)
            {
                cout << "Initial config: "<< vecLength << endl;
                eval.getRatefrom2DVectorCross((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
            else if (bFull==true)
            {
                cout << "Initial config: "<< vecLength << endl;
                eval.getRateFullfrom2DVectorCross((*finalDmrts),(*finalCounts),(*finalUpts),vec);
            }
        }
        else if (bTftp==true)
        {
            if (bBins==true)
            {
                cout << "Initial config: "<< vecLength << endl;
                //eval.getFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
                eval.getTFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
                cout << (*finalDmrts)[0][0] << " " <<  (*finalDmrts)[10][0] << " " << (*finalDmrts)[20][0] << endl;
            }
            else if (bCross==true)
            {
                cout << "Initial config: "<< vecLength << endl;
                cout << "TFTPCROSS not implemented (in fact incorrect!): Calling TFTPBINS"<< endl;
                //eval.getTFPTfrom2DVectorCross((*finalDmrts)[0][0],(*finalDmrts)[0][1],(*finalCounts),vec);
                eval.getTFPTfrom2DVectorBins((*finalDmrts),(*finalCounts),vec);
                cout << (*finalDmrts)[0][0] << " " <<  (*finalDmrts)[10][0] << " " << (*finalDmrts)[20][0] << endl;
            }
        }
    }

    if(this->bVerb){cout << "Finished calculation!" << endl;}
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
    vector<vector<int> >* upts = new vector< vector<int> >;
    execute2(tes,count,upts,argv[1],argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]));
    //dmrtReader reader = dmrtReader();
    return 0;
}
