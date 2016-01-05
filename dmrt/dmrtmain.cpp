#include "dmrtmain.h"
#include "dmrtalg.h"
#include "dmrtalg2.h"
#include "dmrtreader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>

using namespace std;


dmrtMain::dmrtMain()
{
    this->mVerb=1;
    this->mMode="single";
}

dmrtMain::dmrtMain(const char *mode, bool verb)
{
    this->mVerb=verb;
    this->mMode=mode;
}

vector<vector<double> > * dmrtMain::execute(const char *input, const char *output, const double start, const double interval, const double end)
{

    if(this->mVerb){ cout << "Starting in" << input << endl;}

    ifstream myfile;
    myfile.open(input);
    if(!myfile)
    {
        cout << "ERROR cannot open input file!" << endl;
        return 0;
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
        return 0;
    }
    else
    {
        if(this->mVerb){ cout << "Opened ouput:" << output << endl;}
    }

    dmrtReader reader = dmrtReader(&myfile,this->mVerb);
    vector< vector<double> >* vec = reader.read2DvectorSpace(start-interval,end+interval);

    if (!vec || vec->size()== 0)
    {
        cout << "ERROR no data was read!" << endl;
        return 0;
    }


    dmrtAlg eval = dmrtAlg(this->mMode,this->mVerb);
    vector<vector<double> > * dmrt = eval.getAllDMRTfrom2DVector(vec,end,start,interval);
    //reader.displayLines(vec,1,5);
    //reader.displayLines(dmrt,5,5);
    if(this->mVerb){cout << "Finished calculation!" << endl;}
    if (!dmrt || dmrt->size()== 0)
    {
        cout << "ERROR algorithm did not produce results!" << endl;
        return 0;
    }
    else
    {
        size_t c = 0;
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
    //vector<vector<double>> * fpts = eval.getFPTfrom2DVector(vec,10.0);
    delete vec;
    outfile.close();
    myfile.close();
    return dmrt;
}

vector<vector<double> > dmrtMain::executeLongFile(const char *input, const char *output, const double start, const double interval, const double end)
{

    if(this->mVerb){ cout << "Starting in" << input << endl;}

    ifstream myfile;
    myfile.open(input);
    if(!myfile)
    {
        cout << "ERROR cannot open input file!" << endl;
        return vector<vector<double> >(0);
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
        return vector<vector<double> >(0);
    }
    else
    {
        if(this->mVerb){ cout << "Opened ouput:" << output << endl;}
    }

    dmrtReader reader = dmrtReader(&myfile,this->mVerb);
    dmrtAlg eval = dmrtAlg(this->mMode,this->mVerb);
    vector< vector<double> >* radii = new vector< vector<double> >;
    int vecLength = eval.getRadiiVec(radii,end,start,interval);
    vector< vector<double> > finalDmrts(vecLength+2,vector<double>(vecLength+1,0.0));
    vector< vector<int> > finalCounts(vecLength+1,vector<int>(vecLength+1,0));

    bool success = true;
    int part = 0;
    while(success==true)
    {
        vector< vector<double> >* vec= new vector<vector<double> > ;
        vec = reader.read2DvectorSpace4gb(start-interval,end+interval);
        part++;
        if (!vec || vec->size()!= MAXDOUBLEVEC)
        {
            cout << "Reached end of file in this run" << endl;
            success = false;
        }


        vector< vector<double> > dmrts(vecLength+1,vector<double>(vecLength+1));
        vector< vector<int> > counts(vecLength+1,vector<int>(vecLength+1));

        if (mVerb){cout << "Read part "<< part << " of file. Evaluating..." << endl;}
        eval.getAllMFPTfrom2DVector(dmrts,counts,radii,vec);
        for (int i=0;i<vecLength;i++)
        {
            for (int j=i;j<vecLength;j++)
            {
                finalDmrts[i][j+1] += dmrts[i][j+1];
                finalDmrts[j+1][i] += dmrts[j+1][i];
                finalCounts[i][j+1] += counts[i][j+1];
                finalCounts[j+1][i] += counts[j+1][i];
            }
        }
        delete vec;
    }
    for (int i=0;i<vecLength;i++)
    {
        for (int j=i;j<vecLength;j++)
        {
            finalDmrts[i][j+1] = finalDmrts[i][j+1]/finalCounts[i][j+1];
            if (finalCounts[j+1][i] != 0)
            {
                finalDmrts[j+1][i] = finalDmrts[j+1][i]/finalCounts[j+1][i];
            }
        }
        finalDmrts[vecLength+1][i+1]=(*radii)[i][2];
    }
    finalDmrts[vecLength+1][0]=finalDmrts[vecLength+1][1]-(finalDmrts[vecLength+1][2]-finalDmrts[vecLength+1][1]);
    finalDmrts[vecLength+1][vecLength]=finalDmrts[vecLength+1][vecLength-1]+(finalDmrts[vecLength+1][2]-finalDmrts[vecLength+1][1]);

    if(this->mVerb){cout << "Finished calculation!" << endl;}
    if (finalDmrts.size()== 0)
    {
        cout << "ERROR algorithm did not produce results!" << endl;
        return vector<vector<double> >(0);
    }
    else
    {
        size_t c = 0;
        for(size_t i = 0; i< finalDmrts.size(); i++)
        {
            if(finalDmrts[i][1]!=finalDmrts[i][1]) c++;
        }
        if (c == finalDmrts.size())
        {
            cout << "ERROR algorithm did not produce results!" << endl;
            return vector<vector<double> >(0);
        }
    }
    reader.print2DVectorToXVG(&finalDmrts,&outfile);
    //vector<vector<double>> * fpts = eval.getFPTfrom2DVector(vec,10.0);

    outfile.close();
    myfile.close();
    return finalDmrts;
}

void dmrtMain::executeLongFile(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, const char *input, const char *output, const double start, const double interval, const double end)
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
    dmrtAlg eval = dmrtAlg(this->mMode,this->mVerb);
    vector< vector<double> >* radii = new vector< vector<double> >;
    int vecLength = eval.getRadiiVec(radii,end,start,interval);
    (*finalDmrts) = vector< vector<double> >(vecLength+2,vector<double>(vecLength+1,0.0));
    (*finalCounts) = vector< vector<int> >(vecLength+1,vector<int>(vecLength+1,0));

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


        vector< vector<double> > dmrts(vecLength+1,vector<double>(vecLength+1,0.0));
        vector< vector<int> > counts(vecLength+1,vector<int>(vecLength+1,0));

        if (mVerb){cout << "Read part "<< part << " of file. Evaluating..." << endl;}
        eval.getAllMFPTfrom2DVector(dmrts,counts,radii,vec);
        for (int i=0;i<vecLength;i++)
        {
            for (int j=i;j<vecLength;j++)
            {
                (*finalDmrts)[i][j+1] += dmrts[i][j+1];
                (*finalDmrts)[j+1][i] += dmrts[j+1][i];
                (*finalCounts)[i][j+1] += counts[i][j+1];
                (*finalCounts)[j+1][i] += counts[j+1][i];
            }
        }
        delete vec;
    }
    for (int i=0;i<vecLength;i++)
    {   /*
        for (int j=i;j<vecLength;j++)
        {
            (*finalDmrts)[i][j+1] = (*finalDmrts)[i][j+1]/(*finalCounts)[i][j+1];
            if ((*finalCounts)[j+1][i] != 0)
            {
                (*finalDmrts)[j+1][i] = (*finalDmrts)[j+1][i]/(*finalCounts)[j+1][i];
            }
        }
        */
        (*finalDmrts)[vecLength+1][i+1]=(*radii)[i][2];
    }
    (*finalDmrts)[vecLength+1][0]=(*finalDmrts)[vecLength+1][1]-((*finalDmrts)[vecLength+1][2]-(*finalDmrts)[vecLength+1][1]);
    (*finalDmrts)[vecLength+1][vecLength]=(*finalDmrts)[vecLength+1][vecLength-1]+((*finalDmrts)[vecLength+1][2]-(*finalDmrts)[vecLength+1][1]);

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



void dmrtMain::execute2(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, const char *input, const char *output, const double start, const double interval, const double end)
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
    dmrtalg2 eval = dmrtalg2(this->mMode,this->mVerb,end,start,interval);
    int vecLength = eval.getVecLength()-1;
    (*finalDmrts) = vector< vector<double> >(vecLength+2,vector<double>(vecLength+1,0.0));
    (*finalCounts) = vector< vector<int> >(vecLength+1,vector<int>(vecLength+1,0));

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
        delete vec;
    }
    cout << "run complete!"<< endl;
    vector<double> radii = eval.getRadii();


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

    //(*finalDmrts)[vecLength+1][0]=(*finalDmrts)[vecLength+1][1]-((*finalDmrts)[vecLength+1][2]-(*finalDmrts)[vecLength+1][1]);
    //(*finalDmrts)[vecLength+1][vecLength-2]=(*finalDmrts)[vecLength+1][vecLength-3]+((*finalDmrts)[vecLength+1][2]-(*finalDmrts)[vecLength+1][1]);
    //(*finalDmrts)[vecLength+1][vecLength-1]=(*finalDmrts)[vecLength+1][vecLength-2]+((*finalDmrts)[vecLength+1][2]-(*finalDmrts)[vecLength+1][1]);
    (*finalDmrts)[vecLength+1][vecLength]=(*finalDmrts)[vecLength+1][vecLength-1]+((*finalDmrts)[vecLength+1][2]-(*finalDmrts)[vecLength+1][1]);


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

void dmrtMain::executeFly(vector< vector<double> >* finalDmrts, vector< vector<int> >* finalCounts, const vector< vector<double> >* vec, const double start, const double interval, const double end)
{
    dmrtalg2 eval = dmrtalg2(this->mMode,this->mVerb,end,start,interval);
    int vecLength = eval.getVecLength()-1;
    (*finalDmrts) = vector< vector<double> >(vecLength+2,vector<double>(vecLength+1,0.0));
    (*finalCounts) = vector< vector<int> >(vecLength+1,vector<int>(vecLength+1,0));

    bool success = true;
    int part = 0;
    while(success==true)
    {
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
    }
    cout << "run complete!"<< endl;
    vector<double> radii = eval.getRadii();


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


int dmrtMain::execute(int argc, const char *input, const char *output, const double start, const double interval, const double end)
{
    if(argc < 6)
    {
        cout << "Not enough arguments given, please enter as follows:\n\ndmrt inputfile.txt outputfile.txt lowCutoff[nm] deltaR[nm] highCutoff[nm]\n" << endl;
        return 0;
    }

    cout << "Starting..." << endl;

    ifstream myfile;
    myfile.open(input);
    if(!myfile)
    {
        cout << "ERROR cannot open input file!" << endl;
        return 0;
    }
    else
    {
        cout << "Opened input:" << input << endl;
    }
    ofstream outfile;
    outfile.open(output);
    if(!outfile)
    {
        cout << "ERROR cannot open output file!" << endl;
        return 0;
    }
    else
    {
        cout << "Opened ouput:" << output << endl;
    }

    dmrtReader reader = dmrtReader(&myfile);
    vector< vector<double> >* vec = reader.read2DvectorSpace();

    if (!vec || vec->size()== 0)
    {
        cout << "ERROR no data was read!" << endl;
        return 0;
    }


    dmrtAlg eval = dmrtAlg();
    vector<vector<double> > * dmrt = eval.getAllDMRTfrom2DVector(vec,end,start,interval);

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

    //vector<vector<double>> * fpts = eval.getFPTfrom2DVector(vec,10.0);

    outfile.close();
    myfile.close();
    return 0;
}

int dmrtMain::execute(int argc, const char *argv[])
{
    if(argc < 6)
    {
        cout << "Not enough arguments given, please enter as follows:\n\ndmrt inputfile.txt outputfile.txt lowCutoff[nm] deltaR[nm] highCutoff[nm]\n" << endl;
        return 0;
    }
    vector<vector<double> > * result = this->execute(argv[1],argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]));
    dmrtReader reader = dmrtReader();
    reader.displayLines(result,10,10);
    return 0;
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
