#include "dmrtreader.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

dmrtReader::dmrtReader(ifstream *handle)
{
    this->mFilehandle=handle;
    this->mVerb=1;
}

dmrtReader::dmrtReader(ifstream *handle, bool verb)
{
    this->mFilehandle=handle;
    this->mVerb = verb;
}

vector<vector<double> > *dmrtReader::read2Dvector()
{
    string line = string("");
    mData = new vector< vector<double>  >;
    int nLines = 0;
    while(getline(*mFilehandle,line))
    {

        stringstream line_ss(line);
        string column = "";
        vector<double> col = vector<double>(2);
        unsigned int index = 0;
        while(getline(line_ss,column,'\t'))
        {
            if(index < col.size())
            {
                col[index]=strtof(column.c_str(),NULL) ;
                index++;
            }
        }
        mData->push_back(col);
        nLines++;
    }
    if(this->mVerb){cout <<  "read " << nLines << " lines" << endl;}
    return this->mData;
}

vector<vector<double> > *dmrtReader::read2DvectorSpace()
{
    string line = string("");
    mData = new vector< vector<double>  >;
    mData->reserve(MAXDOUBLEVEC);
    int nLines = 0;
    while(getline(*mFilehandle,line))
    {

        stringstream line_ss(line);
        string column = "";
        vector<double> col = vector<double>(2);
        unsigned int index = 0;
        while(getline(line_ss,column,'\ '))
        {
            if(!column.empty() && index<2)
            {
                //if (nLines > 135000000) {cout << nLines <<" "<<mData->max_size() <<" "<< column << " " << line_ss.str() << endl;}
                col[index]=strtod(column.c_str(),NULL) ;
                index++;
            }
        }
        nLines++;
        if (col.size() == 2)
        {
            mData->push_back(col);
        }
        else
        {
            if(this->mVerb == 1){cout <<"Could not read line " << nLines << endl;}
        }
    }
    if(this->mVerb){cout <<  "read " << nLines << " lines" << endl;}
    return this->mData;
}


vector<vector<double> > *dmrtReader::read2DvectorSpace(const double rmin, const double rmax)
{
    string line = string("");
    mData = new vector< vector<double>  >;
    mData->reserve(MAXDOUBLEVEC);
    int nLines = 0;
    while(getline(*mFilehandle,line))
    {

        stringstream line_ss(line);
        string column = "";
        vector<double> col = vector<double>(2);
        unsigned int index = 0;
        while(getline(line_ss,column,'\ '))
        {
            if(!column.empty() && index<2)
            {
                //if (nLines > 135000000) {cout << nLines <<" "<<mData->max_size() <<" "<< column << " " << line_ss.str() << endl;}
                col[index]=strtod(column.c_str(),NULL) ;
                index++;
            }
        }
        nLines++;
        if (col.size() == 2)
        {
            if (col[1]< rmax && col[1] > rmin)
            {
                mData->push_back(col);
            }
        }
        else
        {
            if(this->mVerb == 1){cout <<"Could not read line " << nLines << endl;}
        }
    }
    if(this->mVerb){cout <<  "read " << nLines << " lines" << endl;}
    return this->mData;
}

vector<vector<double> > *dmrtReader::read2DvectorSpace4gb(const double rmin, const double rmax)
{
    string line = string("");
    mData = new vector< vector<double>  >;
    int nLines = 0;
    //ifstream is = *mFilehandle;
    //is.seekg(mFPosition);
    while(getline(*mFilehandle,line) && mData->size() < MAXDOUBLEVEC)
    {

        stringstream line_ss(line);
        string column = "";
        vector<double> col = vector<double>(2);
        unsigned int index = 0;
        while(getline(line_ss,column,'\ '))
        {
            if(!column.empty() && index<2)
            {
                //if (nLines > 135000000) {cout << nLines <<" "<<mData->max_size() <<" "<< column << " " << line_ss.str() << endl;}
                col[index]=strtod(column.c_str(),NULL) ;
                index++;
            }
        }
        nLines++;
        if (col.size() == 2)
        {
            if (col[1]< rmax && col[1] > rmin)
            {
                mData->push_back(col);
            }
        }
        else
        {
            if(this->mVerb == 1){cout <<"Could not read line " << nLines << endl;}
        }
    }
    //this->mFPosition = is.tellg();
    if(this->mVerb){cout <<  "read " << nLines << " lines" << endl;}
    return this->mData;
}

void dmrtReader::print2DVectorToXVG(vector<vector<double> > *vec, ofstream *handle)
{
    (*handle) << "@TYPE xy" << endl;
    for(size_t i=0; i< vec->size() ; i++)
    {
        // check is value is NaN
        if ((*vec)[i][0] == (*vec)[i][0])
        {
            (*handle) << (*vec)[i][0] << "\t" << (*vec)[i][1];
            if ((*vec)[i][2])
            {
                (*handle) << "\t" << (*vec)[i][2];
            }
            (*handle) << endl;
        }
    }
}


void dmrtReader::display(int i , int j)
{
    if((*mData)[i][j])
    {
        cout <<  (*mData)[i][j]<< endl;
    }
    else
    {
        cout << "WARNING: could not display data at given coordinates!" << endl;
    }
}

void dmrtReader::display(vector<vector<double> > *vec, int i, int j)
{
    if((*vec)[i][j])
    {
        cout <<  (*vec)[i][j]<< endl;
    }
    else
    {
        cout << "WARNING: could not display data at given coordinates!" << endl;
    }
}

void dmrtReader::displayLines(vector<vector<double> > *vec, int columns, int lines)
{
    if(vec!=NULL)
    {
        for (int i = 0; i<lines;i++)
        {
            for (int j = 0; j<columns;j++)
            {
                cout <<  (*vec)[i][j]<< "\t";
            }
            cout << endl;
        }
    }
    else
    {
        cout << "WARNING: could not display data at given coordinates!" << endl;
    }
}

void dmrtReader::display(vector<double> *vec, size_t i)
{
    if((*vec).size()>i)
    {
        cout <<  (*vec)[i]<< endl;
    }
    else
    {
        cout << "WARNING: could not display data at given coordinates!" << endl;
    }
}

