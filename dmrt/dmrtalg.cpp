#include "dmrtalg.h"
#include <vector>
#include <iostream>
#include <string.h>
#include <cmath>

using namespace std;

#define DQ 0.005

dmrtAlg::dmrtAlg()
{
    this->mMode="single";
    this->mVerb=1;
    this->mDq=DQ;
}

dmrtAlg::dmrtAlg(const char * mode, bool verb)
{
    this->mMode=mode;
    this->mVerb=verb;
    this->mDq=DQ;
}

dmrtAlg::dmrtAlg(const char * mode, bool verb,const double dq)
{
    this->mMode=mode;
    this->mVerb=verb;
    this->mDq = dq;
}

double dmrtAlg::getMFPTfrom2DVectorCross(const vector<vector<double> > *vec, const double absorbD, const double escapeD)
{
    bool started = false;
    double start = 0.0;
    bool isBelowAbsorb = false;
    double thisFpt = 0.0;
    size_t counts = 0;
    size_t totCounts = 0;
    double mfpt = 0.0;
    int direction = 1;
    if (absorbD==escapeD)
    {
        return 0.0;
    }
    if (escapeD<absorbD)
    {
        direction = -1;
    }
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            started = false;
            thisFpt = 0.0;
            start = 0.0;
            counts = 0;
        }
        if(started == true)
        {
            if(direction*(*vec)[i][1]>direction*escapeD)
            {
                double relFin = (*vec)[i][0]-start;
                thisFpt += (relFin*counts);
                //cout << thisFpt/counts << endl;
                if (thisFpt > this->mEpsilon)
                {
                    mfpt +=thisFpt;
                    totCounts +=counts;
                }
                thisFpt = 0.0;
                counts = 0;
                start = 0.0;
                isBelowAbsorb = false;
                started = false;
            }
            else
            {
                if(isBelowAbsorb == true)
                {
                    if(direction*(*vec)[i][1]>direction*absorbD)
                    {
                        thisFpt -= (*vec)[i][0];
                        counts++;
                        isBelowAbsorb =false;
                    }
                }
                else
                {
                    if(direction*(*vec)[i][1]<direction*absorbD)
                    {
                        thisFpt -= (*vec)[i][0];
                        counts++;
                        isBelowAbsorb = true;
                    }
                }
            }
        }
        else
        {
            if(direction*(*vec)[i][1]<direction*absorbD)
            {
                isBelowAbsorb = true;
                started = true;
                thisFpt -= (*vec)[i][0];
                counts++;
            }
        }
    }
    return mfpt/totCounts;
}

double dmrtAlg::getMFPTfrom2DVectorBins(const vector<vector<double> > *vec, const double absorbD, const double escapeD)
{
    const double dq = this->mDq;
    double thisFpt = 0.0;
    double mfpt = 0.0;
    double start = 0.0;
    bool started = false;
    int counts = 0;
    int totCounts = 0;
    int direction = 1;
    if (escapeD<absorbD)
    {
        direction = -1;
    }
    if (absorbD==escapeD)
    {
        return 0.0;
    }
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            thisFpt = 0.0;
            counts = 0;
            start = 0.0;
        }
        if(started == true && (direction*(*vec)[i][1])>direction*(escapeD-dq) )
        {
            double relFin = (*vec)[i][0]-start;
            thisFpt += (relFin*counts);
            //cout << thisFpt/counts << endl;
            if (thisFpt > this->mEpsilon)
            {
                mfpt +=thisFpt;
                totCounts +=counts;
            }
            started = false;
            thisFpt = 0.0;
            counts = 0;
            start = 0.0;
        }
        else
        {
            if((*vec)[i][1]<absorbD+dq && (*vec)[i][1]>absorbD-dq)
            {
                if(started == false)
                {
                    started = true;
                    start = (*vec)[i][0];
                }
                thisFpt += start;
                thisFpt -= (*vec)[i][0];
                counts++;
            }
        }
    }
    return mfpt/totCounts;
}

void dmrtAlg::getMFPTfrom2DVectorCross(double& retMftp, int& retCounts, const vector<vector<double> > *vec, double absorbD, double escapeD,int direction)
{
    bool started = false;
    double start = 0.0;
    bool isBelowAbsorb = false;
    double thisFpt = 0.0;
    size_t counts = 0;
    size_t totCounts = 0;
    double mfpt = 0.0;
    if (absorbD==escapeD)
    {
        return;
    }
    if (direction == -1)
    {
        double tmp = absorbD;
        absorbD = escapeD;
        escapeD = tmp;
    }
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            started = false;
            thisFpt = 0.0;
            start = 0.0;
            counts = 0;
        }
        if(started == true)
        {
            if(direction*(*vec)[i][1]>direction*escapeD)
            {
                double relFin = (*vec)[i][0];
                thisFpt += (relFin*counts);
                if (thisFpt > this->mEpsilon)
                {
                    mfpt +=thisFpt;
                    totCounts +=counts;
                }
                thisFpt = 0.0;
                counts = 0;
                start = 0.0;
                isBelowAbsorb = false;
                started = false;
            }
            else
            {
                if(isBelowAbsorb == true)
                {
                    if(direction*(*vec)[i][1]>direction*absorbD)
                    {
                        thisFpt -= (*vec)[i][0];
                        counts++;
                        isBelowAbsorb =false;
                    }
                }
                else
                {
                    if(direction*(*vec)[i][1]<direction*absorbD)
                    {
                        thisFpt -= (*vec)[i][0];
                        counts++;
                        isBelowAbsorb = true;
                    }
                }
            }
        }
        else
        {
            if(direction*(*vec)[i][1]<direction*absorbD)
            {
                isBelowAbsorb = true;
                started = true;
                thisFpt -= (*vec)[i][0];
                counts++;
            }
        }
    }
    retMftp += mfpt;
    retCounts += totCounts;
}

void dmrtAlg::getMFPTfrom2DVectorBins(double& retMftp, int& retCounts, const vector<vector<double> > *vec, double absorbD, double escapeD,int direction)
{
    const double dq = this->mDq;
    double thisFpt = 0.0;
    double mfpt = 0.0;
    double start = 0.0;
    bool started = false;
    int counts = 0;
    int totCounts = 0;
    if (absorbD==escapeD)
    {
        return;
    }
    if (direction == -1)
    {
        double tmp = absorbD;
        absorbD = escapeD;
        escapeD = tmp;
    }
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            thisFpt = 0.0;
            counts = 0;
            start = 0.0;
        }
        if(started == true && (direction*(*vec)[i][1])>direction*(escapeD-dq) )
        {
            double relFin = (*vec)[i][0]-start;
            thisFpt += (relFin*counts);
            if (thisFpt > this->mEpsilon)
            {
                mfpt +=thisFpt;
                totCounts +=counts;
            }
            started = false;
            thisFpt = 0.0;
            counts = 0;
            start = 0.0;
        }
        else
        {
            if((*vec)[i][1]<absorbD+dq && (*vec)[i][1]>absorbD-dq)
            {
                if(started == false)
                {
                    started = true;
                    start = (*vec)[i][0];
                }
                thisFpt += start;
                thisFpt -= (*vec)[i][0];
                counts++;
            }
        }
    }
    retMftp += mfpt;
    retCounts += totCounts;
}

int dmrtAlg::getRadiiVec(vector< vector<double> >* dmrt, const double escapeD, const double minD,const double dR)
{
    size_t vecLength;
    if (escapeD>minD)
    {
        vecLength = (size_t)((escapeD-minD)/dR)+1;
        for(size_t i =0 ; i < vecLength; i++)
        {
            for(size_t j =i+1 ; j < vecLength; j++)
            {
                dmrt->push_back(vector<double>({0,minD+dR*(i),minD+dR*(j)}));
            }
        }
    }
    else
    {
        if(this->mVerb){cout << "WARNING: escape distance larger than minimal distance, reverting!" << endl;}
        vecLength = (size_t)((minD-escapeD)/dR)+1;
        for(size_t i =0 ; i < vecLength; i++)
        {
            for(size_t j =i+1 ; j < vecLength; j++)
            {
                dmrt->push_back(vector<double>({0,escapeD+dR*(j),escapeD+dR*(i)}));
            }
        }
    }
    return vecLength-1;
}

void dmrtAlg::getAllMFPTLoop(vector< double >& dmrt,vector< int >& counts, const vector<vector<double> > * radii, const vector<vector<double> > *vec, int dir)
{
    if (strncmp(this->mMode+4,"bins",4)==0 || strncmp(this->mMode+2,"bins",4)==0)
    {
        #pragma omp parallel for
        for(size_t k =0 ; k < radii->size(); k++)
        {
            getMFPTfrom2DVectorBins(dmrt[k],counts[k],vec,radii->at(k).at(1),radii->at(k).at(2),dir);
            //if(this->mVerb){cout << "Rmin: \t" << radii->at(k).at(1) << " nm" << "\tRmax: \t" << radii->at(k).at(2) << " nm" << dmrt[k] << "\t" << counts[k]  << endl;}
        }
    }
    else if (strncmp(this->mMode+4,"cross",5)==0 || strncmp(this->mMode+2,"cross",5)==0)
    {
        #pragma omp parallel for
        for(size_t k =0 ; k < radii->size(); k++)
        {
            this->getMFPTfrom2DVectorCross(dmrt[k],counts[k],vec,radii->at(k).at(1),radii->at(k).at(2),dir);
            //if(this->mVerb){cout << "Rmin: \t" << radii->at(k).at(1) << " nm" << "\tRmax: \t" << radii->at(k).at(2) << " nm" << dmrt[k] << "\t" << counts[k]  << endl;}
        }
    }
}


void dmrtAlg::getAllMFPTfrom2DVector(vector< vector<double> >& dmrt,vector< vector<int> >& counts, const vector<vector<double> > * radii, const vector<vector<double> > *vec)
{
    if (strncmp(this->mMode,"rt",2)==0)
    {
        vector< double >  forwardTimes(radii->size(),0.0);
        vector< int >    forwardCounts(radii->size(),0.0);
        vector< double >  backwardTimes(radii->size(),0.0);
        vector< int >    backwardCounts(radii->size(),0.0);
        this->getAllMFPTLoop(forwardTimes,forwardCounts,radii,vec,1);
        this->getAllMFPTLoop(backwardTimes,backwardCounts,radii,vec,-1);

        int vecLength = dmrt.size()-1;

        int cumsum = 0;
        for (int i=0;i<vecLength;i++)
        {
            for (int j=i;j<vecLength;j++)
            {
                dmrt[i][j+1]+=forwardTimes[cumsum+j-i];
                counts[i][j+1]+=forwardCounts[cumsum+j-i];
            }
            cumsum +=vecLength-i;
        }
        cumsum = 0;
        for (int i=0;i<vecLength;i++)
        {
            for (int j=i;j<vecLength;j++)
            {
                dmrt[j+1][i]+=backwardTimes[cumsum+j-i];
                counts[j+1][i]+=backwardCounts[cumsum+j-i];
            }
            cumsum +=vecLength-i;
        }
    }
    else if (strncmp(this->mMode,"mftp",4)==0)
    {
        vector< double >  forwardTimes(radii->size(),0.0);
        vector< int >    forwardCounts(radii->size(),0.0);
        this->getAllMFPTLoop(forwardTimes,forwardCounts,radii,vec,1);

        int vecLength = dmrt.size()-1;

        int cumsum = 0;
        for (int i=0;i<vecLength;i++)
        {
            for (int j=i;j<vecLength;j++)
            {
                dmrt[i][j+1]+=forwardTimes[cumsum+j-i];
                dmrt[j+1][i]=0.0;
                counts[i][j+1]+=forwardCounts[cumsum+j-i];
                counts[j+1][i]=0.0;
            }
            cumsum +=vecLength-i;
        }
    }
}


vector< vector<double> > *dmrtAlg::getAllDMRTfrom2DVector(const vector<vector<double> > *vec, const double escapeD, const double minD,const double dR)
{
    vector< vector<double> >* dmrt = new vector<vector<double> >;
    if (strncmp(this->mMode,"full",4)==0)
    {
        if (escapeD>minD)
        {
            size_t vecLength = (size_t)((escapeD-minD)/dR)+1;
            for(size_t i =0 ; i < vecLength; i++)
            {
                for(size_t j =i+1 ; j < vecLength; j++)
                {
                    dmrt->push_back(vector<double>({0,minD+dR*(i),minD+dR*(j)}));
                }
            }
        }
        else
        {
            if(this->mVerb){cout << "WARNING: escape distance larger than minimal distance, reverting!" << endl;}
            size_t vecLength = (size_t)((minD-escapeD)/dR)+1;
            for(size_t i =0 ; i < vecLength; i++)
            {
                for(size_t j =i+1 ; j < vecLength; j++)
                {
                    dmrt->push_back(vector<double>({0,escapeD+dR*(j),escapeD+dR*(i)}));
                }
            }
        }
        if (strncmp(this->mMode+4,"bins",4)==0)
        {
            #pragma omp parallel for
            for(size_t k =0 ; k < dmrt->size(); k++)
            {
                (*dmrt)[k][0]= this->getMFPTfrom2DVectorBins(vec,dmrt->at(k).at(1),dmrt->at(k).at(2));
                if(this->mVerb){cout << "Rmin: \t" << dmrt->at(k).at(1) << " nm" << "\tRmax: \t" << dmrt->at(k).at(2) << " nm" << (*dmrt)[k][0] << endl;}
            }
        }
        else if (strncmp(this->mMode+4,"cross",5)==0)
        {
            #pragma omp parallel for
            for(size_t k =0 ; k < dmrt->size(); k++)
            {
                (*dmrt)[k][0]= this->getMFPTfrom2DVectorCross(vec,dmrt->at(k).at(1),dmrt->at(k).at(2));
                if(this->mVerb){cout << "Rmin: \t" << dmrt->at(k).at(1) << " nm" << "\tRmax: \t" << dmrt->at(k).at(2) << " nm" << (*dmrt)[k][0]  << endl;}
            }
        }
    }
    else if (strncmp(this->mMode,"rt",2)==0)
    {

        int vecLength = (size_t)abs((escapeD-minD)/dR);
        dmrtAlg alg;
        if (strncmp(this->mMode+2,"bins",4)==0)
        {
            //dmrtAlg alg = dmrtAlg("full",this->mVerb,dR);
            alg = dmrtAlg("fullbins",this->mVerb);
        }
        else if (strncmp(this->mMode+2,"cross",5)==0)
        {
            alg = dmrtAlg("fullcross",this->mVerb);
        }
        vector< vector<double> > * forward;
        vector< vector<double> > * backward;
        if (escapeD>minD)
        {
            forward = alg.getAllDMRTfrom2DVector(vec,escapeD,minD,dR);
            backward = alg.getAllDMRTfrom2DVector(vec,minD,escapeD,dR);
        }
        else
        {
            forward = alg.getAllDMRTfrom2DVector(vec,minD,escapeD,dR);
            backward = alg.getAllDMRTfrom2DVector(vec,escapeD,minD,dR);
        }
        //vector< vector<double> > dmrtv = vector<vector<double> >(vecLength+1,vector<double>(vecLength,0.0));
        dmrt = new vector<vector<double> >(vecLength+2,vector<double>(vecLength+1,0.0));
        cout << vecLength << endl;
        int cumsum = 0;
        for (int i=0;i<vecLength;i++)
        {
            for (int j=i;j<vecLength;j++)
            {
                (*dmrt)[i][j+1]=(*forward)[cumsum+j-i][0];
                (*dmrt)[j+1][i]=(*backward)[cumsum+j-i][0];
                cout << (*forward)[cumsum+j-i][0] << "\t" << (*backward)[cumsum+j-i][0] << "\t"<< cumsum+j-i << endl;
            }
            //cumsum +=vecLength-i-1;
            cout << i << "\t"<< cumsum<< endl;
            cumsum +=vecLength-i;
            (*dmrt)[vecLength+1][i+1]=(*forward)[i][2];
        }
        (*dmrt)[vecLength+1][0]=(*dmrt)[vecLength+1][1]-((*dmrt)[vecLength+1][2]-(*dmrt)[vecLength+1][1]);

    }
    else if (strncmp(this->mMode,"mftp",4)==0)
    {

        int vecLength = (size_t)abs((escapeD-minD)/dR);
        dmrtAlg alg;
        if (strncmp(this->mMode+4,"bins",4)==0)
        {
            //dmrtAlg alg = dmrtAlg("full",this->mVerb,dR);
            alg = dmrtAlg("fullbins",this->mVerb);
        }
        else if (strncmp(this->mMode+4,"cross",5)==0)
        {
            alg = dmrtAlg("fullcross",this->mVerb);
        }
        vector< vector<double> > * forward;
        vector< vector<double> > * backward;
        if (escapeD>minD)
        {
            forward = alg.getAllDMRTfrom2DVector(vec,escapeD,minD,dR);
            backward = new vector< vector<double> >(vec->size(),vector<double>(3,0.0));
        }
        else
        {
            forward = alg.getAllDMRTfrom2DVector(vec,minD,escapeD,dR);
            backward = new vector< vector<double> >(vec->size(),vector<double>(3,0.0));
        }
        //vector< vector<double> > dmrtv = vector<vector<double> >(vecLength+1,vector<double>(vecLength,0.0));
        dmrt = new vector<vector<double> >(vecLength+2,vector<double>(vecLength+1,0.0));
        cout << vecLength << endl;
        int cumsum = 0;
        for (int i=0;i<vecLength;i++)
        {
            for (int j=i;j<vecLength;j++)
            {
                (*dmrt)[i][j+1]=(*forward)[cumsum+j-i][0];
                (*dmrt)[j+1][i]=(*backward)[cumsum+j-i][0];
                cout << (*forward)[cumsum+j-i][0] << "\t" << (*backward)[cumsum+j-i][0] << "\t"<< cumsum+j-i << endl;
            }
            //cumsum +=vecLength-i-1;
            cout << i << "\t"<< cumsum<< endl;
            cumsum +=vecLength-i;
            (*dmrt)[vecLength+1][i+1]=(*forward)[i][2];
        }
        (*dmrt)[vecLength+1][0]=(*dmrt)[vecLength+1][1]-((*dmrt)[vecLength+1][2]-(*dmrt)[vecLength+1][1]);

    }
    return dmrt;
}

