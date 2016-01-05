#include "dmrtalg2.h"
#include <vector>
#include <iostream>

using namespace std;

dmrtalg2::dmrtalg2()
{
}

dmrtalg2::dmrtalg2(const char *mode, bool verb, const double escapeD, const double minD, const double dR)
{
    this->mMode=mode;
    this->mVerb=verb;
    this->mRadii = new vector<double>;
    this->mInd = 0;
    getRadiiVec(mRadii,escapeD,minD,dR);
    this->mVecLength = mRadii->size();
    initializeLocalVectors();
}

void dmrtalg2::initializeLocalVectors()
{
    this->locDmrt = vector<vector<double> >(mVecLength,vector<double>(mVecLength,0.0));
    this->locStart = vector<vector<double> >(mVecLength,vector<double>(mVecLength,0.0));
    this->locCounts = vector<vector<int> >(mVecLength,vector<int>(mVecLength,0));
}

double dmrtalg2::interpolate(const double t1, const double t2, const double r1, const double r2)
{
    const double dt = t2-t1;
    const double dr = r2-r1;
    double tf = (*mRadii)[mInd]/dr;
    tf *= dt;
    tf +=t1;
    return tf;
}

void dmrtalg2::findStart(bool& started, const double d)
{
    if(d>(*mRadii)[0] && d<(*mRadii)[mVecLength-1])
    {
        started = true;
        if(d>(*mRadii)[0] &&d<(*mRadii)[1])
        {
            mInd = 0;
        }
        else if(d>(*mRadii)[mVecLength-2] &&d<(*mRadii)[mVecLength-1])
        {
            mInd = mVecLength-2;
        }
        else
        {
            mInd = 0;
            while(!(d>(*mRadii)[mInd] && d<(*mRadii)[mInd+1]))
            {
                mInd++;
            }
        }
    }
}

void dmrtalg2::updateDMRTatQf(const int i, vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time)
{
    if (locCounts[i][mInd]!=0)
    {
        double relFin   = time-locStart[i][mInd];
        dmrt[i][mInd]  += locDmrt[i][mInd] + (relFin*locCounts[i][mInd]);
        counts[i][mInd] += locCounts[i][mInd];
        locCounts[i][mInd]=0;
    }
}

void dmrtalg2::updateQfatQ(const int i,const double time)
{
    if (locCounts[mInd][i]==0)
    {
        locStart[mInd][i]=time;
        locDmrt[mInd][i]  =0.0;
    }
    locDmrt[mInd][i]+=  locStart[mInd][i] -time;
    locCounts[mInd][i]++;
}

void dmrtalg2::updateVectorsMFPT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time)
{
    for (int i=mInd+1;i< mVecLength;i++)
    {
        // update forward Qfs for given Q at mInd
        updateQfatQ(i,time);
    }
    for (int i=0;i< int(mInd);i++)
    {
        // update forward dmrts for given Qf at mInd
        updateDMRTatQf(i,dmrt,counts,time);
    }
}

void dmrtalg2::updateVectorsMFPTCont(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time)
{
    for (int i=mInd+1;i< mVecLength;i++)
    {
        // update forward Qfs for given Q at mInd
        if (locCounts[mInd][i]==0)
        {
            locStart[mInd][i]=time;
            locDmrt[mInd][i]  =0.0;
        }
        locDmrt[mInd][i]+=  locStart[mInd][i] -time;
        locCounts[mInd][i]++;
    }
    for (int i=0;i< int(mInd);i++)
    {
        // update forward dmrts for given Qf at mInd
        if (locCounts[i][mInd]!=0)
        {
            double relFin   = time-locStart[i][mInd];
            double relCounts = counts[i][mInd];
            counts[i][mInd] += locCounts[i][mInd];
            relCounts /= counts[i][mInd];
            dmrt[i][mInd]  =  relCounts* dmrt[i][mInd];
            dmrt[i][mInd]  += (locDmrt[i][mInd] + (relFin*locCounts[i][mInd]))/counts[i][mInd];
            locCounts[i][mInd]=0;
            locStart[i][mInd] =0.0;
            locDmrt[i][mInd]  =0.0;
        }
    }
}

void dmrtalg2::updateVectorsRTT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time)
{
    for (int i=mInd+1;i< mVecLength;i++)
    {
        // update forward Qfs for given Q at mInd
        updateQfatQ(i,time);


        // update return dmrts for given Qf at mInd
        updateDMRTatQf(i,dmrt,counts,time);
    }
    for (int i=0;i< int(mInd);i++)
    {
        // update forward dmrts for given Qf at mInd
        updateDMRTatQf(i,dmrt,counts,time);

        // update return Qfs for given Q at mInd
        updateQfatQ(i,time);
    }
}


void dmrtalg2::getRadiiVec(vector<double> *dmrt, const double escapeD, const double minD, const double dR)
{
    size_t vecLength;
    if (escapeD>minD)
    {
        vecLength = (size_t)((escapeD-minD)/dR)+1;
        for(size_t i =0 ; i < vecLength; i++)
        {
            dmrt->push_back(minD+dR*(int(i)-1));
        }
    }
}

void dmrtalg2::getMFPTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts,const vector<vector<double> > *vec)
{
    bool started = false;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            initializeLocalVectors();
            started = false;
        }
        if(started == true)
        {
            if((*vec)[i][1]>(*mRadii)[mInd+1])
            {
                while((*vec)[i][1]>(*mRadii)[mInd+1] && (int)mInd < mVecLength-1)
                {
                    for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        updateDMRTatQf(j,dmrt,counts,(*vec)[i][0]);
                    }
                    mInd++;
                }
                if (mInd == mVecLength-1)
                {
                    started = false;
                }
            }
            else if((*vec)[i][1]<(*mRadii)[mInd])
            {
                while((*vec)[i][1]<(*mRadii)[mInd] && mInd > 0)
                {
                    mInd--;
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
            updateVectorsMFPT(dmrt,counts,(*vec)[i][0]);
        }
        else
        {
            findStart(started,(*vec)[i][1]);
        }
    }
}

void dmrtalg2::getRTTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts,const vector<vector<double> > *vec)
{
    bool started = false;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            initializeLocalVectors();
            started = false;
        }
        if(started == true)
        {
            if((*vec)[i][1]>(*mRadii)[mInd+1])
            {
                while((*vec)[i][1]>(*mRadii)[mInd+1] && (int)mInd < mVecLength-1)
                {
                    for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        updateDMRTatQf(j,dmrt,counts,(*vec)[i][0]);
                    }
                    mInd++;
                }
                if (mInd == mVecLength-1)
                {
                    started = false;
                }
            }
            else if((*vec)[i][1]<(*mRadii)[mInd])
            {
                while((*vec)[i][1]<(*mRadii)[mInd] && mInd > 0)
                {
                    mInd--;
                    for (int j=0;j< int(mInd);j++)
                    {
                        // update retrun dmrts for given Qf at mInd
                        updateDMRTatQf(j,dmrt,counts,(*vec)[i][0]);
                    }
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
            updateVectorsRTT(dmrt,counts,(*vec)[i][0]);
        }
        else
        {
            findStart(started,(*vec)[i][1]);
        }
    }
}

void dmrtalg2::getMFPTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts,const vector<vector<double> > *vec)
{
    bool started = false;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            initializeLocalVectors();
            started = false;
        }
        if(started == true)
        {
            if((*vec)[i][1]>(*mRadii)[mInd+1])
            {
                while((*vec)[i][1]>(*mRadii)[mInd+1] && (int)mInd < mVecLength-1)
                {
                    //MFPT:
                    double interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][1],(*vec)[i][1]);
                    updateVectorsMFPT(dmrt,counts,interTime);
                    //updateVectorsMFPT(dmrt,counts,(*vec)[i][0]);
                    mInd++;
                }
                if (mInd == mVecLength-1)
                {
                    started = false;
                }
            }
            else if((*vec)[i][1]<(*mRadii)[mInd])
            {
                while((*vec)[i][1]<(*mRadii)[mInd] && mInd > 0)
                {
                    mInd--;
                    //MFPT:
                    double interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][1],(*vec)[i][1]);
                    updateVectorsMFPT(dmrt,counts,interTime);
                    //updateVectorsMFPT(dmrt,counts,(*vec)[i][0]);
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
        }
        else
        {
            findStart(started,(*vec)[i][1]);
        }
    }
}


void dmrtalg2::getRTTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts,const vector<vector<double> > *vec)
{
    bool started = false;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            initializeLocalVectors();
            started = false;
        }
        if(started == true)
        {
            if((*vec)[i][1]>(*mRadii)[mInd+1])
            {
                while((*vec)[i][1]>(*mRadii)[mInd+1] && (int)mInd < mVecLength-2)
                {
                    // RTT:
                    updateVectorsRTT(dmrt,counts,(*vec)[i][0]);
                    mInd++;
                }
                if (mInd == mVecLength-2)
                {
                    started = false;
                }
            }
            else if((*vec)[i][1]<(*mRadii)[mInd])
            {
                while((*vec)[i][1]<(*mRadii)[mInd] && mInd > 0)
                {
                    mInd--;
                    // RTT:
                    updateVectorsRTT(dmrt,counts,(*vec)[i][0]);
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
        }
        else
        {
            findStart(started,(*vec)[i][1]);
        }
    }
}
