#include "dmrtalg2.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string.h>


/*
In MFPT mode the matrices contain q_i in the first index and q_f in the second index

Recent changes:
* IMPORTANT: radii vec in rtcross is shifted by one entry to the right,
* i.e. when the   result contains -1.0
* with interval 0.1 in first column, matrix values correspond to 0.9, in forward triangle
* this is inversed in the backward triangle, see illustration below:

Current Matrix in RTCROSS mode (for the input A=-1.0 and B=1.0 and interval 0.1)
pairs correspond to (Qi,Qf) coordinates

0   #   #   #   #   (0.9,1.0)
#   0   #   #   #   #
#   #   0   #   #   #
#   #   #   0   #   #
#   #   #   #   0   #
(0.9,1.0) ###   #   0

*/


using namespace std;

dmrtalg2::dmrtalg2()
{
}

dmrtalg2::~dmrtalg2()
{
}


dmrtalg2::dmrtalg2(const char *mode, bool verb, const double escapeD, const double minD, const double dR, const int dataColumn)
{
    mMode=mode;
    mVerb=verb;
    mRadii = vector<double>(0.0);
    mInd = 0;
    mDataColumn = dataColumn;
    getRadiiVec(mRadii,escapeD,minD,dR);
    mVecLength = mRadii.size();
    if (strncmp(this->mMode+6,"dist",4)==0 || strncmp(this->mMode+8,"dist",4)==0 || strncmp(this->mMode+7,"dist",4)==0 ||strncmp(this->mMode+9,"dist",4)==0){this->recordMFPTdistribution=true;}
    else{recordMFPTdistribution=false;}
    // THESE vectors need to be initialized!!!
    //this->mfptDistribution = vector<vector<vector<double> > > (mVecLength,vector<vector<double>>(mVecLength,vector<double>(0)));
    //this->fptDistribution = vector<vector<double> > (2,vector<double>(0));
}

dmrtalg2::dmrtalg2(const char *mode, bool verb, const vector<double>& radii, const int dataColumn)
{
    mMode=mode;
    mVerb=verb;
    mRadii = vector<double>(radii);
    mInd = 0;
    mDataColumn = dataColumn;
    mVecLength = mRadii.size();
    if (strncmp(this->mMode+6,"dist",4)==0 || strncmp(this->mMode+8,"dist",4)==0 || strncmp(this->mMode+7,"dist",4)==0 ||strncmp(this->mMode+9,"dist",4)==0){this->recordMFPTdistribution=true;}
    else{recordMFPTdistribution=false;}
    // THESE vectors need to be initialized!!!
    //this->mfptDistribution = vector<vector<vector<double> > > (mVecLength,vector<vector<double>>(mVecLength,vector<double>(0)));
    //this->fptDistribution = vector<vector<double> > (2,vector<double>(0));

}


void dmrtalg2::initializeLocalVectors()
{
    /*
    locDmrt.clear();
    locStart.clear();
    locCounts.clear();
    fptDistribution.clear();
    locDmrt = vector<vector<double> >(mVecLength,vector<double>(mVecLength,0.0));
    locStart = vector<vector<double> >(mVecLength,vector<double>(mVecLength,0.0));
    locCounts = vector<vector<int> >(mVecLength,vector<int>(mVecLength,0));
    fptDistribution = vector<vector<vector<double> > > (mVecLength,vector<vector<double> >(mVecLength,vector<double>(0)));
    */

    locDmrt.assign(mVecLength, vector<double>(mVecLength,0.0));
    locStart.assign(mVecLength, vector<double>(mVecLength,0.0));
    locCounts.assign(mVecLength, vector<int>(mVecLength,0));
    fptDistribution.assign(mVecLength, vector<vector<double> >(mVecLength,vector<double>(0)));
}

double dmrtalg2::interpolate(const double t1, const double t2, const double r1, const double r2)
{
    const double dt = t2-t1;
    const double dr = r2-r1;
    double tf = (mRadii[mInd]-r1)/dr;
    tf *= dt;
    tf +=t1;
    return tf;
}


void dmrtalg2::findStart2(bool& started, const double d)
{
    if((d>=mRadii[0]) && (d<mRadii[mVecLength-1]))
    {
        started = true;
        if(d>=mRadii[0] &&d<mRadii[1])
        {
            mInd = 1;
        }
        else if(d>=mRadii[mVecLength-2] &&d<mRadii[mVecLength-1])
        {
            mInd = mVecLength-1;
        }
        else
        {
            mInd = 1;
            while(!(d>=mRadii[mInd-1] && d<mRadii[mInd]))
            {
                mInd++;
            }
        }
    }
}


void dmrtalg2::updateDMRTatQf(const int i, vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const double time)
{
    /*if (i==mInd)
    {
        locCounts[i][mInd] =0;
    }*/
    if (locCounts[i][mInd]!=0)
    {
        double relFin   = time-locStart[i][mInd];
        double mfpt = relFin*locCounts[i][mInd];
        dmrt[i][mInd]  += locDmrt[i][mInd] + mfpt;
        counts[i][mInd] += locCounts[i][mInd];
        locCounts[i][mInd]=0;
        upts[i][mInd]++;
        if (this->recordMFPTdistribution==true && mfpt>0.0)
        {
            (*mfptDistribution)[i][mInd].push_back(relFin);
            for(int j=0; j<fptDistribution[i][mInd].size(); j++)
            {
                (*mfptDistribution)[i][mInd].push_back(relFin+fptDistribution[i][mInd][j]);
            }
            (*tptDistribution)[i][mInd].push_back((*mfptDistribution)[i][mInd][(*mfptDistribution)[i][mInd].size()-1]);
            fptDistribution[i][mInd].assign(0,0.0);
            //cout << "MFPT " << mfpt << " " << i << " " << mInd << endl;
            //(*mfptDistribution)[i][mInd].push_back(mfpt);
        }
    }
}

//DEPRECATED
void dmrtalg2::updateDMRTatQfWithDistribution(const int i, vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const double time)
{
    if (locCounts[i][mInd]!=0)
    {
        double relFin   = time-locStart[i][mInd];
        dmrt[i][mInd]  += locDmrt[i][mInd] + (relFin*locCounts[i][mInd]);
        counts[i][mInd] += locCounts[i][mInd];
        locCounts[i][mInd]=0;
        upts[i][mInd]++;
        // shitty implementation here!!
        /*for(int j=0; j<fptDistribution[i].size(); j++)
        {
            upts[i].push_back(fptDistribution[i][j]);
        }*/
    }
}


void dmrtalg2::updateQfatQ(const int i,const double time)
{
    //cout << mInd << " " << i << endl;
    if (locCounts[mInd][i]==0)
    {
        locStart[mInd][i]=time;
        locDmrt[mInd][i]  =0.0;
        locCounts[mInd][i] = 1;
    }
    else
    {
        locDmrt[mInd][i]+=  locStart[mInd][i]-time;
        locCounts[mInd][i]++;
        if (this->recordMFPTdistribution==true)
        {
            fptDistribution[mInd][i].push_back(locStart[mInd][i]-time);
        }
    }
}

//DEPRECATED
void dmrtalg2::updateQfatQWithDistribution(const int i,const double time)
{
    if (locCounts[mInd][i]==0)
    {
        locStart[mInd][i]=time;
        locDmrt[mInd][i]  =0.0;
        locCounts[mInd][i] = 1;
    }
    else
    {
        locDmrt[mInd][i]+=  locStart[mInd][i]-time;
        locCounts[mInd][i]++;
        //fptDistribution[i].push_back(locStart[mInd][i]-time);
    }
}


void dmrtalg2::updateCMatrixTFPT(vector<vector<int> > &counts)
{
    for (int i=0;i< mVecLength;i++)
    {
       counts[i][0]+=locCounts[i][0];
    }

}

void dmrtalg2::updateCMatrixTFPT(vector<vector<double> > &counts)
{
    for (int i=0;i< mVecLength;i++)
    {
       counts[i][0]+=locCounts[i][0];
    }

}

void dmrtalg2::updateVectorsMFPT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const double time)
{
    for (int i=mInd+1;i< mVecLength;i++)
    {
        // update forward Qfs for given Q at mInd
        updateQfatQ(i,time);
    }
    for (int i=0;i< int(mInd);i++)
    {
        // update forward dmrts for given Qf at mInd
        updateDMRTatQf(i,dmrt,counts ,upts,time);
    }
}

void dmrtalg2::updateVectorsRTT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const double time)
{
    for (int i=mInd;i< mVecLength;i++)
    {

        // update forward Qfs for given Q at mInd
        updateQfatQ(i,time);
        // update return dmrts for given Qf at mInd
        updateDMRTatQf(i,dmrt,counts ,upts,time);

    }
    for (int i=0;i< int(mInd-1);i++)
    {
        // update forward dmrts for given Qf at mInd
        updateDMRTatQf(i,dmrt,counts ,upts,time);
        // update return Qfs for given Q at mInd
        updateQfatQ(i,time);
    }
}

void dmrtalg2::getRadiiVec(vector<double> &dmrt, const double escapeD, const double minD, const double dR)
{
    if(strncmp(this->mMode,"rate",4)==0)
    {
        dmrt.push_back(minD);
        dmrt.push_back(escapeD);
        return;
    }

    size_t vecLength;
    if (escapeD>minD)
    {
        vecLength = (size_t)((escapeD-minD)/dR)+1;
        for(size_t i =0 ; i < vecLength; i++)
        {
            //dmrt->push_back(minD+dR*(int(i)-1));
            dmrt.push_back(minD+dR*(int(i)));
        }
    }
}

void dmrtalg2::getMFPTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const vector<vector<double> > *vec)
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
            if((*vec)[i][mDataColumn]>mRadii[mInd])
            {
                while((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength-1)
                {
                    /*for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        updateDMRTatQf(j,dmrt,counts ,upts,(*vec)[i][0]);
                    }*/
                    mInd++;
                }
                if (mInd == mVecLength-1)
                {
                    started = false;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1])
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
            updateVectorsMFPT(dmrt,counts ,upts,(*vec)[i][0]);
        }
        else
        {
            initializeLocalVectors();
            findStart2(started,(*vec)[i][mDataColumn]);
        }
    }
}

void dmrtalg2::getRTTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const vector<vector<double> > *vec)
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
            if((*vec)[i][mDataColumn]>mRadii[mInd])
            {
                while((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
                {
                    /*for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        updateDMRTatQf(j,dmrt,counts ,upts,(*vec)[i][0]);
                    }*/
                    mInd++;
                }
                if (mInd == mVecLength)
                {
                    started = false;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1])
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                    /*for (int j=0;j< int(mInd);j++)
                    {
                        // update retrun dmrts for given Qf at mInd
                        updateDMRTatQf(j,dmrt,counts ,upts,(*vec)[i][0]);
                    }*/
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
            updateVectorsRTT(dmrt,counts ,upts,(*vec)[i][0]);
        }
        else
        {
            findStart2(started,(*vec)[i][mDataColumn]);
        }
    }
}

void dmrtalg2::getTFPTfrom2DVectorBins(vector<vector<double> > &normal, vector<vector<int> > &counts,const vector<vector<double> > *vec)
{
    bool started = false;
    int timer = 0;
    bool lowstart = false;
    bool highstart = false;
    int forcross = 0;
    int backcross = 0;
    int back = 0;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if(((*vec)[i][0]<(*vec)[i-1][0]) && (i<(*vec).size()-1))
        {
            locCounts = vector<vector<int> >(mVecLength,vector<int>(mVecLength,0));
            bool newStart=false;
            findStart2(newStart,(*vec)[i][mDataColumn]);
            while(newStart==true && i<(*vec).size()-3)
            {
                i++;
                newStart=false;
                findStart2(newStart,(*vec)[i][mDataColumn]);
            }
            started = false;
            if(i<(*vec).size()-1)
            {
                i++;
            }
        }
        if(started == true)
        {
            //cout << mInd << endl;
            timer++;
            if((*vec)[i][mDataColumn]>=mRadii[mInd])
            {
                while((*vec)[i][mDataColumn]>=mRadii[mInd] && (int)mInd < mVecLength)
                {
                    mInd++;
                }

                if (mInd == mVecLength)
                {
                    if(lowstart == true)
                    {
                        updateCMatrixTFPT(counts);
                        backcross++;
                    }
                    else
                    {
                        back++;
                        updateCMatrixTFPT(normal);
                    }
                    timer =0;
                    lowstart = false;
                    highstart = false;
                    locCounts = vector<vector<int> >(mVecLength,vector<int>(mVecLength,0));
                    started = false;
                }
                else
                {
                    locCounts[mInd][0]++;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1])
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                }

                if (mInd == 0)
                {
                    if(highstart == true)
                    {
                        updateCMatrixTFPT(counts);
                        forcross++;
                    }
                    else
                    {
                        updateCMatrixTFPT(normal);
                        back++;
                    }
                    timer =0;
                    lowstart = false;
                    highstart = false;
                    locCounts = vector<vector<int> >(mVecLength,vector<int>(mVecLength,0));
                    started = false;
                }
                else
                {
                    locCounts[mInd][0]++;
                }
            }
            else
            {
                locCounts[mInd][0]++;
            }
        }
        else
        {
            findStart2(started,(*vec)[i][mDataColumn]);
            if (started == true)
            {
                if((*vec)[i][1]>(*vec)[i-1][mDataColumn])
                {
                    lowstart = true;
                }
                else
                {
                    highstart = true;
                }
            }
            timer = 0;

        }
    }
    cout << "Crossing vs Back:" << forcross << " "<< backcross << " " << back << endl;
}

void dmrtalg2::getMFPTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const vector<vector<double> > *vec)
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
            if((*vec)[i][mDataColumn]>mRadii[mInd])
            {
                while((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
                {
                    //MFPT:
                    double interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    updateVectorsMFPT(dmrt,counts ,upts,interTime);
                    mInd++;
                }
                if (mInd == mVecLength)
                {
                    started = false;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1])
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                    //MFPT:
                    double interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    updateVectorsMFPT(dmrt,counts ,upts,interTime);
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
        }
        else
        {
            findStart2(started,(*vec)[i][mDataColumn]);
        }
    }
}

void dmrtalg2::getRTTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const vector<vector<double> > *vec)
{
    bool started = false;
    double interTime = 0.0;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            if(mVerb)
            {
                cout << "Restarting at index:" << i << endl;
            }
            initializeLocalVectors();
            started = false;
        }
        if(started == true)
        {
            if((*vec)[i][mDataColumn]>mRadii[mInd])
            {
                while((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
                {
                    // RTT:
                    interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    updateVectorsRTT(dmrt,counts ,upts,interTime);
                    mInd++;
                }
                if (mInd == mVecLength)
                {
                    started = false;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1])
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                    // RTT:
                    interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    updateVectorsRTT(dmrt,counts ,upts,interTime);
                }
                if (mInd == 0)
                {
                    started = false;
                }
            }
        }
        else
        {
            findStart2(started,(*vec)[i][mDataColumn]);
        }
    }
}

void dmrtalg2::getRatefrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const vector<vector<double> > *vec)
{
    bool started = false;
    double interTime =0.0;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            initializeLocalVectors();
            started = false;
        }
        if(started == true)
        {
            if((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
            {
                while((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
                {
                    //Rate:
                    interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    //updateVectorsRTT(dmrt,counts ,upts,interTime);
                    for (int j=mInd+1;j< mVecLength;j++)
                    {

                        // update forward Qfs for given Q at mInd
                        updateQfatQ(j,interTime);

                        // update return dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,intertime);

                    }
                    for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,interTime);

                        // update return Qfs for given Q at mInd
                        updateQfatQ(j,interTime);
                    }
                    mInd++;
                }
                if (mInd == mVecLength)
                {
                    mInd--;
                    updateDMRTatQf(0,dmrt,counts ,upts,interTime);
                    mInd++;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1]&& mInd > 0)
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                    //Rate:
                    interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    for (int j=mInd+1;j< mVecLength;j++)
                    {

                        // update forward Qfs for given Q at mInd
                        updateQfatQ(j,interTime);

                        // update return dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,interTime);

                    }
                    for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,time);

                        // update return Qfs for given Q at mInd
                        updateQfatQ(j,interTime);
                    }
                }
                if (mInd == 0)
                {
                    updateDMRTatQf(1,dmrt,counts ,upts,interTime);

                }
            }
        }
        else
        {
            findStart2(started,(*vec)[i][mDataColumn]);
        }
    }
}

void dmrtalg2::getRateFullfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, vector<vector<int> > &upts, const vector<vector<double> > *vec)
{
    bool started = false;
    double interTime =0.0;
    for(size_t i = 1; i<(*vec).size(); i++)
    {
        if((*vec)[i][0]<(*vec)[i-1][0])
        {
            initializeLocalVectors();
            started = false;
        }
        if(started == true)
        {
            if((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
            {
                while((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
                {
                    //Rate:
                    interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    //updateVectorsRTT(dmrt,counts ,upts,interTime);
                    for (int j=mInd+1;j< mVecLength;j++)
                    {

                        // update forward Qfs for given Q at mInd
                        updateQfatQWithDistribution(j,interTime);

                        // update return dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,intertime);

                    }
                    for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,interTime);

                        // update return Qfs for given Q at mInd
                        updateQfatQWithDistribution(j,interTime);
                    }
                    mInd++;
                }
                if (mInd == mVecLength)
                {
                    mInd--;
                    updateDMRTatQfWithDistribution(0,dmrt,counts ,upts,interTime);
                    mInd++;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1]&& mInd > 0)
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                    //Rate:
                    interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][mDataColumn],(*vec)[i][mDataColumn]);
                    for (int j=mInd+1;j< mVecLength;j++)
                    {

                        // update forward Qfs for given Q at mInd
                        updateQfatQWithDistribution(j,interTime);

                        // update return dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,interTime);

                    }
                    for (int j=0;j< int(mInd);j++)
                    {
                        // update forward dmrts for given Qf at mInd
                        //updateDMRTatQf(j,dmrt,counts ,upts,time);

                        // update return Qfs for given Q at mInd
                        updateQfatQWithDistribution(j,interTime);
                    }
                }
                if (mInd == 0)
                {
                    updateDMRTatQfWithDistribution(1,dmrt,counts ,upts,interTime);

                }
            }
        }
        else
        {
            findStart2(started,(*vec)[i][mDataColumn]);
        }
    }
}


void dmrtalg2::makeHist(vector<vector<double> > &counts, const vector<vector<double> > *vec)
{
    for(size_t i = 0; i<(*vec).size(); i++)
    {
        for(size_t j = 0; j<mRadii.size(); j++)
        {
            if((*vec)[i][mDataColumn]>mRadii[j])
            {
                for(size_t k = j+1; k<mRadii.size(); k++)
                {
                    if((*vec)[i][mDataColumn]<mRadii[k])
                    {
                        counts[k-1][1]=counts[k-1][1]+1;
                        break;
                    }
                }
                break;
            }
        }
    }
}

//// COMMITOR CALCULATION FUNCTION (NOT FINISHED):

void dmrtalg2::updateCountsatQ(const int i)
{
    if (locCounts[mInd][i]==0)
    {
        locCounts[mInd][i]=1;
    }
}

void dmrtalg2::updateCountsatQf(const int i, vector<vector<int> > &counts)
{
    if (locCounts[i][mInd]==1)
    {
        counts[i][mInd]++;
        locCounts[i][mInd]=0;
    }
}

void dmrtalg2::updateCountsatQReturn(const int i)
{
    if (locCounts[i][mInd]==0)
    {
        locCounts[i][mInd]=1;
    }
}

void dmrtalg2::updateCountsatQfReturn(const int i, vector<vector<int> > &counts)
{
    if (locCounts[mInd][i]==1)
    {
        counts[mInd][i]++;
        locCounts[mInd][i]=0;
    }
}

void dmrtalg2::updateVectorsFPT(vector<vector<int> > &counts)
{
    for (int i=mInd+1;i< mVecLength;i++)
    {
        // update forward Counts for given Q at mInd
        updateCountsatQ(i);


        // update return Counts for given Qf at mInd
        updateCountsatQf(i,counts);
    }
    for (int i=0;i< int(mInd);i++)
    {
        // update forward Counts for given Qf at mInd
        updateCountsatQf(i,counts);

        // update return Counts for given Q at mInd
        updateCountsatQ(i);
    }
    if(mInd == mVecLength-2 || mInd == 0)
    {
        locCounts = vector<vector<int> >(mVecLength,vector<int>(mVecLength,0));
    }
}

void dmrtalg2::getFPTfrom2DVectorCross(vector<vector<int> > &counts,const vector<vector<double> > *vec)
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
            if((*vec)[i][mDataColumn]>mRadii[mInd])
            {
                while((*vec)[i][mDataColumn]>mRadii[mInd] && (int)mInd < mVecLength)
                {
                    //FPT:
                    //double interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][1],(*vec)[i][1]);
                    updateVectorsFPT(counts);
                    //updateVectorsMFPT(dmrt,counts,(*vec)[i][0]);
                    mInd++;
                }
                if (mInd == mVecLength)
                {
                    started = false;
                }
            }
            else if((*vec)[i][mDataColumn]<mRadii[mInd-1])
            {
                while((*vec)[i][mDataColumn]<mRadii[mInd-1] && mInd > 0)
                {
                    mInd--;
                    //FPT:
                    //double interTime = interpolate((*vec)[i-1][0],(*vec)[i][0],(*vec)[i-1][1],(*vec)[i][1]);
                    updateVectorsFPT(counts);
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
            findStart2(started,(*vec)[i][mDataColumn]);
        }
    }
}
/////////////////////////////
