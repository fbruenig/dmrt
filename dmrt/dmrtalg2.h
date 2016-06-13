#ifndef DMRTALG2_H
#define DMRTALG2_H

#include <vector>

using namespace std;

class dmrtalg2
{
public:
    dmrtalg2();
    dmrtalg2(const char *mode, bool verb, const double escapeD, const double minD, const double dR, const int dataColumn = 1);

    //Current implementation is purely based on binning as setup in the member function findstart2
    //Accordingly for radii r[i],r[i+1],r[i+2]... updates are done a follows:
    //
    //CROSS: when crossing r[i] from low to high updates are performed at i
    //BINS:  when crossing r[i] from low to high updates are performed at i+1



    //Helper functions:

    void getRadiiVec(vector<double> *dmrt, const double escapeD, const double minD, const double dR);
    void initializeLocalVectors();
    int getVecLength(){return mVecLength;}
    vector<double> getRadii(){return (*mRadii);}
    //void findStart(bool &started, const double d);
    void findStart2(bool &started, const double d);
    double interpolate(const double t1, const double t2, const double r1, const double r2);
    void makeHist(vector<vector<double> > &counts, const vector<vector<double> > *vec);


    // MFTP/RTT extraction (option == mftp oder rt)

    void getAllMFPTLoop(vector<double> &dmrt, vector<int> &counts, const vector<vector<double> > *radii, const vector<vector<double> > *vec, int dir);
    void getAllMFPTfrom2DVector(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *radii, const vector<vector<double> > *vec);

    void getMFPTfrom2DVectorCross(double &retMftp, int &retCounts, const vector<vector<double> > *vec, double absorbD, double escapeD, int direction);

    void getMFPTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);
    void getMFPTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);
    void getRTTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);
    void getRTTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);
    void updateVectorsMFPT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);
    void updateVectorsRTT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);

    void updateVectorsMFPTCont(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);
    void updateDMRTatQf(const int i, vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);
    void updateQfatQ(const int i, const double time);


    // p(tp|r) extraction: (option== tftp)

    void updateCMatrixTFPT(vector<vector<double> > &counts);
    void updateCMatrixTFPT(vector<vector<int> > &counts);
    void getTFPTfrom2DVectorBins(vector<vector<double> > &normal, vector<vector<int> > &counts, const vector<vector<double> > *vec);



    // Commitor extraction: (option== cftp)

    void getFPTfrom2DVectorCross(vector<vector<int> > &counts, const vector<vector<double> > *vec);
    void updateVectorsFPT(vector<vector<int> > &counts);
    void updateCountsatQ(const int i);
    void updateCountsatQf(const int i, vector<vector<int> > &counts);
    void updateCountsatQReturn(const int i);
    void updateCountsatQfReturn(const int i, vector<vector<int> > &counts);

    void updateDMRTatQfReturn(const int i, vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);
    void updateQfatQReturn(const int i, const double time);
private:

    double mEpsilon = 0.00000001;
    bool mVerb;
    const char* mMode;
    double mDq;
    int mVecLength ;
    size_t mInd;
    vector<double> *mRadii;
    vector<vector<double> > locDmrt;
    vector<vector<double> > locStart;
    vector<vector<int> > locCounts;

    int mDataColumn;

};

#endif // DMRTALG2_H
