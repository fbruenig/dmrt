#ifndef DMRTALG2_H
#define DMRTALG2_H

#include <vector>

using namespace std;

class dmrtalg2
{
public:
    dmrtalg2();
    dmrtalg2(const char *mode, bool verb, const double escapeD, const double minD, const double dR);

    void getRadiiVec(vector<double> *dmrt, const double escapeD, const double minD, const double dR);
    void getAllMFPTLoop(vector<double> &dmrt, vector<int> &counts, const vector<vector<double> > *radii, const vector<vector<double> > *vec, int dir);
    void getAllMFPTfrom2DVector(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *radii, const vector<vector<double> > *vec);
    void getMFPTfrom2DVectorCross(double &retMftp, int &retCounts, const vector<vector<double> > *vec, double absorbD, double escapeD, int direction);

    void getMFPTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);
    void getMFPTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);
    void updateVectorsMFPT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);

    void updateVectorsRTT(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);
    void getRTTfrom2DVectorCross(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);

    void initializeLocalVectors();

    int getVecLength(){return mVecLength;}
    vector<double> getRadii(){return (*mRadii);}


    void updateVectorsMFPTCont(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);
    void getRTTfrom2DVectorBins(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *vec);

    void updateDMRTatQf(const int i, vector<vector<double> > &dmrt, vector<vector<int> > &counts, const double time);
    void updateQfatQ(const int i, const double time);
    void findStart(bool &started, const double d);
    double interpolate(const double t1, const double t2, const double r1, const double r2);
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

};

#endif // DMRTALG2_H
