#ifndef DMRTALG_H
#define DMRTALG_H

#include <vector>

using namespace std;

class dmrtAlg
{

public:
    dmrtAlg();
    dmrtAlg(const char*mode, bool verb);
    dmrtAlg(const char *mode, bool verb, const double dq);

    vector<vector<double> > *getAllDMRTfrom2DVector(const vector<vector<double> > *vec, const double escapeD, const double minD, const double dR);

    vector<double> *getRTTfrom2DVectorReverse(const vector<vector<double> > *vec, const double absorbD, const double escapeD);
    vector<double> *getRTTTfrom2DVector(const vector<vector<double> > *vec, const double absorbD, const double escapeD);

    double getMFPTfrom2DVectorBins(const vector<vector<double> > *vec, const double absorbD, const double escapeD);
    double getMFPTfrom2DVectorCross(const vector<vector<double> > *vec, const double absorbD, const double escapeD);

    int getRadiiVec(vector<vector<double> > *dmrt, const double escapeD, const double minD, const double dR);
    void getAllMFPTLoop(vector<double> &dmrt, vector<int> &counts, const vector<vector<double> > *radii, const vector<vector<double> > *vec, int dir);
    void getAllMFPTfrom2DVector(vector<vector<double> > &dmrt, vector<vector<int> > &counts, const vector<vector<double> > *radii, const vector<vector<double> > *vec);
    void getMFPTfrom2DVectorBins(double &retMftp, int &retCounts, const vector<vector<double> > *vec, double absorbD, double escapeD, int direction);
    void getMFPTfrom2DVectorCross(double &retMftp, int &retCounts, const vector<vector<double> > *vec, double absorbD, double escapeD, int direction);


private:

    double mEpsilon = 0.00000001;
    bool mVerb;
    const char* mMode;
    double mDq;
};

#endif // DMRTALG_H
