#ifndef DMRTMAIN_H
#define DMRTMAIN_H

#include <vector>
#include "dmrtalg2.h"

using namespace std;

class dmrtMain
{
public:
    dmrtMain(){};
    dmrtMain(const char *mode, bool verb);

    void decodeMode();
    void initLocalVectors(const double start, const double interval, const double end, const int dataColumn = 1);
    void initLocalVectors(const double start, const double interval, const double end, const int dataColumn, vector<vector<vector<double> > > *finalDist, vector<vector<vector<double> > > *finalTPDist);
    void initLocalVectors(const vector<double> &radii, const int dataColumn, vector<vector<vector<double> > > *finalDist, vector<vector<vector<double> > > *finalTPDist);

    void execute2(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, vector<vector<int> > *finalUpts, const char *input, const char *output, const double start, const double interval, const double end, const int dataColumn = 1);
    int execute2(int argc, const char *argv[]);

    void executeFly(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, vector<vector<int> > *finalUpts, vector<vector<vector<double> > > *finalDist, vector<vector<vector<double> > > *finalTPDist, const vector<vector<double> > *vec, const double start, const double interval, const double end, const int dataColumn = 1);
    void executeFly(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, vector<vector<int> > *finalUpts, vector<vector<vector<double> > > *finalDist, vector<vector<vector<double> > > *finalTPDist, const vector<vector<double> > *vec, const vector<double> radii, const int dataColumn = 1);

    void executeFly_continue(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, vector<vector<int> > *finalUpts, const vector<vector<double> > *vec);

private:
    bool bVerb;
    const char* mMode;

    bool bRt,bMfpt,bCftp,bTftp,bRate,bBins,bCross,bFull;

    dmrtalg2 eval;
};

#endif // DMRTMAIN_H
