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

    void initLocalVectors(const double start, const double interval, const double end, const int dataColumn = 1);

    void execute2(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, vector<vector<int> > *finalUpts, const char *input, const char *output, const double start, const double interval, const double end, const int dataColumn = 1);
    int execute2(int argc, const char *argv[]);

    void executeFly(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, vector<vector<int> > *finalUpts, const vector<vector<double> > *vec, const double start, const double interval, const double end, const int dataColumn = 1);

    void executeFly(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, vector<vector<int> > *finalUpts, const vector<vector<double> > *vec);

private:
    bool mVerb;
    const char* mMode;

    dmrtalg2 eval;
};

#endif // DMRTMAIN_H
