#ifndef DMRTMAIN_H
#define DMRTMAIN_H

#include <vector>

using namespace std;

class dmrtMain
{
public:
    dmrtMain(){};
    dmrtMain(const char *mode, bool verb);

    void execute2(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, const char *input, const char *output, const double start, const double interval, const double end);
    int execute2(int argc, const char *argv[]);

    void executeFly(vector<vector<double> > *finalDmrts, vector<vector<int> > *finalCounts, const vector<vector<double> > *vec, const double start, const double interval, const double end);

private:
    bool mVerb;
    const char* mMode;
};

#endif // DMRTMAIN_H
