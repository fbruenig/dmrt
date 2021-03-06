#ifndef DMRTREADER_H
#define DMRTREADER_H

#include <fstream>
#include <vector>

using namespace std;

#define MAXDOUBLEVEC 25000000

class dmrtReader
{

// functions:

public:

    dmrtReader(){};

    dmrtReader(ifstream *handle, bool verb=1);

    void print2DVectorToXVG(vector<vector<double> > *vec, ofstream *handle);

    void display(int i, int j);
    void display(vector< vector<double> > *vector, int i, int j);
    void display(vector<double> *vector, size_t i);
    void displayLines(vector<vector<double> > *vec, int column, int lines);

    vector<vector<double> >* read2Dvector();
    vector<vector<double> >* read2Dvector(const vector<int> columnsOfInterest);

    vector<vector<double> >* read2DvectorSpace();
    vector<vector<double> >* read2DvectorSpace(const vector<int> columnsOfInterest);

    vector<vector<double> >* read2DvectorSpace(const double rmin, const double rmax);
    vector<vector<double> >* read2DvectorSpace4gb(const double rmin, const double rmax);

// members:

    vector< vector<double> >* mData;

    ifstream *mFilehandle;
    int mFPosition;

private:

    bool mVerb;

}
;

#endif // DMRTREADER_H
