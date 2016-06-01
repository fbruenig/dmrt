#include <iostream>
#include "../dmrt/dmrtmain.h"
#include <string.h>
using namespace std;

int main(int argc, char** argv)
{

    if(argc>6)
    {
        int verb = atoi(argv[7]);
        if (verb==1){cout << "Entering " << argv[6] << " mode" << endl;}
        dmrtMain prog = dmrtMain(argv[6],verb);
        prog.execute2(argc,argv);
    }
    else
    {
        cout << "Entering" << " single" << " mode" << endl;
        dmrtMain prog = dmrtMain();
        prog.execute2(argc,argv);
    }

    return 0;
}

