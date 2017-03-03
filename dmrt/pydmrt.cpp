#include "Python.h"
#include "numpy/arrayobject.h"
#include "dmrtmain.h"
#include <iostream>
#include <vector>
#include <string.h>

using namespace std;

static PyObject* py_dmrtMain(PyObject* self, PyObject* args)
{
        const char *infile;
        const char *outfile;
        bool verb;
        const char *mode;
        float start;
        float interval;
        float end;
        const char *name = "dmrtlib";

        if (!PyArg_ParseTuple(args, "ssfffsh", &infile, &outfile, &start, &interval,&end,&mode,&verb))
            return NULL;

        dmrtMain prog = dmrtMain(mode,verb);

        if (verb==1){cout << "Entering in " << mode << " mode" << endl;}

        vector<vector<double> >* tes = new vector< vector<double> >;
        vector<vector<int> >* count = new vector< vector<int> >;
        vector<vector<int> >* upts = new vector< vector<int> >;
        prog.execute2(tes,count,upts,infile,outfile,start,interval,end);

        PyObject * TwoDList =NULL;
        PyObject * TwoDListCounts =NULL;

        if (strncmp(mode,"rt",2)==0 || strncmp(mode,"mftp",4)==0 || strncmp(mode,"cftp",4)==0)
        {
            size_t veclength = tes->size()-1;
            TwoDList = PyList_New(veclength+1);
            PyObject * PListIt = NULL;
            PListIt = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList4 = NULL;
                PList4 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList4,j, Py_BuildValue("f", (*tes)[i][j]));
                }
                PyList_SetItem(TwoDList,i,Py_BuildValue("O",PList4));
                PyList_SetItem(PListIt,i,Py_BuildValue("f",(*tes)[veclength][i]));
            }
            PyList_SetItem(TwoDList,veclength,Py_BuildValue("O",PListIt));
            TwoDListCounts = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList4 = NULL;
                PList4 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList4,j, Py_BuildValue("i", (*count)[i][j]));
                }
                PyList_SetItem(TwoDListCounts,i,Py_BuildValue("O",PList4));
            }
        }
        return Py_BuildValue("OO",TwoDList,TwoDListCounts);
}


static PyObject* py_dmrtMainInp(PyObject* self, PyObject* args)
{
        PyObject *inputVec1;
        vector<vector<double> >* vec;
        bool verb;
        const char *mode;
        float start;
        float interval;
        float end;
        const char *name = "dmrtlib";

        if (!PyArg_ParseTuple(args, "Offfsh",&inputVec1, &start, &interval,&end,&mode,&verb))
            return NULL;

        /* Interpret the input objects as numpy arrays. */
        PyObject *input = PyArray_FROM_OTF(inputVec1, NPY_DOUBLE, NPY_IN_ARRAY);

        /* If that didn't work, throw an exception. */
        if (input == NULL) {
            Py_XDECREF(input);
            return NULL;
        }

        /* How many data points are there? */
        int N = (int)PyArray_DIM(input, 0)/2;

        /* Get pointers to the data as C-types. */
        double *in    = (double*)PyArray_DATA(input);

        vector< vector<double> > data = vector< vector<double> >(N,vector<double>(2));
        cout << "Start parsing data!" << endl;
        for (int i =0; i<N; i++)
        {
            //cout << "parsed: " << *(in+i*2) << " " << *(in+i*2+1) << endl;
            data[i].assign(in+i*2,in+i*2+2);
        }

        dmrtMain prog = dmrtMain(mode,verb);
        if (verb==1){cout << "Entering in " << mode << " mode" << endl;}

        vector<vector<double> >* tes = new vector< vector<double> >;
        vector<vector<int> >* count = new vector< vector<int> >;
        vector<vector<int> >* upts = new vector< vector<int> >;
        prog.executeFly(tes,count,upts,&data,start,interval,end);

        /* Clean up. */
        Py_DECREF(input);

        PyObject * TwoDList =NULL;
        PyObject * TwoDListCounts =NULL;
        PyObject * TwoDListUpts =NULL;

        if (strncmp(mode,"rt",2)==0 || strncmp(mode+1,"ftp",3)==0)
        {
            size_t veclength = tes->size()-1;
            TwoDList = PyList_New(veclength+1);
            PyObject * PListIt = NULL;
            PListIt = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList4 = NULL;
                PList4 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList4,j, Py_BuildValue("f", (*tes)[i][j]));
                }
                PyList_SetItem(TwoDList,i,Py_BuildValue("O",PList4));
                PyList_SetItem(PListIt,i,Py_BuildValue("f",(*tes)[veclength][i]));
            }
            PyList_SetItem(TwoDList,veclength,Py_BuildValue("O",PListIt));
            TwoDListCounts = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList4 = NULL;
                PList4 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList4,j, Py_BuildValue("i", (*count)[i][j]));
                    //cout << (*count)[i][j] << endl;
                }
                PyList_SetItem(TwoDListCounts,i,Py_BuildValue("O",PList4));
            }
            TwoDListUpts = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList5 = NULL;
                PList5 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList5,j, Py_BuildValue("i", (*upts)[i][j]));
                    //cout << (*count)[i][j] << endl;
                }
                PyList_SetItem(TwoDListUpts,i,Py_BuildValue("O",PList5));
            }
        }
        return Py_BuildValue("OOO",TwoDList,TwoDListCounts,TwoDListUpts);
}


/*
 * Bind Python function names to C functions
 */
static struct PyMethodDef myModule_methods[] = {
        {"dmrt",(PyCFunction) py_dmrtMain, METH_VARARGS},
        {"dmrtInp",(PyCFunction) py_dmrtMainInp, METH_VARARGS},
        {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

struct module_state {
    PyObject *error;
};

static int module_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int module_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "pydmrt",
        NULL,
        sizeof(struct module_state),
        myModule_methods,
        NULL,
        module_traverse,
        module_clear,
        NULL
};
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_pydmrt (void)
#else
initpydmrt (void)
#endif
{
    #if PY_MAJOR_VERSION >= 3
    import_array();
    return PyModule_Create(&moduledef);
    #else
    (void)Py_InitModule("pydmrt", myModule_methods);
    import_array();
    #endif
}
