#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <vector>
#include <string.h>

#include "dmrtmain.h"

using namespace std;

template<typename T>
void empty_swap(std::vector<T>& vec) {
   vector<T>().swap(vec);
}

PyObject* parseCArraysToNumpyArrays(const vector<vector<double> >* tes,const vector<vector<int> >* count,const vector<vector<int> >* upts, const vector<vector<double> >* vars,const vector<vector<vector<double> >>* dist,const vector<vector<vector<double> >>* tpDist,const char* mode)
{
        PyObject * TwoDList =NULL;
        PyObject * TwoDListCounts =NULL;
        PyObject * TwoDListUpts =NULL;
        PyObject * TwoDListVars =NULL;
        PyObject * ThreeDListDist =NULL;
        PyObject * ThreeDListTPDist =NULL;

        if (strncmp(mode,"rt",2)==0 || strncmp(mode+1,"ftp",3)==0 || strncmp(mode+1,"fpt",3)==0 || strncmp(mode,"ptp",3)==0)
        {
            int veclength = tes->size()-1;
            npy_intp tmsDims[2] = {veclength+1,veclength};
            npy_intp countDims[2] = {veclength,veclength};
            TwoDList = PyArray_SimpleNewFromData(2,tmsDims, NPY_DOUBLE, (double*)tes->data() );
            TwoDListCounts = PyArray_SimpleNewFromData(2,countDims, NPY_INT, (double*)count->data() );
            TwoDListUpts = PyArray_SimpleNewFromData(2,countDims, NPY_INT, (int*)upts->data() );
            TwoDListVars = PyArray_SimpleNewFromData(2,countDims, NPY_INT, (int*)vars->data() );
            ThreeDListDist = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList6 = NULL;
                PList6 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    size_t distVeclength = (*dist)[i][j].size();
                    PyObject * PList7 = NULL;
                    PList7 = PyList_New(distVeclength);
                    for(size_t k = 0; k < distVeclength ; k++ )
                    {
                        PyList_SetItem(PList7,k, Py_BuildValue("f", (*dist)[i][j][k]));
                        //cout <<"Writing: "<< (*dist)[i][j][k] << endl;
                    }
                    PyList_SetItem(PList6,j,PList7);
		        }
                PyList_SetItem(ThreeDListDist,i,PList6);
            }
            ThreeDListTPDist = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList8 = NULL;
                PList8 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    size_t distVeclength = (*tpDist)[i][j].size();
                    PyObject * PList9 = NULL;
                    PList9 = PyList_New(distVeclength);
                    for(size_t k = 0; k < distVeclength ; k++ )
                    {
                        PyList_SetItem(PList9,k, Py_BuildValue("f", (*tpDist)[i][j][k]));
                        //cout <<"Writing: "<< (*dist)[i][j][k] << endl;
                    }
                    PyList_SetItem(PList8,j,PList9);
		        }
                PyList_SetItem(ThreeDListTPDist,i,PList8);
            }
        PyObject* returnList = Py_BuildValue("OOOOOO",TwoDList,TwoDListCounts,TwoDListUpts,TwoDListVars,ThreeDListDist,ThreeDListTPDist);
        Py_XDECREF(TwoDList);
        Py_XDECREF(TwoDListCounts);
        Py_XDECREF(TwoDListUpts);
        Py_XDECREF(TwoDListVars);
        Py_XDECREF(ThreeDListDist);
        Py_XDECREF(ThreeDListTPDist);
        return returnList;
	}
	else
	{
		return NULL;
	}
}

PyObject* convertCArraysToPythonLists(const vector<vector<double> >* tes,const vector<vector<int> >* count,const vector<vector<int> >* upts, const vector<vector<double> >* vars,const vector<vector<vector<double> >>* dist,const vector<vector<vector<double> >>* tpDist,const char* mode)
{
        PyObject * TwoDList =NULL;
        PyObject * TwoDListCounts =NULL;
        PyObject * TwoDListUpts =NULL;
        PyObject * TwoDListVars =NULL;
        PyObject * ThreeDListDist =NULL;
        PyObject * ThreeDListTPDist =NULL;

        if (strncmp(mode,"rt",2)==0 || strncmp(mode+1,"ftp",3)==0|| strncmp(mode+1,"fpt",3)==0 || strncmp(mode,"ptp",3)==0)
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
                PyList_SetItem(TwoDList,i,PList4);
                PyList_SetItem(PListIt,i,Py_BuildValue("f",(*tes)[veclength][i]));
            }
            PyList_SetItem(TwoDList,veclength,PListIt);


            TwoDListCounts = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList4 = NULL;
                PList4 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList4,j, Py_BuildValue("i", (*count)[i][j]));
                }
                PyList_SetItem(TwoDListCounts,i,PList4);
            }
            TwoDListUpts = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList5 = NULL;
                PList5 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList5,j, Py_BuildValue("i", (*upts)[i][j]));
                }
                PyList_SetItem(TwoDListUpts,i,PList5);
            }
            TwoDListVars = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList5 = NULL;
                PList5 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    PyList_SetItem(PList5,j, Py_BuildValue("f", (*vars)[i][j]));
                }
                PyList_SetItem(TwoDListVars,i,PList5);
            }
            ThreeDListDist = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList6 = NULL;
                PList6 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    size_t distVeclength = (*dist)[i][j].size();
                    PyObject * PList7 = NULL;
                    PList7 = PyList_New(distVeclength);
                    for(size_t k = 0; k < distVeclength ; k++ )
                    {
                        PyList_SetItem(PList7,k, Py_BuildValue("f", (*dist)[i][j][k]));
                        //cout <<"Writing: "<< (*dist)[i][j][k] << endl;
                    }
                    PyList_SetItem(PList6,j,PList7);
		}
                PyList_SetItem(ThreeDListDist,i,PList6);
            }
            ThreeDListTPDist = PyList_New(veclength);
            for(size_t i = 0; i < veclength ; i++ )
            {
                PyObject * PList8 = NULL;
                PList8 = PyList_New(veclength);
                for(size_t j = 0; j < veclength ; j++ )
                {
                    size_t distVeclength = (*tpDist)[i][j].size();
                    PyObject * PList9 = NULL;
                    PList9 = PyList_New(distVeclength);
                    for(size_t k = 0; k < distVeclength ; k++ )
                    {
                        PyList_SetItem(PList9,k, Py_BuildValue("f", (*tpDist)[i][j][k]));
                        //cout <<"Writing: "<< (*dist)[i][j][k] << endl;
                    }
                    PyList_SetItem(PList8,j,PList9);
		}
                PyList_SetItem(ThreeDListTPDist,i,PList8);
            }
        PyObject* returnList = Py_BuildValue("OOOOOO",TwoDList,TwoDListCounts,TwoDListUpts,TwoDListVars,ThreeDListDist,ThreeDListTPDist);
        Py_XDECREF(TwoDList);
        Py_XDECREF(TwoDListCounts);
        Py_XDECREF(TwoDListUpts);
        Py_XDECREF(TwoDListVars);
        Py_XDECREF(ThreeDListDist);
        Py_XDECREF(ThreeDListTPDist);
        return returnList;
	}
	else
	{
        Py_XDECREF(TwoDList);
        Py_XDECREF(TwoDListCounts);
        Py_XDECREF(TwoDListUpts);
        Py_XDECREF(TwoDListVars);
        Py_XDECREF(ThreeDListDist);
        Py_XDECREF(ThreeDListTPDist);
		return NULL;
	}
}

static PyObject* py_dmrtMain(PyObject* self, PyObject* args)
{
        const char *infile;
        const char *outfile;
        int verb;
        const char *mode;
        float start;
        float interval;
        float end;

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

        if (strncmp(mode,"rt",2)==0 || strncmp(mode,"mfpt",4)==0 || strncmp(mode,"cftp",4)==0)
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
        delete tes; delete count; delete upts;
        return Py_BuildValue("OO",TwoDList,TwoDListCounts);
}


static PyObject* py_dmrtMainInp(PyObject* self, PyObject* args)
{
        PyObject *inputVec1;
        int verb;
        const char *mode;
        float start;
        float interval;
        float end;

        if (!PyArg_ParseTuple(args, "Offfsh",&inputVec1, &start, &interval,&end,&mode,&verb)){return NULL;}

        /* Interpret the input objects as numpy arrays. */
        PyArrayObject *input1 = reinterpret_cast<PyArrayObject*>(PyArray_FROM_OTF(inputVec1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));

        /* If that didn't work, throw an exception. */
        if (input1 == NULL) {
            Py_XDECREF(input1);
            return NULL;
        }

        /* How many data points are there? */
        int N = (int)PyArray_DIM(input1, 0);

        /* Get pointers to the data as C-types. */
        double *in    = (double*)PyArray_DATA(input1);

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
        vector<vector<double> >* vars = new vector< vector<double> >;
        vector<vector<vector<double> >>* dist = new vector<vector< vector<double> > >;
        vector<vector<vector<double> >>* tpDist = new vector<vector< vector<double> > >;
        //prog.executeFly(tes,count,upts,dist,&data,start,interval,end);
        prog.executeFly(tes,count,upts,vars,dist,tpDist,&data,start,interval,end);

        /* Clean up. */
        Py_DECREF(input1);

        /* Parse output */
        PyObject* result =  convertCArraysToPythonLists(tes, count, upts, vars, dist, tpDist, mode);
        //PyObject* result =  parseCArraysToNumpyArrays(tes, count, upts, dist, mode);

        empty_swap(*tes); empty_swap(*count); empty_swap(*upts); empty_swap(*vars); empty_swap(*dist); empty_swap(*tpDist);
        delete tes; delete count; delete upts; delete vars, delete dist; delete tpDist;
        return result;
}


static PyObject* py_dmrtMainInpRadii(PyObject* self, PyObject* args)
{
        PyObject *inputVec1=NULL;
        PyObject *inputRadii=NULL;
        PyArrayObject* input1 = NULL;
        PyArrayObject* radii = NULL;

        int verb = false;
        const char *mode = NULL;

        if (!PyArg_ParseTuple(args, "shOO",&mode,&verb,&inputVec1, &inputRadii)){return NULL;}
        if (verb==1){cout << "Entering in " << mode << " mode" << endl;}

        input1 = reinterpret_cast<PyArrayObject*>(PyArray_FROM_OTF(inputVec1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));

        if (input1 == NULL) {
            Py_XDECREF(input1);
            return NULL;
        }

        radii = reinterpret_cast<PyArrayObject*>(PyArray_FROM_OTF(inputRadii, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
        if (radii == NULL) {
            Py_XDECREF(radii);
            return NULL;
        }

        /* How many data points are there? */
        int N = (int)PyArray_DIM(input1, 0);
        int rN = (int)PyArray_DIM(radii, 0);

        /* Get pointers to the data as C-types. */
        double *in    = (double*)PyArray_DATA(input1);

        vector< vector<double> > data = vector< vector<double> >(N,vector<double>(2));
        cout << "Start parsing data!" << endl;
        for (int i =0; i<N; i++)
        {
            //cout << "parsed: " << *(in+i*2) << " " << *(in+i*2+1) << endl;
            data[i].assign(in+i*2,in+i*2+2);
        }


        /* Get pointers to the data as C-types. */
        double *rin    = (double*)PyArray_DATA(radii);

        vector<double> rad = vector<double>(rN);
        cout << "Start parsing radii data!" << endl;
        for (int i =0; i<rN; i++)
        {
            rad[i] = *(rin+i);
        }

        dmrtMain prog = dmrtMain(mode,verb);
        if (verb==1){cout << "Entering in " << mode << " mode" << endl;}

        vector<vector<double> >* tes = new vector< vector<double> >;
        vector<vector<int> >* count = new vector< vector<int> >;
        vector<vector<int> >* upts = new vector< vector<int> >;
        vector<vector<double> >* vars = new vector< vector<double> >;
        vector<vector<vector<double> >>* dist = new vector<vector< vector<double> > >;
        vector<vector<vector<double> >>* tpDist = new vector<vector< vector<double> > >;
        prog.executeFly(tes,count,upts,vars,dist,tpDist,&data,rad);

        /* Clean up. */
        Py_XDECREF(input1);
        Py_XDECREF(radii);

        if (verb==1){cout << "Start parsing output." << endl;}
        /* Parse output */
        PyObject* result =  convertCArraysToPythonLists(tes, count, upts, vars, dist, tpDist, mode);
        //PyObject* result =  parseCArraysToNumpyArrays(tes, count, upts, vars, dist, tpDist, mode);

        //delete in; delete rin;
        empty_swap(*tes); empty_swap(*count); empty_swap(*upts); empty_swap(*vars); empty_swap(*dist); empty_swap(*tpDist);
        delete tes; delete count; delete upts; delete vars, delete dist; delete tpDist;
        return result;
}

/*
static PyObject* py_dmrtDEBUG(PyObject* self, PyObject* args)
{
        PyObject *inputVec1=NULL;
        PyObject *inputRadii=NULL;
        PyObject* input1 = NULL;
        PyObject* radii = NULL;
	bool verb = false;
        const char *mode = NULL;

       if (!PyArg_ParseTuple(args, "OO", &inputVec1, &inputRadii)){return NULL;}

        input1 = PyArray_FROM_OTF(inputVec1, NPY_DOUBLE, NPY_IN_ARRAY);
        if (input1 == NULL) {
            Py_XDECREF(input1);
            return NULL;
        }

  	radii  = PyArray_FROM_OTF(inputRadii, NPY_DOUBLE, NPY_IN_ARRAY);
       	if (radii == NULL) {
         	Py_XDECREF(radii);
            	return NULL;
        }

        int N = (int)PyArray_DIM(input1, 0)/2;

        double *in    = (double*)PyArray_DATA(input1);

        vector< vector<double> > data = vector< vector<double> >(N,vector<double>(2));
        cout << "Start parsing data!" << endl;
        for (int i =0; i<N; i++)
        {
            //cout << "parsed: " << *(in+i*2) << " " << *(in+i*2+1) << endl;
            data[i].assign(in+i*2,in+i*2+2);
        }

        int rN = (int)PyArray_DIM(radii, 0);

        double *rin    = (double*)PyArray_DATA(radii);

        vector<double> rad = vector<double>(rN);
        cout << "Start parsing radii data!" << endl;
        for (int i =0; i<rN; i++)
       	{
           	rad[i] = *(rin+i);
        }


        dmrtMain prog = dmrtMain("rtcross",verb);

        vector<vector<double> >* tes = new vector< vector<double> >;
        vector<vector<int> >* count = new vector< vector<int> >;
        vector<vector<int> >* upts = new vector< vector<int> >;
        vector<vector<vector<double> >>* dist = new vector<vector< vector<double> > >;
        prog.executeFly(tes,count,upts,dist,&data,rad);

	cout << tes->size() << endl;

        Py_XDECREF(input1);
        Py_XDECREF(radii);
	return convertCArraysToPythonLists(tes, count, upts, dist, "rtcross");
}
*/

/*
 * Bind Python function names to C functions
 */
static struct PyMethodDef myModule_methods[] = {
        {"dmrt",(PyCFunction) py_dmrtMain, METH_VARARGS},
        {"dmrtInp",(PyCFunction) py_dmrtMainInp, METH_VARARGS},
        {"dmrtInpRadii",(PyCFunction) py_dmrtMainInpRadii, METH_VARARGS},
//        {"dmrtDEBUG",(PyCFunction) py_dmrtDEBUG, METH_VARARGS},
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
