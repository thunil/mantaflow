/******************************************************************************
*
* MantaFlow fluid solver framework
* Copyright 2017 Steffen Wiewel, Moritz Baecher, Rachel Chu
*
* This program is free software, distributed under the terms of the
* GNU General Public License (GPL) 
* http://www.gnu.org/licenses
*
* Convert mantaflow grids to/from numpy arrays
*
******************************************************************************/

#include "manta.h"
#include "pythonInclude.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

namespace Manta
{
    
#if 1
    void initNumpy()// for python 2.7
    {
        import_array()
    }
#else
    int  initNumpy() // for python 3.4
    {
        import_array()
        return 0;
    }
#endif
    
    template<> PyArrayContainer fromPy<PyArrayContainer>(PyObject* obj) {
        if (PyArray_API == NULL){
            initNumpy();
        }
        
        if (!PyArray_Check(obj)) {
            errMsg("argument is not an numpy array");
        }
        
        PyArrayContainer abuf;
        PyArrayObject* obj_p = reinterpret_cast<PyArrayObject*>(obj);
        
        abuf.TotalSize = PyArray_SIZE(obj_p);
        int source_typ = PyArray_TYPE(obj_p);
        abuf.pData = PyArray_DATA(obj_p);
        
        switch (source_typ) {
            case NPY_FLOAT:
                abuf.DataType = N_FLOAT;
                break;
            case NPY_DOUBLE:
                abuf.DataType = N_DOUBLE;
                break;
            default:
                errMsg("unknown type of Numpy array");
                break;
        }
        
        return abuf;
    }
    
    template<> PyArrayContainer* fromPyPtr<PyArrayContainer>(PyObject* obj, std::vector<void*>* tmp)
    {
        if (!tmp) throw Error("dynamic de-ref not supported for this type");
        void* ptr = malloc(sizeof(PyArrayContainer));
        tmp->push_back(ptr);
        
        *((PyArrayContainer*) ptr) = fromPy<PyArrayContainer>(obj);
        return (PyArrayContainer*) ptr;
    }
}
