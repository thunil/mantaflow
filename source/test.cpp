/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Use this file to test new functionality
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include <cmath>

using namespace std;

namespace Manta {

// two simple example kernels

KERNEL(idx, reduce=+) returns (double sum=0)
double reductionTest(const Grid<Real>& v)
{
	sum += v[idx];
}

KERNEL(idx, reduce=min) returns (double sum=0)
double minReduction(const Grid<Real>& v)
{
	if (sum < v[idx])
		sum = v[idx];
}

// ... add own test code here if necessary ...


// test function and kernel with python array

KERNEL(bnd=0)
void knNumpyTest(Grid<Real>& data, PyArrayContainer pyar) 
{
	data(i,j,k) += (reinterpret_cast<float*>(pyar.pData))[0]; // calc access into numpy array!
}

PYTHON() void numpyTest(Grid<Real>& data, PyArrayContainer pyar) {
	knNumpyTest(data,pyar);
}

add scene file and test...


} //namespace

