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

KERNEL() 
void testkkk(const Grid<Real>& av, int a=1)
{
	Real sum = av(i,j,k);
}

PYTHON()
double testttt(const Grid<Real>& v)
{
	testkkk ttt(v);
	return v(0,0,0);
	//check error boris???
}



// ... add own test code ...


} //namespace

