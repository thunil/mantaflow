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

#include "grid.h"
#include "commonkernels.h"
#include <cmath>

using namespace std;

namespace Manta {

//! Kernel: get component (not shifted)
KERNEL(idx) returns(Grid<Real> ret(parent))
Grid<Real> GetComponent2(const Grid<Vec3>& grid, int dim) {
    ret[idx] = grid[idx][dim];
};

PYTHON void testp(Grid<Vec3>& b) {
    Grid<Real> d(parent);
    b(20,20,20) = Vec3(21,22,23); 
    {
        cout <<"middle" << endl;        
        Grid<Real> a = GetComponent2(b,0);
        cout << a(20,20,20) << endl;        
        cout <<"middle" << endl;        
    }
    cout << "end" << endl;errMsg("f");
}


KERNEL(idx, reduce=+) returns (double sum=0)
double ddtest(const Grid<Real>& v)
{
    sum += v[idx];
}

KERNEL(idx, reduce=min) returns (double sum=0)
double detest(const Grid<Real>& v)
{
    if (sum < v[idx])
        sum = v[idx];
}



} //namespace