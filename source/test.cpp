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

using namespace std;

namespace Manta {

//! Kernel: get component (not shifted)
KERNEL(idx) Grid<Real> GetComponent2(const Grid<Vec3>& grid, int dim) {
    _ret[idx] = grid[idx][dim];
};

PYTHON void testp(Grid<Vec3>& b) {
    Grid<Real> d(parent);
    b(20,20,20) = Vec3(21,22,23); 
    {
        cout <<"middle" << endl;        
        Grid<Real> a = GetComponent2(b,0);
        cout << a(20,20,20) << endl;        
        Grid<Real> c = GetComponent2(b,0,d);
        cout << c(20,20,20) << endl;        
        cout <<"middle" << endl;        
    }
    cout << "end" << endl;errMsg("f");
}

KERNEL(idx) REDUCE(+: double sum=0)
void test(const Grid<Real>& v)
{
    sum += v[idx];
}


} //namespace