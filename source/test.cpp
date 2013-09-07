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

PYTHON void checkGrids(Grid<int>& flags1, Grid<int>& flags2, Grid<Real>& phi1, Grid<Real>& phi2, Grid<Vec3>& vel1, Grid<Vec3>& vel2) {
    FOR_IJK(flags1) {
        assertMsg(flags1(i,j,k) == flags2(i,j,k), "flags mismatch");
        assertMsg(norm(vel1(i,j,k)-vel2(i,j,k)) < 1e-1, "vel mismatch");
        assertMsg( fabs(phi1(i,j,k)-phi2(i,j,k)) < 1e-4, "phi mismatch");
    }
}


struct myvec {
    myvec(int n) : x(n) { cout << "constructor" << endl; };
    myvec(const myvec& a) : x(a.x) { cout << "copy constructor" << endl; }
    myvec& operator=(const myvec& a) { x=a.x; cout << "copy operator" << endl; return *this;}
    int& operator[](int idx) { return x[idx]; }
    
    vector<int> x;
};

KERNEL(pts) returns(myvec vec(size)) 
myvec testy(vector<int>& a) {
    vec[i] = a[i];
}

PYTHON void kernelTest() {
    cout << "kernel test" << endl;
    vector<int> a(10);
    for (int i=0;i<10;i++) a[i]=i;
    
    //testy xx(a);
    myvec b = testy(a);
    for (int i=0;i<10;i++) cout << b[i] << endl;
    cout << "kernel end" << endl;
}

PYTHON void getCurl(MACGrid& vel, Grid<Real>& vort, int comp) {
    Grid<Vec3> velCenter(parent), curl(parent);
    
    GetCentered(velCenter, vel);
    CurlOp(velCenter, curl);
    GetComponent(curl, vort, comp);
}

PYTHON void setinflow(FlagGrid& flags, MACGrid& vel, LevelsetGrid& phi, Real h) {
    FOR_IJK(vel) {
        if (i<=2) {
            if (j < h*flags.getSizeY()) {
                vel(i,j,k).x = 1;            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 1;        
                    phi(i,j,k) = -1;
                }                
            } else {
                vel(i,j,k).x = 0;                            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 4;
                    phi(i,j,k) = 1;
                }
            }
        }
        else if (i>=flags.getSizeX()-2) {
            vel(i,j,k).x = 1;            
            /*if (j < 30-12) {
                vel(i,j,k).x = 1;            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 1;        
                    phi(i,j,k) = -1;
                }                
            } else {
                vel(i,j,k).x = 0;                            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 4;
                    phi(i,j,k) = 1;
                }
            }*/
        }
    }
}
    
} //namespace

