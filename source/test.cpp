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


PYTHON() void getCurl(MACGrid& vel, Grid<Real>& vort, int comp) {
	Grid<Vec3> velCenter(vel.getParent()), curl(vel.getParent());
	
	GetCentered(velCenter, vel);
	CurlOp(velCenter, curl);
	GetComponent(curl, vort, comp);
}

PYTHON() void setinflow(FlagGrid& flags, MACGrid& vel, LevelsetGrid& phi, Real h) {
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
		}
	}
}
	
PYTHON() void testDiscardNth (BasicParticleSystem& parts,  int skip=1) { 
	//knSetPdataConst<Real>(pd,value); 
	for(int i=0; i<parts.size(); ++i) {
		if(i%(skip+1) == skip) { // keep 
		} else {
			parts.setPos(i, Vec3(-100000) );
		}
	}
}



} //namespace

