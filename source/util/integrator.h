/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Helper functions for simple integration
 *
 ******************************************************************************/

#ifndef _INTEGRATE_H
#define _INTEGRATE_H

namespace Manta {

enum IntegrationMode { IntEuler=0, IntRK2, IntRK4 };

template<class EVAL>
void IntegratePointSet(vector<Vec3>& pos, EVAL& evalulator, IntegrationMode mode) {
    const int N = pos.size();
    
    if (mode == IntEuler) 
    {
        vector<Vec3> k(N);
        
        evalulator.eval(pos, k);
        for(int i=0; i<N; i++) pos[i] += k[i];
        
    } 
    else if (mode == IntRK2) 
    {
        vector<Vec3> k(N),x(N);
        
        evalulator.eval(pos, k);
        for(int i=0; i<N; i++) x[i] = pos[i] + 0.5*k[i];
        
        evalulator.eval(x, k);
        for(int i=0; i<N; i++) pos[i] += k[i];
    } 
    else if (mode == IntRK4) 
    {
        vector<Vec3> k1(N),k2(N),k3(N),k4(N),x(N);
        
        evalulator.eval(pos, k1);
        for(int i=0; i<N; i++) x[i] = pos[i] + 0.5*k1[i];
        
        evalulator.eval(x, k2);
        for(int i=0; i<N; i++) x[i] = pos[i] + 0.5*k2[i];
        
        evalulator.eval(x, k3);
        for(int i=0; i<N; i++) x[i] = pos[i] + k3[i];
        
        evalulator.eval(x, k4);
        for(int i=0; i<N; i++) pos[i] += (_1/6) * (k1[i] + 2*(k2[i] + k3[i]) + k4[i]);
    } 
    else 
        errMsg("unknown integration type");
}


} // namespace

#endif