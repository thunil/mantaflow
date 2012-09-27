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

#include <vector>
#include "vectorbase.h"
#include "kernel.h"

namespace Manta {
    
// access either X or X.pos
template<class C> inline Vec3& posResolve(C& v) { return v.pos; }
template<> inline Vec3& posResolve<Vec3>(Vec3& v) { return v; }

// check for fixed, deleted flags

//! Parallelize point set integration
KERNEL(particle) template<class POS, class POS2, class DATA, class EVAL>
void pointSetKernel(const std::vector<POS>& x, std::vector<POS2>& orig, const DATA& data, std::vector<Vec3>& u, const EVAL& ev) {
    u[i] = ev.eval(posResolve(x[i]), orig[i], data);
}

enum IntegrationMode { IntEuler=0, IntRK2, IntRK4 };
    
//! helper function for integration
template<class POS, class DATA, class EVAL>
void __integratePointSet(std::vector<POS>& pos, std::vector<Vec3>& x, const DATA& data, const EVAL& eval, int mode) {
    const int N = pos.size();
    
    if (mode == IntEuler) {
        std::vector<Vec3> k(N);
        
        pointSetKernel<POS,POS,DATA,EVAL>(pos, pos, data, k, eval);
        for(int i=0; i<N; i++) posResolve(pos[i]) += k[i];
    } 
    else if (mode == IntRK2) {
        std::vector<Vec3> k(N);
        
        pointSetKernel<POS,POS,DATA,EVAL>(pos, pos, data, k, eval);
        for(int i=0; i<N; i++) x[i] = posResolve(pos[i]) + 0.5*k[i];
        
        pointSetKernel<Vec3,POS,DATA,EVAL>(x, pos, data, k, eval);
        for(int i=0; i<N; i++) posResolve(pos[i]) += k[i];
    } 
    else if (mode == IntRK4) {
        std::vector<Vec3> k1(N),k2(N),k3(N),k4(N);
        
        pointSetKernel<POS,POS,DATA,EVAL>(pos, pos, data, k1, eval);
        for(int i=0; i<N; i++) x[i] = posResolve(pos[i]) + 0.5*k1[i];
        
        pointSetKernel<Vec3,POS,DATA,EVAL>(x, pos, data, k2, eval);
        for(int i=0; i<N; i++) x[i] = posResolve(pos[i]) + 0.5*k2[i];
        
        pointSetKernel<Vec3,POS,DATA,EVAL>(x, pos, data, k3, eval);
        for(int i=0; i<N; i++) x[i] = posResolve(pos[i]) + k3[i];
        
        pointSetKernel<Vec3,POS,DATA,EVAL>(x, pos, data, k4, eval);
        for(int i=0; i<N; i++) posResolve(pos[i]) += (_1/6) * (k1[i] + 2*(k2[i] + k3[i]) + k4[i]);
    } 
    else 
        errMsg("unknown integration type");
}

//! Integrate a point set with position-dependant forces
template<class POS, class DATA, class EVAL>
void integratePointsSelf(std::vector<POS>& pos, const EVAL& eval, int mode) {
    std::vector<Vec3> x(pos.size());
    for (size_t i=0; i<x.size(); i++) x[i] = posResolve(pos[i]);
    
    __integratePointSet(pos, x, x, eval, mode);
}

//! Integrate a point set with a passive force field
template<class POS, class DATA, class EVAL>
void integratePointsPassive(std::vector<POS>& pos, const DATA& data, const EVAL& eval, int mode) {
    std::vector<Vec3> x(pos.size());

    __integratePointSet(pos, x, data, eval, mode);    
}


} // namespace

#endif