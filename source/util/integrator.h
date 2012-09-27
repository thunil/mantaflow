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

namespace Manta {

// access either vec<X> or vec<X>.pos
template<class C> inline Vec3& posElement(C& v, int idx) {
    return v[idx].pos;
}
template<> Vec3& posElement<std::vector<Vec3> >(std::vector<Vec3>& v, int idx) {
    return v[idx];
}


enum IntegrationMode { IntEuler=0, IntRK2, IntRK4 };

template<class POSVEC, class DATA>
class PointSetIntegrator {
public:
    
    void integrate (POSVEC& pos, const DATA& y0, IntegrationMode mode) {
        const int N = pos.size();
        
        if (mode == IntEuler) 
        {
            std::vector<Vec3> k(N);
            
            doEval(pos, pos, y0, k);
            for(int i=0; i<N; i++) pos[i] += k[i];
        } 
        else if (mode == IntRK2) 
        {
            std::vector<Vec3> k(N),x(N);
            
            doEval(pos, pos, y0, k);
            for(int i=0; i<N; i++) x[i] = pos[i] + 0.5*k[i];
            
            doEval(pos, x, y0,k);
            for(int i=0; i<N; i++) pos[i] += k[i];
        } 
        else if (mode == IntRK4) 
        {
            std::vector<Vec3> k1(N),k2(N),k3(N),k4(N),x(N);
            
            doEval(pos, pos, y0, k1);
            for(int i=0; i<N; i++) x[i] = pos[i] + 0.5*k1[i];
            
            doEval(pos, x, y0, k2);
            for(int i=0; i<N; i++) x[i] = pos[i] + 0.5*k2[i];
            
            doEval(pos, y0, x, k3);
            for(int i=0; i<N; i++) x[i] = pos[i] + k3[i];
            
            doEval(pos, x, y0, k4);
            for(int i=0; i<N; i++) pos[i] += (_1/6) * (k1[i] + 2*(k2[i] + k3[i]) + k4[i]);
        } 
        else 
            errMsg("unknown integration type");
    }
    
    virtual void eval(const std::vector<Vec3>& x, const DATA& y0, std::vector<Vec3>& u) = 0;

protected:
    
    inline void doEval(const std::vector<Vec3>& pos, const std::vector<Vec3>& x, const DATA& y0, std::vector<Vec3>& u);
};

// passive evaluation, y0 constant
template<class DATA>
void PointSetIntegrator<DATA>::doEval(const std::vector<Vec3>& pos, const std::vector<Vec3>& x, const DATA& y0, std::vector<Vec3>& u) {
    eval(x,y0,u);
    std::cout << "passive" << std::endl;
}

// possibly active evaluation, y0=x
template<>
void PointSetIntegrator<std::vector<Vec3> >::doEval(const std::vector<Vec3>& pos, const std::vector<Vec3>& x, const std::vector<Vec3>& y0, std::vector<Vec3>& u) {
    if (&pos == & y0)
        eval(x,x,u);
    else
        eval(x,y0,u);
    std::cout << "active" << std::endl;
}

} // namespace

#endif