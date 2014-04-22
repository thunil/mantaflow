/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Helper functions for higher order interpolation 
 *
 ******************************************************************************/

#ifndef _INTERPOLHIGH_H
#define _INTERPOLHIGH_H

#include "vectorbase.h"

namespace Manta {

// copied from interpol.h 
#define BUILD_INDEX \
    Real px=pos.x-0.5f, py=pos.y-0.5f, pz=pos.z-0.5f; \
    int xi = (int)px; \
    int yi = (int)py; \
    int zi = (int)pz; \
    Real s1 = px-(Real)xi, s0 = 1.-s1; \
    Real t1 = py-(Real)yi, t0 = 1.-t1; \
    Real f1 = pz-(Real)zi, f0 = 1.-f1; \
    /* clamp to border */ \
    if (px < 0.) { xi = 0; s0 = 1.0; s1 = 0.0; } \
    if (py < 0.) { yi = 0; t0 = 1.0; t1 = 0.0; } \
    if (pz < 0.) { zi = 0; f0 = 1.0; f1 = 0.0; } \
    if (xi >= size.x-1) { xi = size.x-2; s0 = 0.0; s1 = 1.0; } \
    if (yi >= size.y-1) { yi = size.y-2; t0 = 0.0; t1 = 1.0; } \
    if (size.z>1) { if (zi >= size.z-1) { zi = size.z-2; f0 = 0.0; f1 = 1.0; } } \
    const int X = 1; \
    const int Y = size.x;    
        
template <class T>
inline T interpolCubic(const T* data, const Vec3i& size, const int Z, const Vec3& pos) {
    BUILD_INDEX
    int idx = xi + Y * yi + Z * zi;    
    DEBUG_ONLY(checkIndexInterpol(size,idx)); DEBUG_ONLY(checkIndexInterpol(size,idx+X+Y+Z));
    
    return  ((data[idx]        *t0 + data[idx+Y]        *t1) * s0
           + (data[idx+X]*t0 + data[idx+X+Y]*t1) * s1) * f0
           +((data[idx+Z]*t0 + data[idx+Y+Z]*t1) * s0
           + (data[idx+X+Z]*t0 + data[idx+X+Y+Z]*t1) * s1) * f1;
}

#undef BUILD_INDEX

} //namespace

#endif


