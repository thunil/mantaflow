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

enum IntegrationMode { EULER=0, RK2 , RK4, _MAXMODES };

#define DefineIntegrator(Name, VelocityClass, Func) \
    template<IntegrationMode mode> inline Vec3 Name(const Vec3& pos, const VelocityClass& obj, Real dt); \
    template<> inline Vec3 Name<EULER>(const Vec3& pos, const VelocityClass& obj, Real dt) { \
        return dt * obj.Func(pos); \
    } \
    template<> inline Vec3 Name<RK2>(const Vec3& pos, const VelocityClass& obj, Real dt) { \
        Vec3 v1 = obj.Func(pos); \
        Vec3 pos2 = pos + v1*dt; \
        Vec3 v2 = obj.Func(pos2); \
        return (v1+v2) * 0.5 * dt; \
    } \
    template<> inline Vec3 Name<RK4>(const Vec3& pos, const VelocityClass& obj, Real dt) { \
        Vec3 v1 = obj.Func(pos); \
        Vec3 pos2 = pos + v1*(0.5*dt); \
        Vec3 v2 = obj.Func(pos2); \
        Vec3 pos3 = pos + v2*(0.5*dt); \
        Vec3 v3 = obj.Func(pos3); \
        Vec3 pos4 = pos + v3*dt; \
        Vec3 v4 = obj.Func(pos4); \
        return (v1 + (v2+v3)*2.0 + v4) * (dt/6.0); \
    }

} // namespace

#endif