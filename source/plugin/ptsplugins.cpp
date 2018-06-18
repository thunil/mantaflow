// ----------------------------------------------------------------------------
//
// MantaFlow fluid solver framework
// Copyright 2018 Kiwon Um, Nils Thuerey
//
// This program is free software, distributed under the terms of the
// GNU General Public License (GPL)
// http://www.gnu.org/licenses
//
// Particle system helper
//
// ----------------------------------------------------------------------------

#include "particle.h"


namespace Manta {

KERNEL(pts)
void KnAddForcePvel(ParticleDataImpl<Vec3> &v, const Vec3 &da, const ParticleDataImpl<int> *ptype, const int exclude)
{
	if(ptype && ((*ptype)[idx] & exclude)) return;
	v[idx] += da;
}
//! add force to vec3 particle data; a: acceleration
PYTHON() void addForcePvel(ParticleDataImpl<Vec3> &vel, const Vec3 &a, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude)
{
	KnAddForcePvel(vel, a*dt, ptype, exclude);
}

KERNEL(pts)
void KnUpdateVelocityFromDeltaPos(const BasicParticleSystem &p, ParticleDataImpl<Vec3> &v, const ParticleDataImpl<Vec3> &x_prev, const Real over_dt, const ParticleDataImpl<int> *ptype, const int exclude)
{
	if(ptype && ((*ptype)[idx] & exclude)) return;
	v[idx] = (p[idx].pos - x_prev[idx])*over_dt;
}
//! retrieve velocity from position change
PYTHON() void updateVelocityFromDeltaPos(const BasicParticleSystem& parts, ParticleDataImpl<Vec3> &vel, const ParticleDataImpl<Vec3> &x_prev, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude)
{
	KnUpdateVelocityFromDeltaPos(parts, vel, x_prev, 1.0/dt, ptype, exclude);
}

KERNEL(pts)
void KnStepEuler(BasicParticleSystem &p, const ParticleDataImpl<Vec3> &v, const Real dt, const ParticleDataImpl<int> *ptype, const int exclude)
{
	if(ptype && ((*ptype)[idx] & exclude)) return;
	p[idx].pos += v[idx]*dt;
}
//! simple foward Euler integration for particle system
PYTHON() void eulerStep(BasicParticleSystem& parts, const ParticleDataImpl<Vec3> &vel, const ParticleDataImpl<int> *ptype, const int exclude)
{
	KnStepEuler(parts, vel, parts.getParent()->getDt(), ptype, exclude);
}


KERNEL(pts)
void KnSetPartType(ParticleDataImpl<int> &ptype, const BasicParticleSystem &part, const int mark, const int stype, const FlagGrid &flags, const int cflag)
{
	if(flags.isInBounds(part.getPos(idx), 0) && (flags.getAt(part.getPos(idx))&cflag) && (ptype[idx]&stype)) ptype[idx] = mark;
}
//! if particle is stype and in cflag cell, set ptype as mark
PYTHON() void setPartType(const BasicParticleSystem &parts, ParticleDataImpl<int> &ptype, const int mark, const int stype, const FlagGrid &flags, const int cflag)
{
	KnSetPartType(ptype, parts, mark, stype, flags, cflag);
}

} // namespace
