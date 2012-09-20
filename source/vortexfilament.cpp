/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex filament
 *
 ******************************************************************************/

#include "vortexfilament.h"
#include "integrator.h"
#include "mesh.h"

using namespace std;
namespace Manta {

// vortex filament effect
struct FilamentKernel {
    FilamentKernel(const Vec3& g0, const Vec3& g1, Real circ, Real scale, Real cutoff, Real a) : 
        gamma0(g0), gamma1(g1), cutoff2(cutoff*cutoff), a2(a*a) {
            gammaMid = 0.5*(g0+g1);
            gammaDir = g1-g0;
            strength = scale*circ;
            regFact = a*a*normSquare(gammaDir);
        }

    Vec3 gamma0, gamma1, gammaMid, gammaDir;
    Real strength;
    Real cutoff2;
    Real a2, regFact;
    
    inline bool isValid(const Vec3& p) const {
        const Real rlen2 = normSquare(p - gammaMid);
        const Real r0_2 = normSquare(p - gamma0);
        const Real r1_2 = normSquare(p - gamma1);
        return rlen2 < cutoff2 && r0_2 > 1e-8 && r1_2 > 1e-8;
    }
    
    inline Vec3 evaluate(const Vec3& p) const {
        // vortex line integral
        const Vec3 r0 = gamma0-p, r1 = gamma1-p;
        const Real r0n = 1.0f/sqrt(a2+normSquare(r0));
        const Real r1n = 1.0f/sqrt(a2+normSquare(r1));
        const Vec3 cp = cross(r0,r1);
        const Real upper = dot(r1,gammaDir)*r1n - dot(r0,gammaDir)*r0n;
        const Real lower = regFact + normSquare(cp);
        return upper/lower*cp;
    }
};


void VortexFilamentSystem::addRing(const Vec3& position, Real circulation, Real radius, Vec3 normal, int number) {
	normalize(normal);
    Vec3 worldup (0,1,0);
	if (norm(normal - worldup) < 1e-5) worldup = Vec3(1,0,0);
	
	Vec3 u = cross(normal, worldup); normalize(u);
	Vec3 v = cross(normal, u); normalize(v);
	
	int firstNum = size();
	for (int i=0; i<number; i++) {
		Real phi = (Real)i/(Real)number * M_PI * 2.0;
		Vec3 p = position + radius * (u*cos(phi) + v*sin(phi));
		
		int num = add(BasicParticleData(p));
        mSegments.push_back(VortexFilamentData(num, num+1, circulation));
	}
	mSegments.rbegin()->idx1 = firstNum;
}

void VortexFilamentSystem::addLine(const Vec3& p0, const Vec3& p1, Real circulation) {
	mSegments.push_back(VortexFilamentData(add(BasicParticleData(p0)), add(BasicParticleData(p1)), circulation));
}

VortexFilamentSystem::VortexFilamentSystem(FluidSolver* parent) :
    ConnectedParticleSystem<BasicParticleData, VortexFilamentData>(parent)
{ 
}

DefineIntegrator(integrateVortexKernel, FilamentKernel, evaluate);

KERNEL(particle) template<IntegrationMode mode> 
void advectNodes(vector<Vec3>& nodesNew, const vector<Vec3>& nodesOld, const FilamentKernel& kernel, Real dt) {
    const Vec3 p = nodesOld[i];
    if (kernel.isValid(p))
        nodesNew[i] += integrateVortexKernel<mode>(p, kernel, dt);    
}

void VortexFilamentSystem::advectSelf(Real scale, Real regularization, int integrationMode) {
    const Real dt = getParent()->getDt();
    const Real cutoff = 1e7;
    
    // backup
    vector<Vec3> nodesOld(size()), nodesNew(size());
    for (int i=0; i<size(); i++)
        nodesOld[i] = nodesNew[i] = mData[i].pos;
    
    // loop over the vortices
    for (size_t i=0; i<mSegments.size(); i++) {
        if (!isSegActive(i)) continue;
        FilamentKernel kernel = FilamentKernel(mData[mSegments[i].idx0].pos, mData[mSegments[i].idx1].pos, mSegments[i].circulation, scale, cutoff, regularization);
        if (integrationMode==EULER) advectNodes<EULER>(nodesNew, nodesOld, kernel, dt);
        else if (integrationMode==RK2) advectNodes<RK2>(nodesNew, nodesOld, kernel, dt);
        else if (integrationMode==RK4) advectNodes<RK4>(nodesNew, nodesOld, kernel, dt);
        else errMsg("unknown integration type");
    }
    
    for (int i=0; i<size(); i++)
        mData[i].pos = nodesNew[i];
}

void VortexFilamentSystem::applyToMesh(Mesh& mesh, Real scale, Real regularization, int integrationMode) {
    const Real dt = getParent()->getDt();
    const Real cutoff = 1e7;
    
    // copy node array
    const int nodes = mesh.numNodes();
    vector<Vec3> nodesOld(nodes), nodesNew(nodes);
    for (int i=0; i<nodes; i++)
        nodesOld[i] = nodesNew[i] = mesh.nodes(i).pos;
    
    // loop over the vortices
    for (size_t i=0; i<mSegments.size(); i++) {
        if (!isSegActive(i)) continue;
        FilamentKernel kernel = FilamentKernel(mData[mSegments[i].idx0].pos, mData[mSegments[i].idx1].pos, mSegments[i].circulation, scale, cutoff, regularization);
        if (integrationMode==EULER) advectNodes<EULER>(nodesNew, nodesOld, kernel, dt);
        else if (integrationMode==RK2) advectNodes<RK2>(nodesNew, nodesOld, kernel, dt);
        else if (integrationMode==RK4) advectNodes<RK4>(nodesNew, nodesOld, kernel, dt);
        else errMsg("unknown integration type");
    }
    
    // copy back
    for (int i=0; i<nodes; i++) {
        if (!mesh.isNodeFixed(i))
            mesh.nodes(i).pos = nodesNew[i];
    }    
}

ParticleBase* VortexFilamentSystem::clone() {
    VortexFilamentSystem* nm = new VortexFilamentSystem(getParent());
    //compress();
    
    nm->mData = mData;
    nm->mSegments = mSegments;
    nm->setName(getName());
    return nm;
}


} // namespace
