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
#include "quaternion.h"

using namespace std;
namespace Manta {

void VortexRing::renumber(int *_renumber) {
    for (size_t i=0; i<indices.size(); i++)
        indices[i] = _renumber[indices[i]];
}
    
// vortex filament effect
struct FilamentKernel {
    FilamentKernel() : gamma0(0.), gamma1(0.), gammaMid(0.), gammaDir(0.), strength(0.), cutoff2(0.), a2(0.) {}
    FilamentKernel(const Vec3& g0, const Vec3& g1, Real circ, Real scale, Real cutoff, Real a) : 
        gamma0(g0), gamma1(g1), cutoff2(cutoff*cutoff), a2(a*a) {
            gammaMid = 0.5*(g0+g1);
            gammaDir = getNormalized(g1-g0);
            strength = 0.25/M_PI*scale*circ;
        }

    Vec3 gamma0, gamma1, gammaMid, gammaDir;
    Real strength;
    Real cutoff2;
    Real a2;
    
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
        const Vec3 cp = cross(r0,gammaDir);
        const Real upper = dot(r1,gammaDir)*r1n - dot(r0,gammaDir)*r0n;
        const Real lower = a2 + normSquare(cp);
        return (strength*upper/lower)*cp;
    }
};

DefineIntegrator(integrateVortexKernel, FilamentKernel, evaluate);

KERNEL(particle) template<IntegrationMode mode> 
void advectNodes(vector<Vec3>& nodesNew, const vector<Vec3>& nodesOld, const FilamentKernel& kernel, Real dt) {
    const Vec3 p = nodesOld[i];
    if (kernel.isValid(p))
        nodesNew[i] += integrateVortexKernel<mode>(p, kernel, dt);    
}

void VortexFilamentSystem::integrate(const vector<Vec3>& nodesOld, vector<Vec3>& nodesNew, Real scale, Real reg, int integrationMode) {
    const Real dt = getParent()->getDt();
    const Real cutoff = 1e7;
    
    // loop over the vortices
    for (size_t i=0; i<mSegments.size(); i++) {
        if (!isSegActive(i)) continue;
        const VortexRing& r = mSegments[i];
        for (int j=0; j<r.size(); j++) {
            FilamentKernel kernel = FilamentKernel(mData[r.idx0(j)].pos, mData[r.idx1(j)].pos, mSegments[i].circulation, scale, cutoff, reg);
            if (integrationMode==EULER) advectNodes<EULER>(nodesNew, nodesOld, kernel, dt);
            else if (integrationMode==RK2) advectNodes<RK2>(nodesNew, nodesOld, kernel, dt);
            else if (integrationMode==RK4) advectNodes<RK4>(nodesNew, nodesOld, kernel, dt);
            else errMsg("unknown integration type");
        }
    }    
}

void VortexFilamentSystem::advectSelf(Real scale, Real regularization, int integrationMode) {    
    // backup
    vector<Vec3> nodesOld(size()), nodesNew(size());
    for (int i=0; i<size(); i++)
        nodesOld[i] = nodesNew[i] = mData[i].pos;
    
    integrate(nodesOld, nodesNew, scale, regularization, integrationMode);
    
    for (int i=0; i<size(); i++)
        mData[i].pos = nodesNew[i];
}

void VortexFilamentSystem::applyToMesh(Mesh& mesh, Real scale, Real regularization, int integrationMode) {
    // copy node array
    const int nodes = mesh.numNodes();
    vector<Vec3> nodesOld(nodes), nodesNew(nodes);
    for (int i=0; i<nodes; i++)
        nodesOld[i] = nodesNew[i] = mesh.nodes(i).pos;
    
    integrate(nodesOld, nodesNew, scale, regularization, integrationMode);
    
    // copy back
    for (int i=0; i<nodes; i++) {
        if (!mesh.isNodeFixed(i))
            mesh.nodes(i).pos = nodesNew[i];
    }    
}

VortexFilamentSystem::VortexFilamentSystem(FluidSolver* parent) :
    ConnectedParticleSystem<BasicParticleData, VortexRing>(parent)
{     
}

ParticleBase* VortexFilamentSystem::clone() {
    VortexFilamentSystem* nm = new VortexFilamentSystem(getParent());
    compress();
    
    nm->mData = mData;
    nm->mSegments = mSegments;
    nm->setName(getName());
    return nm;
}

// ------------------------------------------------------------------------------
// Functions needed for doubly-discrete smoke flow using Darboux transforms
// see [Weissmann,Pinkall 2009]
// ------------------------------------------------------------------------------

Real evaluateRefU(int N, Real L, Real circ, Real reg) {
    // construct regular n-polygon
    const Real l = L/(Real)N;
    const Real r = 0.5*l/sin(M_PI/(Real)N);
    
    vector<Vec3> pos(N);
    for(int i=0; i<N; i++) {
        pos[i] = Vec3( r*cos(2.0*M_PI*(Real)i/N), r*sin(2.0*M_PI*(Real)i/N), 0);
    }
    
    // evaluate impact on pos[0]
    Vec3 sum(0.);
    for(int i=1; i<N-1; i++) {
        FilamentKernel kernel (pos[i], pos[i+1], circ, 1.0, 1e10, reg);
        sum += kernel.evaluate(pos[0]);
    }
    return norm(sum);
}

Vec3 darbouxStep(const Vec3& Si, const Vec3& lTi, Real r) {
    Quaternion rlTS (lTi - Si, -r);
    Quaternion lT (lTi, 0);
    Quaternion lTnext = rlTS * lT * rlTS.inverse();
    return lTnext.imag();
}

Vec3 monodromy(const vector<Vec3>& gamma, const Vec3& lT_1, Real r) {
    const int N = gamma.size();
    Vec3 lT (lT_1);
    
    for (int i=0; i<N; i++) {
        Vec3 Si = gamma[(i+1)%N]-gamma[i];
        lT = darbouxStep(Si, lT, r);
    }
    return lT;
}

bool powerMethod(const vector<Vec3>& gamma, Real l, Real r, Vec3& lT) {
    const int maxIter = 100;
    const Real epsilon = 1e-4;
    
    lT = Vec3(0,l,0);
    for (int i=0; i<maxIter; i++) {
        Vec3 lastLT (lT);
        lT = monodromy(gamma, lT, r);
        //if ((i%1) == 0) cout << "iteration " << i << " residual: " << norm(lT-lastLT) << endl;
        if (norm(lT-lastLT) < epsilon) 
            return true;
    }   
    return false;
}

bool darboux(const vector<Vec3>& from, vector<Vec3>& to, Real l, Real r) {
    const int N = from.size();
    Vec3 lT(0.);
    if (!powerMethod(from, l, r, lT))
        return false;
    
    for (int i=0; i<N; i++) {
        to[i] = from[i] + lT;
        Vec3 Si = from[(i+1)%N] - from[i];
        lT = darbouxStep(Si, lT, r);
    }
    return true;
}
        

void VortexFilamentSystem::doublyDiscreteUpdate(Real reg) {
    const Real dt = getParent()->getDt();
    
    for (int rc=0; rc<segSize(); rc++) {
        if (!isSegActive(rc)) continue;
        
         VortexRing& r = mSegments[rc];
         int N = r.size();
        
        // compute arc length
        Real L=0;
        for (int i=0; i<N; i++)
            L += norm(mData[r.idx0(i)].pos - mData[r.idx1(i)].pos);
        
        // build gamma
        vector<Vec3> gamma(N);
        for (int i=0; i<N; i++) gamma[i] = mData[r.indices[i]].pos;
        
        //N=1000; L=2.0*M_PI; reg=0.1; r.circulation=1;
        
        // compute reference parameters 
        const Real U = 0.5*r.circulation/L * (log(4.0*L/(M_PI*reg)) - 1.0);
        const Real Ur = evaluateRefU(N, L, r.circulation, reg);
        const Real d = 0.5*dt*(U-Ur);
        const Real l = sqrt( square(L/N) + square(d) );
        const Real ra = d*tan(M_PI * (0.5 - 1.0/N)); // d*cot(pi/n)
        
        /*cout << "L   : " << L << endl;
        cout << "N   : " << N << endl;
        cout << "U   : " << U << endl;
        cout << "Ur  : " << Ur << endl;
        cout << "d   : " << d << endl;
        cout << "l   : " << l << endl;
        cout << "ra  : " << ra << endl;*/
        
        // fwd darboux transform
        vector<Vec3> eta(N);
        if (!darboux(gamma, eta, l, ra)) {
            cout << "Fwd Darboux correction failed, skipped." << endl; 
            continue;
        }
        
        // bwd darboux transform
        if (!darboux(eta, gamma, l, -ra)) {
            cout << "Bwd Darboux correction failed, skipped." << endl; 
            continue;
        }
        
        // copy back
        for (int i=0; i<N; i++) {
            mData[r.indices[i]].pos = gamma[i];
        }        
    }
}


void VortexFilamentSystem::addRing(const Vec3& position, Real circulation, Real radius, Vec3 normal, int number) {
    normalize(normal);
    Vec3 worldup (0,1,0);
    if (norm(normal - worldup) < 1e-5) worldup = Vec3(1,0,0);
    
    Vec3 u = cross(normal, worldup); normalize(u);
    Vec3 v = cross(normal, u); normalize(v);
    
    VortexRing ring(circulation);
    
    for (int i=0; i<number; i++) {
        Real phi = (Real)i/(Real)number * M_PI * 2.0;
        Vec3 p = position + radius * (u*cos(phi) + v*sin(phi));
        
        int num = add(BasicParticleData(p));
        ring.indices.push_back(num);
    }
    mSegments.push_back(ring);
}

    
} // namespace
