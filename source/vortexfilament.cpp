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
#include "interpol.h"
#include "mesh.h"
#include "quaternion.h"

using namespace std;
namespace Manta {

void VortexRing::renumber(int *_renumber) {
	for (size_t i=0; i<indices.size(); i++)
		indices[i] = _renumber[indices[i]];
}

inline Vec3 FilamentKernel(const Vec3& pos, const vector<VortexRing>& rings, const vector<BasicParticleData>& fp, Real reg, Real cutoff, Real scale) {
	const Real strength = 0.25 / M_PI * scale;
	const Real a2 = square(reg);
	const Real cutoff2 = square(cutoff);
	const Real mindist = 1e-6;
	Vec3 u(_0);
	
	for (size_t i=0; i<rings.size(); i++) {
		const VortexRing& r = rings[i];
		if (r.flag & ParticleBase::PDELETE) continue;
		
		const int N = r.isClosed ? (r.size()) : (r.size()-1);
		const Real str = strength * r.circulation;
		for (int j=0; j<N; j++) {
			const Vec3 r0 = fp[r.idx0(j)].pos - pos;
			const Vec3 r1 = fp[r.idx1(j)].pos - pos;
			const Real r0_2 = normSquare(r0), r1_2 = normSquare(r1);
			if (r0_2 > cutoff2 || r1_2 > cutoff2 || r0_2 < mindist || r1_2 < mindist)
				continue;
			
			const Vec3 e = getNormalized(r1-r0);
			const Real r0n = 1.0f/sqrt(a2+r0_2);
			const Real r1n = 1.0f/sqrt(a2+r1_2);
			const Vec3 cp = cross(r0,e);
			const Real A = str * (dot(r1,e)*r1n - dot(r0,e)*r0n) / (a2 + normSquare(cp));
			u += A * cp;
		}
	}
	return u;
}

KERNEL(pts) returns(vector<Vec3> u(size))
vector<Vec3> KnFilamentAdvectParts(vector<BasicParticleData>& nodes, vector<BasicParticleData>& fp, const vector<VortexRing>& rings, Real reg, Real cutoff, Real scale) {
	if (nodes[idx].flag & ParticleBase::PDELETE)
		u[idx] = _0;
	else
		u[idx] = FilamentKernel(nodes[idx].pos, rings, fp, reg, cutoff, scale);
}

KERNEL(pts) returns(vector<Vec3> u(size))
vector<Vec3> KnFilamentAdvectMesh(vector<Node>& nodes, const vector<VortexRing>& rings, const vector<BasicParticleData>& fp, Real reg, Real cutoff, Real scale) {
	if (nodes[idx].flags & Mesh::NfFixed)
		u[idx] = _0;
	else
		u[idx] = FilamentKernel(nodes[idx].pos, rings, fp, reg, cutoff, scale);
}

void VortexFilamentSystem::advectSelf(Real scale, Real regularization, int integrationMode) {
	KnFilamentAdvectParts kernel(mData, mData, mSegments , regularization, 1e10, scale * getParent()->getDt());
	integratePointSet( kernel, integrationMode);
}

void VortexFilamentSystem::advectMesh(Mesh& mesh, Real scale, Real regularization, int integrationMode) {
	KnFilamentAdvectMesh kernel(mesh.getNodeData(), mSegments, mData, regularization, 1e10, scale * getParent()->getDt());
	integratePointSet( kernel, integrationMode);   
}

void VortexFilamentSystem::advectParticles(BasicParticleSystem& sys, Real scale, Real regularization, int integrationMode) {
	KnFilamentAdvectParts kernel(sys.getData(), mData, mSegments, regularization, 1e10, scale * getParent()->getDt());
	integratePointSet( kernel, integrationMode);   
}

void VortexFilamentSystem::remesh(Real maxLen, Real minLen) {
	const Real maxLen2 = maxLen*maxLen, minLen2 = minLen*minLen;
	
	for (int i=0; i < segSize(); i++) {
		VortexRing& r = mSegments[i];
			
		// insert edges
		for(;;) {
			const int oldLen = r.size();
			map<int,int> insert;        
			int offset = 1;
			
			for (int j=0; j<oldLen; j++) {
				const Vec3 p0 = mData[r.idx0(j)].pos;
				const Vec3 p1 = mData[r.idx1(j)].pos;
				const Real l2 = normSquare(p1-p0);
				
				if (l2 > maxLen2) {
					// insert midpoint
					const Vec3 p_1 = mData[r.idx(j-1)].pos;
					const Vec3 p2 = mData[r.idx(j+2)].pos;
					const Vec3 mp = hermiteSpline(p0,p1,crTangent(p_1,p0,p1),crTangent(p0,p1,p2), 0.5);
					insert.insert(pair<int,int>(j+offset, add(mp)));
					offset++;
				}
			}
			if (insert.empty()) 
				break;
			
			// renumber indices
			const int newLen = oldLen + insert.size();
			int num=oldLen-1;
			r.indices.resize(newLen);
			for (int j=newLen-1; j>=0; j--) {
				map<int,int>::const_iterator f = insert.find(j);
				if (f==insert.end())
					r.indices[j] = r.indices[num--];
				else
					r.indices[j] = f->second;
			}
		}
		
		// remove edges
		for(;;) {
			const int oldLen = r.size();
			const int N = r.isClosed ? oldLen : (oldLen-1);
			std::vector<bool> deleted(r.size());
			
			int newLen=oldLen;
			for (int j=0; j<N; j++) {
				if (mData[r.idx0(j)].flag & PDELETE || mData[r.idx1(j)].flag & PDELETE) continue;
				const Vec3 p0 = mData[r.idx0(j)].pos;
				const Vec3 p1 = mData[r.idx1(j)].pos;
				const Real l2 = normSquare(p1-p0);
				
				if (l2 < minLen2) {
					// kill edge
					mData[r.idx0(j)].flag |= PDELETE;
					mData[r.idx1(j)].pos = 0.5*(p0+p1);
					deleted[j] = true;
					newLen--;
					j++;
				}
			}
			if (newLen == oldLen)
				break;
			
			// renumber indices
			for (int j=0, copyFrom=0; j<newLen; j++,copyFrom++) {
				while (deleted[copyFrom]) 
					copyFrom++;
				if (j!=copyFrom)
					r.indices[j] = r.indices[copyFrom];
			}            
			r.indices.resize(newLen);
		}
	}
	
	// remove deleted particles
	compress();
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
// doesn't really work yet (can't reverse rotation dir)
// ------------------------------------------------------------------------------

Real evaluateRefU(int N, Real L, Real circ, Real reg) {
	// construct regular n-polygon
	const Real l = L/(Real)N;
	const Real r = 0.5*l/sin(M_PI/(Real)N);
	cout << r << " " << l << endl;
	// build vortex ring
	VortexRing ring (circ);
	vector<BasicParticleData> pos(N);
	for(int i=0; i<N; i++) {
		pos[i].pos = Vec3( r*cos(2.0*M_PI*(Real)i/N), r*sin(2.0*M_PI*(Real)i/N), 0);
		pos[i].flag =0;
		ring.indices.push_back(i);
	}
	
	// Build kernel
	vector<VortexRing> rings;
	rings.push_back(ring);
	
	// evaluate impact on pos[0]
	return norm(FilamentKernel(pos[0].pos, rings, pos, reg, 1e10, 1.0));    
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
	Vec3 lT(0,0,l);
	if (!powerMethod(from, l, r, lT))
		return false;
	cout << "iniLT " << lT << " norm " << lT/l<< endl;
	
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
		if (!isSegActive(rc) || !mSegments[rc].isClosed) continue;
		
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
		cout << U << " <-< " << Ur << endl;
		
		// fwd darboux transform
		vector<Vec3> eta(N);
		if (!darboux(gamma, eta, l, ra)) {
			cout << "Fwd Darboux correction failed, skipped." << endl; 
			continue;
		}
		
		// bwd darboux transform
		if (!darboux(eta, gamma, l, ra)) {
			cout << "Bwd Darboux correction failed, skipped." << endl; 
			continue;
		}
		
		// copy back
		for (int i=0; i<N; i++) {
			mData[r.indices[i]].pos = gamma[i];
		}
	}
}

void VortexFilamentSystem::addLine(const Vec3& p0, const Vec3& p1, Real circulation) {
	VortexRing ring(circulation, false);
	
	ring.indices.push_back(add(BasicParticleData(p0)));
	ring.indices.push_back(add(BasicParticleData(p1)));
	mSegments.push_back(ring);
}

void VortexFilamentSystem::addRing(const Vec3& position, Real circulation, Real radius, Vec3 normal, int number) {
	normalize(normal);
	Vec3 worldup (0,1,0);
	if (norm(normal - worldup) < 1e-5) worldup = Vec3(1,0,0);
	
	Vec3 u = cross(normal, worldup); normalize(u);
	Vec3 v = cross(normal, u); normalize(v);
	
	VortexRing ring(circulation, true);
	
	for (int i=0; i<number; i++) {
		Real phi = (Real)i/(Real)number * M_PI * 2.0;
		Vec3 p = position + radius * (u*cos(phi) + v*sin(phi));
		
		int num = add(BasicParticleData(p));
		ring.indices.push_back(num);
	}
	mSegments.push_back(ring);
}

	
} // namespace
