/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Conjugate gradient solver
 *
 ******************************************************************************/

#ifndef _CONJUGATEGRADIENT_H
#define _CONJUGATEGRADIENT_H

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"

namespace Manta { 

static const bool CG_DEBUG = false;

//! Basic CG interface 
class GridCgInterface {
	public:
		enum PreconditionType { PC_None=0, PC_ICP, PC_mICP };
        
        GridCgInterface() : mUseResNorm(true) {};
		virtual ~GridCgInterface() {};

		// solving functions
		virtual bool iterate() = 0;
		virtual void solve(int maxIter) = 0;

		// precond
		virtual void setPreconditioner(PreconditionType method, Grid3<Real> *A0, Grid3<Real> *Ai, Grid3<Real> *Aj, Grid3<Real> *Ak) = 0;

		// access
		virtual Real getSigma() const = 0;
		virtual Real getIterations() const = 0;
		virtual Real getResNorm() const = 0;
		virtual void setAccuracy(Real set) = 0;
		virtual Real getAccuracy() const = 0;

		void setUseResNorm(bool set) { mUseResNorm = set; }

	protected:

		// use norm of residual, or max value for threshold?
		bool mUseResNorm; 
};


//! Run single iteration of the cg solver
/*! the template argument determines the type of matrix multiplication,
    typically a ApplyMatrix kernel, another one is needed e.g. for the
    mesh-based wave equation solver */
template<class APPLYMAT>
class GridCg : public GridCgInterface {
	public:
        //! constructor
		GridCg(Grid3<Real>& dst, Grid3<Real>& rhs, Grid3<Real>& residual, Grid3<Real>& search, FlagGrid3& flags, Grid3<Real>& tmp, 
				Grid3<Real>* A0, Grid3<Real>* pAi, Grid3<Real>* pAj, Grid3<Real>* pAk);
        ~GridCg() {}
        
        void doInit();
        bool iterate();
        void solve(int maxIter);
        //! init pointers, and copy values from "normal" matrix
        void setPreconditioner(PreconditionType method, Grid3<Real> *A0, Grid3<Real> *Ai, Grid3<Real> *Aj, Grid3<Real> *Ak);
        
        // Accessors        
        Real getSigma() const { return mSigma; }
        Real getIterations() const { return mIterations; }

        Real getResNorm() const { return mResNorm; }

        void setAccuracy(Real set) { mAccuracy=set; }
        Real getAccuracy() const { return mAccuracy; }

	protected:
		bool mInited;
		int mIterations;
		// grids
		Grid3<Real>& mDst;
		Grid3<Real>& mRhs;
		Grid3<Real>& mResidual;
		Grid3<Real>& mSearch;
		FlagGrid3& mFlags;
		Grid3<Real>& mTmp;

		Grid3<Real> *mpA0, *mpAi, *mpAj, *mpAk;

		PreconditionType mPcMethod;
		//! preconditioning grids
		Grid3<Real> *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk;

		//! sigma / residual
		Real mSigma;
		//! accuracy of solver (max. residuum)
		Real mAccuracy;
		//! norm of the residual
		Real mResNorm;
}; // GridCg


//! Kernel: Apply symmetric stored Matrix
KERNEL(idx) 
ApplyMatrix (FlagGrid3& flags, Grid3<Real>& dst, Grid3<Real>& src, 
             Grid3<Real>& A0, Grid3<Real>& Ai, Grid3<Real>& Aj, Grid3<Real>& Ak)
{
    if (!flags.isFluid(idx)) {
        dst[idx] = src[idx];
        return;
    }    
    dst[idx] =  src[idx] * A0[idx]
                + src[idx-X] * Ai[idx-X]
                + src[idx+X] * Ai[idx]
                + src[idx-Y] * Aj[idx-Y]
                + src[idx+Y] * Aj[idx]
                + src[idx-Z] * Ak[idx-Z] 
                + src[idx+Z] * Ak[idx];
}

//! Kernel: Construct the matrix for the poisson equation
KERNEL (bnd=1) MakeLaplaceMatrix(FlagGrid3& flags, Grid3<Real>& A0, Grid3<Real>& Ai, Grid3<Real>& Aj, Grid3<Real>& Ak) {
    if (!flags.isFluid(i,j,k))
        return;
    
    // center
    if (!flags.isObstacle(i-1,j,k)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i+1,j,k)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i,j-1,k)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i,j+1,k)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i,j,k-1)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i,j,k+1)) A0(i,j,k) += 1.;
    
    if (flags.isFluid(i+1,j,k)) Ai(i,j,k) = -1.;
    if (flags.isFluid(i,j+1,k)) Aj(i,j,k) = -1.;
    if (flags.isFluid(i,j,k+1)) Ak(i,j,k) = -1.;    
}




} // namespace

#endif 