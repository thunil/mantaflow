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

#include "conjugategrad.h"
#include "commonkernels.h"

using namespace std;
namespace Manta {

const int CG_DEBUGLEVEL = 5;
	
//*****************************************************************************
//  Precondition helpers

//! Preconditioning a la Wavelet Turbulence (needs 4 add. grids)
void InitPreconditionIncompCholesky(FlagGrid& flags,
				Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak,
				Grid<Real>& orgA0, Grid<Real>& orgAi, Grid<Real>& orgAj, Grid<Real>& orgAk) 
{
	// compute IC according to Golub and Van Loan
	A0.copyFrom( orgA0 );
	Ai.copyFrom( orgAi );
	Aj.copyFrom( orgAj );
	Ak.copyFrom( orgAk );
	
	FOR_IJK(A0) {
		if (flags.isFluid(i,j,k)) {
			const IndexInt idx = A0.index(i,j,k);
			A0[idx] = sqrt(A0[idx]);
			
			// correct left and top stencil in other entries
			// for i = k+1:n
			//  if (A(i,k) != 0)
			//    A(i,k) = A(i,k) / A(k,k)
			Real invDiagonal = 1.0f / A0[idx];
			Ai[idx] *= invDiagonal;
			Aj[idx] *= invDiagonal;
			Ak[idx] *= invDiagonal;
			
			// correct the right and bottom stencil in other entries
			// for j = k+1:n 
			//   for i = j:n
			//      if (A(i,j) != 0)
			//        A(i,j) = A(i,j) - A(i,k) * A(j,k)
			A0(i+1,j,k) -= square(Ai[idx]);
			A0(i,j+1,k) -= square(Aj[idx]);
			A0(i,j,k+1) -= square(Ak[idx]);
		}
	}
	
	// invert A0 for faster computation later
	InvertCheckFluid (flags, A0);
};

//! Preconditioning using modified IC ala Bridson (needs 1 add. grid)
void InitPreconditionModifiedIncompCholesky2(FlagGrid& flags,
				Grid<Real>&Aprecond, 
				Grid<Real>&A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak) 
{
	// compute IC according to Golub and Van Loan
	Aprecond.clear();
	
	FOR_IJK(flags) {
		if (!flags.isFluid(i,j,k)) continue;

		const Real tau = 0.97;
		const Real sigma = 0.25;
			
		// compute modified incomplete cholesky
		Real e = 0.;
		e = A0(i,j,k) 
			- square(Ai(i-1,j,k) * Aprecond(i-1,j,k) )
			- square(Aj(i,j-1,k) * Aprecond(i,j-1,k) )
			- square(Ak(i,j,k-1) * Aprecond(i,j,k-1) ) ;
		e -= tau * (
				Ai(i-1,j,k) * ( Aj(i-1,j,k) + Ak(i-1,j,k) )* square( Aprecond(i-1,j,k) ) +
				Aj(i,j-1,k) * ( Ai(i,j-1,k) + Ak(i,j-1,k) )* square( Aprecond(i,j-1,k) ) +
				Ak(i,j,k-1) * ( Ai(i,j,k-1) + Aj(i,j,k-1) )* square( Aprecond(i,j,k-1) ) +
				0. );

		// stability cutoff
		if(e < sigma * A0(i,j,k))
			e = A0(i,j,k);

		Aprecond(i,j,k) = 1. / sqrt( e );
	}
};

//! Preconditioning using multigrid ala Dick et al.
void InitPreconditionMultigrid(GridMg* MG, Grid<Real>&A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak, Real mAccuracy) 
{
	// build multigrid hierarchy if necessary
	if (!MG->isASet()) MG->setA(&A0, &Ai, &Aj, &Ak);
	MG->setCoarsestLevelAccuracy(mAccuracy * 1E-4);
	MG->setSmoothing(1,1);
};

//! Apply WT-style ICP
void ApplyPreconditionIncompCholesky(Grid<Real>& dst, Grid<Real>& Var1, FlagGrid& flags,
				Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak,
				Grid<Real>& orgA0, Grid<Real>& orgAi, Grid<Real>& orgAj, Grid<Real>& orgAk)
{
	
	// forward substitution        
	FOR_IJK(dst) {
		if (!flags.isFluid(i,j,k)) continue;
		dst(i,j,k) = A0(i,j,k) * (Var1(i,j,k)
				 - dst(i-1,j,k) * Ai(i-1,j,k)
				 - dst(i,j-1,k) * Aj(i,j-1,k)
				 - dst(i,j,k-1) * Ak(i,j,k-1));    
	}
	
	// backward substitution
	FOR_IJK_REVERSE(dst) {
		const IndexInt idx = A0.index(i,j,k);
		if (!flags.isFluid(idx)) continue;
		dst[idx] = A0[idx] * ( dst[idx] 
			   - dst(i+1,j,k) * Ai[idx]
			   - dst(i,j+1,k) * Aj[idx]
			   - dst(i,j,k+1) * Ak[idx]);
	}
}

//! Apply Bridson-style mICP
void ApplyPreconditionModifiedIncompCholesky2(Grid<Real>& dst, Grid<Real>& Var1, FlagGrid& flags,
				Grid<Real>& Aprecond, 
				Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak) 
{
	// forward substitution        
	FOR_IJK(dst) {
		if (!flags.isFluid(i,j,k)) continue;
		const Real p = Aprecond(i,j,k);
		dst(i,j,k) = p * (Var1(i,j,k)
				 - dst(i-1,j,k) * Ai(i-1,j,k) * Aprecond(i-1,j,k)
				 - dst(i,j-1,k) * Aj(i,j-1,k) * Aprecond(i,j-1,k)
				 - dst(i,j,k-1) * Ak(i,j,k-1) * Aprecond(i,j,k-1) );
	}
	
	// backward substitution
	FOR_IJK_REVERSE(dst) {            
		const IndexInt idx = A0.index(i,j,k);
		if (!flags.isFluid(idx)) continue;
		const Real p = Aprecond[idx];
		dst[idx] = p * ( dst[idx] 
			   - dst(i+1,j,k) * Ai[idx] * p
			   - dst(i,j+1,k) * Aj[idx] * p
			   - dst(i,j,k+1) * Ak[idx] * p);
	}
}

//! Perform one Multigrid VCycle
void ApplyPreconditionMultigrid(GridMg* pMG, Grid<Real>& dst, Grid<Real>& Var1) 
{
	// one VCycle on "A*dst = Var1" with initial guess dst=0
	pMG->setRhs(Var1);
	pMG->doVCycle(dst); 
}


//*****************************************************************************
// Kernels    

//! Kernel: Compute the dot product between two Real grids
/*! Uses double precision internally */
KERNEL(idx, reduce=+) returns(double result=0.0)
double GridDotProduct (const Grid<Real>& a, const Grid<Real>& b) {
	result += (a[idx] * b[idx]);    
};

//! Kernel: compute residual (init) and add to sigma
KERNEL(idx, reduce=+) returns(double sigma=0)
double InitSigma (FlagGrid& flags, Grid<Real>& dst, Grid<Real>& rhs, Grid<Real>& temp) 
{    
	const double res = rhs[idx] - temp[idx]; 
	dst[idx] = (Real)res;

	// only compute residual in fluid region
	if(flags.isFluid(idx)) 
		sigma += res*res;
};

//! Kernel: update search vector
KERNEL(idx) void UpdateSearchVec (Grid<Real>& dst, Grid<Real>& src, Real factor)
{
	dst[idx] = src[idx] + factor * dst[idx];
}

//*****************************************************************************
//  CG class

template<class APPLYMAT>
GridCg<APPLYMAT>::GridCg(Grid<Real>& dst, Grid<Real>& rhs, Grid<Real>& residual, Grid<Real>& search, FlagGrid& flags, Grid<Real>& tmp, 
			   Grid<Real>* pA0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk) :
	GridCgInterface(), mInited(false), mIterations(0), mDst(dst), mRhs(rhs), mResidual(residual),
	mSearch(search), mFlags(flags), mTmp(tmp), mpA0(pA0), mpAi(pAi), mpAj(pAj), mpAk(pAk),
	mPcMethod(PC_None), mpPCA0(nullptr), mpPCAi(nullptr), mpPCAj(nullptr), mpPCAk(nullptr), mMG(nullptr), mSigma(0.), mAccuracy(VECTOR_EPSILON), mResNorm(1e20) 
{
	dst.clear();
	residual.clear();
	search.clear();
	tmp.clear();            
}

template<class APPLYMAT>
void GridCg<APPLYMAT>::doInit() {
	mInited = true;

	mResidual.copyFrom( mRhs ); // p=0, residual = b
	
	if (mPcMethod == PC_ICP) {
		assertMsg(mDst.is3D(), "ICP only supports 3D grids so far");
		InitPreconditionIncompCholesky(mFlags, *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk, *mpA0, *mpAi, *mpAj, *mpAk);
		ApplyPreconditionIncompCholesky(mTmp, mResidual, mFlags, *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk, *mpA0, *mpAi, *mpAj, *mpAk);
	} else if (mPcMethod == PC_mICP) {
		assertMsg(mDst.is3D(), "mICP only supports 3D grids so far");
		InitPreconditionModifiedIncompCholesky2(mFlags, *mpPCA0, *mpA0, *mpAi, *mpAj, *mpAk);
		ApplyPreconditionModifiedIncompCholesky2(mTmp, mResidual, mFlags, *mpPCA0, *mpA0, *mpAi, *mpAj, *mpAk);
	} else if (mPcMethod == PC_MGP) {
		InitPreconditionMultigrid(mMG, *mpA0, *mpAi, *mpAj, *mpAk, mAccuracy);
		ApplyPreconditionMultigrid(mMG, mTmp, mResidual);
	} else {
		mTmp.copyFrom( mResidual );
	}
	
	mSearch.copyFrom( mTmp );
	
	mSigma = GridDotProduct(mTmp, mResidual);    
}

template<class APPLYMAT>
bool GridCg<APPLYMAT>::iterate() {
	if(!mInited) doInit();

	mIterations++;

	// create matrix application operator passed as template argument,
	// this could reinterpret the mpA pointers (not so clean right now)
	// tmp = applyMat(search)
	
	APPLYMAT (mFlags, mTmp, mSearch, *mpA0, *mpAi, *mpAj, *mpAk);
	
	// alpha = sigma/dot(tmp, search)
	Real dp = GridDotProduct(mTmp, mSearch);
	Real alpha = 0.;
	if(fabs(dp)>0.) alpha = mSigma / (Real)dp;
	
	gridScaledAdd<Real,Real>(mDst, mSearch, alpha);    // dst += search * alpha
	gridScaledAdd<Real,Real>(mResidual, mTmp, -alpha); // residual += tmp * -alpha
	
	if (mPcMethod == PC_ICP)
		ApplyPreconditionIncompCholesky(mTmp, mResidual, mFlags, *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk, *mpA0, *mpAi, *mpAj, *mpAk);
	else if (mPcMethod == PC_mICP)
		ApplyPreconditionModifiedIncompCholesky2(mTmp, mResidual, mFlags, *mpPCA0, *mpA0, *mpAi, *mpAj, *mpAk);
	else if (mPcMethod == PC_MGP)
		ApplyPreconditionMultigrid(mMG, mTmp, mResidual);
	else
		mTmp.copyFrom( mResidual );
		
	// use the l2 norm of the residual for convergence check? (usually max norm is recommended instead)
	if(this->mUseL2Norm) { 
		mResNorm = GridSumSqr(mResidual).sum; 
	} else {
		mResNorm = mResidual.getMaxAbsValue();        
	}
	//if(mIterations % 10 == 9) debMsg("GridCg::Iteration i="<<mIterations<<", resNorm="<<mResNorm<<" accuracy="<<mAccuracy, 1);

	// abort here to safe some work...
	if(mResNorm<mAccuracy) {
		mSigma = mResNorm; // this will be returned later on to the caller...
		return false;
	}

	Real sigmaNew = GridDotProduct(mTmp, mResidual);
	Real beta = sigmaNew / mSigma;
	
	// search =  tmp + beta * search
	UpdateSearchVec (mSearch, mTmp, beta);

	debMsg("GridCg::iterate i="<<mIterations<<" sigmaNew="<<sigmaNew<<" sigmaLast="<<mSigma<<" alpha="<<alpha<<" beta="<<beta<<" ", CG_DEBUGLEVEL);
	mSigma = sigmaNew;
	
	//debMsg("PB-CG-Norms::p"<<sqrt( GridOpNormNosqrt(mpDst, mpFlags).getValue() ) <<" search"<<sqrt( GridOpNormNosqrt(mpSearch, mpFlags).getValue(), CG_DEBUGLEVEL ) <<" res"<<sqrt( GridOpNormNosqrt(mpResidual, mpFlags).getValue() ) <<" tmp"<<sqrt( GridOpNormNosqrt(mpTmp, mpFlags).getValue() ), CG_DEBUGLEVEL ); // debug
	return true;
}

template<class APPLYMAT>
void GridCg<APPLYMAT>::solve(int maxIter) {
	for (int iter=0; iter<maxIter; iter++) {
		if (!iterate()) iter=maxIter;
	} 
	return;
}

static bool gPrint2dWarning = true;
template<class APPLYMAT>
void GridCg<APPLYMAT>::setICPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) {
	assertMsg(method==PC_ICP || method==PC_mICP, "GridCg<APPLYMAT>::setICPreconditioner: Invalid method specified.");

	mPcMethod = method;
	if( (!A0->is3D())) {
		if(gPrint2dWarning) {
			debMsg("ICP/mICP pre-conditioning only supported in 3D for now, disabling it.", 1);
			gPrint2dWarning = false;
		}
		mPcMethod=PC_None;
	}
	mpPCA0 = A0;
	mpPCAi = Ai;
	mpPCAj = Aj;
	mpPCAk = Ak;
}

template<class APPLYMAT>
void GridCg<APPLYMAT>::setMGPreconditioner(PreconditionType method, GridMg* MG) {
    assertMsg(method==PC_MGP, "GridCg<APPLYMAT>::setMGPreconditioner: Invalid method specified.");

	mPcMethod = method;
	
	mMG = MG;
}

// explicit instantiation
template class GridCg<ApplyMatrix>;
template class GridCg<ApplyMatrix2D>;

}; // DDF
