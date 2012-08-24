/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for pressure correction:
 * - solve_pressure
 *
 ******************************************************************************/
#include "vectorbase.h"
#include "kernel.h"
#include "conjugategrad.h"

using namespace std;
namespace Manta {

//! Kernel: Construct the right-hand side of the poisson equation
KERNEL(bnd=1, reduce) struct MakeRhs (FlagGrid3& flags, Grid3<Real>& rhs, MACGrid3& vel, 
                                      Grid3<Real>* perCellCorr, Real corr) 
{
    double redSum;
    int redCnt;
    
    void operator()(int i, int j, int k) {
        if (!flags.isFluid(i,j,k)) {
            rhs(i,j,k) = 0;
            return;
        }
            
        Real set;
        set  = vel(i,j,k).x - vel(i+1,j,k).x;
        set += vel(i,j,k).y - vel(i,j+1,k).y;
        set += vel(i,j,k).z - vel(i,j,k+1).z;
        set += corr;
        
        // per cell divergence correction
        if(perCellCorr) 
            set += perCellCorr->get(i,j,k);
        
        // obtain sum, cell count
        redSum += set;
        redCnt++;
        
        rhs(i,j,k) = set;
    }
    
    void setup() {
        redSum = 0;
        redCnt = 0;
    }
    void join(const MakeRhs& other) {
        redSum += other.redSum;
        redCnt += other.redCnt;
    }
};

//! Kernel: Apply velocity update from poisson equation
KERNEL(bnd=1) CorrectVelocity(FlagGrid3& flags, MACGrid3& vel, Grid3<Real>& pressure) 
{
	// correct all faces between fluid-fluid and fluid-empty cells
	// skip everything with obstacles...
	if (flags.isObstacle(i,j,k))
        return;

	// skip faces between two empty cells
	const bool myFlagIsEmpty = flags.isEmpty(i,j,k);

	if(! ((flags.isEmpty(i-1,j,k)) && (myFlagIsEmpty)) )
	{
        vel(i,j,k).x -= (pressure(i,j,k) - pressure(i-1,j,k) );
    }
	if(! ((flags.isEmpty(i,j-1,k)) && (myFlagIsEmpty)) )
	{
        vel(i,j,k).y -= (pressure(i,j,k) - pressure(i,j-1,k) );
    }
	if(! ((flags.isEmpty(i,j,k-1)) && (myFlagIsEmpty)) )
	{
        vel(i,j,k).z -= (pressure(i,j,k) - pressure(i,j,k-1) );
    }    
}

//! Kernel: Set matrix stencils and velocities to enable open boundaries
KERNEL SetOpenBound(Grid3<Real>& A0, Grid3<Real>& Ai, Grid3<Real>& Aj, Grid3<Real>& Ak, MACGrid3& vel,
                    Vector3D<bool> lowerBound, Vector3D<bool> upperBound)
{    
    // set velocity boundary conditions
    if (lowerBound.x && i == 0) vel(0,j,k) = vel(1,j,k);
    if (lowerBound.y && j == 0) vel(i,0,k) = vel(i,1,k);
    if (lowerBound.z && k == 0) vel(i,j,0) = vel(i,j,1);
    if (upperBound.x && i == maxX-1) vel(maxX-1,j,k) = vel(maxX-2,j,k);
    if (upperBound.y && j == maxY-1) vel(i,maxY-1,k) = vel(i,maxY-2,k);
    if (upperBound.z && k == maxZ-1) vel(i,j,maxZ-1) = vel(i,j,maxZ-2);
    
    // set matrix stencils at boundary
    if ((lowerBound.x && i<=1) || (upperBound.x && i>=maxX-2) ||
        (lowerBound.y && j<=1) || (upperBound.y && j>=maxY-2) ||
        (lowerBound.z && k<=1) || (upperBound.z && k>=maxZ-2)) {
        A0(i,j,k) = 6.0;
        Ai(i,j,k) = -1.0;
        Aj(i,j,k) = -1.0;
        Ak(i,j,k) = -1.0;
    }
}

//! Kernel: Set matrix rhs for outflow
KERNEL SetOutflow (Grid3<Real>& rhs, Vector3D<bool> lowerBound, Vector3D<bool> upperBound, int height)
{
    if ((lowerBound.x && i < height) || (upperBound.x && i >= maxX-1-height) ||
        (lowerBound.y && j < height) || (upperBound.y && j >= maxY-1-height) ||
        (lowerBound.z && k < height) || (upperBound.z && k >= maxZ-1-height))
        rhs(i,j,k) = 0;
}


// *****************************************************************************
// Ghost fluid helpers
// TODO set sides individually?

// iso surface level, usually zero
static const int LEVELSET_ISOSURFACE = 0.;

static inline Real getGhostMatrixAddition(Real a, Real b, const Real accuracy) {
	Real ret = 0.f;

	if(a < 0 && b < 0)
		ret = 1;
	else if( (a >= 0) && (b < 0))
		ret = b / (b - a);
	else if ( (a < 0) && (b >= 0))
		ret = a / (a - b);
	else
		ret = 0;

	if(ret < accuracy)
		ret = accuracy;

	Real invret = 1./ret;
	return invret;
}

//! Kernel: Adapt A0 for ghost fluid
KERNEL(bnd=1) ApplyGhostFluid (FlagGrid3& flags, Grid3<Real>& phi, Grid3<Real>& A0, Real accuracy) 
{
    if (!flags.isFluid(i,j,k))
        return;
    
    const Real curPhi = phi(i,j,k);
    
    if (flags.isEmpty(i-1,j,k)) A0(i,j,k) += getGhostMatrixAddition( phi(i-1,j,k), curPhi, accuracy);
    if (flags.isEmpty(i,j-1,k)) A0(i,j,k) += getGhostMatrixAddition( phi(i,j-1,k), curPhi, accuracy);
    if (flags.isEmpty(i,j,k-1)) A0(i,j,k) += getGhostMatrixAddition( phi(i,j,k-1), curPhi, accuracy);
    if (flags.isEmpty(i+1,j,k)) A0(i,j,k) += getGhostMatrixAddition( curPhi, phi(i+1,j,k), accuracy);
    if (flags.isEmpty(i,j+1,k)) A0(i,j,k) += getGhostMatrixAddition( curPhi, phi(i,j+1,k), accuracy);
    if (flags.isEmpty(i,j,k+1)) A0(i,j,k) += getGhostMatrixAddition( curPhi, phi(i,j,k+1), accuracy);
}

//! Kernel: Correct velocities for ghost fluids
KERNEL(bnd=1) CorrectVelGhostFluid (FlagGrid3& flags, MACGrid3& vel, Grid3<Real>& pressure) //, Grid3<Real>& phi)
{
    bool curFluid = flags.isFluid(i,j,k);
    if (!curFluid && !flags.isEmpty(i,j,k))
        return;
    
    const Real curPress = pressure(i,j,k);

    //const Real curPhi = phi(i,j,k);
	// TODO - include ghost fluid factor  NT_DEBUG
    
    // in contrast to old implementation:
    // make sure to add gradient for all fluid-empty or fluid-fluid combinations
    // of neighbors...

    if (!flags.isObstacle(i-1,j,k) && (curFluid || flags.isFluid(i-1,j,k)))
        vel(i,j,k).x -= curPress - pressure(i-1,j,k);
    
    if (!flags.isObstacle(i,j-1,k) && (curFluid || flags.isFluid(i,j-1,k)))
        vel(i,j,k).y -= curPress - pressure(i,j-1,k);
    
    if (!flags.isObstacle(i,j,k-1) && (curFluid || flags.isFluid(i,j,k-1)))
        vel(i,j,k).z -= curPress - pressure(i,j,k-1);
}

inline void convertDescToVec(const string& desc, Vector3D<bool>& lo, Vector3D<bool>& up) {
    for(size_t i=0; i<desc.size(); i++) {
        if (desc[i] == 'x') lo.x = true;
        else if (desc[i] == 'y') lo.y = true;
        else if (desc[i] == 'z') lo.z = true;
        else if (desc[i] == 'X') up.x = true;
        else if (desc[i] == 'Y') up.y = true;
        else if (desc[i] == 'Z') up.z = true;
        else throw Error("invalid character in boundary description string. Only [xyzXYZ] allowed.");
    }
}

//! Perform pressure projection of the velocity grid
PLUGIN void solvePressure(MACGrid3& vel, Grid3<Real>& pressure, FlagGrid3& flags,
                     Grid3<Real>* phi = 0, 
                     Grid3<Real>* perCellCorr = 0, 
                     Real divCorr = 0,
                     Real ghostAccuracy = 0, 
                     Real cgMaxIterFac = 1.5,
                     Real cgAccuracy = 1e-3,
                     string openBound = "",
                     string outflow = "",
                     int outflowHeight = 1,
                     int precondition = 0,
                     bool enforceCompatibility = false,
                     bool useResNorm = true )
{
    // parse strings
    Vector3D<bool> loOpenBound, upOpenBound, loOutflow, upOutflow;
    convertDescToVec(openBound, loOpenBound, upOpenBound);
    convertDescToVec(outflow, loOutflow, upOutflow);
    
    // reserve temp grids
    Grid3<Real> rhs(parent);
    Grid3<Real> residual(parent);
    Grid3<Real> search(parent);
    Grid3<Real> A0(parent);
    Grid3<Real> Ai(parent);
    Grid3<Real> Aj(parent);
    Grid3<Real> Ak(parent);
    Grid3<Real> tmp(parent);
    Grid3<Real> pca0(parent);
    Grid3<Real> pca1(parent);
    Grid3<Real> pca2(parent);
    Grid3<Real> pca3(parent);
        
    // setup matrix and boundaries
    MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak);
    SetOpenBound (A0, Ai, Aj, Ak, vel, loOpenBound, upOpenBound);
    
    if (ghostAccuracy > 0) {
        if (!phi) throw("solve_pressure: if ghostAccuracy>0, need to specify levelset phi=xxx");
        ApplyGhostFluid (flags, A0, *phi, ghostAccuracy);
    }
    
    // compute divergence and init right hand side
    MakeRhs kernMakeRhs (flags, rhs, vel, perCellCorr, divCorr);
    
    if (!outflow.empty())
        SetOutflow (rhs, loOutflow, upOutflow, outflowHeight);
    
    if (enforceCompatibility)
        rhs += (Real)(-kernMakeRhs.redSum / (Real)kernMakeRhs.redCnt);
    
    // CG
    const int maxIter = (int)(cgMaxIterFac * flags.getSize().max());
    GridCgInterface *gcg = new GridCg<ApplyMatrix>(pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
    
    gcg->setAccuracy( cgAccuracy ); 
    gcg->setUseResNorm( useResNorm );

    // optional preconditioning
    gcg->setPreconditioner( (GridCgInterface::PreconditionType)precondition, &pca0, &pca1, &pca2, &pca3);

    for (int iter=0; iter<maxIter; iter++) {
        if (!gcg->iterate()) iter=maxIter;
    } 
    debMsg("FluidSolver::solvePressure iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
    delete gcg;
    
    if(ghostAccuracy<=0.) {
        // ghost fluid off, normal correction
        CorrectVelocity (flags, vel, pressure );
    } else {        
        CorrectVelGhostFluid (flags, vel, pressure);
    }    
}

} // end namespace

