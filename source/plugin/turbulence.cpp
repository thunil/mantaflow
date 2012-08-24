/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for using vortex sheet meshes 
 *
 ******************************************************************************/
 
#include "grid.h"
#include "commonkernels.h"
#include "vortexsheet.h"

using namespace std;

namespace Manta {

// k-epsilon model constants
const Real keCmu = 0.09;
const Real keC1 = 1.44;
const Real keC2 = 1.92;
const Real keS1 = 1.0;
const Real keS2 = 1.3;

// k-epsilon limiters
const Real keU0 = 1.0;
const Real keImin = 2e-3;
const Real keImax = 1.0;
const Real keNuMin = 1e-3;
const Real keNuMax = 4.0;

//! clamp k and epsilon to limits    
KERNEL(idx) KnTurbulenceClamp(Grid3<Real>& kgrid, Grid3<Real>& egrid, Real minK, Real maxK, Real minNu, Real maxNu) {
    Real eps = egrid[idx];
    Real ke = clamp(kgrid[idx],minK,maxK);
    Real nu = keCmu*square(ke)/eps;
    if (nu > maxNu) 
        eps = keCmu*square(ke)/maxNu;
    if (nu < minNu) 
        eps = keCmu*square(ke)/minNu;

    kgrid[idx] = ke;
    egrid[idx] = eps;
}

//! clamp k and epsilon to limits    
void TurbulenceClampMesh(VortexSheetMesh& mesh, Real minK, Real maxK, Real minNu, Real maxNu) {
    for (int idx=0; idx<mesh.numNodes(); idx++) {
        Real eps = mesh.turb(idx).epsilon;
        Real ke = clamp(mesh.turb(idx).k,minK,maxK);
        Real nu = keCmu*square(ke)/eps;
        if (nu > maxNu) 
            eps = keCmu*square(ke)/maxNu;
        if (nu < minNu) 
            eps = keCmu*square(ke)/minNu;

        mesh.turb(idx).k = ke;
        mesh.turb(idx).epsilon = eps;
    }
}

//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Sij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps
KERNEL(bnd=1) KnComputeProductionStrain(const MACGrid3& vel, const Grid3<Vec3>& velCenter, const Grid3<Real>& ke, const Grid3<Real>& eps, 
                                  Grid3<Real>& prod, Grid3<Real>& nuT, Real pscale = 1.0f) 
{
    Real curEps = eps(i,j,k);
    if (curEps > 0) {
        // turbulent viscosity: nu_T = C_mu * k^2/eps
        Real curNu = keCmu * square(ke(i,j,k)) / curEps;
        
        // compute Sij = 1/2 * (dU_i/dx_j + dU_j/dx_i)
        Vec3 diag = Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, vel(i,j,k+1).z) - vel(i,j,k);
        Vec3 ux = 0.5*(velCenter(i+1,j,k)-velCenter(i-1,j,k));
        Vec3 uy = 0.5*(velCenter(i,j+1,k)-velCenter(i,j-1,k));
        Vec3 uz = 0.5*(velCenter(i,j,k+1)-velCenter(i,j,k-1));
        Real S12 = 0.5*(ux.y+uy.x);
        Real S13 = 0.5*(ux.z+uz.x);
        Real S23 = 0.5*(uy.z+uz.y);
        Real S2 = square(diag.x) + square(diag.y) + square(diag.z) +
                  2.0*square(S12) + 2.0*square(S13) + 2.0*square(S23);
        
        // P = 2*nu_T*sum_ij(Sij^2)
        prod(i,j,k) = 2.0 * curNu * S2 * pscale;
        nuT(i,j,k) = curNu;        
    } 
    else {
        prod(i,j,k) = 0;
        nuT(i,j,k) = 0;
    }
}

//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Omegaij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps
KERNEL(bnd=1) KnComputeProductionCurl(const Grid3<Vec3>& curl, const Grid3<Vec3>& bcurl, const Grid3<Real>& ke, const Grid3<Real>& eps, 
                                  Grid3<Real>& prod, Grid3<Real>& nuT, Real pscale =1.0f, Grid3<Vec3>* debug=NULL) 
{
    Real curEps = eps(i,j,k);
    if (curEps > 0) {
        // turbulent viscosity: nu_T = C_mu * k^2/eps
        Real curNu = keCmu * square(ke(i,j,k)) / curEps;
        
        // diff with clamping
        Vec3 actual = curl(i,j,k), buoyant=bcurl(i,j,k);
        Vec3 diff = actual - buoyant;
        for (int c=0; c<3; c++) {
            if (actual[c]*diff[c] < 0) 
                diff[c] = 0; // clamp to 0 if buoyant is bigger than actual
            if (actual[c]*buoyant[c] < 0) 
                diff[c] = actual[c]; // ignore if initial dirs point in different direction            
            if (fabs(diff[c]) > fabs(actual[c]))
                diff[c] = actual[c];
        }
        if (debug)
            (*debug)(i,j,k) = diff;
        
        prod(i,j,k) = 2.0 * curNu * 2.0 * normSquare(diff) * pscale;
        nuT(i,j,k) = curNu;                
    } 
    else {
        prod(i,j,k) = 0;
        nuT(i,j,k) = 0;
    }
}
    
//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Sij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps
PLUGIN void KEpsilonComputeProduction(MACGrid3& vel, Grid3<Real>& k, Grid3<Real>& eps, Grid3<Real>& prod, Grid3<Real>& nuT, Real pscale = 1.0f, Grid3<Vec3>* bcurl=NULL, Grid3<Vec3>* debug=NULL) 
{
    // get centered velocity grid
    Grid3<Vec3> vcenter(parent);
    GetCentered(vcenter, vel);
    
    // compute limits
    const Real minK = 1.5*square(keU0)*square(keImin);
    const Real maxK = 1.5*square(keU0)*square(keImax);    
    KnTurbulenceClamp(k, eps, minK, maxK, keNuMin, keNuMax);
    
    if (bcurl) {
        Grid3<Vec3> curl(parent);
        CurlOp(vcenter, curl);
        
        // compute production field
        KnComputeProductionCurl(curl,*bcurl, k, eps, prod, nuT, pscale, debug);
    } else {
        KnComputeProductionStrain(vel, vcenter, k, eps, prod, nuT, pscale);
    }
}

//! Integrate source terms of k-epsilon equation
KERNEL(idx) KnAddTurbulenceSource(Grid3<Real>& kgrid, Grid3<Real>& egrid, const Grid3<Real>& pgrid, Real dt) {
    Real eps = egrid[idx], prod = pgrid[idx], ke = kgrid[idx];
    if (ke <= 0) ke = 1e-3; // pre-clamp to avoid nan
    
    Real newK = ke + dt*(prod - eps);
    Real newEps = eps + dt*(prod * keC1 - eps * keC2) * (eps / ke);
    if (newEps <= 0) newEps = 1e-4; // pre-clamp to avoid nan

    kgrid[idx] = newK;
    egrid[idx] = newEps;
}


//! Integrate source terms of k-epsilon equation
PLUGIN void KEpsilonSources(Grid3<Real>& k, Grid3<Real>& eps, Grid3<Real>& prod) {
    Real dt = parent->getDt();
        
    KnAddTurbulenceSource(k, eps, prod, dt);
    
    // compute limits
    const Real minK = 1.5*square(keU0)*square(keImin);
    const Real maxK = 1.5*square(keU0)*square(keImax);
    KnTurbulenceClamp(k, eps, minK, maxK, keNuMin, keNuMax);    
}

//! Integrate source terms of k-epsilon equation
void AddTurbulenceSourceMesh(VortexSheetMesh& mesh, Grid3<Vec3>& curl, Grid3<Vec3>& bcurl, Real dt) {
    for (int idx=0; idx<mesh.numNodes(); idx++) {
        Real eps = mesh.turb(idx).epsilon;
        Real k = mesh.turb(idx).k;
        const Vec3& pos = mesh.nodes(idx).pos;
        
        // turbulent viscosity: nu_T = C_mu * k^2/eps
        Real curNu = keCmu * square(k) / eps;
        
        // diff with clamping
        Vec3 actual = curl.getInterpolated(pos), buoyant=bcurl.getInterpolated(pos);
        Vec3 diff = actual - buoyant;
        for (int c=0; c<3; c++) {
            if (actual[c]*diff[c] < 0) 
                diff[c] = 0; // clamp to 0 if buoyant is bigger than actual
            if (actual[c]*buoyant[c] < 0) 
                diff[c] = actual[c]; // ignore if initial dirs point in different direction            
            if (fabs(diff[c]) > fabs(actual[c]))
                diff[c] = actual[c];
        }
        
        // production
        Real prod = 2.0 * curNu * 2.0 * normSquare(diff);
        Real newK = k + dt*(prod - eps);
        Real newEps = eps + dt*(prod * keC1 - eps * keC2) * (eps / k);
        
        mesh.turb(idx).k = newK;
        mesh.turb(idx).epsilon = newEps;
    }    
}

//! Integrate source terms of k-epsilon equation
PLUGIN void KEpsilonSourcesMesh(VortexSheetMesh& mesh, MACGrid3& vel, Grid3<Vec3>& bcurl, Grid3<Real>& kgrid, Grid3<Real>& epsGrid) {
    Real dt = parent->getDt();
    
    // get centered velocity grid
    Grid3<Vec3> vcenter(parent);
    GetCentered(vcenter, vel);
    
    // compute limits
    const Real minK = 1.5*square(keU0)*square(keImin);
    const Real maxK = 1.5*square(keU0)*square(keImax);    
    
    Grid3<Vec3> curl(parent);
    CurlOp(vcenter, curl);
    
    TurbulenceClampMesh(mesh, minK, maxK, keNuMin, keNuMax);
    AddTurbulenceSourceMesh(mesh, curl, bcurl, dt);
    TurbulenceClampMesh(mesh, minK, maxK, keNuMin, keNuMax);
    
    Grid3<Real> sum(parent);
    kgrid.clear();
    epsGrid.clear();
    for (int i=0; i<mesh.numNodes(); i++) {
        const Vec3& p = mesh.nodes(i).pos;
        kgrid.setInterpolated(p, mesh.turb(i).k, sum);
        epsGrid.setInterpolated(p, mesh.turb(i).epsilon, sum);
    }
    sum *= 0.5;
    kgrid.safeDivide(sum);
    epsGrid.safeDivide(sum);
}

PLUGIN void KEpsilonInit(Grid3<Real>& k, Grid3<Real>& eps, Real intensity, Real nu) {
    // compute limits
    const Real vk = 1.5*square(keU0)*square(intensity);
    const Real ve = keCmu*square(vk) / nu;
    
    FOR_IDX(k) {
        k[idx] = vk;
        eps[idx] = ve;
    }
}

PLUGIN void KEpsilonInitMesh(VortexSheetMesh& mesh, Real intensity, Real nu) {
    // compute limits
    const Real vk = 1.5*square(keU0)*square(intensity);
    const Real ve = keCmu*square(vk) / nu;
    
    for(int i=0; i<mesh.numNodes(); i++) {
        mesh.turb(i).k = vk;
        mesh.turb(i).epsilon = ve;
    }
}

//! Compute k-epsilon turbulent viscosity
PLUGIN void KEpsilonGradientDiffusion(Grid3<Real>& k, Grid3<Real>& eps, Grid3<Real>& nuT) {
    Real dt = parent->getDt();
    MACGrid3 grad(parent);
    Grid3<Real> div(parent);
    
    // gradient diffusion of k
    GradientOpMAC(grad, k);
    grad *= nuT;
    DivergenceOpMAC(div, grad);
    div *= dt/keC1;
    k += div;

    // gradient diffusion of epsilon
    GradientOpMAC(grad, eps);
    grad *= nuT;
    DivergenceOpMAC(div, grad);
    div *= dt/keC2;
    eps += div;
}

} // namespace