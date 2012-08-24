/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Levelset
 *
 ******************************************************************************/

#include "levelset.h"
#include "fastmarch.h"
#include "kernel.h"
#include "mcubes.h"
#include "mesh.h"

using namespace std;
namespace Manta {

//************************************************************************
// Helper functions and kernels for marching

static const int FlagInited = FastMarch<3,FmHeapComparatorOut, +1>::FlagInited;

// neighbor lookup vectors
static const Vec3i neighbors[6] = { Vec3i(-1,0,0), Vec3i(1,0,0), Vec3i(0,-1,0), Vec3i(0,1,0), Vec3i(0,0,-1), Vec3i(0,0,1) };

KERNEL(bnd=1) template<int DIM>
InitFmIn (FlagGrid<DIM>& flags, Grid<DIM,int>& fmFlags, LevelsetGrid<DIM>& phi, bool ignoreWalls) {
    const int idx = flags.index(i,j,k);
    const Real v = phi[idx];
    if (v>=0 && (!ignoreWalls || !flags.isObstacle(idx)))
        fmFlags[idx] = FlagInited;
    else
        fmFlags[idx] = 0;
}

KERNEL(bnd=1) template<int DIM>
InitFmOut (FlagGrid<DIM>& flags, Grid<DIM,int>& fmFlags, LevelsetGrid<DIM>& phi, bool ignoreWalls) {
    const int idx = flags.index(i,j,k);
    const Real v = phi[idx];
    if (ignoreWalls) {
        fmFlags[idx] = (v<0) ? FlagInited : 0;
        if (flags.isObstacle(idx)) {
            fmFlags[idx] = 0;
            phi[idx] = 0;
        }
    }
    else
        fmFlags[idx] = (v<0) ? FlagInited : 0;
}

KERNEL(bnd=1) template<int DIM>
SetUninitialized (Grid<DIM,int>& fmFlags, LevelsetGrid<DIM>& phi, const Real val) {
    if (fmFlags(i,j,k) != FlagInited)
        phi(i,j,k) = val;
}

template<int DIM, bool inward>
inline bool isAtInterface(Grid<DIM,int>& fmFlags, LevelsetGrid<DIM>& phi, const Vec3i& p) {
    // check for interface
    for (int nb=0; nb<4+2*DIM; nb++) {
        const Vec3i pn(p + neighbors[nb]);
        if (!fmFlags.isInBounds(pn)) continue;
        
        if (fmFlags(pn) != FlagInited) continue;
        if ((inward && phi(pn) >= 0) || 
            (!inward && phi(pn) < 0)) return true;
    }
    return false;
}

//************************************************************************
// Levelset class def

template <int DIM>
LevelsetGrid<DIM>::LevelsetGrid(FluidSolver* parent, bool show) 
    : Grid<DIM,Real>(parent, show) 
{ 
    this->mType = (GridBase::GridType)(this->TypeLevelset | this->TypeReal | ((DIM==2) ? this->Type2D : 0));
}    

//extern void updateQtGui(bool full, int frame); // HACK

template<int DIM>
Real LevelsetGrid<DIM>::invalidTimeValue() {
    return InvalidTime;
}

template<int DIM>
void LevelsetGrid<DIM>::reinitMarching(FlagGrid<DIM>& flags, Real maxTime, MACGrid<DIM>* velTransport, bool ignoreWalls, bool correctOuterLayer)
{
    Grid<DIM,int> fmFlags(this->getParent());
    LevelsetGrid<DIM>& phi = *this;
    
    FastMarch<DIM, FmHeapComparatorOut, +1> marchOut(flags, fmFlags, phi, maxTime, velTransport);
    FastMarch<DIM, FmHeapComparatorIn, -1> marchIn(flags, fmFlags, phi, maxTime, NULL);
    
    // march inside
    InitFmIn<DIM> (flags, fmFlags, phi, ignoreWalls);
    
    FOR_IJK_BND(flags, 1) {
        if (fmFlags(i,j,k) == FlagInited) continue;
        if (flags.isObstacle(i,j,k)) continue;
        const Vec3i p(i,j,k);
                
        if(isAtInterface<DIM,true>(fmFlags, phi, p)) {
            // set value
            fmFlags(p) = FlagInited;
            
            // add neighbors that are not at the interface
            for (int nb=0; nb<6; nb++) {
                const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                if (flags.isObstacle(pn)) continue;
                
                // check neighbors of neighbor
                if (phi(pn) < 0 && !isAtInterface<DIM,true>(fmFlags, phi, pn)) {
                    marchIn.addToList(pn, p);
                }
            }            
        }
    }
    marchIn.performMarching(); 
        
    // set un initialized regions
    SetUninitialized<DIM> (fmFlags, phi, -maxTime); 
    
    // done with inwards marching, now march out...    
    InitFmOut<DIM> (flags, fmFlags, phi, ignoreWalls);
    
    // by default, correctOuterLayer is on
    if (correctOuterLayer) {    
        // normal version, inwards march is done, now add all outside values (0..2] to list
        // note, this might move the interface a bit! but keeps a nice signed distance field...        
        FOR_IJK_BND(flags, 1) {
            if (flags.isObstacle(i,j,k)) continue;
            const Vec3i p(i,j,k);
            
            // check nbs
            for (int nb=0; nb<6; nb++) {
                const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                
                if (fmFlags(pn) != FlagInited) continue;
                if (flags.isObstacle(pn)) continue;
                
                const Real nbPhi = phi(pn);
                
                // only add nodes near interface, not e.g. outer boundary vs. invalid region                
                if (nbPhi < 0 && nbPhi >= -2)
                    marchOut.addToList(p, pn);
            }
        }         
    } else {
        // alternative version, keep interface, do not distort outer cells
        // add all ouside values, but not those at the IF layer
        FOR_IJK_BND(flags, 1) {
            if (flags.isObstacle(i,j,k)) continue;
            
            // only look at ouside values
            const Vec3i p(i,j,k);
            if (phi(p) < 0) continue;
            
            if (isAtInterface<DIM,false>(fmFlags, phi, p)) {
                // now add all non, interface neighbors
                fmFlags(p) = FlagInited;
                
                // add neighbors that are not at the interface
                for (int nb=0; nb<6; nb++) {
                    const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                    if (flags.isObstacle(pn)) continue;
                
                    // check neighbors of neighbor
                    if (phi(pn) > 0 && !isAtInterface<DIM,false>(fmFlags, phi, pn))
                        marchOut.addToList(pn, p);
                }            
            }
        }
    }    
    marchOut.performMarching();
    
    // set un initialized regions
    SetUninitialized<DIM> (fmFlags, phi, +maxTime);    
}

// helper function
inline Vec3 getNormal(const Grid<3,Real>& data, int i, int j, int k) {
    if (i > data.getSizeX()-2) i= data.getSizeX()-2;
    if (j > data.getSizeY()-2) j= data.getSizeY()-2;
    if (k > data.getSizeZ()-2) k= data.getSizeZ()-2;
    return Vec3 (data(i-1,j  ,k  ) - data(i+1,j  ,k  ),
                 data(i  ,j-1,k  ) - data(i  ,j+1,k  ),
                 data(i  ,j  ,k-1) - data(i  ,j  ,k+1));
}

template<> void LevelsetGrid<2>::createMesh(Mesh& mesh) {
    throw Error("createMesh only defined for 3D Levelset grids.");
}

template<> void LevelsetGrid<3>::createMesh(Mesh& mesh) {
    mesh.clear();
        
    const Real invalidTime = InvalidTime;
    const Real isoValue = 1e-4;
    
    // create some temp grids
    Grid3<int> edgeVX(this->getParent());
    Grid3<int> edgeVY(this->getParent());
    Grid3<int> edgeVZ(this->getParent());
    
    FOR_IJK_BND(*this,1) {
         Real value[8] = { get(i,j,k),   get(i+1,j,k),   get(i+1,j+1,k),   get(i,j+1,k),
                           get(i,j,k+1), get(i+1,j,k+1), get(i+1,j+1,k+1), get(i,j+1,k+1) };
        
        // build lookup index, check for invalid times
        bool skip = false;
        int cubeIdx = 0;
        for (int l=0;l<8;l++) {
            value[l] *= -1;
            if (-value[l] <= invalidTime)
                skip = true;
            if (value[l] < isoValue) 
                cubeIdx |= 1<<l;
        }
        if (skip || (mcEdgeTable[cubeIdx] == 0)) continue;
        
        // where to look up if this point already exists
        int triIndices[12];
        int *eVert[12] = { &edgeVX(i,j,k),   &edgeVY(i+1,j,k),   &edgeVX(i,j+1,k),   &edgeVY(i,j,k), 
                           &edgeVX(i,j,k+1), &edgeVY(i+1,j,k+1), &edgeVX(i,j+1,k+1), &edgeVY(i,j,k+1), 
                           &edgeVZ(i,j,k),   &edgeVZ(i+1,j,k),   &edgeVZ(i+1,j+1,k), &edgeVZ(i,j+1,k) };
        
        const Vec3 pos[9] = { Vec3(i,j,k),   Vec3(i+1,j,k),   Vec3(i+1,j+1,k),   Vec3(i,j+1,k),
                        Vec3(i,j,k+1), Vec3(i+1,j,k+1), Vec3(i+1,j+1,k+1), Vec3(i,j+1,k+1) };
        
        for (int e=0; e<12; e++) {
            if (mcEdgeTable[cubeIdx] & (1<<e)) {
                // vertex already calculated ?
                if (*eVert[e] == 0) {
                    // interpolate edge
                    const int e1 = mcEdges[e*2  ];
                    const int e2 = mcEdges[e*2+1];
                    const Vec3 p1 = pos[ e1  ];    // scalar field pos 1
                    const Vec3 p2 = pos[ e2  ];    // scalar field pos 2
                    const float valp1  = value[ e1  ];  // scalar field val 1
                    const float valp2  = value[ e2  ];  // scalar field val 2
                    const float mu = (isoValue - valp1) / (valp2 - valp1);

                    // init isolevel vertex
                    Node vertex;
                    vertex.pos = p1 + (p2-p1)*mu;
                    vertex.normal = getNormalized( 
                                        getNormal( *this, i+cubieOffsetX[e1], j+cubieOffsetY[e1], k+cubieOffsetZ[e1]) * (1.0-mu) +
                                        getNormal( *this, i+cubieOffsetX[e2], j+cubieOffsetY[e2], k+cubieOffsetZ[e2]) * (    mu)) ;
                    
                    triIndices[e] = mesh.addNode(vertex) + 1;
                    
                    // store vertex 
                    *eVert[e] = triIndices[e];
                } else {
                    // retrieve  from vert array
                    triIndices[e] = *eVert[e];
                }
            }
        }
        
        // Create the triangles... 
        for(int e=0; mcTriTable[cubeIdx][e]!=-1; e+=3) {
            mesh.addTri( Triangle( triIndices[ mcTriTable[cubeIdx][e+0]] - 1,
                                        triIndices[ mcTriTable[cubeIdx][e+1]] - 1,
                                        triIndices[ mcTriTable[cubeIdx][e+2]] - 1));
        }
    }
    
    //mesh.rebuildCorners();
    //mesh.rebuildLookup();
}

template class LevelsetGrid<2>;
template class LevelsetGrid<3>;

} //namespace