/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * shapes classes
 *
 ******************************************************************************/

#ifndef _SHAPES_H
#define _SHAPES_H

#include "pclass.h"
#include "vectorbase.h"
#include "levelset.h"

namespace Manta {

// forward declaration
class Mesh;
    
//! Base class for all shapes
PYTHON class Shape : public PbClass {
public:
    enum GridType { TypeNone = 0, TypeBox = 1, TypeSphere = 2, TypeCylinder };
    
    PYTHON Shape(FluidSolver* parent);
    
    //! Get the type of grid
    inline GridType getType() const { return mType; }
    
    //! Apply shape to flag grid, set inside cells to <value>
    PYTHON void applyToGrid(GridBase* grid);
    PYTHON void applyToGridSmooth(GridBase* grid, Real sigma=1.0, Real shift=0);
    PYTHON LevelsetGrid computeLevelset();
    PYTHON void collideMesh(Mesh& mesh);
    
    //! Inside test of the shape
    virtual bool isInside(const Vec3& pos) const;
    inline bool isInsideGrid(int i, int j, int k) const { return isInside(Vec3(i+0.5,j+0.5,k+0.5)); };
    
    virtual void generateMesh(Mesh* mesh) {} ;    
protected:
    virtual void generateLevelset(Grid<Real>& phi) {};    
    
    GridType mType;
};

//! Cylindrical shape
PYTHON class NullShape : public Shape {    
public:
    PYTHON NullShape (FluidSolver* parent) : Shape(parent) {}
    
    virtual bool isInside(const Vec3& pos) const { return false; }
    virtual void generateMesh(Mesh* mesh) {}
       
protected:
    virtual void generateLevelset(Grid<Real>& phi) { phi = 1000.0f; }
};

//! Box shape
PYTHON class Box : public Shape {    
public:
    PYTHON Box(FluidSolver* parent, Vec3 center = Vec3::Invalid, Vec3 p0 = Vec3::Invalid, Vec3 p1 = Vec3::Invalid, Vec3 size = Vec3::Invalid);
    
    inline Vec3 getCenter() const { return 0.5*(mP0+mP1); }
    inline Vec3 getSize() const { return mP1-mP0; }
    inline Vec3 getP0() const { return mP0; }
    inline Vec3 getP1() const { return mP1; }
    virtual bool isInside(const Vec3& pos) const;
    virtual void generateMesh(Mesh* mesh);
        
protected:
    virtual void generateLevelset(Grid<Real>& phi);
    
    Vec3 mP0, mP1;
};

//! Spherical shape
PYTHON class Sphere : public Shape {    
public:
    PYTHON Sphere (FluidSolver* parent, Vec3 center, Real radius, Vec3 scale=Vec3(1,1,1));
    
    inline Vec3 getCenter() const { return mCenter; }
    inline Real getRadius() const { return mRadius; }
    virtual bool isInside(const Vec3& pos) const;
    virtual void generateMesh(Mesh* mesh);
       
protected:
    virtual void generateLevelset(Grid<Real>& phi);
    
    Vec3 mCenter, mScale;
    Real mRadius;
};

//! Cylindrical shape
PYTHON class Cylinder : public Shape {    
public:
    PYTHON Cylinder (FluidSolver* parent, Vec3 center, Real radius, Vec3 z);
    
    PYTHON void setCenter(Vec3 center) { mCenter=center; }
    PYTHON void setRadius(Real r) { mRadius = r; }
    PYTHON void setZ(Vec3 z) { mZDir=z; mZ=normalize(mZDir); }
    
    inline Vec3 getCenter() const { return mCenter; }
    inline Real getRadius() const { return mRadius; }
    inline Vec3 getZ() const { return mZ*mZDir; }
    virtual bool isInside(const Vec3& pos) const;
    virtual void generateMesh(Mesh* mesh);
        
protected:
    virtual void generateLevelset(Grid<Real>& phi);
    
    Vec3 mCenter, mZDir;
    Real mRadius, mZ;
};

    

} //namespace
#endif