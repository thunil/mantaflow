/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Shape classes
 *
 ******************************************************************************/

#include "shapes.h"
#include "commonkernels.h"
#include "mesh.h"

using namespace std;
namespace Manta {

//******************************************************************************
// Shape class members

Shape::Shape (FluidSolver* parent) 
    : PbClass(parent), mType(TypeNone)
{
}

bool Shape::isInside(const Vec3& pos) const {
    return false;
}

//! Kernel: Apply a shape to a grid, setting value inside
KERNEL template<class T> ApplyShapeToGrid (Grid<T>* grid, Shape* shape, T value) {
    if (shape->isInsideGrid(i,j,k))
        (*grid)(i,j,k) = value;
}

//! Kernel: Apply a shape to a grid, setting value inside (scaling by SDF value)
KERNEL template<class T> ApplyShapeToGridSmooth (Grid<T>* grid, Grid<Real>& phi, Real sigma, Real shift, T value) {
    const Real p = phi(i,j,k) - shift;
    if (p < -sigma)
        (*grid)(i,j,k) = value;
    else if (p < sigma)
        (*grid)(i,j,k) = value*(0.5f*(1.0f-p/sigma));
}

//! Kernel: Apply a shape to a MAC grid, setting value inside
KERNEL ApplyShapeToMACGrid (MACGrid* grid, Shape* shape, Vec3 value) 
{
    if (shape->isInside(Vec3(i,j+0.5,k+0.5))) (*grid)(i,j,k).x = value.x;
    if (shape->isInside(Vec3(i+0.5,j,k+0.5))) (*grid)(i,j,k).y = value.y;
    if (shape->isInside(Vec3(i+0.5,j+0.5,k))) (*grid)(i,j,k).z = value.z;
}

void Shape::applyToGrid(GridBase* grid) {
    if (grid->getType() & GridBase::TypeInt)
        ApplyShapeToGrid<int> ((Grid<int>*)grid, this, _args.get<int>("value"));
    else if (grid->getType() & GridBase::TypeReal)
        ApplyShapeToGrid<Real> ((Grid<Real>*)grid, this, _args.get<Real>("value"));
    else if (grid->getType() & GridBase::TypeMAC)
        ApplyShapeToMACGrid ((MACGrid*)grid, this, _args.get<Vec3>("value"));
    else if (grid->getType() & GridBase::TypeVec3)
        ApplyShapeToGrid<Vec3> ((Grid<Vec3>*)grid, this, _args.get<Vec3>("value"));
    else
        errMsg("Shape::applyToGrid(): unknown grid type");
}

void Shape::applyToGridSmooth(GridBase* grid, Real sigma, Real shift) {
    Grid<Real> phi(grid->getParent());
    generateLevelset(phi);

    if (grid->getType() & GridBase::TypeInt)
        ApplyShapeToGridSmooth<int> ((Grid<int>*)grid, phi, sigma, shift, _args.get<int>("value"));
    else if (grid->getType() & GridBase::TypeReal)
        ApplyShapeToGridSmooth<Real> ((Grid<Real>*)grid, phi, sigma, shift, _args.get<Real>("value"));
    else if (grid->getType() & GridBase::TypeVec3)
        ApplyShapeToGridSmooth<Vec3> ((Grid<Vec3>*)grid, phi, sigma, shift, _args.get<Vec3>("value"));
    else
        errMsg("Shape::applyToGridSmooth(): unknown grid type");
}

void Shape::collideMesh(Mesh& mesh) {
    const Real margin = 0.2;
    
    Grid<Real> phi(getParent());
    Grid<Vec3> grad(getParent());
    generateLevelset(phi);
    GradientOp(grad, phi);
    
    const int num=mesh.numNodes();
    for(int i=0; i<num; i++) {
        const Vec3& p = mesh.nodes(i).pos;
        mesh.nodes(i).flags &= ~(Mesh::NfCollide | Mesh::NfMarked);
        if (!phi.isInBounds(p,1)) continue;        
        
        for (int iter=0; iter<10; iter++) {
            const Real dist= phi.getInterpolated(p);
            if (dist<margin) {
                Vec3 n = grad.getInterpolated(p);
                normalize(n);
                mesh.nodes(i).pos += (margin-dist) * n;
                mesh.nodes(i).flags |= Mesh::NfCollide | Mesh::NfMarked;
            }
            else break;
        }
    }
}

//******************************************************************************
// Derived shape class members

Box::Box(FluidSolver* parent, Vec3 center, Vec3 p0, Vec3 p1, Vec3 size)
    : Shape(parent)
{
    mType = TypeBox;
    if (center.isValid() && size.isValid()) {
        mP0 = center - size;
        mP1 = center + size;
    } else if (p0.isValid() && p1.isValid()) {
        mP0 = p0;
        mP1 = p1;
    } else 
		errMsg("Box: specify either p0,p1 or size,center");
	
}

bool Box::isInside(const Vec3& pos) const {
    return (pos.x >= mP0.x && pos.y >= mP0.y && pos.z >= mP0.z &&
            pos.x <= mP1.x && pos.y <= mP1.y && pos.z <= mP1.z);
}

void Box::generateMesh(Mesh* mesh) {
    const int quadidx[24] = { 0,4,6,2, 3,7,5,1, 0,1,5,4, 6,7,3,2, 0,2,3,1, 5,7,6,4 };
    const int nodebase = mesh->numNodes();
    int oldtri = mesh->numTris();
    for (int i=0; i<8; i++) {
        Node p; 
        p.flags = 0;
        p.pos = mP0;
        if (i&1) p.pos.x=mP1.x;
        if (i&2) p.pos.y=mP1.y;
        if (i&4) p.pos.z=mP1.z;
        mesh->addNode(p);
    }
    for (int i=0; i<6; i++) {
        mesh->addTri(Triangle(nodebase + quadidx[i*4+0], nodebase + quadidx[i*4+1], nodebase + quadidx[i*4+3]));
        mesh->addTri(Triangle(nodebase + quadidx[i*4+1], nodebase + quadidx[i*4+2], nodebase + quadidx[i*4+3]));        
    }
    mesh->rebuildCorners(oldtri,-1);
    mesh->rebuildLookup(oldtri,-1);
}

//! Kernel: Analytic SDF for box shape
KERNEL BoxSDF(Grid<Real>& phi, const Vec3& p1, const Vec3& p2) {
    const Vec3 p(i+0.5, j+0.5, k+0.5);
    if (p.x <= p2.x && p.x >= p1.x && p.y <= p2.y && p.y >= p1.y && p.z <= p2.z && p.z >= p1.z) {
        // inside: minimal surface distance
        Real mx = max(p.x-p2.x, p1.x-p.x);
        Real my = max(p.y-p2.y, p1.y-p.y);
        Real mz = max(p.z-p2.z, p1.z-p.z);
        phi(i,j,k) = max(mx,max(my,mz));
    } else if (p.y <= p2.y && p.y >= p1.y && p.z <= p2.z && p.z >= p1.z) {
        // outside plane X
        phi(i,j,k) = max(p.x-p2.x, p1.x-p.x);
    } else if (p.x <= p2.x && p.x >= p1.x && p.z <= p2.z && p.z >= p1.z) {
        // outside plane Y
        phi(i,j,k) = max(p.y-p2.y, p1.y-p.y);
    } else if (p.x <= p2.x && p.x >= p1.x && p.y <= p2.y && p.y >= p1.y) {
        // outside plane Z
        phi(i,j,k) = max(p.z-p2.z, p1.z-p.z);
    } else if (p.x > p1.x && p.x < p2.x) {
        // lines X
        Real m1 = sqrt(square(p1.y-p.y)+square(p1.z-p.z));
        Real m2 = sqrt(square(p2.y-p.y)+square(p1.z-p.z));
        Real m3 = sqrt(square(p1.y-p.y)+square(p2.z-p.z));
        Real m4 = sqrt(square(p2.y-p.y)+square(p2.z-p.z));
        phi(i,j,k) = min(m1,min(m2,min(m3,m4)));
    } else if (p.y > p1.y && p.y < p2.y) {
        // lines Y
        Real m1 = sqrt(square(p1.x-p.x)+square(p1.z-p.z));
        Real m2 = sqrt(square(p2.x-p.x)+square(p1.z-p.z));
        Real m3 = sqrt(square(p1.x-p.x)+square(p2.z-p.z));
        Real m4 = sqrt(square(p2.x-p.x)+square(p2.z-p.z));
        phi(i,j,k) = min(m1,min(m2,min(m3,m4)));
    } else if (p.z > p1.x && p.z < p2.z) {
        // lines Z
        Real m1 = sqrt(square(p1.y-p.y)+square(p1.x-p.x));
        Real m2 = sqrt(square(p2.y-p.y)+square(p1.x-p.x));
        Real m3 = sqrt(square(p1.y-p.y)+square(p2.x-p.x));
        Real m4 = sqrt(square(p2.y-p.y)+square(p2.x-p.x));
        phi(i,j,k) = min(m1,min(m2,min(m3,m4)));
    } else {
        // points
        Real m =   norm(p-Vec3(p1.x,p1.y,p1.z));
        m = min(m, norm(p-Vec3(p1.x,p1.y,p2.z)));
        m = min(m, norm(p-Vec3(p1.x,p2.y,p1.z)));
        m = min(m, norm(p-Vec3(p1.x,p2.y,p2.z)));
        m = min(m, norm(p-Vec3(p2.x,p1.y,p1.z)));
        m = min(m, norm(p-Vec3(p2.x,p1.y,p2.z)));
        m = min(m, norm(p-Vec3(p2.x,p2.y,p1.z)));
        m = min(m, norm(p-Vec3(p2.x,p2.y,p2.z)));
        phi(i,j,k) = m;
    }
}
void Box::generateLevelset(Grid<Real>& phi) {
    BoxSDF(phi, mP0, mP1);
}

Sphere::Sphere (FluidSolver* parent, Vec3 center, Real radius, Vec3 scale) 
    : Shape(parent), mCenter(center), mRadius(radius), mScale(scale)
{
    mType = TypeSphere;
}

bool Sphere::isInside(const Vec3& pos) const {
    return normSquare((pos - mCenter) / mScale) <= mRadius * mRadius;
}

struct Tri { Vec3 t[3]; int i[3]; Tri(Vec3 a,Vec3 b, Vec3 c) {t[0]=a;t[1]=b;t[2]=c;}};        
void Sphere::generateMesh(Mesh* mesh) {
    vector<Tri> tris;
    const int iterations = 3;
    int oldtri = mesh->numTris();
    
    // start with octahedron
    const Real d = sqrt(0.5);
    Vec3 p[6] = {Vec3(0,1,0), Vec3(0,-1,0), Vec3(-d,0,-d), Vec3(d,0,-d), Vec3(d,0,d), Vec3(-d,0,d)};
    tris.push_back(Tri(p[0],p[4],p[3]));
    tris.push_back(Tri(p[0],p[5],p[4]));
    tris.push_back(Tri(p[0],p[2],p[5]));
    tris.push_back(Tri(p[0],p[3],p[2]));
    tris.push_back(Tri(p[1],p[3],p[4]));
    tris.push_back(Tri(p[1],p[4],p[5]));
    tris.push_back(Tri(p[1],p[5],p[2]));
    tris.push_back(Tri(p[1],p[2],p[3]));
    
    // Bisect each edge and move to the surface of a unit sphere
    for (int it=0; it<iterations; it++) {
        int ntold = tris.size();
        for (int i=0; i<ntold; i++) {
            Vec3 pa = 0.5 * (tris[i].t[0] + tris[i].t[1]);
            Vec3 pb = 0.5 * (tris[i].t[1] + tris[i].t[2]);
            Vec3 pc = 0.5 * (tris[i].t[2] + tris[i].t[0]);
            normalize(pa); normalize(pb); normalize(pc);
            
            tris.push_back(Tri(tris[i].t[0], pa, pc));
            tris.push_back(Tri(pa, tris[i].t[1], pb));
            tris.push_back(Tri(pb, tris[i].t[2], pc));
            tris[i].t[0] = pa;
            tris[i].t[1] = pb;
            tris[i].t[2] = pc;         
        }
    }
    
    // index + scale
    vector<Vec3> nodes;
    for (size_t i=0; i<tris.size(); i++) {
        for (int t=0; t<3; t++) {
            Vec3 p = mCenter + tris[i].t[t] * mRadius * mScale;
            // vector already there ?
            int idx=nodes.size();
            for (size_t j=0; j<nodes.size(); j++) {                
                if (p==nodes[j]) {
                    idx = j; break;
                }
            }
            if (idx == (int)nodes.size())
                nodes.push_back(p);
            tris[i].i[t] = idx;
        }
    }
   
    // add the to mesh
    const int ni = mesh->numNodes();
    for (size_t i=0; i<nodes.size(); i++) {
        mesh->addNode(Node(nodes[i]));}
    for (size_t t=0; t<tris.size(); t++)
        mesh->addTri(Triangle(tris[t].i[0]+ni, tris[t].i[1]+ni, tris[t].i[2]+ni));

    mesh->rebuildCorners(oldtri,-1);
    mesh->rebuildLookup(oldtri,-1);
}
    
KERNEL SphereSDF(Grid<Real>& phi, Vec3 center, Real radius, Vec3 scale) {
    phi(i,j,k) = norm((Vec3(i+0.5,j+0.5,k+0.5)-center)/scale)-radius;
}
void Sphere::generateLevelset(Grid<Real>& phi) {
    SphereSDF(phi, mCenter, mRadius, mScale);
} 

Cylinder::Cylinder(FluidSolver* parent, Vec3 center, Real radius, Vec3 z)
    : Shape(parent), mCenter(center), mRadius(radius)
{
    mType = TypeCylinder;
    mZDir = z;
    mZ = normalize(mZDir);
}

bool Cylinder::isInside(const Vec3& pos) const {
    Real z = dot(pos-mCenter, mZDir);
    if (fabs(z) > mZ) return false;
    Real r2 = normSquare(pos-mCenter)-z*z;
    return r2 < (mRadius*mRadius);
}

void Cylinder::generateMesh(Mesh* mesh) {
    // generate coordinate system
    Vec3 x = getOrthogonalVector(mZDir)*mRadius;
    Vec3 y = cross(x, mZDir);
    Vec3 z = mZDir*mZ;
    int oldtri = mesh->numTris();
    
    // construct node ring
    const int N = 20;
    const int base = mesh->numNodes();
    for (int i=0;i<N;i++) {
        const Real phi = 2.0*M_PI*(Real)i/(Real)N;
        Vec3 r = x*cos(phi) + y*sin(phi) + mCenter;
        mesh->addNode(Node(r+z));
        mesh->addNode(Node(r-z));
    }
    // top/bottom center
    mesh->addNode(Node(mCenter+z));
    mesh->addNode(Node(mCenter-z));
        
    // connect with tris
    for (int i=0;i<N;i++) {
        int cur = base+2*i;
        int next = base+2*((i+1)%N);
        // outside
        mesh->addTri(Triangle(cur, next, cur+1));
        mesh->addTri(Triangle(next, next+1, cur+1));        
        // upper / lower
        mesh->addTri(Triangle(cur,base+2*N,next));
        mesh->addTri(Triangle(cur+1,next+1,base+2*N+1));
    }
    
    mesh->rebuildCorners(oldtri, -1);
    mesh->rebuildLookup(oldtri,-1);
}
    
KERNEL CylinderSDF(Grid<Real>& phi, Vec3 center, Real radius, Vec3 zaxis, Real maxz) {
    Vec3 p=Vec3(i+0.5,j+0.5,k+0.5)-center;
    Real z = fabs(dot(p, zaxis));
    Real r = sqrt(normSquare(p)-z*z);
    if (z < maxz) {        
        // cylinder z area
        if (r < radius) 
            phi(i,j,k) = max(r-radius,z-maxz);
        else
            phi(i,j,k) = r-radius;        
    } else if (r < radius) {
        // cylinder top area
        phi(i,j,k) = fabs(z-maxz);
    } else {
        // edge
        phi(i,j,k) = sqrt(square(z-maxz)+square(r-radius));
    }
}
void Cylinder::generateLevelset(Grid<Real>& phi) {
    CylinderSDF(phi, mCenter, mRadius, mZDir, mZ);
}

} //namespace