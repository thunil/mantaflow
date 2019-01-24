/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Meshes
 * 
 *  note: this is only a temporary solution, details are bound to change
 *        long term goal is integration with Split&Merge code by Wojtan et al.
 *
 ******************************************************************************/

#ifndef _MESH_H
#define _MESH_H

#include <vector>
#include "manta.h"
#include "vectorbase.h"
#include <set>
#include "levelset.h"

namespace Manta {

// fwd decl
class GridBase;
//class LevelsetGrid;
class FlagGrid;
class MACGrid;
class Shape;
class MeshDataBase;
template<class T> class MeshDataImpl;


//! Node position and flags
struct Node {
    Node() : flags(0), pos(Vec3::Zero), normal(Vec3::Zero) {}
    Node(const Vec3& p) : flags(0), pos(p) {}
    int flags;
    Vec3 pos, normal;
};

//! Carries indices of its nodes
struct Triangle {    
    Triangle() : flags(0) { c[0] = c[1] = c[2] = 0; }
    Triangle(int n0, int n1, int n2) : flags(0) { c[0]=n0; c[1]=n1; c[2]=n2; }
    
    int c[3];
    int flags;
};

//! For fast access to nodes and neighboring triangles
struct Corner {
    Corner() : tri(-1), node(-1), opposite(-1), next(-1), prev(-1) {};
    Corner(int t, int n) : tri(t), node(n), opposite(-1), next(-1), prev(-1) {}    
    
    int tri;
    int node;
    int opposite;
    int next;
    int prev;
};

//! Base class for mesh data channels (texture coords, vorticity, ...)
struct NodeChannel {
    virtual ~NodeChannel() {};
    virtual void resize(int num) = 0;
    virtual int size() = 0;
    virtual NodeChannel* clone() = 0;
    
    virtual void addInterpol(int a, int b, Real alpha) = 0;
    virtual void mergeWith(int node, int delnode, Real alpha) = 0;
    virtual void renumber(const std::vector<int>& newIndex, int newsize) = 0;
};

//! Node channel using only a vector
template<class T>
struct SimpleNodeChannel : public NodeChannel {
    SimpleNodeChannel() {};
    SimpleNodeChannel(const SimpleNodeChannel<T>& a) : data(a.data) {}
    void resize(int num) { data.resize(num); }
    virtual int size() { return data.size(); }
    virtual void renumber(const std::vector<int>& newIndex, int newsize);
    
    //virtual void addSplit(int from, Real alpha) { data.push_back(data[from]); }    

    std::vector<T> data;
};

//! Base class for mesh data channels (texture coords, vorticity, ...)
struct TriChannel {
    virtual ~TriChannel() {};
    virtual void resize(int num) = 0;
    virtual TriChannel* clone() = 0;
    virtual int size() = 0;
    
    virtual void addNew() = 0;
    virtual void addSplit(int from, Real alpha) = 0;
    virtual void remove(int tri) = 0;
};

//! Tri channel using only a vector
template<class T>
struct SimpleTriChannel : public TriChannel {
    SimpleTriChannel() {};
    SimpleTriChannel(const SimpleTriChannel<T>& a) : data(a.data) {}
    void resize(int num) { data.resize(num); }
    void remove(int tri) { if (tri!=(int)data.size()-1) data[tri] = *data.rbegin(); data.pop_back(); }
    virtual int size() { return data.size(); }
        
    virtual void addSplit(int from, Real alpha) { data.push_back(data[from]); }    
    virtual void addNew() { data.push_back(T()); }
    
    std::vector<T> data;
};

struct OneRing {
    OneRing() {}
    std::set<int> nodes;
    std::set<int> tris;
};

//! Triangle mesh class
/*! note: this is only a temporary solution, details are bound to change
          long term goal is integration with Split&Merge code by Wojtan et al.*/
PYTHON() class Mesh : public PbClass {
public:
    PYTHON() Mesh(FluidSolver* parent);
    virtual ~Mesh();
    virtual Mesh* clone();
    
    enum NodeFlags { NfNone = 0, NfFixed = 1, NfMarked = 2, NfKillme = 4, NfCollide = 8 };
    enum FaceFlags { FfNone = 0, FfDoubled = 1, FfMarked = 2 };
    enum MeshType { TypeNormal = 0, TypeVortexSheet };
    
    virtual MeshType getType() { return TypeNormal; }
        
    Real computeCenterOfMass(Vec3& cm) const;
    void computeVertexNormals();
    
    // plugins
    PYTHON() void clear();
    PYTHON() void load (std::string name, bool append = false);
    PYTHON() void fromShape (Shape& shape, bool append = false);
    PYTHON() void save (std::string name);
    PYTHON() void advectInGrid(FlagGrid& flags, MACGrid& vel, int integrationMode);
    PYTHON() void scale(Vec3 s);
    PYTHON() void offset(Vec3 o);
	PYTHON() void rotate(Vec3 thetas);
    PYTHON() void computeVelocity(Mesh& oldMesh, MACGrid& vel);

	PYTHON() void computeLevelset(LevelsetGrid& levelset, Real sigma, Real cutoff=-1.);
	PYTHON() LevelsetGrid getLevelset(Real sigma, Real cutoff = -1.);

	//! map mesh to grid with sdf
	PYTHON() void applyMeshToGrid(GridBase* grid, FlagGrid* respectFlags=0, Real cutoff=-1.);

	//! get data pointer of nodes
	PYTHON() std::string getNodesDataPointer();
	//! get data pointer of tris
	PYTHON() std::string getTrisDataPointer();

    // ops
    Mesh& operator=(const Mesh& o);
    
    // accessors
    inline int numTris() const { return mTris.size(); }
    inline int numNodes() const { return mNodes.size(); }
    inline int numTriChannels() const { return mTriChannels.size(); }
    inline int numNodeChannels() const { return mNodeChannels.size(); }

	//! return size of container
	//! note , python binding disabled for now! cannot yet deal with long-long types
	inline IndexInt size() const { return mNodes.size(); }
	//! slow virtual function of base class, also returns size
	virtual IndexInt getSizeSlow() const { return size(); }
    
    inline Triangle& tris(int i) { return mTris[i]; }
    inline Node& nodes(int i) { return mNodes[i]; }    
    inline Corner& corners(int tri, int c) { return mCorners[tri*3+c]; }
    inline Corner& corners(int c) { return mCorners[c]; }
    inline NodeChannel* nodeChannel(int i) { return mNodeChannels[i]; }
    inline TriChannel* triChannel(int i) { return mTriChannels[i]; }

	// allocate memory (eg upon load)
	void resizeTris(int numTris);
	void resizeNodes(int numNodes);
    
    inline bool isNodeFixed(int n) { return mNodes[n].flags & NfFixed; }
    inline bool isTriangleFixed(int t) { return (mNodes[mTris[t].c[0]].flags & NfFixed) || (mNodes[mTris[t].c[1]].flags & NfFixed) || (mNodes[mTris[t].c[2]].flags & NfFixed); }
    
    inline const Vec3 getNode(int tri, int c) const { return mNodes[mTris[tri].c[c]].pos; }
    inline Vec3& getNode(int tri, int c) { return mNodes[mTris[tri].c[c]].pos; }
    inline const Vec3 getEdge(int tri, int e) const { return getNode(tri,(e+1)%3) - getNode(tri,e); }
    inline OneRing& get1Ring(int node) { return m1RingLookup[node]; }
    inline Real getFaceArea(int t) const { Vec3 c0 = mNodes[mTris[t].c[0]].pos; return 0.5*norm(cross(mNodes[mTris[t].c[1]].pos - c0, mNodes[mTris[t].c[2]].pos - c0)); }
    inline Vec3 getFaceNormal(int t) { Vec3 c0 = mNodes[mTris[t].c[0]].pos; return getNormalized(cross(mNodes[mTris[t].c[1]].pos - c0, mNodes[mTris[t].c[2]].pos - c0)); }
    inline Vec3 getFaceCenter(int t) const { return (mNodes[mTris[t].c[0]].pos + mNodes[mTris[t].c[1]].pos + mNodes[mTris[t].c[2]].pos) / 3.0; }
    inline std::vector<Node>& getNodeData() { return mNodes; }
    
    void mergeNode(int node, int delnode);
    int addNode(Node a);
    int addTri(Triangle a);
    void addCorner(Corner a);
    void removeTri(int tri);
    void removeTriFromLookup(int tri);
    void removeNodes(const std::vector<int>& deletedNodes);
    void rebuildCorners(int from=0, int to=-1);
    void rebuildLookup(int from=0, int to=-1);
    void rebuildQuickCheck();
    void fastNodeLookupRebuild(int corner);
    void sanityCheck(bool strict=true, std::vector<int>* deletedNodes=0, std::map<int,bool>* taintedTris=0);
    
    void addTriChannel(TriChannel* c) { mTriChannels.push_back(c); rebuildChannels(); }
    void addNodeChannel(NodeChannel* c) { mNodeChannels.push_back(c); rebuildChannels(); }

	//! mesh data functions

	//! create a mesh data object
	PYTHON() PbClass* create(PbType type, PbTypeVec T=PbTypeVec(), const std::string& name = "");
	//! add a mesh data field, set its parent mesh pointer
	void registerMdata(MeshDataBase* mdata);
	void registerMdataReal(MeshDataImpl<Real>* mdata);
	void registerMdataVec3(MeshDataImpl<Vec3>* mdata);
	void registerMdataInt (MeshDataImpl<int >* mdata);
	//! remove a mesh data entry
	void deregister(MeshDataBase* mdata);
	//! add one zero entry to all data fields
	void addAllMdata();
	// note - deletion of mdata is handled in compress function

	//! how many are there?
	IndexInt getNumMdata() const { return mMeshData.size(); }
	//! access one of the fields
	MeshDataBase* getMdata(int i) { return mMeshData[i]; }

	//! update data fields
	void updateDataFields();

protected:    
    void rebuildChannels();
    
    std::vector<Node> mNodes;
    std::vector<Triangle> mTris;
    std::vector<Corner> mCorners;
    std::vector<NodeChannel*> mNodeChannels;
    std::vector<TriChannel*> mTriChannels;
    std::vector<OneRing> m1RingLookup;

	//! store mesh data , each pointer has its own storage vector of a certain type (int, real, vec3)
	std::vector<MeshDataBase*> mMeshData;
	//! lists of different types, for fast operations w/o virtual function calls
	std::vector< MeshDataImpl<Real> *> mMdataReal;
	std::vector< MeshDataImpl<Vec3> *> mMdataVec3;
	std::vector< MeshDataImpl<int> *>  mMdataInt;
	//! indicate that mdata of this mesh is copied, and needs to be freed
	bool mFreeMdata;
};

//******************************************************************************

//! abstract interface for mesh data
PYTHON() class MeshDataBase : public PbClass {
public:
	PYTHON() MeshDataBase(FluidSolver* parent);
	virtual ~MeshDataBase();

	//! data type IDs, in line with those for grids
	enum MdataType { TypeNone = 0, TypeReal = 1, TypeInt = 2, TypeVec3 = 4 };

	//! interface functions, using assert instead of pure virtual for python compatibility
	virtual IndexInt  getSizeSlow() const { assertMsg( false , "Dont use, override..."); return 0; }
	virtual void addEntry()   { assertMsg( false , "Dont use, override..."); return;   }
	virtual MeshDataBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; }
	virtual MdataType getType() const { assertMsg( false , "Dont use, override..."); return TypeNone; }
	virtual void resize(IndexInt size)     { assertMsg( false , "Dont use, override..."); return;  }
	virtual void copyValueSlow(IndexInt from, IndexInt to) { assertMsg( false , "Dont use, override..."); return;  }

	//! set base pointer
	void setMesh(Mesh* set) { mMesh = set; }

	//! debugging
	inline void checkNodeIndex(IndexInt idx) const;

protected:
	Mesh* mMesh;
};


//! abstract interface for mesh data
PYTHON() template<class T>
class MeshDataImpl : public MeshDataBase {
public:
	PYTHON() MeshDataImpl(FluidSolver* parent);
	MeshDataImpl(FluidSolver* parent, MeshDataImpl<T>* other);
	virtual ~MeshDataImpl();

	//! access data
	inline       T& get(IndexInt idx)              { DEBUG_ONLY(checkNodeIndex(idx)); return mData[idx]; }
	inline const T& get(IndexInt idx) const        { DEBUG_ONLY(checkNodeIndex(idx)); return mData[idx]; }
	inline       T& operator[](IndexInt idx)       { DEBUG_ONLY(checkNodeIndex(idx)); return mData[idx]; }
	inline const T& operator[](IndexInt idx) const { DEBUG_ONLY(checkNodeIndex(idx)); return mData[idx]; }

	//! set all values to 0, note - different from meshSystem::clear! doesnt modify size of array (has to stay in sync with parent system)
	PYTHON() void clear();

	//! set grid from which to get data...
	PYTHON() void setSource(Grid<T>* grid, bool isMAC=false );

	//! mesh data base interface
	virtual IndexInt  getSizeSlow() const;
	virtual void addEntry();
	virtual MeshDataBase* clone();
	virtual MdataType getType() const;
	virtual void resize(IndexInt s);
	virtual void copyValueSlow(IndexInt from, IndexInt to);

	IndexInt  size() const { return mData.size(); }

	//! fast inlined functions for per mesh operations
	inline void copyValue(IndexInt from, IndexInt to) { get(to) = get(from); }
	void initNewValue(IndexInt idx, Vec3 pos);

	//! python interface (similar to grid data)
	PYTHON() void setConst(T s);
	PYTHON() void setConstRange(T s, const int begin, const int end);
	PYTHON() MeshDataImpl<T>& copyFrom(const MeshDataImpl<T>& a);
	PYTHON() void add(const MeshDataImpl<T>& a);
	PYTHON() void sub(const MeshDataImpl<T>& a);
	PYTHON() void addConst(T s);
	PYTHON() void addScaled(const MeshDataImpl<T>& a, const T& factor);
	PYTHON() void mult( const MeshDataImpl<T>& a);
	PYTHON() void multConst(T s);
	PYTHON() void safeDiv(const MeshDataImpl<T>& a);
	PYTHON() void clamp(Real min, Real max);
	PYTHON() void clampMin(Real vmin);
	PYTHON() void clampMax(Real vmax);

	PYTHON() Real getMaxAbs();
	PYTHON() Real getMax();
	PYTHON() Real getMin();

	PYTHON() T    sum(const MeshDataImpl<int> *t=NULL, const int itype=0) const;
	PYTHON() Real sumSquare() const;
	PYTHON() Real sumMagnitude() const;

	//! special, set if int flag in t has "flag"
	PYTHON() void setConstIntFlag(T s, const MeshDataImpl<int>& t, const int flag);

	PYTHON() void printMdata(IndexInt start=-1, IndexInt stop=-1, bool printIndex=false);

	//! file io
	PYTHON() void save(const std::string name);
	PYTHON() void load(const std::string name);

	//! get data pointer of mesh data
	PYTHON() std::string getDataPointer();
protected:
	//! data storage
	std::vector<T> mData;

	//! optionally , we might have an associated grid from which to grab new data
	Grid<T>* mpGridSource;
	//! unfortunately , we need to distinguish mac vs regular vec3
	bool mGridSourceMAC;
};

PYTHON() alias MeshDataImpl<int>  MdataInt;
PYTHON() alias MeshDataImpl<Real> MdataReal;
PYTHON() alias MeshDataImpl<Vec3> MdataVec3;


// ***************************************************************************************************************
// Implementation

template<class T>
void SimpleNodeChannel<T>::renumber(const std::vector<int>& newIndex, int newsize) {
    for(size_t i=0; i<newIndex.size(); i++) {
        if(newIndex[i]!=-1)
            data[newIndex[i]] = data[newsize+i];
    }
    data.resize(newsize);
}

inline void MeshDataBase::checkNodeIndex(IndexInt idx) const {
	IndexInt mySize = this->getSizeSlow();
	if (idx<0 || idx > mySize ) {
		errMsg( "MeshData " << " size " << mySize << " : index " << idx << " out of bound " );
	}
	if ( mMesh && mMesh->getSizeSlow()!=mySize ) {
		errMsg( "MeshData " << " size " << mySize << " does not match parent! (" << mMesh->getSizeSlow() << ") " );
	}
}

template<class T>
void MeshDataImpl<T>::clear() {
	for(IndexInt i=0; i<(IndexInt)mData.size(); ++i) mData[i] = 0.;
}


} //namespace
#endif
