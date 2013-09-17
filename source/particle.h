/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Base class for particle systems
 *
 ******************************************************************************/

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <vector>
#include "grid.h"
#include "vectorbase.h"
#include "integrator.h"
#include "randomstream.h"
namespace Manta {

// fwd decl
template<class T> class Grid;
class ParticleDataBase;
template<class T> class ParticleDataImpl;

//! Baseclass for particle systems. Does not implement any data
PYTHON class ParticleBase : public PbClass {
public:
    enum SystemType { BASE=0, PARTICLE, VORTEX, FILAMENT, FLIP, TURBULENCE };
    
    enum ParticleStatus {
        PNONE         = 0,
        PNEW          = (1<<1),  // particles newly created in this step
        PDELETE       = (1<<10), // mark as deleted, will be deleted in next compress() step
        PINVALID      = (1<<30), // unused
    };

    PYTHON ParticleBase(FluidSolver* parent);
    virtual ~ParticleBase();

	//! copy all the particle data thats registered with the other particle system to this one
    virtual void cloneParticleData(ParticleBase* nm);

    virtual SystemType getType() const { return BASE; }
    virtual std::string infoString() const; 
    virtual ParticleBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; } 

	// slow virtual function to query size, do not use in kernels! use size() instead
	virtual int getSizeSlow() const { assertMsg( false , "Dont use, override..."); return 0; } 

	//! add a position as potential candidate for new particle (todo, make usable from parallel threads)
	inline void addBuffered(const Vec3& pos);

	//! debug info about pdata
	std::string debugInfoPdata();

	// particle data functions

    //! create a particle data object
    PYTHON PbClass* create(PbType type, const std::string& name = "");
	//! add a particle data field, set its parent particle-system pointer
	void registerPdata(ParticleDataBase* pdata);
	void registerPdataReal(ParticleDataImpl<Real>* pdata);
	void registerPdataVec3(ParticleDataImpl<Vec3>* pdata);
	void registerPdataInt (ParticleDataImpl<int >* pdata);
	//! remove a particle data entry
	void deregister(ParticleDataBase* pdata);
	//! add one zero entry to all data fields
	void addAllPdata();
	// note - deletion of pdata is handled in compress function

	//! how many are there?
	int getNumPdata() const { return mPartData.size(); }
	//! access one of the fields
	ParticleDataBase* getPdata(int i) { return mPartData[i]; }

protected:  
	//! new particle candidates
	std::vector<Vec3> mNewBuffer;

	//! allow automatic compression / resize? disallowed for, eg, flip particle systems
	bool mAllowCompress;

	//! store particle data , each pointer has its own storage vector of a certain type (int, real, vec3)
	std::vector<ParticleDataBase*> mPartData;
	//! lists of different types, for fast operations w/o virtual function calls (all calls necessary per particle)
	std::vector< ParticleDataImpl<Real> *> mPdataReal;
	std::vector< ParticleDataImpl<Vec3> *> mPdataVec3;
	std::vector< ParticleDataImpl<int> *>  mPdataInt;
	//! indicate that pdata of this particle system is copied, and needs to be freed
	bool mFreePdata;
};


//! Main class for particle systems
/*! Basetype S must at least contain flag, pos fields */
PYTHON template<class S> class ParticleSystem : public ParticleBase {
public:    
    PYTHON ParticleSystem(FluidSolver* parent) : ParticleBase(parent), mDeletes(0), mDeleteChunk(0) {}
    virtual ~ParticleSystem() {};
    
    virtual SystemType getType() const { return S::getType(); };
    
    // accessors
    inline S& operator[](int i) { return mData[i]; }
    inline const S& operator[](int i) const { return mData[i]; }
    PYTHON inline int size() const { return mData.size(); }
	// slow virtual function of base class
	virtual int getSizeSlow() const { return size(); }
    std::vector<S>& getData() { return mData; }

	//! explicitly trigger compression from outside
	void doCompress() { if ( mDeletes > mDeleteChunk) compress(); }
	//! insert buffered positions as new particles, update additional particle data
	void insertBufferedParticles();
    
    // adding and deleting 
    inline void kill(int i);
    int add(const S& data);
    void clear();
	// query status
    inline bool isActive(int i);
    inline int  getStatus(int i);
    
    //! safe accessor for python
    PYTHON void setPos(int idx, const Vec3& pos);
    //! safe accessor for python
    PYTHON Vec3 getPos(int idx);
            
    //! Advect particle in grid velocity field
    PYTHON void advectInGrid(FlagGrid& flaggrid, MACGrid& vel, int integrationMode, bool deleteInObstacle=true );
    
    //! Project particles outside obstacles
    PYTHON void projectOutside(Grid<Vec3>& gradient);
    
    virtual ParticleBase* clone();
    virtual std::string infoString() const;
    
protected:  
   	//! deletion count , and interval for re-compressing 
    int mDeletes, mDeleteChunk;    
	//! the particle data
    std::vector<S> mData;    

	//! reduce storage , called by doCompress
    virtual void compress(); 
};

//! Simplest data class for particle systems
struct BasicParticleData {
public:
    BasicParticleData() : pos(0.), flag(0) {}
    BasicParticleData(const Vec3& p) : pos(p), flag(0) {}
    static ParticleBase::SystemType getType() { return ParticleBase::PARTICLE; }

	//! data
    Vec3 pos;
    int flag;
};

PYTHON class BasicParticleSystem : public ParticleSystem<BasicParticleData> {
public:
    PYTHON BasicParticleSystem(FluidSolver* parent);
    
	//! file io
    PYTHON void save(std::string name);
    PYTHON void load(std::string name);

	// save to text file
	void writeParticlesText(std::string name);

    PYTHON void addParticle(Vec3 pos) { add(BasicParticleData(pos)); }
};




//! Particle set with connectivity
PYTHON template<class DATA, class CON> 
class ConnectedParticleSystem : public ParticleSystem<DATA> {
public:
    PYTHON ConnectedParticleSystem(FluidSolver* parent) : ParticleSystem<DATA>(parent) {}
    
    // accessors
    inline bool isSegActive(int i) { return (mSegments[i].flag & ParticleBase::PDELETE) == 0; }    
    inline int segSize() const { return mSegments.size(); }    
    inline CON& seg(int i) { return mSegments[i]; }
    inline const CON& seg(int i) const { return mSegments[i]; }
        
    virtual ParticleBase* clone();
    
protected:
    std::vector<CON> mSegments;
    virtual void compress();    
};


//! abstract interface for particle data
PYTHON class ParticleDataBase : public PbClass {
public:
    PYTHON ParticleDataBase(FluidSolver* parent);
	virtual ~ParticleDataBase(); 

    enum PdataType { UNKNOWN=0, DATA_INT, DATA_REAL, DATA_VEC3 };

	// interface functions, using assert instead of pure virtual for python compatibility
	virtual int  size() const { assertMsg( false , "Dont use, override..."); return 0; } 
	virtual void add()        { assertMsg( false , "Dont use, override..."); return;   }
    virtual ParticleDataBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; }
	virtual PdataType getType() const { assertMsg( false , "Dont use, override..."); return UNKNOWN; } 
	virtual void resize(int i)        { assertMsg( false , "Dont use, override..."); return;  }
	virtual void copyValueSlow(int from, int to) { assertMsg( false , "Dont use, override..."); return;  }

	//! set base pointer
	void setParticleSys(ParticleBase* set) { mpParticleSys = set; }

	//! debugging
	inline void checkPartIndex(int idx) const;

protected:
	ParticleBase* mpParticleSys;
};

//! abstract interface for particle data
PYTHON template<class T>
class ParticleDataImpl : public ParticleDataBase {
public:
	PYTHON ParticleDataImpl(FluidSolver* parent);
	ParticleDataImpl(FluidSolver* parent, ParticleDataImpl<T>* other);
	virtual ~ParticleDataImpl();

    //! access data
    inline T& get(int idx)            { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
    inline const T get(int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
    inline T& operator[](int idx)            { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
    inline const T operator[](int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }

	//! set grid from which to get data...
	PYTHON void setSource(Grid<T>* grid, bool isMAC=false );

	// particle data base interface
	virtual int  size() const;
	virtual void add();
    virtual ParticleDataBase* clone();
	virtual PdataType getType() const;
	virtual void resize(int s);
	virtual void copyValueSlow(int from, int to);

	// fast inlined functions for per particle operations
	inline void copyValue(int from, int to) { get(to) = get(from); } 
	void initNewValue(int idx, Vec3 pos);
protected:
	//! data storage
	std::vector<T> mData; 

	//! optionally , we might have an associated grid from which to grab new data
	Grid<T>* mpGridSource;
	//! unfortunately , we need to distinguish mac vs regular vec3
	bool mGridSourceMAC;
};

PYTHON alias ParticleDataImpl<int>  PdataInt;
PYTHON alias ParticleDataImpl<Real> PdataReal;
PYTHON alias ParticleDataImpl<Vec3> PdataVec3;


//******************************************************************************
// Implementation
//******************************************************************************

const int DELETE_PART = 20; // chunk size for compression

void ParticleBase::addBuffered(const Vec3& pos) {
	mNewBuffer.push_back(pos);
}
   
template<class S>
void ParticleSystem<S>::clear() {
    mDeleteChunk = mDeletes = 0;
    mData.clear();
}

template<class S>
int ParticleSystem<S>::add(const S& data) {
    mData.push_back(data); 
    mDeleteChunk = mData.size() / DELETE_PART;
	this->addAllPdata();
    return mData.size()-1;
}

template<class S>
inline void ParticleSystem<S>::kill(int idx)     { 
    assertMsg(idx>=0 && idx<size(), "Index out of bounds");
	mData[idx].flag |= PDELETE; 
	if ( (++mDeletes > mDeleteChunk) && (mAllowCompress) ) compress(); 
}

template<class S>
inline bool ParticleSystem<S>::isActive(int idx) { 
    assertMsg(idx>=0 && idx<size(), "Index out of bounds");
	return (mData[idx].flag & PDELETE) == 0; 
}  
template<class S>
inline int ParticleSystem<S>::getStatus(int idx) { 
    assertMsg(idx>=0 && idx<size(), "Index out of bounds");
	return mData[idx].flag;
}  

template<class S> Vec3 ParticleSystem<S>::getPos(int idx) {
    assertMsg(idx>=0 && idx<size(), "Index out of bounds");
    return mData[idx].pos;
}

template<class S> void ParticleSystem<S>::setPos(int idx, const Vec3& pos) {
    assertMsg(idx>=0 && idx<size(), "Index out of bounds");
    mData[idx].pos = pos;
}
KERNEL(pts) template<class S> returns(std::vector<Vec3> u(size))
std::vector<Vec3> GridAdvectKernel (std::vector<S>& p, const MACGrid& vel, const FlagGrid& flaggrid, Real dt,
		bool deleteInObstacle )
{
    if (p[i].flag & ParticleBase::PDELETE) {
        u[i] =_0;
	} else if (!flaggrid.isInBounds(p[i].pos,1) || flaggrid.isObstacle(p[i].pos)) {
        u[i] = _0;

		// for simple tracer particles, its convenient to delete particles right away
		// eg for flip sim, its not a good idea
		if(deleteInObstacle) 
			p[i].flag |= ParticleBase::PDELETE;
    } else {
        u[i] = vel.getInterpolated(p[i].pos) * dt;
	}
};

// advection plugin
template<class S>
void ParticleSystem<S>::advectInGrid(FlagGrid& flaggrid, MACGrid& vel, int integrationMode, bool deleteInObstacle ) {
    GridAdvectKernel<S> kernel(mData, vel, flaggrid, getParent()->getDt(), deleteInObstacle );
    integratePointSet(kernel, integrationMode);
}

KERNEL(pts, single) // no thread-safe random gen yet
template<class S>
void KnProjectParticles(ParticleSystem<S>& part, Grid<Vec3>& gradient) {
    static RandomStream rand (3123984);
    const double jlen = 0.1;
    
    if (part.isActive(i)) {
        // project along levelset gradient
        Vec3 p = part[i].pos;
        if (gradient.isInBounds(p)) {
            Vec3 n = gradient.getInterpolated(p);
            Real dist = normalize(n);
            Vec3 dx = n * (-dist + jlen * (1 + rand.getReal()));
            p += dx;            
        }
        // clamp to outer boundaries (+jitter)
        const double jlen = 0.1;
        Vec3 jitter = jlen * rand.getVec3();
        part[i].pos = clamp(p, Vec3(1,1,1)+jitter, toVec3(gradient.getSize()-1)-jitter);
    }
}

template<class S>
void ParticleSystem<S>::projectOutside(Grid<Vec3>& gradient) {
    KnProjectParticles<S>(*this, gradient);
}

template<class S>
void ParticleSystem<S>::compress() {
    int nextRead = mData.size();
    for (int i=0; i<(int)mData.size(); i++) {
        while ((mData[i].flag & PDELETE) != 0) {
            nextRead--;
            mData[i] = mData[nextRead];
			// ugly, but prevent virtual function calls here:
			for(int pd=0; pd<(int)mPdataReal.size(); ++pd) mPdataReal[pd]->copyValue(nextRead, i);
			for(int pd=0; pd<(int)mPdataVec3.size(); ++pd) mPdataVec3[pd]->copyValue(nextRead, i);
			for(int pd=0; pd<(int)mPdataInt .size(); ++pd) mPdataInt [pd]->copyValue(nextRead, i);
            mData[nextRead].flag = PINVALID;
        }
    }
	if(nextRead<(int)mData.size()) debMsg("Deleted "<<((int)mData.size() - nextRead)<<" particles", 1); // debug info

    mData.resize(nextRead);
	for(int i=0; i<(int)mPartData.size(); ++i)
		mPartData[i]->resize(nextRead);

    mDeletes = 0;
    mDeleteChunk = mData.size() / DELETE_PART;
}

//! insert buffered positions as new particles, update additional particle data
template<class S>
void ParticleSystem<S>::insertBufferedParticles() {
	if(mNewBuffer.size()==0) return;
	int newCnt = mData.size();

	// resize all buffers to target size in 1 go
    mData.resize(newCnt + mNewBuffer.size());
	for(int i=0; i<(int)mPartData.size(); ++i)
		mPartData[i]->resize(newCnt + mNewBuffer.size());

	// clear new flag everywhere
	for(int i=0; i<(int)mData.size(); ++i) mData[i].flag &= ~PNEW;

	for(int i=0; i<(int)mNewBuffer.size(); ++i) {
		// note, other fields are not initialized here...
		mData[newCnt].pos  = mNewBuffer[i];
		mData[newCnt].flag = PNEW;
		// now init pdata fields from associated grids...
		for(int pd=0; pd<(int)mPdataReal.size(); ++pd) 
			mPdataReal[pd]->initNewValue(newCnt, mNewBuffer[i] );
		for(int pd=0; pd<(int)mPdataVec3.size(); ++pd) 
			mPdataVec3[pd]->initNewValue(newCnt, mNewBuffer[i] );
		for(int pd=0; pd<(int)mPdataInt.size(); ++pd) 
			mPdataInt[pd]->initNewValue(newCnt, mNewBuffer[i] );
		newCnt++;
	}
	if(mNewBuffer.size()>0) debMsg("Added & initialized "<<(int)mNewBuffer.size()<<" particles", 1); // debug info
	mNewBuffer.clear();
}


template<class DATA, class CON>
void ConnectedParticleSystem<DATA,CON>::compress() {
    const int sz = ParticleSystem<DATA>::size();
    int *renumber_back = new int[sz];
    int *renumber = new int[sz];
    for (int i=0; i<sz; i++)
        renumber[i] = renumber_back[i] = -1;
        
    // reorder elements
    std::vector<DATA>& data = ParticleSystem<DATA>::mData;
    int nextRead = sz;
    for (int i=0; i<nextRead; i++) {
        if ((data[i].flag & ParticleBase::PDELETE) != 0) {
            nextRead--;
            data[i] = data[nextRead];
            data[nextRead].flag = 0;           
            renumber_back[i] = nextRead;
        } else 
            renumber_back[i] = i;
    }
    
    // acceleration structure
    for (int i=0; i<nextRead; i++)
        renumber[renumber_back[i]] = i;
    
    // rename indices in filaments
    for (int i=0; i<(int)mSegments.size(); i++)
        mSegments[i].renumber(renumber);
        
    ParticleSystem<DATA>::mData.resize(nextRead);
    ParticleSystem<DATA>::mDeletes = 0;
    ParticleSystem<DATA>::mDeleteChunk = ParticleSystem<DATA>::size() / DELETE_PART;
    
    delete[] renumber;
    delete[] renumber_back;
}

template<class S>
ParticleBase* ParticleSystem<S>::clone() {
    ParticleSystem<S>* nm = new ParticleSystem<S>(getParent());
    if(this->mAllowCompress) compress();
    
    nm->mData = mData;
    nm->setName(getName());
	this->cloneParticleData(nm);
    return nm;
}

template<class DATA,class CON>
ParticleBase* ConnectedParticleSystem<DATA,CON>::clone() {
    ConnectedParticleSystem<DATA,CON>* nm = new ConnectedParticleSystem<DATA,CON>(this->getParent());
    if(this->mAllowCompress) compress();
    
    nm->mData = this->mData;
    nm->mSegments = mSegments;
    nm->setName(this->getName());
	this->cloneParticleData(nm);
    return nm;
}

template<class S>  
std::string ParticleSystem<S>::infoString() const { 
    std::stringstream s;
    s << "ParticleSys '" << getName() << "' [" << size() << " parts]";
	if(this->getNumPdata()>0) s<< " Pdata: "<< this->getNumPdata();
    return s.str();
}
    
    
inline void ParticleDataBase::checkPartIndex(int idx) const {
	int mySize = this->size();
    if (idx<0 || idx > mySize ) {
        errMsg( "ParticleData " << " size " << mySize << " : index " << idx << " out of bound " );
    }
    if ( mpParticleSys && mpParticleSys->getSizeSlow()!=mySize ) {
        errMsg( "ParticleData " << " size " << mySize << " does not match parent! (" << mpParticleSys->getSizeSlow() << ") " );
    }
}

} // namespace

#endif

