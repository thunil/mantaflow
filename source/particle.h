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
PYTHON() class ParticleBase : public PbClass {
public:
	enum SystemType { BASE=0, PARTICLE, VORTEX, FILAMENT, FLIP, TURBULENCE, INDEX };
	
	enum ParticleStatus {
		PNONE         = 0,
		PNEW          = (1<<1),  // particles newly created in this step
		PDELETE       = (1<<10), // mark as deleted, will be deleted in next compress() step
		PINVALID      = (1<<30), // unused
	};

	PYTHON() ParticleBase(FluidSolver* parent);
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

	// particle data functions

	//! create a particle data object
	PYTHON() PbClass* create(PbType type, PbTypeVec T=PbTypeVec(), const std::string& name = "");
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
PYTHON() template<class S> class ParticleSystem : public ParticleBase {
public:    
	PYTHON() ParticleSystem(FluidSolver* parent) : ParticleBase(parent), mDeletes(0), mDeleteChunk(0) {}
	virtual ~ParticleSystem() {};
	
	virtual SystemType getType() const { return S::getType(); };
	
	// accessors
	inline S& operator[](int idx)             { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const S& operator[](int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	// return size of container
	PYTHON() inline int size() const { return mData.size(); }
	// slow virtual function of base class, also returns size
	virtual int getSizeSlow() const { return size(); }

	// query status
	inline int  getStatus(int idx) { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx].flag; }
	inline bool isActive(int idx)  { DEBUG_ONLY(checkPartIndex(idx)); return (mData[idx].flag & PDELETE) == 0; }
	
	//! safe accessor for python
	PYTHON() void setPos(int idx, const Vec3& pos) { DEBUG_ONLY(checkPartIndex(idx)); mData[idx].pos = pos; }
	PYTHON() Vec3 getPos(int idx)                  { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx].pos; }
	//! copy all positions into pdata vec3 field
	PYTHON() void getPosPdata(ParticleDataImpl<Vec3>& target);
	PYTHON() void setPosPdata(ParticleDataImpl<Vec3>& source);
	//! transform coordinate system from one grid size to another (usually upon load)
	void transformPositions( Vec3i dimOld, Vec3i dimNew );

	//! explicitly trigger compression from outside
	void doCompress() { if ( mDeletes > mDeleteChunk) compress(); }
	//! insert buffered positions as new particles, update additional particle data
	void insertBufferedParticles();
	//! resize data vector, and all pdata fields
	void resizeAll(int newsize);
	
	// adding and deleting 
	inline void kill(int idx);
	int add(const S& data);
	// remove all particles, init 0 length arrays (also pdata)
	PYTHON() void clear();
			
	//! Advect particle in grid velocity field
	PYTHON() void advectInGrid(FlagGrid& flags, MACGrid& vel, int integrationMode, bool deleteInObstacle=true, bool stopInObstacle=true );
	
	//! Project particles outside obstacles
	PYTHON() void projectOutside(Grid<Vec3>& gradient);
	
	virtual ParticleBase* clone();
	virtual std::string infoString() const;

	//! debugging
	inline void checkPartIndex(int idx) const;
	
protected:  
	//! deletion count , and interval for re-compressing 
	int mDeletes, mDeleteChunk;    
	//! the particle data
	std::vector<S> mData;    

	//! reduce storage , called by doCompress
	virtual void compress(); 
};

//******************************************************************************

//! Simplest data class for particle systems
struct BasicParticleData {
public:
	BasicParticleData() : pos(0.), flag(0) {}
	BasicParticleData(const Vec3& p) : pos(p), flag(0) {}
	static ParticleBase::SystemType getType() { return ParticleBase::PARTICLE; }

	//! data (note, this size is currently hard coded for uni i/o)
	Vec3 pos;
	int  flag;
};

PYTHON() class BasicParticleSystem : public ParticleSystem<BasicParticleData> {
public:
	PYTHON() BasicParticleSystem(FluidSolver* parent);
	
	//! file io
	PYTHON() void save(std::string name);
	PYTHON() void load(std::string name);

	// save to text file
	void writeParticlesText(std::string name);
	// other output formats
	void writeParticlesRawPositionsGz(std::string name);
	void writeParticlesRawVelocityGz(std::string name);

	// add particles in python
	PYTHON() void addParticle(Vec3 pos) { add(BasicParticleData(pos)); }

	// dangerous, get low level access - avoid usage, only used in vortex filament advection for now
	std::vector<BasicParticleData>& getData() { return mData; }

	PYTHON() void printParts(int start=-1, int stop=-1, bool printIndex=false); 
};


//******************************************************************************

//! Index into other particle system
//  used for grid based neighborhood searches on generic particle systems (stores
//  only active particles, and reduces copied data)
//  note - pos & flag are disabled here, do not use!
struct ParticleIndexData {
public:
	ParticleIndexData() : sourceIndex(0) {}
	static ParticleBase::SystemType getType() { return ParticleBase::INDEX; }

	int  sourceIndex; // index of this particle in the original particle system
	// note - the following two are needed for template instantiation, but not used
	// for the particle index system (use values from original one!)
	static Vec3 pos;  // do not use... 
	static int  flag; // not needed usally 
	//Vec3 pos; // enable for debugging
};

PYTHON() class ParticleIndexSystem : public ParticleSystem<ParticleIndexData> {
public:
	PYTHON() ParticleIndexSystem(FluidSolver* parent) : ParticleSystem<ParticleIndexData>(parent) {};
	
	//! we only need a resize function...
	void resize(int size) { mData.resize(size); }
};



//******************************************************************************

//! Particle set with connectivity
PYTHON() template<class DATA, class CON> 
class ConnectedParticleSystem : public ParticleSystem<DATA> {
public:
	PYTHON() ConnectedParticleSystem(FluidSolver* parent) : ParticleSystem<DATA>(parent) {}
	
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

//******************************************************************************

//! abstract interface for particle data
PYTHON() class ParticleDataBase : public PbClass {
public:
	PYTHON() ParticleDataBase(FluidSolver* parent);
	virtual ~ParticleDataBase(); 

	// data type IDs, in line with those for grids
	enum PdataType { TypeNone = 0, TypeReal = 1, TypeInt = 2, TypeVec3 = 4 };

	// interface functions, using assert instead of pure virtual for python compatibility
	virtual int  getSizeSlow() const { assertMsg( false , "Dont use, override..."); return 0; } 
	virtual void addEntry()   { assertMsg( false , "Dont use, override..."); return;   }
	virtual ParticleDataBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; }
	virtual PdataType getType() const { assertMsg( false , "Dont use, override..."); return TypeNone; } 
	virtual void resize(int size)     { assertMsg( false , "Dont use, override..."); return;  }
	virtual void copyValueSlow(int from, int to) { assertMsg( false , "Dont use, override..."); return;  }

	//! set base pointer
	void setParticleSys(ParticleBase* set) { mpParticleSys = set; }

	//! debugging
	inline void checkPartIndex(int idx) const;

protected:
	ParticleBase* mpParticleSys;
};


//! abstract interface for particle data
PYTHON() template<class T>
class ParticleDataImpl : public ParticleDataBase {
public:
	PYTHON() ParticleDataImpl(FluidSolver* parent);
	ParticleDataImpl(FluidSolver* parent, ParticleDataImpl<T>* other);
	virtual ~ParticleDataImpl();

	//! access data
	inline T& get(int idx)            { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const T get(int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline T& operator[](int idx)            { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const T operator[](int idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }

	// set all values to 0, note - different from particleSystem::clear! doesnt modify size of array (has to stay in sync with parent system)
	PYTHON() void clear();

	//! set grid from which to get data...
	PYTHON() void setSource(Grid<T>* grid, bool isMAC=false );

	// particle data base interface
	virtual int  getSizeSlow() const;
	virtual void addEntry();
	virtual ParticleDataBase* clone();
	virtual PdataType getType() const;
	virtual void resize(int s);
	virtual void copyValueSlow(int from, int to);

	int  size() const { return mData.size(); }

	// fast inlined functions for per particle operations
	inline void copyValue(int from, int to) { get(to) = get(from); } 
	void initNewValue(int idx, Vec3 pos);

	// python interface (similar to grid data)
	PYTHON() void setConst(T s);
	PYTHON() ParticleDataImpl<T>& copyFrom(const ParticleDataImpl<T>& a);
	PYTHON() void add(const ParticleDataImpl<T>& a);
	PYTHON() void sub(const ParticleDataImpl<T>& a);
	PYTHON() void addConst(T s);
	PYTHON() void addScaled(const ParticleDataImpl<T>& a, const T& factor); 
	PYTHON() void mult( const ParticleDataImpl<T>& a);
	PYTHON() void multConst(T s);
	PYTHON() void clamp(Real min, Real max);
	PYTHON() Real getMaxAbsValue();
	PYTHON() Real getMaxValue();
	PYTHON() Real getMinValue();    

	PYTHON() void printPdata(int start=-1, int stop=-1, bool printIndex=false); 
	
	//! file io
	PYTHON() void save(std::string name);
	PYTHON() void load(std::string name);
protected:
	//! data storage
	std::vector<T> mData; 

	//! optionally , we might have an associated grid from which to grab new data
	Grid<T>* mpGridSource;
	//! unfortunately , we need to distinguish mac vs regular vec3
	bool mGridSourceMAC;
};

PYTHON() alias ParticleDataImpl<int>  PdataInt;
PYTHON() alias ParticleDataImpl<Real> PdataReal;
PYTHON() alias ParticleDataImpl<Vec3> PdataVec3;


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
	this->resizeAll(0); // instead of mData.clear
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
void ParticleSystem<S>::getPosPdata(ParticleDataImpl<Vec3>& target) {
	for(int i=0; i<(int)this->size(); ++i) {
		target[i] = this->getPos(i);
	}
}
template<class S>
void ParticleSystem<S>::setPosPdata(ParticleDataImpl<Vec3>& target) {
	for(int i=0; i<(int)this->size(); ++i) {
		this->getPos(i) = target[i];
	}
}

template<class S>
void ParticleSystem<S>::transformPositions( Vec3i dimOld, Vec3i dimNew )
{
	Vec3 factor = calcGridSizeFactor( dimNew, dimOld );
	for(int i=0; i<(int)this->size(); ++i) {
		this->setPos(i, this->getPos(i) * factor );
	}
}

// check for deletion/invalid position, otherwise return velocity
KERNEL(pts) returns(std::vector<Vec3> u(size)) template<class S>
std::vector<Vec3> GridAdvectKernel (std::vector<S>& p, const MACGrid& vel, const FlagGrid& flags, Real dt,
		bool deleteInObstacle, bool stopInObstacle )
{
	if (p[idx].flag & ParticleBase::PDELETE) {
		u[idx] =_0; return;
	} 
	// special handling
	/*NT_DEBUG if(deleteInObstacle || stopInObstacle) {
		if (!flags.isInBounds(p[idx].pos, 1) || flags.isObstacle(p[idx].pos) ) {
			if(stopInObstacle)
				u[idx] = _0; 
			// for simple tracer particles, its convenient to delete particles right away
			// for other sim types, eg flip, we can try to fix positions later on
			if(deleteInObstacle) 
				p[idx].flag |= ParticleBase::PDELETE; 
			return;
		} 
	}*/
	u[idx] = vel.getInterpolated(p[idx].pos) * dt;
};

// final check after advection to make sure particles haven't escaped
// (similar to particle advection kernel)
KERNEL(pts) template<class S>
void KnDeleteInObstacle(std::vector<S>& p, const FlagGrid& flags) {
	if (p[idx].flag & ParticleBase::PDELETE) return;
	if (!flags.isInBounds(p[idx].pos,1) || flags.isObstacle(p[idx].pos)) {
		p[idx].flag |= ParticleBase::PDELETE;
	} 
}

// try to get closer to actual obstacle boundary
static inline Vec3 bisectBacktracePos(const FlagGrid& flags, const Vec3& oldp, const Vec3& newp)
{
	Real s = 0.;
	for(int i=1; i<5; ++i) {
		Real ds = 1./(Real)(1<<i);
		if (!flags.isObstacle( oldp*(1.-(s+ds)) + newp*(s+ds) )) {
			s += ds;
		}
	}
	return( oldp*(1.-(s)) + newp*(s) );
}

// at least make sure all particles are inside domain
KERNEL(pts) template<class S>
void KnClampPositions(std::vector<S>& p, const FlagGrid& flags, ParticleDataImpl<Vec3> *posOld = NULL, bool stopInObstacle=true)
{
	if (p[idx].flag & ParticleBase::PDELETE) return;
	if (!flags.isInBounds(p[idx].pos,0) ) {
		p[idx].pos = clamp( p[idx].pos, Vec3(0.), toVec3(flags.getSize())-Vec3(1.) );
	} 
	if (stopInObstacle && (flags.isObstacle(p[idx].pos)) ) {
		p[idx].pos = bisectBacktracePos(flags, (*posOld)[idx], p[idx].pos);
	}
}

// advection plugin
template<class S>
void ParticleSystem<S>::advectInGrid(FlagGrid& flags, MACGrid& vel, int integrationMode, bool deleteInObstacle, bool stopInObstacle ) {
	// position clamp requires old positions, backup
	ParticleDataImpl<Vec3> *posOld = NULL;
	if(!deleteInObstacle) {
		posOld = new ParticleDataImpl<Vec3>(this->getParent());
		posOld->resize(mData.size());
		for(int i=0; i<mData.size();++i) (*posOld)[i] = mData[i].pos;
	}

	// update positions
	GridAdvectKernel<S> kernel(mData, vel, flags, getParent()->getDt(), deleteInObstacle, stopInObstacle );
	integratePointSet(kernel, integrationMode);

	if(!deleteInObstacle) {
		KnClampPositions<S>  ( mData, flags, posOld , stopInObstacle );
		delete posOld;
	} else {
		KnDeleteInObstacle<S>( mData, flags);
	}
}

KERNEL(pts, single) // no thread-safe random gen yet
template<class S>
void KnProjectParticles(ParticleSystem<S>& part, Grid<Vec3>& gradient) {
	static RandomStream rand (3123984);
	const double jlen = 0.1;
	
	if (part.isActive(idx)) {
		// project along levelset gradient
		Vec3 p = part[idx].pos;
		if (gradient.isInBounds(p)) {
			Vec3 n = gradient.getInterpolated(p);
			Real dist = normalize(n);
			Vec3 dx = n * (-dist + jlen * (1 + rand.getReal()));
			p += dx;            
		}
		// clamp to outer boundaries (+jitter)
		const double jlen = 0.1;
		Vec3 jitter = jlen * rand.getVec3();
		part[idx].pos = clamp(p, Vec3(1,1,1)+jitter, toVec3(gradient.getSize()-1)-jitter);
	}
}

template<class S>
void ParticleSystem<S>::projectOutside(Grid<Vec3>& gradient) {
	KnProjectParticles<S>(*this, gradient);
}

template<class S>
void ParticleSystem<S>::resizeAll(int size) {
	// resize all buffers to target size in 1 go
	mData.resize(size);
	for(int i=0; i<(int)mPartData.size(); ++i)
		mPartData[i]->resize(size);
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

	resizeAll(nextRead);
	mDeletes = 0;
	mDeleteChunk = mData.size() / DELETE_PART;
}

//! insert buffered positions as new particles, update additional particle data
template<class S>
void ParticleSystem<S>::insertBufferedParticles() {
	if(mNewBuffer.size()==0) return;
	int newCnt = mData.size();
	resizeAll(newCnt + mNewBuffer.size());

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
	s << "ParticleSys '" << getName() << "' [" << size() << " parts";
	if(this->getNumPdata()>0) s<< " "<< this->getNumPdata()<<" pd";
	s << "]";
	//for(int i=0; i<(int)mPartData.size(); ++i) { sstr << i<<":" << mPartData[i]->size() <<" "; } 
	return s.str();
}
	
template<class S>  
inline void ParticleSystem<S>::checkPartIndex(int idx) const {
	int mySize = this->size();
	if (idx<0 || idx > mySize ) {
		errMsg( "ParticleBase " << " size " << mySize << " : index " << idx << " out of bound " );
	}
}
	
inline void ParticleDataBase::checkPartIndex(int idx) const {
	int mySize = this->getSizeSlow();
	if (idx<0 || idx > mySize ) {
		errMsg( "ParticleData " << " size " << mySize << " : index " << idx << " out of bound " );
	}
	if ( mpParticleSys && mpParticleSys->getSizeSlow()!=mySize ) {
		errMsg( "ParticleData " << " size " << mySize << " does not match parent! (" << mpParticleSys->getSizeSlow() << ") " );
	}
}

// set contents to zero, as for a grid
template<class T>
void ParticleDataImpl<T>::clear() {
	for(int i=0; i<(int)mData.size(); ++i) mData[i] = 0.;
}


} // namespace

#endif

