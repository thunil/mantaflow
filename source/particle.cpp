/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2013 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Particle data functionality
 *
 ******************************************************************************/

#include <fstream>
#include  <cstring>
#include  <limits>
#if NO_ZLIB!=1
#include <zlib.h>
#endif
#include "particle.h"
#include "levelset.h"
#include "fileio.h"

using namespace std;
namespace Manta {


ParticleBase::ParticleBase(FluidSolver* parent) : 
	PbClass(parent), mAllowCompress(true), mFreePdata(false) {
}

ParticleBase::~ParticleBase()
{
	// make sure data fields now parent system is deleted
	for(int i=0; i<(int)mPartData.size(); ++i)
		mPartData[i]->setParticleSys(NULL);
	
	if(mFreePdata) {
		for(int i=0; i<(int)mPartData.size(); ++i)
			delete mPartData[i];
	}
	
}

std::string ParticleBase::infoString() const { 
	return "ParticleSystem " + mName + " <no info>"; 
}

void ParticleBase::cloneParticleData(ParticleBase* nm) {
	// clone additional data , and make sure the copied particle system deletes it
	nm->mFreePdata = true;
	for(int i=0; i<(int)mPartData.size(); ++i) {
		ParticleDataBase* pdata = mPartData[i]->clone();
		nm->registerPdata(pdata);
	} 
}

void ParticleBase::deregister(ParticleDataBase* pdata) {
	bool done = false;
	// remove pointer from particle data list
	for(int i=0; i<(int)mPartData.size(); ++i) {
		if(mPartData[i] == pdata) {
			if(i<(int)mPartData.size()-1)
				mPartData[i] = mPartData[mPartData.size()-1];
			mPartData.pop_back();
			done = true;
		}
	} 
	if(!done)
		errMsg("Invalid pointer given, not registered!");
}

// create and attach a new pdata field to this particle system
PbClass* ParticleBase::create(PbType t, PbTypeVec T, const string& name) {        
#	if NOPYTHON!=1
	_args.add("nocheck",true);
	if (t.str() == "")
		errMsg("Specify particle data type to create");
	//debMsg( "Pdata creating '"<< t.str <<" with size "<< this->getSizeSlow(), 5 );
	
	PbClass* pyObj = PbClass::createPyObject(t.str() + T.str(), name, _args, this->getParent() );

	ParticleDataBase* pdata = dynamic_cast<ParticleDataBase*>(pyObj);
	if(!pdata) {
		errMsg("Unable to get particle data pointer from newly created object. Only create ParticleData type with a ParticleSys.creat() call, eg, PdataReal, PdataVec3 etc.");
		delete pyObj;
		return NULL;
	} else {
		this->registerPdata(pdata);
	}

	// directly init size of new pdata field:
	pdata->resize( this->getSizeSlow() );
#	else
	PbClass* pyObj = NULL;
#	endif
	return pyObj;
}

void ParticleBase::registerPdata(ParticleDataBase* pdata) {
	pdata->setParticleSys(this);
	mPartData.push_back(pdata);

	if( pdata->getType() == ParticleDataBase::TypeReal ) {
		ParticleDataImpl<Real>* pd = dynamic_cast< ParticleDataImpl<Real>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as real!");
		this->registerPdataReal(pd);
	}
	else if( pdata->getType() == ParticleDataBase::TypeInt ) {
		ParticleDataImpl<int>* pd = dynamic_cast< ParticleDataImpl<int>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as int!");
		this->registerPdataInt(pd);
	}
	else if( pdata->getType() == ParticleDataBase::TypeVec3 ) {
		ParticleDataImpl<Vec3>* pd = dynamic_cast< ParticleDataImpl<Vec3>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as vec3!");
		this->registerPdataVec3(pd);
	}
}
void ParticleBase::registerPdataReal(ParticleDataImpl<Real>* pd) { mPdataReal.push_back(pd); }
void ParticleBase::registerPdataVec3(ParticleDataImpl<Vec3>* pd) { mPdataVec3.push_back(pd); }
void ParticleBase::registerPdataInt (ParticleDataImpl<int >* pd) { mPdataInt .push_back(pd); }

void ParticleBase::addAllPdata() {
	for(int i=0; i<(int)mPartData.size(); ++i) {
		mPartData[i]->addEntry();
	} 
}

 
BasicParticleSystem::BasicParticleSystem(FluidSolver* parent)
	   : ParticleSystem<BasicParticleData>(parent) {
	this->mAllowCompress = false;
}

// file io

void BasicParticleSystem::writeParticlesText(string name) {
	ofstream ofs(name.c_str());
	if (!ofs.good())
		errMsg("can't open file!");
	ofs << this->size()<<", pdata: "<< mPartData.size()<<" ("<<mPdataInt.size()<<","<<mPdataReal.size()<<","<<mPdataVec3.size()<<") \n";
	for(int i=0; i<this->size(); ++i) {
		ofs << i<<": "<< this->getPos(i) <<" , "<< this->getStatus(i) <<". "; 
		for(int pd=0; pd<(int)mPdataInt.size() ; ++pd) ofs << mPdataInt [pd]->get(i)<<" ";
		for(int pd=0; pd<(int)mPdataReal.size(); ++pd) ofs << mPdataReal[pd]->get(i)<<" ";
		for(int pd=0; pd<(int)mPdataVec3.size(); ++pd) ofs << mPdataVec3[pd]->get(i)<<" ";
		ofs << "\n"; 
	}
	ofs.close();
}

void BasicParticleSystem::writeParticlesRawPositionsGz(string name) {
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "wb1");
	if (!gzf) errMsg("can't open file "<<name);
	for(int i=0; i<this->size(); ++i) {
		Vector3D<float> p = toVec3f( this->getPos(i) );
		gzwrite(gzf, &p, sizeof(float)*3);
	}
	gzclose(gzf);
#	else
	cout << "file format not supported without zlib" << endl;
#	endif
}

void BasicParticleSystem::writeParticlesRawVelocityGz(string name) {
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "wb1");
	if (!gzf) errMsg("can't open file "<<name);
	if( mPdataVec3.size() < 1 ) errMsg("no vec3 particle data channel found!");
	// note , assuming particle data vec3 0 is velocity! make optional...
	for(int i=0; i<this->size(); ++i) {		
		Vector3D<float> p = toVec3f( mPdataVec3[0]->get(i) );
		gzwrite(gzf, &p, sizeof(float)*3);
	}
	gzclose(gzf);
#	else
	cout << "file format not supported without zlib" << endl;
#	endif
}


void BasicParticleSystem::load(string name ) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if ( ext == ".uni") 
		readParticlesUni(name, this );
	else if ( ext == ".raw") // raw = uni for now
		readParticlesUni(name, this );
	else 
		errMsg("particle '" + name +"' filetype not supported for loading");
}

void BasicParticleSystem::save(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".txt") 
		this->writeParticlesText(name);
	else if (ext == ".uni") 
		writeParticlesUni(name, this);
	else if (ext == ".raw") // raw = uni for now
		writeParticlesUni(name, this);
	// raw data formats, very basic for simple data transfer to other programs
	else if (ext == ".posgz") 
		this->writeParticlesRawPositionsGz(name);
	else if (ext == ".velgz") 
		this->writeParticlesRawVelocityGz(name);
	else
		errMsg("particle '" + name +"' filetype not supported for saving");
}

void BasicParticleSystem::printParts(int start, int stop, bool printIndex)
{
	std::ostringstream sstr;
	int s = (start>0 ? start : 0                 );
	int e = (stop>0  ? stop  : (int)mData.size() );
	s = Manta::clamp(s, 0, (int)mData.size());
	e = Manta::clamp(e, 0, (int)mData.size());

	for(int i=s; i<e; ++i) {
		if(printIndex) sstr << i<<": ";
		sstr<<mData[i].pos<<" "<<mData[i].flag<<"\n";
	} 
	debMsg( sstr.str() , 1 );
}

// particle data

ParticleDataBase::ParticleDataBase(FluidSolver* parent) : 
		PbClass(parent) , mpParticleSys(NULL) {
}

ParticleDataBase::~ParticleDataBase()
{
	// notify parent of deletion 
	if(mpParticleSys)
		mpParticleSys->deregister(this);
}


// actual data implementation

template<class T>
ParticleDataImpl<T>::ParticleDataImpl(FluidSolver* parent) : 
	ParticleDataBase(parent) , mpGridSource(NULL), mGridSourceMAC(false) {
}

template<class T>
ParticleDataImpl<T>::ParticleDataImpl(FluidSolver* parent, ParticleDataImpl<T>* other) : 
	ParticleDataBase(parent) , mpGridSource(NULL), mGridSourceMAC(false) {
	this->mData = other->mData;
}

template<class T>
ParticleDataImpl<T>::~ParticleDataImpl() {
}

template<class T>
int ParticleDataImpl<T>::getSizeSlow() const {
	return mData.size();
}
template<class T>
void ParticleDataImpl<T>::addEntry() {
	// add zero'ed entry
	T tmp = T(0.);
	// for debugging, force init:
	//tmp = T(0.02 * mData.size()); // increasing
	//tmp = T(1.); // constant 1
	return mData.push_back(tmp);
}
template<class T>
void ParticleDataImpl<T>::resize(int s) {
	mData.resize(s);
}
template<class T>
void ParticleDataImpl<T>::copyValueSlow(int from, int to) {
	this->copyValue(from,to);
}
template<class T>
ParticleDataBase* ParticleDataImpl<T>::clone() {
	ParticleDataImpl<T>* npd = new ParticleDataImpl<T>( getParent(), this );
	return npd;
}

template<class T>
void ParticleDataImpl<T>::setSource(Grid<T>* grid, bool isMAC ) {
	mpGridSource = grid;
	mGridSourceMAC = isMAC;
	if(isMAC) assertMsg( dynamic_cast<MACGrid*>(grid) != NULL , "Given grid is not a valid MAC grid");
}

template<class T>
void ParticleDataImpl<T>::initNewValue(int idx, Vec3 pos) {
	if(!mpGridSource)
		mData[idx] = 0; 
	else {
		mData[idx] = mpGridSource->getInterpolated(pos);
	}
}
// special handling needed for velocities
template<>
void ParticleDataImpl<Vec3>::initNewValue(int idx, Vec3 pos) {
	if(!mpGridSource)
		mData[idx] = 0;
	else {
		if(!mGridSourceMAC)
			mData[idx] = mpGridSource->getInterpolated(pos);
		else
			mData[idx] = ((MACGrid*)mpGridSource)->getInterpolated(pos);
	}
}

template<typename T>
void ParticleDataImpl<T>::load(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if ( ext == ".uni") 
		readPdataUni<T>(name, this);
	else if ( ext == ".raw") // raw = uni for now 
		readPdataUni<T>(name, this);
	else 
		errMsg("particle data '" + name +"' filetype not supported for loading");
}

template<typename T>
void ParticleDataImpl<T>::save(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".uni") 
		writePdataUni<T>(name, this);
	else if (ext == ".raw") // raw = uni for now
		writePdataUni<T>(name, this);
	else
		errMsg("particle data '" + name +"' filetype not supported for saving");
}

// specializations

template<>
ParticleDataBase::PdataType ParticleDataImpl<Real>::getType() const {
	return ParticleDataBase::TypeReal;
} 
template<>
ParticleDataBase::PdataType ParticleDataImpl<int>::getType() const {
	return ParticleDataBase::TypeInt;
} 
template<>
ParticleDataBase::PdataType ParticleDataImpl<Vec3>::getType() const {
	return ParticleDataBase::TypeVec3;
}

// note, we need a flag value for functions such as advection
// ideally, this value should never be modified
int ParticleIndexData::flag = 0; 
Vec3 ParticleIndexData::pos = Vec3(0.,0.,0.); 

KERNEL(pts) template<class T> void knSetPdataConst(ParticleDataImpl<T>& pdata, T value) { pdata[idx] = value; }

KERNEL(pts) template<class T, class S> void knPdataSet  (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] += other[idx]; }
KERNEL(pts) template<class T, class S> void knPdataAdd  (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] += other[idx]; }
KERNEL(pts) template<class T, class S> void knPdataSub  (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] -= other[idx]; }
KERNEL(pts) template<class T, class S> void knPdataMult (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] *= other[idx]; }
KERNEL(pts) template<class T, class S> void knPdataDiv  (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] /= other[idx]; }

KERNEL(pts) template<class T, class S> void knPdataSetScalar (ParticleDataImpl<T>& me, const S& other)  { me[idx]  = other; }
KERNEL(pts) template<class T, class S> void knPdataAddScalar (ParticleDataImpl<T>& me, const S& other)  { me[idx] += other; }
KERNEL(pts) template<class T, class S> void knPdataMultScalar(ParticleDataImpl<T>& me, const S& other)  { me[idx] *= other; }
KERNEL(pts) template<class T, class S> void knPdataScaledAdd (ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other, const S& factor) { me[idx] += factor * other[idx]; }

KERNEL(pts) template<class T> void knPdataSafeDiv (ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other) { me[idx] = safeDivide(me[idx], other[idx]); }
KERNEL(pts) template<class T> void knPdataSetConst(ParticleDataImpl<T>& pdata, T value) { pdata[idx] = value; }

KERNEL(pts) template<class T> void knPdataClamp (ParticleDataImpl<T>& me, T min, T max) { me[idx] = clamp( me[idx], min, max); }

// python operators


template<typename T>
ParticleDataImpl<T>& ParticleDataImpl<T>::copyFrom(const ParticleDataImpl<T>& a) {
	assertMsg (a.mData.size() == mData.size() , "different pdata size "<<a.mData.size()<<" vs "<<this->mData.size() );
	memcpy( &mData[0], &a.mData[0], sizeof(T) * mData.size() );
	return *this; 
}

template<typename T>
void ParticleDataImpl<T>::setConst(T s) {
	knPdataSetScalar<T,T> op( *this, s );
}

template<typename T>
void ParticleDataImpl<T>::add(const ParticleDataImpl<T>& a) {
	knPdataAdd<T,T> op( *this, a );
}
template<typename T>
void ParticleDataImpl<T>::sub(const ParticleDataImpl<T>& a) {
	knPdataSub<T,T> op( *this, a );
}

template<typename T>
void ParticleDataImpl<T>::addConst(T s) {
	knPdataAddScalar<T,T> op( *this, s );
}

template<typename T>
void ParticleDataImpl<T>::addScaled(const ParticleDataImpl<T>& a, const T& factor) {
	knPdataScaledAdd<T,T> op( *this, a, factor );
}

template<typename T>
void ParticleDataImpl<T>::mult( const ParticleDataImpl<T>& a) {
	knPdataMult<T,T> op( *this, a );
}

template<typename T>
void ParticleDataImpl<T>::multConst(T s) {
	knPdataMultScalar<T,T> op( *this, s );
}

template<typename T>
void ParticleDataImpl<T>::clamp(Real min, Real max) {
	knPdataClamp<T> op( *this, min,max );
}

template<typename T>
KERNEL(pts, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real CompPdata_Min(const ParticleDataImpl<T>& val) {
	if (val[idx] < minVal)
		minVal = val[idx];
}

template<typename T>
KERNEL(pts, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
Real CompPdata_Max(const ParticleDataImpl<T>& val) {
	if (val[idx] > maxVal)
		maxVal = val[idx];
}

template<typename T>
Real ParticleDataImpl<T>::getMinValue() {
	return sqrt(CompPdata_Min<T> (*this));
}

template<typename T>
Real ParticleDataImpl<T>::getMaxAbsValue() {
	Real amin = CompPdata_Min<T> (*this);
	Real amax = CompPdata_Max<T> (*this);
	return max( fabs(amin), fabs(amax));
}

template<typename T>
Real ParticleDataImpl<T>::getMaxValue() {
	return sqrt(CompPdata_Max<T> (*this));
} 

template<typename T>
void ParticleDataImpl<T>::printPdata(int start, int stop, bool printIndex)
{
	std::ostringstream sstr;
	int s = (start>0 ? start : 0                 );
	int e = (stop>0  ? stop  : (int)mData.size() );
	s = Manta::clamp(s, 0, (int)mData.size());
	e = Manta::clamp(e, 0, (int)mData.size());

	for(int i=s; i<e; ++i) {
		if(printIndex) sstr << i<<": ";
		sstr<<mData[i]<<" "<<"\n";
	} 
	debMsg( sstr.str() , 1 );
}

// specials for vec3


KERNEL(pts, reduce=min) returns(Real minVal=-std::numeric_limits<Real>::max())
Real CompPdata_MinVec3(const ParticleDataImpl<Vec3>& val) {
	const Real s = normSquare(val[idx]);
	if (s < minVal)
		minVal = s;
}

KERNEL(pts, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::min())
Real CompPdata_MaxVec3(const ParticleDataImpl<Vec3>& val) {
	const Real s = normSquare(val[idx]);
	if (s > maxVal)
		maxVal = s;
}

template<>
Real ParticleDataImpl<Vec3>::getMinValue() {
	return sqrt(CompPdata_MinVec3 (*this));
}

template<>
Real ParticleDataImpl<Vec3>::getMaxAbsValue() {
	Real amin = CompPdata_MinVec3 (*this);
	Real amax = CompPdata_MaxVec3 (*this);
	return max( fabs(amin), fabs(amax));
}

template<>
Real ParticleDataImpl<Vec3>::getMaxValue() {
	return sqrt(CompPdata_MaxVec3 (*this));
}


// explicit instantiation
template class ParticleDataImpl<int>;
template class ParticleDataImpl<Real>;
template class ParticleDataImpl<Vec3>;


} // namespace

