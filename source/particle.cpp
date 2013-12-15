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
#if NO_ZLIB!=1
#include <zlib.h>
#endif
#include "particle.h"
#include "levelset.h"
#include  <cstring>

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

PbClass* ParticleBase::create(PbType t, const string& name) {        
    _args.add("nocheck",true);
    if (t.str == "")
        errMsg("Specify particle data type to create");
	//debMsg( "Pdata creating '"<< t.str , 5 );
    
    PbClass* pyObj = PbClass::createPyObject(t.str, name, _args, this->getParent() );

	ParticleDataBase* pdata = dynamic_cast<ParticleDataBase*>(pyObj);
	if(!pdata) {
		errMsg("Unable to get particle data pointer from newly created object. Only create ParticleData type with a ParticleSys.creat() call, eg, PdataReal, PdataVec3 etc.");
		delete pyObj;
		return NULL;
	} else {
		this->registerPdata(pdata);
	}

	return pyObj;
}

void ParticleBase::registerPdata(ParticleDataBase* pdata) {
	pdata->setParticleSys(this);
	mPartData.push_back(pdata);

	if( pdata->getType() == ParticleDataBase::DATA_REAL ) {
		ParticleDataImpl<Real>* pd = dynamic_cast< ParticleDataImpl<Real>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as real!");
		this->registerPdataReal(pd);
	}
	else if( pdata->getType() == ParticleDataBase::DATA_VEC3 ) {
		ParticleDataImpl<Vec3>* pd = dynamic_cast< ParticleDataImpl<Vec3>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as vec3!");
		this->registerPdataVec3(pd);
	}
	else if( pdata->getType() == ParticleDataBase::DATA_INT ) {
		ParticleDataImpl<int>* pd = dynamic_cast< ParticleDataImpl<int>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as int!");
		this->registerPdataInt(pd);
	}
}
void ParticleBase::registerPdataReal(ParticleDataImpl<Real>* pd) { mPdataReal.push_back(pd); }
void ParticleBase::registerPdataVec3(ParticleDataImpl<Vec3>* pd) { mPdataVec3.push_back(pd); }
void ParticleBase::registerPdataInt (ParticleDataImpl<int >* pd) { mPdataInt .push_back(pd); }

void ParticleBase::addAllPdata() {
	for(int i=0; i<(int)mPartData.size(); ++i) {
		mPartData[i]->add();
	} 
}

std::string ParticleBase::debugInfoPdata()
{
	std::ostringstream sstr;
	sstr << "Particle system "<<mName<<" , size: "<< this->getSizeSlow() <<", data ";
	for(int i=0; i<(int)mPartData.size(); ++i) {
		sstr << i<<":" << mPartData[i]->size() <<" ";
	} 
	sstr << ".";
	return sstr.str();
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

//! in line with grid uni header
typedef struct {
    int dim; // number of partilces
    int elementType, bytesPerElement; // type id and byte size
	char info[256]; // mantaflow build information
    unsigned long timestamp; // creation time
} UniPartHeader;

template <class T>
void writeParticlesUni(const string& name, BasicParticleSystem* parts ) {
    cout << "writing particles " << parts->getName() << " to uni file " << name << endl;
    
#	if NO_ZLIB!=1
    char ID[5] = "PB01";
    UniPartHeader head;
	head.dim      = parts->size();
    head.bytesPerElement = sizeof(T);
    head.elementType = 0; // 0 for base data
	snprintf( head.info, 256, "%s", buildInfoString().c_str() );	
	MuTime stamp; stamp.get();
	head.timestamp = stamp.time;
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf) errMsg("can't open file");
    
    gzwrite(gzf, ID, 4);
    gzwrite(gzf, &head, sizeof(UniPartHeader));
    gzwrite(gzf, &(parts->getData()[0]), sizeof(T)*head.dim);
    gzclose(gzf);
#	else
    cout << "file format not supported without zlib" << endl;
#	endif
};

template <class T>
void readParticlesUni(const string& name, BasicParticleSystem* parts ) {
    cout << "reading particles " << parts->getName() << " from uni file " << name << endl;
    
#	if NO_ZLIB!=1
    gzFile gzf = gzopen(name.c_str(), "rb");
    if (!gzf) errMsg("can't open file");

    char ID[5]={0,0,0,0,0};
	gzread(gzf, ID, 4);
    
    if (!strcmp(ID, "PB01")) {
        // current file format
        UniPartHeader head;
        assertMsg (gzread(gzf, &head, sizeof(UniPartHeader)) == sizeof(UniPartHeader), "can't read file, no header present");
        assertMsg ( ((head.bytesPerElement == sizeof(T)) && (head.elementType==0) ), "particle type doesn't match");

		// re-allocate all data
		parts->resizeAll( head.dim );

        assertMsg (head.dim == parts->size() , "particle size doesn't match");
    	int bytes = sizeof(T)*head.dim;
        int readBytes = gzread(gzf, &(parts->getData()[0]), sizeof(T)*head.dim);
    	assertMsg(bytes==readBytes, "can't read uni file, stream length does not match, "<<bytes<<" vs "<<readBytes );
    }
    gzclose(gzf);
#	else
    cout << "file format not supported without zlib" << endl;
#	endif
};


void BasicParticleSystem::load(string name) {
    if (name.find_last_of('.') == string::npos)
        errMsg("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if ( ext == ".uni") 
        readParticlesUni<BasicParticleData>(name, this);
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
        writeParticlesUni<BasicParticleData>(name, this);
	// raw data formats, very basic for simple data transfer to other programs
    else if (ext == ".posgz") 
        this->writeParticlesRawPositionsGz(name);
    else if (ext == ".velgz") 
        this->writeParticlesRawVelocityGz(name);
    else
        errMsg("particle '" + name +"' filetype not supported for saving");
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
int ParticleDataImpl<T>::size() const {
	return mData.size();
}
template<class T>
void ParticleDataImpl<T>::add() {
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


template <class T>
void writePdataUni(const string& name, ParticleDataImpl<T>* pdata ) {
    cout << "writing particle data " << pdata->getName() << " to uni file " << name << endl;
    
#	if NO_ZLIB!=1
    char ID[5] = "PD01";
    UniPartHeader head;
	head.dim      = pdata->size();
    head.bytesPerElement = sizeof(T);
    head.elementType = 1; // 1 for particle data, todo - add sub types?
	snprintf( head.info, 256, "%s", buildInfoString().c_str() );	
	MuTime stamp; stamp.get();
	head.timestamp = stamp.time;
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf) errMsg("can't open file");
    
    gzwrite(gzf, ID, 4);
    gzwrite(gzf, &head, sizeof(UniPartHeader));
    gzwrite(gzf, &(pdata->get(0)), sizeof(T)*head.dim);
    gzclose(gzf);
#	else
    cout << "file format not supported without zlib" << endl;
#	endif
};

template <class T>
void readPdataUni(const string& name, ParticleDataImpl<T>* pdata ) {
    cout << "reading particle data " << pdata->getName() << " from uni file " << name << endl;
    
#	if NO_ZLIB!=1
    gzFile gzf = gzopen(name.c_str(), "rb");
    if (!gzf) errMsg("can't open file");

    char ID[5]={0,0,0,0,0};
	gzread(gzf, ID, 4);
    
    if (!strcmp(ID, "PD01")) {
        UniPartHeader head;
        assertMsg (gzread(gzf, &head, sizeof(UniPartHeader)) == sizeof(UniPartHeader), "can't read file, no header present");
        assertMsg ( ((head.bytesPerElement == sizeof(T)) && (head.elementType==1) ), "pdata type doesn't match");
        assertMsg (head.dim == pdata->size() , "pdata size doesn't match");
    	int bytes = sizeof(T)*head.dim;
        int readBytes = gzread(gzf, &(pdata->get(0)), sizeof(T)*head.dim);
    	assertMsg(bytes==readBytes, "can't read uni file, stream length does not match, "<<bytes<<" vs "<<readBytes );
    }
    gzclose(gzf);
#	else
    cout << "file format not supported without zlib" << endl;
#	endif
}

template<typename T>
void ParticleDataImpl<T>::load(string name) {
    if (name.find_last_of('.') == string::npos)
        errMsg("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if ( ext == ".uni") 
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
    else
        errMsg("particle data '" + name +"' filetype not supported for saving");
}

// specializations

template<>
ParticleDataBase::PdataType ParticleDataImpl<int>::getType() const {
	return ParticleDataBase::DATA_INT;
} 
template<>
ParticleDataBase::PdataType ParticleDataImpl<Real>::getType() const {
	return ParticleDataBase::DATA_REAL;
} 
template<>
ParticleDataBase::PdataType ParticleDataImpl<Vec3>::getType() const {
	return ParticleDataBase::DATA_VEC3;
}

// note, we need a flag value for functions such as advection
// ideally, this value should never be modified
int ParticleIndexData::flag = 0; 

// explicit instantiation
template class ParticleDataImpl<int>;
template class ParticleDataImpl<Real>;
template class ParticleDataImpl<Vec3>;

KERNEL(pts) template<class T>
void knSetPdataConst(ParticleDataImpl<T>& pdata, T value) {
	pdata[i] = value;
}
PYTHON void setConstPdata    (ParticleDataImpl<Real>& pd, Real value=0.) { knSetPdataConst<Real>(pd,value); }
PYTHON void setConstPdataVec3(ParticleDataImpl<Vec3>& pd, Vec3 value=0.) { knSetPdataConst<Vec3>(pd,value); }
PYTHON void setConstPdataInt (ParticleDataImpl<int >& pd, int  value=0.) { knSetPdataConst<int> (pd,value); }



} // namespace

