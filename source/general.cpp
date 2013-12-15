/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Globally used macros and functions
 *
 ******************************************************************************/

#include "general.h"
#ifdef WIN32
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#  undef WIN32_LEAN_AND_MEAN
#  undef NOMINMAX
#else
#   include <sys/time.h>
#endif

using namespace std;

namespace Manta {

int gDebugLevel = 1;
 
void MuTime::get() {    
#ifdef WIN32
    LARGE_INTEGER liTimerFrequency;
    QueryPerformanceFrequency(&liTimerFrequency);
    LARGE_INTEGER liLastTime;
    QueryPerformanceCounter(&liLastTime);
    time = (INT)( ((double)liLastTime.QuadPart / liTimerFrequency.QuadPart)*1000 );
#else
    struct timeval tv;
    struct timezone tz;
    tz.tz_minuteswest = 0;
    tz.tz_dsttime = 0;
    gettimeofday(&tv,&tz);
    time = (tv.tv_sec*1000)+(tv.tv_usec/1000);
#endif    
}

MuTime MuTime::update() {
    MuTime o = *this;
    get();
    return *this - o;
}

string MuTime::toString() {
    stringstream ss;
    ss << *this;
    return ss.str();
}

ostream& operator<<(ostream& os, const MuTime& t) {
    unsigned long ms = (unsigned long)(   (double)t.time / (60.0*1000.0)  );
    unsigned long ss = (unsigned long)(  ((double)t.time / 1000.0) - ((double)ms*60.0)  );
    int      ps = (int)(       ((double)t.time - (double)ss*1000.0)/1.0 );

    if(ms>0) {
        os << ms<<"m"<< ss<<"s" ;
    } else {
        if(ps>0) {
            os << ss<<".";
            if(ps<10) { os <<"0"; }
            if(ps<100) { os <<"0"; }
            os <<ps<<"s" ;
        } else {
            os << ss<<"s" ;
        }
    }
    return os;
}

std::string buildInfoString() {
	std::ostringstream infoStr;
	infoStr << "mantaflow";
	// TODO , include hg branch info

	// os
#	ifdef WIN32
		infoStr << " win";
#	endif
#	ifdef __APPLE__
		infoStr << " mac";
#	endif
#	ifdef LINUX
		infoStr << " linux";
#	endif

	// 32/64
#	ifdef USE64
		infoStr << " 64bit";
#	else
		infoStr << " 32bit";
#	endif

	// fp precision
#	if FLOATINGPOINT_PRECISION==2
		infoStr << " fp2";
#	else
		infoStr << " fp1";
#	endif

	// other compile switches
#	ifdef DEBUG
		infoStr << " debug";
#	endif 
#	ifdef OPENMP
		infoStr << " omp";
#	endif

	infoStr << " from "<< __DATE__<<", "<<__TIME__;
	return infoStr.str();
}

} // namespace
