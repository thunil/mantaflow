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

} // namespace