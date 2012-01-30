/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Basic vector class
 *
 ******************************************************************************/

#include "vectorbase.h"
#include <limits>
#include <math.h>

using namespace std;

namespace Manta {

template<> const Vector3D<int> Vector3D<int>::Zero( 0, 0, 0 );
template<> const Vector3D<float> Vector3D<float>::Zero( 0.f, 0.f, 0.f );
template<> const Vector3D<double> Vector3D<double>::Zero( 0., 0., 0. );
template<> const Vector3D<float> Vector3D<float>::Invalid( numeric_limits<float>::quiet_NaN(), numeric_limits<float>::quiet_NaN(), numeric_limits<float>::quiet_NaN() );
template<> const Vector3D<double> Vector3D<double>::Invalid( numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN() );
//template<> const Vector3D<int> Vector3D<int>::Invalid( -1, -1, -1 );
template<> bool Vector3D<float>::isValid() const { return !isnan(x) && !isnan(y) && !isnan(z); }
template<> bool Vector3D<double>::isValid() const { return !isnan(x) && !isnan(y) && !isnan(z); }
//template<> bool Vector3D<int>::isValid() const { return x!=-1 || y!=-1 || z!=-1; }




}