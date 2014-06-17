/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Some helper functions for using CUDA
 *
 ******************************************************************************/

#ifndef _CUDATOOLS_H
#define _CUDATOOLS_H

#include <cuda.h>
#include <vector>
#include <thrust/device_vector.h>
#include "vectorbase.h"

// *****************************************************************
// extending float3 type

__device__ inline float3 operator+(const float3 &a, const float3 &b) {
  return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}
__device__ inline float3 operator*(const float3 &a, const float b) {
  return make_float3(a.x*b, a.y*b, a.z*b);
}
__device__ inline float3 operator/(const float3 &a, const float b) {
  return make_float3(a.x/b, a.y/b, a.z/b);
}
__device__ inline float3 operator*(const float b, const float3 &a) {
  return make_float3(a.x*b, a.y*b, a.z*b);
}
__device__ inline float3 operator-(const float3 &a, const float3 &b) {
  return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
}
__device__ inline void operator+=(float3 &a, const float3 &b) {
  a.x += b.x; a.y += b.y; a.z += b.z; 
}
__device__ inline float normSqr(const float3& v) {
	return (v.x*v.x)+(v.y*v.y)+(v.z*v.z);
}
__device__ inline float norm(const float3& v) {
	return sqrtf((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
}
__device__ inline float normalize(float3& v) {
	float n = norm(v);
	if (n==0)
		v = make_float3(0,0,0);
	else {
		float in = 1.0f/n;
		v = make_float3(v.x * in, v.y * in, v.z * in);
	}
	return n;
}
__device__ inline float3 cross(const float3& a, const float3& b) {
	return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
__device__ inline float dot(const float3&a, const float3& b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}


namespace Manta {

// *****************************************************************
// some general macros

inline void cudaAssert(cudaError_t er) {
	if (er != cudaSuccess) {
		throw Error("CUDA assert error: " + std::string(cudaGetErrorString(er)));
	}
}
// *****************************************************************
// define some helpers for unpacking Vec3s etc.

//! Device object to access data provided by Vec3AArray
struct CVec3Ptr {
	float *x, *y, *z; 
	__device__ inline float3 get(int i) const { return make_float3(x[i],y[i],z[i]); };
	__device__ inline void set(int i, const float3& v) { x[i]=v.x; y[i]=v.y; z[i]=v.z; };
};

//! Provide host and data for easy coalescing by the GPU.
//! Has the same methods as CVector
struct CVec3Array {    
	
	CVec3Array(int sz) {
		x.resize(sz);
		y.resize(sz);
		z.resize(sz);        
	}    
	CVec3Array(const std::vector<Vec3>& v) {
		x.resize(v.size());
		y.resize(v.size());
		z.resize(v.size());
		for (size_t i=0; i<v.size(); i++) {
			x[i] = v[i].x;
			y[i] = v[i].y;
			z[i] = v[i].z;
		}
		upload();
	}
	void upload() {
		dx = x; 
		dy = y;
		dz = z;
	}
	void download() {
		x = dx;
		y = dy;
		z = dz;
	}
	void downloadTo(std::vector<Vec3>& v) {
		download();
		if (v.size() != x.size())
			v.resize(x.size());
		for (size_t i=0; i<v.size(); i++)
			v[i] = Vec3(x[i],y[i],z[i]);
	}
	
	
	CVec3Ptr data() {
		CVec3Ptr a = { thrust::raw_pointer_cast(dx.data()), thrust::raw_pointer_cast(dy.data()), thrust::raw_pointer_cast(dz.data()) };
		return a;
	}
	inline const Vec3 operator[](int idx) const { return Vec3((Real)x[idx], (Real)y[idx], (Real)z[idx]); }
	inline void set(int idx, const Vec3& v) { x[idx] = v.x; y[idx] = v.y; z[idx] = v.z; }
		
	inline int size() { return x.size(); }    
	inline int blocks(int blockSize) { return (x.size() - 1) / blockSize + 1; }
	
	thrust::host_vector<float> x, y, z;
	thrust::device_vector<float> dx, dy, dz;
};

//! wrapper around thrust device vector for easier access in CUDA
template<class T>
struct CArray {
	CArray(const std::vector<T>& v) : 
		dev(v.begin(), v.end()) 
	{
	}    
	CArray(int sz) {
		dev.resize(sz);
		host.resize(sz);
	}    
	T* data() {
		return thrust::raw_pointer_cast(dev.data());
	}
	void upload() {
		if (host.size() == dev.size())
			dev = host;
	}
	void download() {
		host = dev;
	}
	void downloadTo(std::vector<T>& v) {
		thrust::copy(dev.begin(), dev.end(), v.begin());
	}
	inline T& operator[](int idx) { return host[idx]; }
	inline void set(int idx, T val) { host[idx] = val; }
	
	int size() { return dev.size(); }
	int blocks(int blockSize) { return (dev.size()-1) / blockSize + 1; }
	
	thrust::host_vector<T> host;
	thrust::device_vector<T> dev;
};
  
} // namespace

#endif