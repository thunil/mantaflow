/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * CUDA curl noise evaluation
 *
 ******************************************************************************/

#ifndef _CURLNOISE_H
#define _CURLNOISE_H

#include "cudatools.h"

namespace Manta {
	
const int NOISE_SIZE = 128;
typedef texture<short4, 3, cudaReadModeNormalizedFloat> CompressedTex;

//! Device interface for CUDA noise texture
struct CudaNoiseDev {
	CudaNoiseDev(float s) : scale(s) {};
	float scale;
	
	__device__ inline float3 eval(CompressedTex tex, const float3& p);
	__device__ inline float3 synthesizeK41(CompressedTex tex, const float3& p, int octaves, float L0);
};

//! Host interface for CUDA noise texture
class CudaNoiseTexture {
public:
	CudaNoiseTexture() { init(); }
	~CudaNoiseTexture();
	void init();
	CudaNoiseDev bind(CompressedTex& tex);
	
//private:
	cudaArray* mNoiseArray;
	float mScale;
};




// *****************************************************************************
// Implementation
// *****************************************************************************

template<char X, int N, int M>
__device__ inline float derivativeWeight(const float t) 
{
	if ((X == 'x' && N == 0) || (X == 'y' && N == 1) || (X == 'z' && N == 2)) {
		if (M == 0)
			return t;
		else if (M == 1)
			return 1.0f - 2.0f * t;
		else
			return t - 1.0f;
	} else {
		if (M == 0)
			return t * t * 0.5f;
		else if (M == 1)
			return 0.5f + t - t*t;
		else
			return 0.5f - t + 0.5f * t * t;        
	}
}

__device__ inline float3 CudaNoiseDev::eval(CompressedTex noise, const float3& p) {
	float cache[12];
	float ext = 1.0f / (float) NOISE_SIZE;

	float t0, t1, t2;
	float midX, midY, midZ;

	midX = ceil(p.x - 0.5f);
	t0   = midX - (p.x - 0.5f);
	midX -= 1.0f;
	midX *= ext;

	midY = ceil(p.y - 0.5f);
	t1   = midY - (p.y - 0.5f);
	midY -= 1.0f;
	midY *= ext;

	midZ = ceil(p.z - 0.5f);
	t2   = midZ - (p.z - 0.5f);
	midZ -= 1.0f;
	midZ *= ext;

	cache[0] = derivativeWeight<'x', 2, 0>(t2);
	cache[1] = derivativeWeight<'x', 2, 1>(t2);
	cache[2] = derivativeWeight<'x', 2, 2>(t2);
	
	cache[3] = derivativeWeight<'z', 2, 0>(t2);
	cache[4] = derivativeWeight<'z', 2, 1>(t2);
	cache[5] = derivativeWeight<'z', 2, 2>(t2);
	
	cache[6] = derivativeWeight<'x', 1, 0>(t1);
	cache[7] = derivativeWeight<'x', 1, 1>(t1);
	cache[8] = derivativeWeight<'x', 1, 2>(t1);

	cache[9] = derivativeWeight<'y', 1, 0>(t1);
	cache[10] = derivativeWeight<'y', 1, 1>(t1);
	cache[11] = derivativeWeight<'y', 1, 2>(t1);

	float3 v;
	v.x = 0.0f;
	v.y = 0.0f;
	v.z = 0.0f;

	///////////////////////////////////////////////////////////////////////////////////////
	// x, y, z derivatives
	///////////////////////////////////////////////////////////////////////////////////////
	for (int k = 0; k < 3; k++, midZ += ext)
	{ 
		for (int j = 0; j < 3; j++, midY += ext)
		{
			for (int i = 0; i < 3; i++, midX += ext)
			{
				// Read the noise texture
				float4 n = tex3D(noise, midX, midY, midZ);
				
				float w2_x = cache[0 + k];
				float w2_y = w2_x;
				float w2_z = cache[3 + k];
				
				float w1_x = cache[6 + j];
				float w1_y = cache[9 + j];
				float w1_z = w1_x;

				float w0_x = (i == 0) ? derivativeWeight<'x', 0, 0>(t0) :
							 (i == 1) ? derivativeWeight<'x', 0, 1>(t0) :
										derivativeWeight<'x', 0, 2>(t0);

				float w0_y = (i == 0) ? derivativeWeight<'y', 0, 0>(t0) :
							 (i == 1) ? derivativeWeight<'y', 0, 1>(t0) :
										derivativeWeight<'y', 0, 2>(t0);
				float w0_z = w0_y;

				// Decompress noise
				n.x *= scale;
				n.y *= scale;
				n.z *= scale;

				// Calculate the final weights
				float w_x = w0_x * w1_x * w2_x;
				float w_y = w0_y * w1_y * w2_y;
				float w_z = w0_z * w1_z * w2_z;

				// Add the weighted noise
				v.z += n.y * w_x;
				v.y -= n.z * w_x;

				v.z -= n.x * w_y;
				v.x += n.z * w_y;

				v.y += n.x * w_z;
				v.x -= n.y * w_z;
			}   
			midX -= 3.0f * ext;
		}   
		midY -= 3.0f * ext;
	}

	return v;
}

const __device__ float PERSISTENCE = 0.56123f;

__device__ inline float3 CudaNoiseDev::synthesizeK41(CompressedTex tex, const float3& p, int octaves, float L0) {
	// Base parameters for octave 0
	float multiplier = L0;
	float amplitude = 1.0f;
	float3 vel = make_float3(0,0,0);

	for (int octave = 0; octave < octaves; octave++) {
		float3 noiseLookup = eval(tex, multiplier * p);        
		vel += noiseLookup * amplitude;
		
		// next scale
		amplitude *= PERSISTENCE; 
		multiplier *= 2.0f;
	}
	return vel;
}



} // namespace

#endif
