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
 
#include "curlnoise.h"
#include "cudatools.h"
#include "noisefield.h"
#include "manta.h"

using namespace std;

namespace Manta {


void CudaNoiseTexture::init() {
	cudaExtent noiseSize = {NOISE_SIZE, NOISE_SIZE, NOISE_SIZE};
	size_t memSize = noiseSize.width*noiseSize.height*noiseSize.depth;

	// register channels
	cudaChannelFormatDesc noiseChannelDesc = cudaCreateChannelDesc(16, 16, 16, 16, cudaChannelFormatKindSigned);
	cudaAssert(cudaMalloc3DArray(&mNoiseArray, &noiseChannelDesc, noiseSize));

	// create noise objects
	WaveletNoiseField nx(NULL), ny(NULL), nz(NULL);
	
	// find scales
	float maxval = -1e20f;
	for (size_t i = 0; i < memSize; i++)
		maxval = std::max(std::max(std::max(fabs(nx.data()[i]), fabs(ny.data()[i])), fabs(nz.data()[i])), maxval);
		
	float scale = 65535.0f / (2.0 * maxval);
	mScale = maxval;
	
	// alloc local buffer
	short4* buffer = new short4[memSize];
	for (size_t i = 0; i < memSize; i++) {
		buffer[i].x = (signed short) (nx.data()[i] * scale);
		buffer[i].y = (signed short) (ny.data()[i] * scale);
		buffer[i].z = (signed short) (nz.data()[i] * scale);
		buffer[i].w = 0;
	}
	
	// copy to page-locked mem
	cudaPitchedPtr page_locked_ptr;
	page_locked_ptr.pitch = noiseSize.width*sizeof(short4);
	page_locked_ptr.xsize = noiseSize.width;
	page_locked_ptr.ysize = noiseSize.height;
	page_locked_ptr.ptr = buffer;
	
	// copy data to 3D array
	cudaMemcpy3DParms copy_params = {0};
	copy_params.srcPtr   = page_locked_ptr;
	copy_params.dstArray = mNoiseArray;
	copy_params.extent   = noiseSize;
	copy_params.kind     = cudaMemcpyHostToDevice;
	cudaAssert(cudaMemcpy3D(&copy_params));

	//FILE* fp=fopen("/tmp/a.bin","wb"); fwrite(nx.data(), sizeof(float), 128*128*128, fp); fclose(fp);
	delete[] buffer;
}

CudaNoiseDev CudaNoiseTexture::bind(CompressedTex& tex) {
	// Bind the texture
	tex.normalized = true;
	tex.filterMode = cudaFilterModeLinear; // this kills high frequencies ! 
	tex.addressMode[0] = cudaAddressModeWrap;
	tex.addressMode[1] = cudaAddressModeWrap;
	tex.addressMode[2] = cudaAddressModeWrap;
	cudaAssert(cudaBindTextureToArray(tex, mNoiseArray));
	
	return CudaNoiseDev(mScale);    
}

CudaNoiseTexture::~CudaNoiseTexture() {
	cudaFreeArray(mNoiseArray);
}

/*
CompressedTex noise;

__global__ void rewriteNoise(CudaNoiseDev nd, float* a) {
	int3 cell = make_int3(threadIdx.x + blockDim.x*blockIdx.x, threadIdx.y + blockDim.y*blockIdx.y, threadIdx.z + blockDim.z*blockIdx.z);
	float scale = 1.0f/128.0f;
	float3 p = make_float3((float)(cell.x+0.5f) * scale, (float)(cell.y+0.5f) * scale, (float)(cell.z+0.5f) * scale);
	float4 n = tex3D(noise, p.x, p.y, p.z);
	a[cell.x + 128*cell.y + 128*128*cell.z] = (float)n.x * nd.scale;
}

PYTHON void testme() {
	CudaNoiseTexture nd;
	CArray<float> a(128*128*128);
	a.upload();
	
	dim3 blocksize(8, 8, 8);
	dim3 blocks(128/8, 128/8, 128/8);
	rewriteNoise<<<blocks, blocksize>>>(nd.bind(noise), a.data());
	a.download();
	FILE* fp=fopen("/tmp/re.bin","wb"); fwrite(&a[0], sizeof(float), 128*128*128, fp); fclose(fp);
}*/
	
} //namespace