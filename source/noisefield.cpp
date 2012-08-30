/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Noise field
 *
 ******************************************************************************/

#include "noisefield.h"
#include "randomstream.h"

using namespace std;

//*****************************************************************************
// Wavelet sampling

static float _aCoeffs[32] = {
        0.000334f,-0.001528f, 0.000410f, 0.003545f,-0.000938f,-0.008233f, 0.002172f, 0.019120f,
        -0.005040f,-0.044412f, 0.011655f, 0.103311f,-0.025936f,-0.243780f, 0.033979f, 0.655340f,
        0.655340f, 0.033979f,-0.243780f,-0.025936f, 0.103311f, 0.011655f,-0.044412f,-0.005040f,
        0.019120f, 0.002172f,-0.008233f,-0.000938f, 0.003546f, 0.000410f,-0.001528f, 0.000334f};

namespace Manta {

int WaveletNoiseField::seed = 13322223;

void WaveletNoiseField::Downsample(float *from, float *to, int n, int stride){
    const float *a = &_aCoeffs[16];
    for (int i = 0; i < n / 2; i++) {
        to[i * stride] = 0;
        for (int k = 2 * i - 16; k <= 2 * i + 16; k++)
            to[i * stride] += a[k - 2 * i] * from[modFast128(k) * stride];
    }
}

static float _pCoeffs[4] = {0.25f, 0.75f, 0.75f, 0.25f};
void WaveletNoiseField::Upsample(float *from, float *to, int n, int stride) {
    const float *p = &_pCoeffs[2];

    for (int i = 0; i < n; i++) {
        to[i * stride] = 0;
        for (int k = i / 2; k <= i / 2 + 1; k++)
            to[i * stride] += p[i - 2 * k] * from[modSlow(k, n / 2) * stride];
    }
}

WaveletNoiseField::WaveletNoiseField(FluidSolver* parent) :
    PbClass(parent), mPosOffset(0.), mPosScale(1.), mValOffset(0.), mValScale(1.), mClamp(false), 
    mClampNeg(0), mClampPos(1), mTimeAnim(0), mGsInvX(0), mGsInvY(0), mGsInvZ(0)
{
    assertMsg(parent->is3D(), "Only works for 3D solvers for now");
    mGsInvX = 1.0/(parent->getGridSize().x);
    mGsInvY = 1.0/(parent->getGridSize().y);
    mGsInvZ = 1.0/(parent->getGridSize().z);
    generateTile();
};

string WaveletNoiseField::toString() {
    std::ostringstream out;
    out <<  "NoiseField: name '"<<mName<<"' "<<
        "  pos off="<<mPosOffset<<" scale="<<mPosScale<<
        "  val off="<<mValOffset<<" scale="<<mValScale<<
        "  clamp ="<<mClamp<<" val="<<mClampNeg<<" to "<<mClampPos<<
        "  timeAni ="<<mTimeAnim<<
        "  gridInv ="<<Vec3(mGsInvX,mGsInvY,mGsInvZ) ;
    return out.str();
}

void WaveletNoiseField::generateTile() {
    // generate tile
    const int n = NOISE_TILE_SIZE;
    const int n3 = n*n*n;
    cout << "generating " << n << "^3 noise tile" << endl;
    
    RandomStream randStream (seed);
    seed += 123;
    
    float *temp13 = new float[n3];
    float *temp23 = new float[n3];
    float *noise3 = new float[n3];

    // initialize
    for (int i = 0; i < n3; i++) {
        temp13[i] = temp23[i] =
            noise3[i] = 0.;
    }

    // Step 1. Fill the tile with random numbers in the range -1 to 1.
    for (int i = 0; i < n3; i++) {
        //noise3[i] = (randStream.getFloat() + randStream2.getFloat()) -1.; // produces repeated values??
        noise3[i] = randStream.getRandNorm(0,1);
    }

    // Steps 2 and 3. Downsample and upsample the tile
    for (int iy = 0; iy < n; iy++) 
        for (int iz = 0; iz < n; iz++) {
            const int i = iy * n + iz*n*n;
            Downsample(&noise3[i], &temp13[i], n, 1);
            Upsample  (&temp13[i], &temp23[i], n, 1);
        }
    for (int ix = 0; ix < n; ix++) 
        for (int iz = 0; iz < n; iz++) {
            const int i = ix + iz*n*n;
            Downsample(&temp23[i], &temp13[i], n, n);
            Upsample  (&temp13[i], &temp23[i], n, n);
        }
    for (int ix = 0; ix < n; ix++) 
        for (int iy = 0; iy < n; iy++) {
            const int i = ix + iy*n;
            Downsample(&temp23[i], &temp13[i], n, n*n);
            Upsample  (&temp13[i], &temp23[i], n, n*n);
        }

    // Step 4. Subtract out the coarse-scale contribution
    for (int i = 0; i < n3; i++) { 
        noise3[i] -= temp23[i];
    }

    // Avoid even/odd variance difference by adding odd-offset version of noise to itself.
    int offset = n / 2;
    if (offset % 2 == 0) offset++;

    if (n != 128) errMsg("WaveletNoise::Fast 128 mod used, change for non-128 resolution");
    
    int icnt=0;
    for (int ix = 0; ix < n; ix++)
        for (int iy = 0; iy < n; iy++)
            for (int iz = 0; iz < n; iz++) { 
                temp13[icnt] = noise3[modFast128(ix+offset) + modFast128(iy+offset)*n + modFast128(iz+offset)*n*n];
                icnt++;
            }

    for (int i = 0; i < n3; i++) {
        noise3[i] += temp13[i];
    }
    
    mNoiseTile = noise3;
    delete[] temp13;
    delete[] temp23;
    
    //FILE* fp = fopen("/tmp/bla.bin","wb"); fwrite(noise3, sizeof(float), n3, fp); fclose(fp);
}


    
}