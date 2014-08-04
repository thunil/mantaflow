/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Simple image IO
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "simpleimage.h"

namespace Manta {

// write rectangle to ppm
bool SimpleImage::writePpm(std::string filename, int minx, int miny, int maxx, int maxy, bool invertXY) {
	int w = maxx-minx;
	int h = maxy-miny;

	if (w<=0 || h<=0 || w>mSize[0] || h>mSize[1]) {
		errMsg("SimpleImage::WritePPM Invalid rect: w="<<w<<", h="<<h<<", size="<<mSize[0]<<","<<mSize[1]<<" min/max: "<<minx<<","<<miny<<" to "<<maxx<<","<<maxy<<", resetting... " );
		minx = miny = 0;
		maxx = mSize[0]-1; 
		maxy = mSize[1]-1;
		w = mSize[0]-1;
		h = mSize[1]-1;
	}

	FILE* fp = fopen(filename.c_str(), "wb");
	if (fp==NULL) { 
		errMsg("SimpleImage::WritePPM Unable to open '"<<filename<<"' for writing");
		return false; 
	}
	fprintf(fp, "P6\n%d %d\n255\n", w, h);

	int pixCnt = 0;
	for(int j=maxy-1; j>=miny; j--) 
		for(int i=minx; i<maxx; i++) {
			unsigned char col[3];
			for(int l=0; l<3; l++) {
				float val;
				if(invertXY) val = (float)get(j,i)[l];
				else val =(float)get(i,j)[l];

				val = clamp(val, (float)0.,(float)1.);
				col[l] = (unsigned char)(255. * val );
			}
			//col[1] = col[2] = col[0];
			//if (fwrite(col,1,3, fp) != 3) errMsg("SimpleImage::writePpm fwrite failed");
			fwrite(col,1,3, fp);
			pixCnt++;
			//fprintf(stderr,"%d %d %d \n",col[0],i,j);
		}

	fclose(fp);
	debMsg("WritePPM Wrote '"<<filename<<"', region="<<minx<<","<<miny<<" to "<<maxx<<","<<maxy<<"; "<<pixCnt, 1);

	return true;
}

bool SimpleImage::writePpm(std::string filename) {
	return writePpm(filename, 0,0, getSize()[0], getSize()[1] );
}


// read in a ppm file, and init the image accordingly
bool SimpleImage::initFromPpm (std::string filename) { 
	// maximum length of a line of text 
	const int MAXLINE=1024;

	int filetype=0;
	enum {PGM, PPM};	// possible file types 

	FILE *fp;
	char line[MAXLINE];
	int size, rowsize;

	// Read in file type  
	fp = fopen(filename.c_str(), "rb");
	if(!fp) {
		if(mAbortOnError) debMsg("SimpleImage Error - unable to open file '"<< filename<<"' for reading", 1);
		return 0;
	}

	// 1st line: PPM or PGM 
	if (fgets (line, MAXLINE, fp) == NULL) {
		if(mAbortOnError) errMsg("SimpleImage::initFromPpm fgets failed");
		return 0;
	}
	
	if (line[1] == '5')
		filetype = PGM;
	else if (line[1] == '6')
		filetype = PPM;
	else {
		if(mAbortOnError) debMsg("SimpleImage Error: need PPM or PGM file as input!", 1);
		return 0;
	}

	// Read in width and height, & allocate space  
	 // 2nd line: width height 
	if (fgets (line, MAXLINE, fp) == NULL) {
		if(mAbortOnError) errMsg("SimpleImage::initFromPpm fgets failed");
		return 0;
	}
	int windW=0, windH=0;	// size of the window on the screen  
	int intsFound = sscanf(line, "%d %d", &windW, &windH);
	if( intsFound == 1) {
		// only X found, search on next line as well for Y...
		if( sscanf(line, "%d", &windH) != 1) {
			if(mAbortOnError) errMsg("initFromPpm Ppm dimensions not found!"<<windW<<","<<windH);
			return 0; 
		} else {
			// ok, found 2 lines
			//debMsg("initFromPpm Ppm dimensions found!"<<windW<<","<<windH, 1);
		}
	} else if( intsFound == 2) {
		// ok!
	} else {
		if(mAbortOnError) errMsg("initFromPpm Ppm dimensions not found at all!"<<windW<<","<<windH);
		return 0;
	}

	if (filetype == PGM) {
		size = windH * windW;		// greymap: 1 byte per pixel
		rowsize = windW;
	} else {
		// filetype == PPM  
		size = windH * windW * 3;	// pixmap: 3 bytes per pixel
		rowsize = windW * 3;
	}

	unsigned char *pic = new unsigned char[size]; // (GLubyte *)malloc (size);

	// Read in maximum value (ignore) , could be scanned with sscanf as well, but this should be 255...
	// 3rd line
	if (fgets (line, MAXLINE, fp) == NULL) {
		if(mAbortOnError) errMsg("SimpleImage::initFromPpm fgets failed");
		return 0;
	}

	// Read in the pixel array row-by-row: 1st row = top scanline */
	unsigned char *ptr = NULL;
	ptr = &pic[(windH-1) * rowsize];
	for (int i = windH; i > 0; i--) {
		fread((void *)ptr, 1, rowsize, fp);
		ptr -= rowsize;
	}

	// init image
	this->init(windW, windH);
	if (filetype == PGM) {
		// grayscale
		for(int i=0; i<windW; i++) {
			for(int j=0; j<windH; j++) {
				double r = (double)pic[(j*windW+i)*1+0] / 255.;
				(*this)(i,j) = Vec3(r,r,r);
			}
		}
	} else {
		// convert grid to RGB vec's
		for(int i=0; i<windW; i++) {
			for(int j=0; j<windH; j++) {
				// return mpData[y*mSize[0]+x];
				double r = (double)pic[(j*windW+i)*3+0] / 255.;
				double g = (double)pic[(j*windW+i)*3+1] / 255.;
				double b = (double)pic[(j*windW+i)*3+2] / 255.;
				
				//(*this)(i,j) = Vec3(r,g,b);

				// RGB values have to be rotated to get the right colors!?
				// this might also be an artifact of photoshop export...?
				(*this)(i,j) = Vec3(g,b,r);
			}
		}
	}

	delete [] pic;
	fclose(fp);
	return 1;
}


// check index is valid
bool SimpleImage::indexIsValid(int i, int j) 
{
	if(i<0) return false;
	if(j<0) return false;
	if(i>=mSize[0]) return false;
	if(j>=mSize[1]) return false;
	return true;
}

}; // manta

