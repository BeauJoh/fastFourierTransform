//
//  global.h
//  fastFourierTransform
//
//  Created by Beau Johnston on 19/09/11.
//  Copyright 2011 University Of New England. All rights reserved.
//

#ifndef fastFourierTransform_global_h
#define fastFourierTransform_global_h

#include <stdio.h>
#include <math.h>

#define uint8 unsigned char

void FFT(short int dir,long m, float *x,float *y);
void ThreeDimensionalFFT(short int dir,long width, long height, long depth, float* real, float* imaginary);
uint8* readImage(char* fileName);
void saveImage(char* fileName, uint8*buffer);

#endif
