//
//  main.c
//  fastFourierTransform
//
//  Created by Beau Johnston on 19/09/11.
//  Copyright 2011 University Of New England. All rights reserved.
//

#include "global.h"
#include "RGBAUtilities.h"
#include "FileHandler.h"

/* Utility functions not needed on Vanilla Essence C kernel */
int width, height;
uint8* readImage(char* fileName){
    readPngFile(fileName);
    
    width = getImageWidth();
    height = getImageLength();
    
    uint8 * buffer = malloc(sizeof(uint8)*getImageSize());    
    memcpy(buffer, getImage(), getImageSize());
    return buffer;
}

void saveImage(char* fileName, uint8*buffer){
    setImage(buffer);
    writePngFile(fileName);
    return;
}

/* Utility functions essential for computation on the Vanilla Essence C kernel */
void FFT(short int dir,long m, float *x,float *y)
{
    long n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
    /* Calculate the number of points */
    n = 1;
    for (i=0;i<m;i++) 
        n *= 2;
    
    /* Do the bit reversal */
    i2 = n >> 1;
    j = 0;
    for (i=0;i<n-1;i++) {
        if (i < j) {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    
    /* Compute the FFT */
    c1 = -1.0; 
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0; 
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<n;i+=l2) {
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1; 
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1) 
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    
    /* Scaling for forward transform */
    if (dir == 1) {
        for (i=0;i<n;i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }
    
    return;
}

void ThreeDimensionalFFT(short int dir,long width, long height, long depth, float* real, float* imaginary){
    //create space for result array
    float resultReal[width*height*depth];
    float resultImag[width*height*depth];
    
    for(int z = 0; z < depth; z++){
        for(int y = 0; y < height; y++){
            float rowWiseReal[width];
            float rowWiseImag[width];
            
            //extract a tmp row
            for(int x = 0; x < width; x++){
                rowWiseReal[x] = real[z*(width*height)+y*(width)+x];
                rowWiseImag[x] = imaginary[z*(width*height)+y*(width)+x];
            }
            
            //apply the Fast Fourier Transform
            FFT(dir,width,rowWiseReal,rowWiseImag);
            
            //store the tmp result into result array
            for(int x = 0; x < width; x++){
                resultReal[z*(width*height)+y*(width)+x] = rowWiseReal[x];
                resultImag[z*(width*height)+y*(width)+x] = rowWiseImag[x];
            }
            
            //done with this row
        }
    }
    
    
}

int main (int argc, const char * argv[])
{
    //char* fileName="../../../../../dataResources/High-Res-Stage-24-Take-4/out.png";
    //char* outputImageFileName = "../../../../../dataResources/output/result.png";
    
    //a copy in executables neighbour used for debugging!
    char* fileName="dataSet/out.png";
    char* outputImageFileName = "dataSetOut/result.png";
    
    generateListOfAssociatedFiles(fileName);
        
    int depth = numberOfFiles();
    
    uint8*bigBuffer;
    bool firstRun = true;
    printf("Set up ... about to read image stack\n");
    system("pwd");
    printf("\n");
    
    //load all images into a buffer
    for (int i = 0; i < numberOfFiles(); i++) {
        char* tmp = getNextFileName();
        //printf("next file name from main is :%s \n", tmp);
        readPngFile(tmp);
        width = getImageWidth();
        height = getImageLength();
        uint8 *buffer = malloc(sizeof(uint8)*getImageSize());
        memcpy(buffer, getImage(), getImageSize());
        if (firstRun) {
            //if its the first run we don't know the dimensions of the image
            //and thus don't know how much memory to statically allocate
            bigBuffer = malloc(sizeof(uint8)*getImageSize()*depth);
            firstRun = false;
        }
        memcpy(bigBuffer+(i*getImageSize()), buffer, getImageSize());
        free(buffer);
    } 
    
//    float filtX[3] = {-1, 0, 1};
//    float filtY[3] = {-1, 0, 1};
//    float filtZ[3] = {-1, 0, 1};
//    
//    float dvfXYZ[3][3][3];
//    float dvfYZX[3][3][3];
//    float dvfZXY[3][3][3];
//    
//    for(int i = 0; i < 3; i++){
//        for(int j = 0; j < 3; j++){
//            for(int k = 0; k < 3; k++){
//                dvfXYZ[i][j][k] = - filtX[i] * exp(-((pow(filtX[i],2)+pow(filtY[j],2)+pow(filtZ[k],2))/2));
//                dvfYZX[i][j][k] = - filtY[j] * exp(-((pow(filtX[i],2)+pow(filtY[j],2)+pow(filtZ[k],2))/2));
//                dvfZXY[i][j][k] = - filtZ[k] * exp(-((pow(filtX[i],2)+pow(filtY[j],2)+pow(filtZ[k],2))/2));
//            }
//        }
//    }

    //save all images from buffer
    for (int i = 0; i < depth; i++) {
        uint8 *buffer = malloc(sizeof(uint8)*getImageSize());
        memcpy(buffer, bigBuffer+(i*getImageSize()), getImageSize());
        
        char* file = substring((int)strrchr(outputImageFileName, '/')+1, (int)strlen(outputImageFileName), outputImageFileName);
        
        char* path = substring(0, (int)strrchr(outputImageFileName, '/')+1, outputImageFileName);
        
        char* cutDownFile = substring(0, (int)strrchr(file, '.'), file);

        char* extension = substring((int)strrchr(file, '.'), (int)strlen(file),file);

        char* newName = cutDownFile;
        char numericalRepresentation[200];
        sprintf(numericalRepresentation, "%d", i);
        newName = strcat(newName, numericalRepresentation);
        newName = strcat(newName, extension);
        newName = strcat(path, newName);
        
        saveImage(newName, buffer);   
    }
    
    return 0;
}

