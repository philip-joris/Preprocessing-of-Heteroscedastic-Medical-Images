
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>

// -------------------------------------------------- //
// -- auxiliary function: calculate labelling cost -- //
// -------------------------------------------------- //

// Full labelling cost
double calculateCost(int* counts, int* labels, int bins, double* D)
{
    double cost = 0;
    for(int i = 0; i < bins; ++i)
    {
        for(int j = 0; j < bins; ++j)
        {
            if(labels[i]!=labels[j]) continue; // Not same cluster
            cost += D[i*bins + j]*counts[i]*counts[j];
        }
    }
    return(cost);
}

// Single labelling cost
double calculateSingleCost(int* counts, int*labels, int bins, double* D, int index)
{
    double cost = 0;
    for(int j = 0; j < bins; ++j)
    {
        if(labels[index]!=labels[j]) continue; // Not same cluster
        cost += D[index*bins + j]*counts[index]*counts[j];
    }
    return(cost);
}

// -------------------------------------------- //
// -- main method: estimates the images ISDs -- //
// -------------------------------------------- //

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check the inputs and outputs
    if(nrhs != 3)   mexErrMsgTxt("Expected 3 inputs (image,distanceMatrix,nrOfClusters)!");
    if(nlhs != 1)    mexErrMsgTxt("Expected exactly one output argument!");
    int imgNDims = mxGetNumberOfDimensions(prhs[0]);
    int DNDims   = mxGetNumberOfDimensions(prhs[1]);
    if(imgNDims!=3) { mexErrMsgTxt("Input image must have 3 dimensions!"); }
    if(DNDims!=2) { mexErrMsgTxt("Distance matrix must have 2 dimensions!"); }
    if( !mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsClass(prhs[0],"int32") )
                    mexErrMsgTxt("Input image must be numeric matrix of type int32!");
    if( !mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || !mxIsClass(prhs[1],"double") )
                    mexErrMsgTxt("Distance matrix must be numeric matrix of type double!");
    if( !mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 )
                    mexErrMsgTxt("Argument 3 expected a scalar value.");
    
    // Parse the method inputs
    const size_t* imgDims   = mxGetDimensions(prhs[0]);
    const size_t* DDims     = mxGetDimensions(prhs[1]);
    int nrOfClusters        = static_cast<int>(mxGetScalar(prhs[2]));
    if(DDims[0]!=DDims[1]) mexErrMsgTxt("Distance matrix must be square!");
    int nrOfBins            = DDims[0];
    int nrOfVoxels          = imgDims[2]*imgDims[1]*imgDims[0];
    int* image              = reinterpret_cast<int*>(mxGetData(prhs[0]));
    double* D               = reinterpret_cast<double*>(mxGetData(prhs[1]));

    // Build histogram
    int* counts = new int[nrOfBins](); // init to zero
    for(int i = 0; i < nrOfVoxels; ++i)
    {
        counts[image[i]-1] += 1; // Convert Matlab index to cpp index
    }
    
    // Initialise labels
    int* labels = new int[nrOfBins];
    for(int i = 0; i < nrOfBins; ++i)
    {
        labels[i] = std::floor(i / (nrOfBins/nrOfClusters));
    }
    
    // Start optimization
    int maxIter = 50;
    double cost = calculateCost(counts, labels, nrOfBins, D);
    double oldCost, binCost, newBinCost;
    int bestLabel;
    for(int iter = 0; iter < maxIter; ++iter)
    {
        oldCost = cost;
        for(int i = 0; i < nrOfBins; ++i)
        {
            bestLabel   = labels[i];
            binCost     = calculateSingleCost(counts, labels, nrOfBins, D, i);
            for(int j = 0; j < nrOfClusters; ++j)
            {
                labels[i]   = j;
                newBinCost  = calculateSingleCost(counts, labels, nrOfBins, D, i);
                if(newBinCost < binCost)
                    bestLabel = j;
            }
            labels[i] = bestLabel;
        }
        
        cost = calculateCost(counts, labels, nrOfBins, D);
        std::cout << "\nCurrent cost: " << cost;
        if(oldCost == cost)
            break;
    }

    plhs[0] = mxCreateNumericArray(imgNDims, imgDims, mxINT32_CLASS, mxREAL);
    int *ptr = (int *)(mxGetData(plhs[0]));
        for(int i = 0; i < nrOfVoxels; ++ i)
            ptr[i] = labels[image[i]];
    
    // Clean up
    delete [] counts;
    delete [] labels;
}

