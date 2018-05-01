
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>

// ---------------------------------------------------------- //
// -- auxiliary function: find image min and max intensity -- //
// ---------------------------------------------------------- //

void imgMinMax(double* image, int nrOfVoxels, double&min, double&max)
{
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();
    for(int i = 0; i < nrOfVoxels; ++i)
    {
        if(image[i] < min) min = image[i];
        if(image[i] > max) max = image[i];
    }
}

// -------------------------------------------- //
// -- main method: estimates the images ISDs -- //
// -------------------------------------------- //

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check the inputs and outputs
    if(nrhs != 2)   mexErrMsgTxt("Expected 2 inputs (image,nrOfBins)!");
    if(nlhs < 1)    mexErrMsgTxt("Expected at least one output argument!");
    if(nlhs > 2)    mexErrMsgTxt("Expected no more than 2 output arguments!");
    int ndims = mxGetNumberOfDimensions(prhs[0]);
    if(ndims!=3) { mexErrMsgTxt("Input image must have 3 dimensions!"); }
    if( !mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsClass(prhs[0],"double") )
                    mexErrMsgTxt("Input image must be numeric matrix of type double!");
    if( !mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1 )
                    mexErrMsgTxt("Invalid argument (2)");
    
    // Parse the method inputs
    const size_t* dims  = mxGetDimensions(prhs[0]);
    int nrOfVoxels      = dims[2]*dims[1]*dims[0];
    double* image       = reinterpret_cast<double*>(mxGetData(prhs[0]));
    int nrOfBins        = static_cast<int>(mxGetScalar(prhs[1]));
    
    // Find image min/max and discretize image
    double min,max,binSize;
    imgMinMax(image,nrOfVoxels,min,max);
    binSize = (max-min)/double(nrOfBins);
    int* imageD = new int[nrOfVoxels];
    for(int i = 0; i < nrOfVoxels; ++i)
    {
        for(int j = 0; j < nrOfBins; ++j)
        {
            if(image[i] >= min+j*binSize && image[i] <= min+(j+1)*binSize)
            {
                imageD[i] = j;
                break;
            }
        }
    }   
    
    // Find most similar neighbour for each voxel
    std::vector<std::vector<double> > values(nrOfBins, std::vector<double>(0));
    int offset[6] = {-1,1,                                                      // top,bottom
                     -int(dims[0]),int(dims[0]),                                // left,right
                     -int(dims[0])*int(dims[1]),int(dims[0])*int(dims[1])};     // front back
    int v, diff,idx;
    double vbest;
    for (int i = 0; i < nrOfVoxels; ++i) {
        v = imageD[i];
        diff = std::numeric_limits<int>::max();
        for (int j = 0; j < 6; ++j) {
            idx = i+offset[j];
            if(idx < 0 || idx >= nrOfVoxels) continue; // outside image
            if(std::abs(v-imageD[idx]) < diff)
            {
                diff = std::abs(v-imageD[idx]);
                vbest = image[idx];
            }
        }
        values.at(v).push_back(vbest);
    }
    
    // Create output arguments
    plhs[0] = mxCreateCellMatrix(nrOfBins,1);
    for(int i = 0; i < values.size(); ++i)
    {
        mxArray* m = mxCreateDoubleMatrix(values.at(i).size(),1,mxREAL);
        double *ptr = (double *)(mxGetData(m));
        for(int j = 0; j < values.at(i).size(); ++j)
            ptr[j] = values.at(i).at(j);
        mxSetCell(plhs[0], i, m);
    }
    if(nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(nrOfBins+1,1,mxREAL);
        double *ptr = (double *)(mxGetData(plhs[1]));
        for(int i = 0; i <= nrOfBins; ++ i)
            ptr[i] = min+i*binSize;
    }

    delete [] imageD;
}

