
#include "mex.h"
#include "tscu_saga_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *x,*y,*k,*z,*w,*s,*sbest;
    int    n;
    double *d;
    double *path1,*path2;
    double *bmat;
    
    if (nrhs != 8) {
        mexErrMsgTxt("Usage: [path1 path2 d]=tscu_saga_register(x y k z w s sbest bmat)");
    } 
    
    x     = mxGetPr(prhs[0]);   
    y     = mxGetPr(prhs[1]);   
    k     = mxGetPr(prhs[2]);   
    z     = mxGetPr(prhs[3]);    
    w     = mxGetPr(prhs[4]); 
    s     = mxGetPr(prhs[5]); 
    sbest = mxGetPr(prhs[6]); 
    bmat  = mxGetPr(prhs[7]);

    
    n = mxGetNumberOfElements(prhs[0]);
        
    plhs[0] = mxCreateDoubleMatrix(1, n,    mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, n,    mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1,    mxREAL);
    
    path1 = mxGetPr(plhs[0]);    
    path2 = mxGetPr(plhs[1]); 
    d     = mxGetPr(plhs[2]); 

    tscu_saga_register(x,y,z,path1,path2,d,n,(int)*k,w,s,sbest,bmat);  
    
}

        
