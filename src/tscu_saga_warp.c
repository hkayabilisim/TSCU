#include "mex.h"
#include "tscu_saga_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *y,*s,*w,*z;
    int    n,k;
    
    if (nrhs != 2) {
        mexErrMsgTxt("Usage: [z w]=tscu_saga_warp(y s)");
    } 
    
    y     = mxGetPr(prhs[0]);   
    s     = mxGetPr(prhs[1]);   

    
    n = mxGetNumberOfElements(prhs[0]);
    k = mxGetNumberOfElements(prhs[1]);
        
    plhs[0] = mxCreateDoubleMatrix(1, n,    mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, n,    mxREAL);
    
    z     = mxGetPr(plhs[0]); 
    w     = mxGetPr(plhs[1]); 

    tscu_saga_base(w,s,n,k);
    tscu_saga_solve(w,n);
    tscu_saga_warp(w,y,z,n);
}

        
