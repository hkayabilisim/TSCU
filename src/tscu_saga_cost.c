#include "mex.h"
#include "tscu_saga_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *x,*y,*s,*w,*z;
    int    n,k;
    double *d;
    
    if (nrhs != 5) {
        mexErrMsgTxt("Usage: x y s w z");
    } 
    
    x = mxGetPr(prhs[0]);   
    y = mxGetPr(prhs[1]);   
    s = mxGetPr(prhs[2]);   
    w = mxGetPr(prhs[3]);   
    z = mxGetPr(prhs[4]);   
    
    n = mxGetNumberOfElements(prhs[0]);
    k = mxGetNumberOfElements(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1,    mxREAL);
    
    d = mxGetPr(plhs[0]); 
    
    tscu_saga_cost(x,y,s,n,k,d,0,w,z);     
}

        
