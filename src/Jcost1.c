#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define min(x, y) ((x)<(y)?(x):(y))
#define max(x, y) ((x)>(y)?(x):(y))
#define dist(x, y) ((x-y)*(x-y))


#define INF 1e20       

#define eps 0.00001
void ramsay(double *w, double *b, int m, int nt, double *u,
        double **A,double *rhs) {
    int i, j;
    double *Awork;
    int *pivot;
    int n,one,lda,ldb,info;
    
    /*
    printf("----- Ramsay ----\n");
    printf("m : %d\n",m);
    printf("n : %d\n",n); */
    
    for(i=0;i<2*m;i++) {
        rhs[i]=0.0;
        for(j=0;j<2*m;j++)
            A[i][j] = 0.0;
    }
    
    
    for (i=0;i<m-1;i++)
    {
        if (fabs(w[i]) < eps ) {
            A[i][i]     = b[i+1];
            A[i+m-1][i] = 1;
        }
        else {
            A[i][i]     = exp(w[i]*b[i+1]);
            A[i+m-1][i] = w[i]*exp(w[i]*b[i+1]);
        }
        if (fabs(w[i+1]) < eps )
        {
            A[i][i+1]     = -b[i+1];
            A[i+m-1][i+1] = -1;
        }
        else
        {
            A[i][i+1]     = -exp(w[i+1]*b[i+1]);
            A[i+m-1][i+1] = -w[i+1]*exp(w[i+1]*b[i+1]);
        }
        
        A[i][m+i] = 1;
        A[i][m+i+1] = -1;
        
    }
    
    if (fabs(w[0]) > eps)  A[2*m-2][0] = 1;
    
    A[2*m-2][m]=1;
    if (fabs(w[m-1]) < eps)
        A[2*m-1][m-1] = 1;
    else
        A[2*m-1][m-1] = exp(w[m-1]);
    
    A[2*m-1][2*m-1]=1;
    
    
    rhs[2*m-1] = 1;
    
    /*
    for(i=0;i<2*m;i++) {
        for(j=0;j<2*m;j++)
        {
            printf("%8.2f ",A[i][j]);
        }
        printf("\n");
    } */

    Awork = (double *)malloc(sizeof(double)*2*m*2*m);
    pivot = (int *)malloc(sizeof(int)*2*m);

    for (j=0;j<2*m;j++)
   	 for (i=0;i<2*m;i++)
	 	Awork[j*2*m+i]=A[i][j];

    /* solution with lapack */
    n = 2*m;
    one = 1;
    lda = 2*m;
    ldb = 2*m;
   dgesv_(&n,&one,Awork,&lda,pivot,rhs,&ldb,&info);
/*
   printf("Solution: ");
   for (i=0;i<2*m;i++)
   	printf("%8.2f ",rhs[i]);
	printf("\n"); */
 
  free(Awork);
  free(pivot);
}

double ramsay_main(double *tstObject,double *trnObject, double *s, int n, int m)
{
    int i,j;
    double *b;
    double distance = 0;
    
    double **A; /* Linear system */
    double *rhs; /* right hand side of linear system */
    double *solution; /* Solution of linear system */    
    double *u;
   /*  double *warpedSignal;
    double *t; */
    
    /* breakpoints */
    b = (double *)malloc(sizeof(double)*(m+1));
    for(i=0;i<m+1;i++) {
        b[i] = i*1.0/m;
    }
       /* breakpoints */
    u = (double *)malloc(sizeof(double)*n);
    for(i=0;i<n;i++) {
        u[i] = 0;
    }
    /* 
    warpedSignal = (double *)malloc(sizeof(double)*n);
    for(i=0;i<n;i++)     warpedSignal[i] = 0;
    t = (double *)malloc(sizeof(double)*n);
    for(i=0;i<n;i++)     t[i] = i/(n-1.0); */
     
    /* Linear system first: A*/
    A = (double **)malloc(sizeof(double *)*2*m);
    for(i=0;i<2*m;i++)
        A[i] = (double *)malloc(sizeof(double)*2*m);
    solution = (double *)malloc(sizeof(double *)*2*m);
    
    ramsay(s, b, m, n, u, A,solution);
    
    for(i=0;i<n;i++){
        for (j=0;j<m-1;j++){
                if ((i*1.0/(n-1)) >= b[j] && (i*1.0/(n-1)) < b[j+1]) {
                    /* printf("%d %d\n",i,j); */
                    if (fabs(s[j]) < eps)
                        u[i] = solution[j]*(i*1.0/(n-1))+solution[j+m];
                    else
                        u[i] = solution[j]*exp(s[j]*(i*1.0/(n-1)))+solution[j+m];
                    
                }
        }
        /* If it is in the last bin */
        if ((i*1.0/(n-1)) >= b[m-1] && (i*1.0/(n-1)) <= b[m]) {
                    /* printf("%d %d\n",i,m); */
                    if (fabs(s[j]) < eps)
                        u[i] = solution[m-1]*(i*1.0/(n-1))+solution[2*m-1];
                    else
                        u[i] = solution[m-1]*exp(s[m-1]*(i*1.0/(n-1)))+solution[2*m-1];
                    break;
                }
    }
    u[0]=0;
    u[n-1]=1;
    
    /* Applying warping function to the testing */
     distance = 0;
    double temp=0.0;
    for (i=0;i<n;i++) {
        for (j=0;j<n-2;j++) {
            if (u[i] >= j/(n-1.0) && u[i] < (j+1)/(n-1.0)) {
                temp = (tstObject[j+1]-tstObject[j])/(((j+1)/(n-1.0))-(j/(n-1.0)))*(u[i]-(j/(n-1.0)))+tstObject[j];
                distance = distance + (temp-trnObject[i])*(temp-trnObject[i]);
              /*  printf("%d %d %8.2f\n",i,j,temp); */
            }
                
        }
            if (u[i] >= j/(n-1.0) && u[i] < (j+1)/(n-1.0)) {
           
                           temp = (tstObject[j+1]-tstObject[j])/(((j+1)/(n-1.0))-(j/(n-1.0)))*(u[i]-(j/(n-1.0)))+tstObject[j];
                distance = distance + (temp-trnObject[i])*(temp-trnObject[i]); 
                  /*       printf("%d %d %8.2f\n",i,j,temp); */
       
        }
            
    }
    
    /*
     printf("----\n");
    for (i=0;i<n;i++) {
        for (j=0;j<n-2;j++) {
            if (u[i] >= t[j] && u[i] < t[j+1]) {

                warpedSignal[i]=(tstObject[j+1]-tstObject[j])/(t[j+1]-t[j])*(u[i]-t[j])+tstObject[j]; 
                               printf("%d %d %8.2f\n",i,j,warpedSignal[i]);

            }
                
        }
        if (u[i] >= t[j] && u[i] <= t[j+1]) {
               warpedSignal[i]=(tstObject[j+1]-tstObject[j])/(t[j+1]-t[j])*(u[i]-t[j])+tstObject[j]; 
         printf("%d %d %8.2f\n",i,j,warpedSignal[i]);

        }
            
    }
    double distance2 = 0.0;
    for (i=0;i<n;i++)
        distance2 = distance2 + (trnObject[i]-warpedSignal[i])*(trnObject[i]-warpedSignal[i]); */
     
     /* 0    1     2  --> i
      * 0   0.5    1  --> t
      * 0   0.5    1
      */
         
    
    /*
   printf("warping: ");
   for (i=0;i<n;i++)
   	printf("%8.2f ",u[i]);
	printf("\n");
  
   printf("Distance = %8.3f\n",distance); */
    
    /* printf("Distance2 = %8.3f\n",distance2); */

    free(b);
    
    for(i=0;i<m;i++)
        free(A[i]);
    free(A);
    free(solution);
    /* free(t);
    free(warpedSignal); */

	return distance;
}


/* J = @(s) norm(ScaleTime(tstObject',ramsay3(par.t,s)*(length(tstObject)-1)+1)-trnObject'); */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *tstObject,*trnObject,*s; /* inputs */
    int n,m;
    double *d; /* output */
    
    /* check number of inputs and outputs */
    if (nrhs != 3) {
        mexErrMsgTxt("Usage: tstObject trnObject s");
    } 
    
    /* retrieve input arguments */
    tstObject = mxGetPr(prhs[0]);    /* pointer to real values of first  argument  */
    trnObject = mxGetPr(prhs[1]);    /* pointer to real values of second argument */
    s         = mxGetPr(prhs[2]);    /* pointer to real value  of third  argument   */
    
    /* check series lengths */
    n = mxGetNumberOfElements(prhs[0]);
    m = mxGetNumberOfElements(prhs[2]);
    
    /* allocate memory for the return value */
    plhs[0] = mxCreateDoubleMatrix(1, 1,    mxREAL);
    
    /* printf("Query Length:%d Path Length:%d\n", ql,ql*(2*((int)r[0])+1)); */
    d = mxGetPr(plhs[0]);    /* pointer to Matlab managed memory for result */
    
    d[0] = ramsay_main(tstObject,trnObject,s,n,m); 
    
}

        
