//
//  main.c
//  SAGA
//
//  Created by Hüseyin Kaya on 19/04/14.
//  Copyright (c) 2014 Hüseyin Kaya. All rights reserved.
//

//
//  main.c
//  SAGA
//
//  Created by Huseyin Kaya on 25/03/14.
//  Copyright (c) 2014 Huseyin Kaya. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tscu_saga_util.h"

#define PI 3.14159265358979323846264338327
#define EPS 0.000001

int main(int argc, const char * argv[])
{
    
	int j,n,k;
	double *x,*y,*z,*w;
    double *path1,*path2;
    double *d;
    double *s,*sbest;
    double *bmat;
    double t;
	n = 29;
	k = 8;
    
	x = (double *)malloc(sizeof(double)*(unsigned long)n);
	y = (double *)malloc(sizeof(double)*(unsigned long)n);
	z = (double *)malloc(sizeof(double)*(unsigned long)n);
    w = (double *)malloc(sizeof(double)*(unsigned long)n);
    path1 = (double *)malloc(sizeof(double)*(unsigned long)n);
    path2 = (double *)malloc(sizeof(double)*(unsigned long)n);
    sbest = (double *)malloc(sizeof(double)*(unsigned long)k);
    s  = (double *)malloc(sizeof(double)*(unsigned long)k);
    bmat = (double *)malloc(sizeof(double)*(unsigned long)(n*k));
    d  = (double *)malloc(sizeof(double)*(unsigned long)1);
    
    for (j=0; j<n; j++) {
        t = j/(n-1.0);
        x[j] = sin(2*PI*t);
        y[j] = sin(2*PI*t*t);
        z[j] = 0;
        w[n] = 0;
        path1[j] = 0;
        path2[j] = 0;
    }
    for (j=0; j<k; j++) {
        s[j]=0.0;
        sbest[j]=0.0;
    }
    for (j=0; j<n*k; j++) {
        bmat[j]=0.0;
    }
    d[0] = 0;
    
    tscu_saga_base(bmat, n,k);
    /* tscu_saga_disp_bmat(bmat, n, k); */
    /* 30 * 1000 */
    for (j=0; j<1; j++) {
        tscu_saga_register(x,y,z,path1,path2,d,n,k,w,s,sbest,bmat);
        
        printf("xx = "); tscu_saga_disp(x, n);
        printf("xy = "); tscu_saga_disp(y, n);
        printf("xz = "); tscu_saga_disp(z, n);
        printf("xpath1 = "); tscu_saga_disp(path1, n);
        printf("xpath2 = "); tscu_saga_disp(path2, n);
        printf("xd = "); tscu_saga_disp(d, 1);

    }
    
    free(bmat);
	free(x);
	free(y);
    free(z);
    free(w);
    free(s);
    free(sbest);
    free(path1);
    free(path2);
    free(d);
    return 0;
}



