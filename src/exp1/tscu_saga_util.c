//
//  saga_util.c
//  SAGA
//
//  Created by Huseyin Kaya on 25/03/14.
//  Copyright (c) 2014 Huseyin Kaya. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "tscu_saga_util.h"

void tscu_saga_rms(double *rms,double *x,double *z,int n)
{
    rms[0] = 0.0;
    int i;
    for (i=0;i<n;i++)
        rms[0] = rms[0] + (x[i]-z[i])*(x[i]-z[i]);
    rms[0] = rms[0]/(n-1);
}

void tscu_saga_cumsum(double *x,int m)
{
	int i;
	for (i=1;i<m;i++)
        x[i] = x[i-1]+x[i];
}

void tscu_saga_disp(double *x, int m)
{
	int i;
	printf("x = [");
	for (i=0;i<m;i++)
		printf("%5.2f ",x[i]);
	printf("];\n");
	return;
}
void tscu_saga_dispi(int *x, int m)
{
	int i;
	for (i=0;i<m;i++)
		printf("%4d ",x[i]);
	printf("\n");
	return;
}

void tscu_saga_base(double *w,double *s,int n,int k)
{
    int i,j;
    int delta;
    int *breakpoints;
    
    /* Breakpoints */
    delta = (n-1)/k;
    breakpoints = (int *)malloc(sizeof(int)*(unsigned long)(k+1));
	breakpoints[0] = 0;
	for (i=1;i<k;i++)
		breakpoints[i]=breakpoints[i-1]+delta;
	breakpoints[k]=n-1;
    
    /* linear combination */
	for (i=0;i<k;i++)
		for(j=breakpoints[i];j<= breakpoints[i+1];j++)
			w[j]=s[i];
    /*
     printf("breakpoints"); saga_dispi(breakpoints,k+1);
     printf("s          "); saga_disp(s,k);
     printf("w          "); saga_disp(w,n); */
    free(breakpoints);
}

void tscu_saga_solve(double *w,int n)
{
    int i;
    
    tscu_saga_cumsum(w,n);
    
	for (i=0;i<n;i++)
    {
        w[i] = expf((w[i] - w[0])/(n-1.0));
        
    }
    
    
	tscu_saga_cumsum(w,n);
    
    for (i=1;i<n;i++)
		w[i] = 1 + (w[i]-w[0])/(w[n-1]-w[0])*(n-1.0);
    w[0]  =1;
    w[n-1]=n;
    
}

void tscu_saga_warp(double *h, double *y, double *z, int n)
{
    int i;
    int l;
    /* y(h(t)) = z(t) */
    z[0]  =y[0];
    z[n-1]=y[n-1];
    
    for (i=1;i<n-1;i++)
    {
        l = (int)floor(h[i]);
        z[i]=(y[l]-y[l-1])*(h[i]-l)+y[l-1];
    }
}
void tscu_saga_register(double *x,double *y,double *z,double *path1, double *path2, double *d,int n,int k,double *w,double *s,double *sbest)
{
    int i;
    
    int PopulationSize = 20;
    int eliteChild     =  2;
    int mutationChild  =  3;
    int crossoverChild = 15;
    int nParents       = 2*crossoverChild + mutationChild;
    double *Fitness,*Expectation;
    double **Population1;
    double **Population2;
    double **swap;
    int *parents;
    int *rank;
    int generation,maxgeneration=100;
    
    
    srand((unsigned int)time(NULL));
    
    /* Allocations */
    Population1   = (double **)malloc(sizeof(double *)*PopulationSize);
    Population2 = (double **)malloc(sizeof(double *)*PopulationSize);
    for (i=0; i<PopulationSize; i++)
    {
        Population1[i] = (double *)malloc(sizeof(double)*k);
        Population2[i] = (double *)malloc(sizeof(double)*k);
    }
    Fitness = (double *)malloc(sizeof(double)*PopulationSize);
    Expectation = (double *)malloc(sizeof(double)*PopulationSize);
    rank = (int *)malloc(sizeof(int)*PopulationSize);
    parents = (int *)malloc(sizeof(int)*nParents);
    for (i=0; i<PopulationSize; i++)    rank[i] = i;
    
    tscu_saga_ga_init_population(Population1,PopulationSize,k);
    tscu_saga_ga_fitness(Population1,PopulationSize,Fitness,x,y,n,k,w,z);
    /* printf("**** Mutation\n"); tscu_saga_disp_pop(Population1, PopulationSize, k);
    for (i=0; i<PopulationSize; i++) printf("%5.2f ",Fitness[i]); printf("\n");
    for (i=0; i<PopulationSize; i++) printf("%5d "  ,rank[i]); printf("\n");*/

    tscu_saga_sort(Fitness, PopulationSize, rank);
    /* for (i=0; i<PopulationSize; i++) printf("%5d "  ,rank[i]); printf("\n");*/

    for (generation=0; generation<maxgeneration; generation++) {

        /* printf("**** Initial Population\n"); tscu_saga_disp_pop(Population, PopulationSize, k);
        for (i=0; i<PopulationSize; i++) printf("%5.2f ",Fitness[i]); printf("\n"); */
        
        tscu_saga_ga_expectation(Fitness,PopulationSize,Expectation);
        tscu_saga_ga_selection(Expectation,PopulationSize,parents,nParents,rank);
        tscu_saga_shuffle(parents,nParents);
        /* printf("**** Initial Population\n"); tscu_saga_disp_pop(Population, PopulationSize, k);
        for (i=0; i<PopulationSize; i++) printf("%5.2f ",Fitness[i]); printf("\n");
        for (i=0; i<PopulationSize; i++) printf("%5d "  ,rank[i]); printf("\n"); */
        
        tscu_saga_elite(Population1,Population2,PopulationSize,eliteChild,rank,k);
        /* printf("**** Elit Population\n"); tscu_saga_disp_pop(NewPopulation, PopulationSize, k); */
        
        tscu_saga_crossover(Population1,Population2,PopulationSize,eliteChild,crossoverChild,parents,k);
        /* printf("**** CrossOver\n"); tscu_saga_disp_pop(NewPopulation, PopulationSize, k); */
        
        tscu_saga_mutation(Population1,Population2,PopulationSize,eliteChild,crossoverChild,mutationChild,parents,k);
        tscu_saga_ga_fitness(Population2,PopulationSize,Fitness,x,y,n,k,w,z);
        /* printf("**** Mutation\n"); tscu_saga_disp_pop(Population2, PopulationSize, k);
        for (i=0;i<PopulationSize;i++) printf("%5.2f ",Fitness[i]); printf("\n");
        for (i=0; i<PopulationSize; i++) printf("%5d "  ,rank[i]); printf("\n"); */

        
        tscu_saga_sort(Fitness, PopulationSize, rank);
        /* for (i=0; i<PopulationSize; i++) printf("%5d "  ,rank[i]); printf("\n"); */
/*
        printf("best score: %8.5f [",Fitness[0]);
        for (i=0; i<k; i++) {
            printf("%5.4f ",Population2[rank[0]][i]);
        }
        printf("]\n"); */

        swap = Population1;
        Population1 = Population2;
        Population2 = swap;
    }
    
    /*
    jmax = 10;
    cmin = 10000;
    for (j=0;j<jmax;j++)
	{
		for (i=0;i<k;i++)
		{
			s[i]= ((double)rand()/(double)(RAND_MAX)) * 4 - 2;
		}
		tscu_saga_cost(x,y,s,n,k,&c,0,w,z);
        if (c < cmin)
        {
            cmin = c;
            for (i=0;i<k;i++) sbest[i] = s[i];
        }
        
        
	} */
    tscu_saga_base(w,Population1[rank[0]],n,k);
    tscu_saga_solve(w,n);
    tscu_saga_warp(w,y,z,n);
    
    for (i=0;i<n;i++)  path1[i] = i+1;
    for (i=0;i<n;i++)  path2[i] = w[i];
    d[0] = Fitness[0];
    
    /* Deallocations */
    for (i=0; i<PopulationSize; i++) {
        free(Population1[i]);
        free(Population2[i]);
    }
    free(Population1);
    free(Population2);
    free(Fitness);
    free(Expectation);
    free(parents);
    free(rank);
}
void   tscu_saga_min(double *y,int n,double *min,int *minindex)
{
    *min=1000000;
    *minindex = 0;
    int i;
    for (i=0; i<n; i++) {
        if (y[i] < *min)
        {
            *min = y[i];
            *minindex = i;
        }
    }
}
void   tscu_saga_disp_pop(double **Population,int PopulationSize,int k)
{
    int i,j;
    for (j=0; j<k;j++) {
        for (i=0; i<PopulationSize; i++) {
            printf("%5.2f ",Population[i][j]);
        }
        printf("\n");
    }
    
}
void   tscu_saga_crossover(double **Population,double **NewPopulation,int PopulationSize, int eliteChild,int crossoverChild,int *parents,int k){
    int i,j;
    for (i=0; i<2*crossoverChild; i=i+2) {
        for (j=0; j<k; j++) {
            NewPopulation[i/2+eliteChild][j] = Population[parents[i+(rand() % 2)]][j];
        }
    }
}
void   tscu_saga_mutation(double **Population,double **NewPopulation,int PopulationSize, int eliteChild,int crossoverChild,int mutationChild,int *parents,int k){
    int i,j;
    for (i=eliteChild+crossoverChild; i<PopulationSize; i++) {
        for (j=0; j<k; j++) {
            /* printf("%d::",parents[i]); */
            NewPopulation[i][j] = Population[parents[i]][j]+rand()/RAND_MAX;
        }
    }
}
void   tscu_saga_elite(double **Population,double **NewPopulation,int PopulationSize,int eliteChild,int *rank,int k)
{
    int i,j;
    for (i=0; i<eliteChild; i++) {
        for (j=0; j<k; j++) {
            NewPopulation[i][j] = Population[rank[i]][j];
        }
    }
}

void   tscu_saga_shuffle(int *array, int n)
{
    int i;
    int t;
    for (i = 0; i < n - 1; i++)
    {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

void tscu_saga_sort(double *a,int n,int *index)
{
    int i,j;
    double temp;
    int l;
    for (i=0; i<n; i++) {
        index[i] = i;
    }
    /*
    int k;
    printf("Unsorted Data:");
    for(k=0;k<n;k++)
        printf("%5.2f",a[k]);
    printf("\nUnsorted Data:");
    for(k=0;k<n;k++)
        printf("%5.2d",index[k]); */
    for(i=1;i< n;i++)
    {
        for(j=0;j< n-1;j++)
            if(a[j]>a[j+1])
            {
                temp=a[j];
                a[j]=a[j+1];
                a[j+1]=temp;
                l = index[j];
                index[j] = index[j+1];
                index[j+1] = l;
            }
    }
    /* printf("\n");
    for(k=0;k< n;k++)
    {
        printf("%5d %5d %5.2f %5.2f\n",k,index[k],a[index[k]],a[k]);
    } */
}

void tscu_saga_ga_expectation(double *Fitness,int PopulationSize,double *Expectation)
{
    int i;
    double sum=0.0;
    for (i=0; i<PopulationSize; i++) {
        Expectation[i] = 1.0/sqrt(i+1.0);
        sum += Expectation[i];
    }
    for (i=0; i<PopulationSize;i++)
    {
        Expectation[i] = PopulationSize*Expectation[i]/sum;
    }
}

void tscu_saga_ga_selection(double *Expectation,int PopulationSize,int *parents,int nParents,int *rank)
{
    int lowest = 0;
    int i,j;
    
    double stepsize = PopulationSize*1.0/nParents;
    double position = stepsize;
    
    tscu_saga_cumsum(Expectation, PopulationSize);
    
    for (i=0; i<nParents; i++) {
        for (j=lowest; j<PopulationSize; j++) {
            /* printf("%5.2f position %5.2f \n",position,Expectation[j]); */
            if (position<Expectation[j]) {
                parents[i] = rank[j];
                lowest = j;
                break;
            }
        }
        position += stepsize;
    }
}

void  tscu_saga_ga_fitness(double **Population,int PopulationSize,double *Fitness,double *x,double *y,int n,int k,double *w,double *z)
{
    int i;
    double c;
    for (i=0; i<PopulationSize; i++) {
        tscu_saga_cost(x,y,Population[i],n,k,&c,0,w,z);
        Fitness[i] = c;
    }
    
}

void tscu_saga_ga_init_population(double **Population,int PopulationSize,int Length)
{
    int i,j;
    
    for (i=0; i<PopulationSize; i++) {
        for (j=0; j<Length; j++) {
            Population[i][j] = ((double)rand()/(double)(RAND_MAX)) * 4 - 2;
        }
    }
}

void tscu_saga_cost(double *x,double *y, double *s, int n, int k,double *d,int debug,double *w,double *z)
{
    
    
    tscu_saga_base(w,s,n,k);
    tscu_saga_solve(w,n);
    tscu_saga_warp(w,y,z,n);
    tscu_saga_rms(d,x,z,n);
    
    if (debug == 1)
    {
        printf("w"); tscu_saga_disp(w,n);
        printf("x"); tscu_saga_disp(x,n);
        printf("y"); tscu_saga_disp(y,n);
        printf("z"); tscu_saga_disp(z,n);
        printf("d"); tscu_saga_disp(d,1);
    }
    
}