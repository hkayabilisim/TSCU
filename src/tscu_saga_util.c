#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "tscu_saga_util.h"
#include <Accelerate/Accelerate.h>

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
    printf("[");
	for (i=0;i<m;i++)
		printf("%8.5f ",x[i]);
    printf("]';\n");
}

void tscu_saga_disp_bmat(double *x,int n, int k)
{
    int i,j;
    for (j=0;j<k;j++)
    {
        printf("xbmat(:,%d)= [",j+1);
        for (i=0; i<n; i++) {
            printf("%5.3f ",x[j+i*k]);
        }
        printf("];\n");
    }
}
void tscu_saga_dispi(int *x, int m)
{
	int i;
	for (i=0;i<m;i++)
		printf("%4d ",x[i]);
	printf("\n");
	return;
}

void tscu_saga_base(double *bmat,int n,int k)
{
    int i,j;
    int delta,deltalast;
    
    delta = n/k;
    deltalast = n-(k-1)*delta;
    
	for (j=0;j<k-1;j++)
    {
		for(i=j*delta;i<j*delta+delta;i++)
        {
			bmat[j+i*k] = 1;
        }
    }
    for (i=n-deltalast; i< n ; i++) {
        bmat[k-1+i*k] = 1;
    }

    /* Double integration */
    tscu_saga_integration(bmat, n, k);
    tscu_saga_integration(bmat, n, k);
}

/* x = n x k [row-wise allocation]
 */
void tscu_saga_integration(double *x,int n, int k)
{
    int i,j;
    double dh = 1.0/(n-1.0);
    for (j=0; j<k; j++) {
        for (i=1; i<n; i++) {
            x[j+i*k] = x[j+i*k-k] + x[j+i*k];
        }
    }

    for (i=0; i<n*k; i++) {
        x[i] = x[i]*dh;
    }

}

void tscu_saga_solve(double *u,double *s,double *bmat,int n, int k)
{
    int i,j;
    double a,b,dh;
    dh = 1.0/(n-1.0);
    
    /* cblas_dgemv(CblasRowMajor, CblasNoTrans, n, k, 1, bmat, n, s, 1, 0, u, 1); */
    
    /* Solution with Taylor approximation */
    for (i=0; i<n; i++) {
        u[i] = 0.0;
        for (j=0; j<k; j++) {
            u[i] = u[i] + bmat[j+i*k]*s[j];
        }
    }
        
    for (i=0; i<n; i++) {
        u[i] = u[i] + i*dh;
    }
    
    /* Without Taylor */
    
    /* Normalize to [1 n] */
    a = (n-1.0)/(u[n-1]-u[0]);
    b = 1-u[0]*a;
    for (i=0;i<n;i++)
        u[i] = b+a*u[i];
    u[0]  =1;
    u[n-1]=n;
    
    /* bmat (nxk)   s (1xk)*/
    
    
    /*
     tscu_saga_cumsum(w,n);
     
     h=1.0/(n-1);
     b=1.0/expf(h*w[0]);
     for (i=0;i<n;i++)
     {
     w[i] = b*expf(h*w[i]);
     }
     tscu_saga_cumsum(w,n);
     a = (n-1.0)/(w[n-1]-w[0]);
     b = 1-w[0]*a;
     for (i=0;i<n;i++)
     w[i] = b+a*w[i];
     w[0]  =1;
     w[n-1]=n; */
    
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
void tscu_saga_register(double *x,double *y,double *z,double *path1, double *path2, double *d,int n,int k,double *u,double *s,double *sbest,double *bmat)
{
    int i;
    
    int PopulationSize =  40;
    int eliteChild     =  4;
    int mutationChild  =  18;
    int crossoverChild =  18;
    int nParents       = 2*crossoverChild + mutationChild;
    double *Fitness,*Expectation;
    double **Population1;
    double **Population2;
    double **swap;
    int *parents;
    int *rank;
    int generation,maxgeneration=10;
    
    
    /* srand((unsigned int)time(NULL)); */
    srand((unsigned int)1);
    
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
    for (i=0; i<nParents; i++) {
        parents[i] = rand() % PopulationSize;
    }
    for (i=0; i<PopulationSize; i++) {
        rank[i] = i;
    }
    
    tscu_saga_ga_init_population(Population1,PopulationSize,k);
    tscu_saga_ga_fitness(Population1,PopulationSize,Fitness,x,y,n,k,u,z,bmat);
    tscu_saga_sort(Fitness,PopulationSize, rank);
    /* tscu_saga_disp_pop(Population1, PopulationSize, k, -1, Fitness, rank); */
    
    for (generation=0; generation<maxgeneration; generation++) {
        tscu_saga_ga_expectation(Fitness,PopulationSize,Expectation);
        tscu_saga_ga_selection(Expectation,PopulationSize,parents,nParents,rank);
        tscu_saga_shuffle(parents,nParents);
        tscu_saga_elite(Population1,Population2,PopulationSize,eliteChild,rank,k);
        tscu_saga_crossover(Population1,Population2,PopulationSize,eliteChild,crossoverChild,parents,k);
        tscu_saga_mutation(Population1,Population2,PopulationSize,eliteChild,crossoverChild,mutationChild,parents,k);
        tscu_saga_ga_fitness(Population2,PopulationSize,Fitness,x,y,n,k,u,z,bmat);
        tscu_saga_sort(Fitness, PopulationSize, rank);
        
        swap = Population1;
        Population1 = Population2;
        Population2 = swap;
        /* tscu_saga_disp_pop(Population2, PopulationSize, k, generation, Fitness, rank); */
        
    }
    
    tscu_saga_solve(u,Population1[rank[0]],bmat,n,k);
    tscu_saga_warp(u,y,z,n);
    
    for (i=0;i<n;i++)  path1[i] = i+1;
    for (i=0;i<n;i++)  path2[i] = u[i];
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
void   tscu_saga_disp_pop(double **Population,int PopulationSize,int k,int generation,double *Fitness,int *rank)
{
    int i,j;
    printf("Generation: %d\n",generation);
    printf("best score: %8.5f :",Fitness[0]);
    for (i=0; i<k; i++) printf("%8.5f ",Population[rank[0]][i]);  printf("\n");
    for (j=1;j<PopulationSize; j++) {
        printf("            %8.5f :",Fitness[j]);
        for (i=0; i<k; i++) printf("%8.5f ",Population[rank[j]][i]);  printf("\n");
    }
    
}
void   tscu_saga_crossover(double **Population,double **NewPopulation,int PopulationSize, int eliteChild,int crossoverChild,int *parents,int k){
    int i,j,p;
    for (i=0; i<2*crossoverChild; i=i+2) {
        for (j=0; j<k; j++) {
            p = i+(rand() % 2);
            NewPopulation[i/2+eliteChild][j] = Population[parents[p]][j];
        }
    }
}
void   tscu_saga_mutation(double **Population,double **NewPopulation,int PopulationSize, int eliteChild,int crossoverChild,int mutationChild,int *parents,int k){
    int i,j;
    for (i=eliteChild+crossoverChild; i<PopulationSize; i++) {
        for (j=0; j<k; j++) {
            NewPopulation[i][j] = Population[parents[i]][j]+0.01*(2.0*rand()/RAND_MAX-1.0);
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

/* Buble-Sort */
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

void  tscu_saga_ga_fitness(double **Population,int PopulationSize,double *Fitness,double *x,double *y,int n,int k,double *u,double *z,double *bmat)
{
    int i;
    double c;
    for (i=0; i<PopulationSize; i++) {
        tscu_saga_cost(x,y,Population[i],n,k,&c,0,u,z,bmat);
        Fitness[i] = c;
    }
    
}

void tscu_saga_ga_init_population(double **Population,int PopulationSize,int Length)
{
    int i,j;
    
    for (j=0;j<Length;j++)
        Population[0][j] = 0.0;
    for (i=1; i<PopulationSize; i++) {
        for (j=0; j<Length; j++) {
            Population[i][j] = ((double)rand()/(double)(RAND_MAX)) * 2.0 - 1.0;
        }
    }
}

void tscu_saga_cost(double *x,double *y, double *s, int n, int k,double *d,int debug,double *u,double *z,double *bmat)
{
    
    tscu_saga_solve(u,s,bmat,n,k);
    tscu_saga_warp(u,y,z,n);
    tscu_saga_rms(d,x,z,n);
    
    if (debug == 1)
    {
        printf("w"); tscu_saga_disp(u,n);
        printf("x"); tscu_saga_disp(x,n);
        printf("y"); tscu_saga_disp(y,n);
        printf("z"); tscu_saga_disp(z,n);
        printf("d"); tscu_saga_disp(d,1);
    }
    
}
void tscu_saga_qsort(double *a, int lo, int hi,int *index)
{
    int h, l;
    int i;
    double p,t;
    if (lo < hi) {
        l = lo;
        h = hi;
        p = a[hi];
        
        do {
            while ((l < h) && (a[l] <= p))
                l = l+1;
            while ((h > l) && (a[h] >= p))
                h = h-1;
            if (l < h) {
                t = a[l];
                a[l] = a[h];
                a[h] = t;
                
                i = index[l];
                index[l]=index[h];
                index[h] = i;
                
            }
        } while (l < h);
        
        a[hi] = a[l];
        a[l] = p;
        
        tscu_saga_qsort( a, lo, l-1,index );
        tscu_saga_qsort( a, l+1, hi,index );
    }
}
