//
//  saga_util.h
//  SAGA
//
//  Created by Huseyin Kaya on 25/03/14.
//  Copyright (c) 2014 Huseyin Kaya. All rights reserved.
//

#ifndef SAGA_saga_util_h
#define SAGA_saga_util_h

void   tscu_saga_rms(double *d,double *x,double *z,int n);
void   tscu_saga_cumsum(double *x,int m);
void   tscu_saga_disp(double *x, int m);
void   tscu_saga_dispi(int *x, int m);
void   tscu_saga_base(double *w,int n,int k);
void   tscu_saga_solve(double *w,double *s,double *bmat,int n, int k);
void   tscu_saga_warp(double *h, double *y, double *z, int n);
void   tscu_saga_register(double *x,double *y,double *z,double *path1,double *path2,double *d,int n,int k,double *w,double *s,double *sbest,double *bmat);
void   tscu_saga_cost(double *x,double *y, double *s, int n, int k,double *d,int debug,double *w, double *z,double *bmat);
void   tscu_saga_ga_init_population(double **Population,int PopulationSize,int Length);
void   tscu_saga_ga_fitness(double **Population,int PopulationSize,double *Fitness,double *x,double *y,int n,int k,double *w,double *z,double *bmat);
void   tscu_saga_ga_selection(double *Fitness,int PopulationSize,int *parents,int nParents,int *rank);
void   tscu_saga_ga_expectation(double *Fitness,int PopulationSize,double *Expectation);
void   tscu_saga_sort(double *a,int n,int *index);
void   tscu_saga_shuffle(int *parents, int nparents);
void   tscu_saga_elite(double **Population,double **NewPopulation,int PopulationSize,int eliteChild,int *rank,int k);
void   tscu_saga_crossover(double **Population,double **NewPopulation,int PopulationSize, int eliteChild,int crossoverChild,int *parents,int k);
void   tscu_saga_mutation(double **Population,double **NewPopulation,int PopulationSize, int eliteChild,int crossoverChild,int mutationChild,int *parents,int k);

void   tscu_saga_disp_pop(double **Population,int PopulationSize,int k,int generation,double *Fitness,int *rank);
void   tscu_saga_min(double *y,int n,double *min,int *minindex);
void   tscu_saga_qsort(double *x, int lo,int hi,int *index);
void   tscu_saga_integration(double *x,int rowsize, int colsize);
void   tscu_saga_disp_bmat(double *x,int rowsize, int colsize);

#endif