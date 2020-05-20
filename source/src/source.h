#ifndef __SOURCE_H__
#define __SOURCE_H__


/* History: Feb 20 2019 Initial coding
            Mar 24 2020 Try LAMBDA_MINITER = LAMBDA_MONOITER = 2,
            Apr 04 2020 In non-parametric EM, change to let stopping criteria
                        be on beta instead of log-likelihood   
            Apr 29 2020 Add LRT_NP function        
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define DEBUG 0
#define EQUAL_EPS 1e-6
#define ALMOST_ZERO 1e-16
#define DBL_MISS -9.0e150
#define DBL_NOTMISS_GT -1.0e100
#define LAMBDA_MINITER 3
#define LAMBDA_MONOITER 3
#define NUMERIC_ZERO 1.0e-100

/* Set to 0 to use loglike, 1 for beta (lambda) */
#define EM_STOP_CRITERIA 0

#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}


/* printing */
void compare_dVec(double *, double *, int);
void compare_iVec(int *, int *, int);
void compare_iMat(int **, int **, int, int);
void print_dVec(double *, int, char *);
void print_iVec(int *, int, char *);
void print_dMat(double **, int, int, char *);
void print_iMat(int **, int, int, char *);

/* memory */
double * dVec_alloc(int, int, double);
int * iVec_alloc(int, int, int);
double ** dMat_alloc(int, int, int, double);
int ** iMat_alloc(int, int, int, int);
void matrix_free(void **, int);

/* matrix alg */
void copy_dVec(double *, double *, int);
void copy_iVec(int *, int *, int);
double sum_dVec(double *, int);
int sum_iVec(int *, int);
void permute_iVec(int *, int, int *);
double dotProd(double *, double *, int);
void matTimesVec(double **, double *, int, int, double *);
double sumLogdVec(double *, int);

/* Parametric */
void C_EM_mainP(int *, double *, int *, int *, double *, double *, double *,
   int *, int *, int *, double *, double *, double *, double *);
void C_ReRandP_new(int *, int *, double *, int *, int *, double *, double *,
   double *, int *, int *, double *, int *, int *, double *, double *, double *);

/* Non-parametric */
void C_EM_mainNP(int *, double *, int *, int *, int *, double *, double *,
    int *, double *, int *, int *, int *, int *, int *,
    int *, double *, int *, double *, int *, double *, double *,
    int *, double *, double *, double *, double *, double *, 
    double *, double *);
void C_ReRandNP(int *, int *, double *, int *, int *, double *, double *,
    int *, double *, double *, int *, int *, int *, double *,
    int *, int *, int *, int *, double *,
    double *, int *, int *, double *, int *, double *);
void C_LRT_NP(int *, double *, int *, int *, int *, double *, double *,
    int *, double *, int *, int *, int *, int *, int *,
    int *, double *, int *, double *, int *, double *, double *,
    int *, double *, double *, double *, double *, double *, 
    double *, double *, double *);







#endif







