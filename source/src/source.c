#include "./source.h"


/* Printing */
void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %g ", vec[i]);
  }
  Rprintf("\n");
}
void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %d ", vec[i]);
  }
  Rprintf("\n");
}

void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %g ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_iMat(mat, nr, nc, name)
int **mat;
int nr, nc;
char name[10];
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %d ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

/* Compare */
void compare_dVec(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double maxdiff, diff;

  maxdiff = -1.0;
  for (i=0; i<n; i++) {
    diff = fabs(v1[i] - v2[i]);
    if (diff > maxdiff) maxdiff = diff;
  }
  Rprintf("ddiff=%g\n", maxdiff);

} 

void compare_iVec(v1, v2, n)
int *v1, *v2;
int n;
{
  int i;
  int maxdiff, diff;

  maxdiff = -1;
  for (i=0; i<n; i++) {
    diff = fabs(v1[i] - v2[i]);
    if (diff > maxdiff) maxdiff = diff;
  }
  Rprintf("idiff=%d\n", maxdiff);

} 

void compare_iMat(v1, v2, nr, nc)
int **v1, **v2;
int nr, nc;
{
  int i, j;
  int maxdiff, diff;

  maxdiff = -1;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      diff = fabs(v1[i][j] - v2[i][j]);
      if (diff > maxdiff) maxdiff = diff;
    }
  }
  Rprintf("idiff=%d\n", maxdiff);

}

/* Memory */
double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: iVec_alloc */

/* Function to allocate a double matrix */
double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  if (ncol > 0) {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);
  } else {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = NULL;
  }

  return(mat);

} /* END: dMat_alloc */

/* Function to allocate a integer matrix */
int ** iMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag, initVal;
{
  int i, **mat, **ptr;

  mat = (int **) malloc(nrow*sizeof(int *));
  CHECK_MEM(mat);
  if (ncol > 0) {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = iVec_alloc(ncol, initFlag, initVal);
  } else {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = NULL;
  }

  return(mat);

} /* END: iMat_alloc */

/* Function to free a matrix */
void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* matrix alg */
void copy_iVec(v1, v2, n)
int *v1, *v2, n;
{
  int i, *p1, *p2;
  
  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) *p1 = *p2;

} /* END: copy_iVec */

void copy_dVec(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double *p1, *p2;
  
  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) *p1 = *p2;

} /* END: copy_dVec */


void permute_iVec(vec, n, ret)
int *vec, *ret, n;
{
  /*
  To shuffle an array a of n elements (indices 0..n-1):
  for i from n-1 downto 1 do
     j = random integer such that 0 <= j <= i
     exchange a[j] and a[i]
  */

  int i, j, tmp;

  copy_iVec(ret, vec, n);
  
  for (i=n-1; i>0; i--) {
    j = floor(runif(0.0, i+1.0));
    if (j > i) j = i;

    tmp    = ret[i];
    ret[i] = ret[j];
    ret[j] = tmp;
  }

} /* END: permute_iVec */

double dotProd(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double *pd1, *pd2, sum=0.0;

  for (i=0, pd1=v1, pd2=v2; i<n; i++, pd1++, pd2++) sum += *pd1 * *pd2;

  return(sum);

} /* END: dotProd */

void matTimesVec(mat, vec, nr, nc, ret)
double **mat, *vec, *ret;
int nr, nc;
{
  int i;
  double **prow, *pret;

  for (i=0, pret=ret, prow=mat; i<nr; i++, pret++, prow++) *pret = dotProd(*prow, vec, nc);

}

double sum_dVec(v, n)
double *v;
int n;
{
  int i;
  double sum=0.0;

  for (i=0; i<n; i++) sum += v[i];
  return(sum);

} 

int sum_iVec(v, n)
int *v;
int n;
{
  int i;
  int sum=0;

  for (i=0; i<n; i++) sum += v[i];
  return(sum);

} 

double sumLogdVec(vec, n)
double *vec;
int n;
{
  int i;
  double sum=0.0, *p;

  for (i=0, p=vec; i<n; i++, p++) sum += log(*p);

  return(sum);

} /* END: sumLogdVec */








