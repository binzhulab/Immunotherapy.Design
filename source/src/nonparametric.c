#include "./source.h"


static void M_fn(t, t_all, N, Nfail, ret)
double *t, *t_all, **ret;
int N, Nfail;
{
  /*
  M=matrix(0, nrow=n.failure, ncol=N)
  for (l in 2:(n.failure+1))
    for (i in 1:N)
    {
      M[(l-1),i]=(min(t.all[l], t[i])-t.all[l-1])*ifelse(t[i]>=t.all[l-1], 1 ,0)
    } 
  */


  int i, k, km1; 
  double tmin, tk, tkm1, ti, test, **prow, *pd1, *pt;

  /* Return by row */
  for (k=1, prow=ret; k<Nfail+1; k++, prow++) {
    km1  = k - 1;
    tk   = t_all[k];
    tkm1 = t_all[km1];
    test = tkm1 - EQUAL_EPS;

    for (i=0, pt=t, pd1=*prow; i<N; i++, pt++, pd1++) {
      ti = *pt;
      if (ti >= test) {
        tmin = MIN(tk, ti);
        *pd1 = tmin - tkm1;
      } else {
        *pd1 = 0.0;
      }
    }
  }

} /* END: M_fn */


static void M1_fn(Mfn, lambda, n, nfail, ret)
double **Mfn, *lambda, *ret;
int n, nfail;
{
  /* as.vector(t(lambda)%*%Mfn)
     Mfn is nfail by n, lambda is length nfail,
     ret length n */
  int i, j;
  double *pd1, *pret, sum;

  for (i=0, pret=ret; i<n; i++, pret++) {
    sum = 0.0;
    for (j=0, pd1=lambda; j<nfail; j++, pd1++) sum += *pd1 * Mfn[j][i]; 
    *pret = sum;
  }

} /* END: M1_fn */

static void Mstar_fn(tall, tstar, tallLEtstar, nfail, ret)
double *tall, tstar, *ret;
int nfail, *tallLEtstar;
{
  /*
  tvec   <- t.all[1:n.failure]
  tmp1   <- tvec <= tstar
  tmp2   <- pmin(rep(tstar, n.failure), t.all[-1]) - tvec
  M      <- tmp1*tmp2
  dim(M) <- c(n.failure, 1)
  */
  int i, *pi1;
  double *pret;

  pi1  = tallLEtstar;
  pret = ret;
  for (i=0; i<nfail; i++, pi1++, pret++) {
    if (*pi1) {
      *pret = MIN(tstar, tall[i+1]) - tall[i];
    } else {
      *pret = 0.0;
    }
  }

} /* END: Mstar_fn */
 
static double Mstar1_fn(Mstarfn, lambda, n)
double *Mstarfn, *lambda;
int n;
{
  /* t(lambda)%*%Mstarfn */

  int i;
  double sum=0.0, *pd1, *pd2;

  for (i=0, pd1=Mstarfn, pd2=lambda; i<n; i++, pd1++, pd2++) sum += *pd1 * *pd2;

  return(sum);

} /* END: Mstar1_fn */ 

static void M2_fn(M1fn, n, Mstar1fn, ret)
double *M1fn, Mstar1fn, *ret;
int n;
{
  /* M1fn - rep(Mstar1fn, N) */
  int i;
  double *pd1, *pd2;

  pd1 = M1fn;
  pd2 = ret;
  for (i=0; i<n; i++, pd1++, pd2++) *pd2 = *pd1 - Mstar1fn;

} /* END: M2_fn */

static void A_fn(Mstarfn, Mfn, Zt, tGEtstar, n, nfail, tmpVec1, ret)
double *Mstarfn, **Mfn, *Zt, *ret, *tmpVec1;
int *tGEtstar, n, nfail;
{
  /*
  tmp0     <- Zt*(tstar <= t)
  tmp      <- tmp0
  dim(tmp) <- c(N, 1)
  mat      <- Mfn%*%(1-tmp)
  tmp2     <- tmp0*matrix(rep(Mstarfn, each=N), nrow=N, ncol=n.failure, byrow=FALSE)
  M        <- colSums(tmp2) + as.vector(mat)

  Mfn is nfail by n
  */

  int i, *pi1;
  double *pd1, *pd2, sum=0.0, d1, *pret;

  /* compute tmp0 and 1-tmp above, all we need from tmp0 is the sum */
  for (i=0, pd1=Zt, pd2=tmpVec1, pi1=tGEtstar; i<n; i++, pd1++, pd2++, pi1++) {
    d1   = *pd1 * *pi1;
    *pd2 = 1.0 - d1;
    sum += d1;
  }

  /* compute Mfn%*%(1-tmp) above and store in ret */
  matTimesVec(Mfn, tmpVec1, nfail, n, ret);

  for (i=0, pret=ret, pd1=Mstarfn; i<nfail; i++, pret++, pd1++) {
    d1    = *pd1 * sum; /* ith element of colSums(tmp2) above */
    *pret = d1 + *pret; /* ith element of colSums(tmp2) + as.vector(mat) */
  }

} /* END: A_fn */

static void B_fn(Mstarfn, Mfn, Zt, tGEtstar, n, nfail, beta, tmpVec1, ret)
double *Mstarfn, **Mfn, *Zt, *ret, *tmpVec1, beta;
int *tGEtstar, n, nfail;
{
  /*
  tmp      <- Zt*(tstar <= t)
  dim(tmp) <- c(N, 1)
  mat      <- Mfn - matrix(rep(Mstarfn, each=N), nrow=n.failure, ncol=N, byrow=TRUE)
  M        <- exp(beta)*as.vector((mat %*% tmp))
  */

  int i, j, *pi1;
  double ebeta, *pd1, *pd2, *pd3, **prow, d2, *pret, sum;

  /* Compute tmp */
  for (i=0, pd1=Zt, pd2=tmpVec1, pi1=tGEtstar; i<n; i++, pd1++, pd2++, pi1++) {
    *pd2 = *pd1 * *pi1;
  }

  ebeta = exp(beta);

  for (i=0, pret=ret, prow=Mfn, pd2=Mstarfn; i<nfail; i++, pret++, prow++, pd2++) {
    sum = 0.0;
    d2  = *pd2;
    for (j=0, pd1=*prow, pd3=tmpVec1; j<n; j++, pd1++, pd3++) sum += (*pd1 - d2) * *pd3;
    *pret = ebeta*sum;
  } 

} /* END: B_fn */

static void getLambdaNP(A1, B1, deltal, nfail, ret)
double *A1, *B1, *deltal, *ret;
int nfail;
{
  int i;
  double *pd1, *pd2, *pd3, *pret;

  for (i=0, pd1=A1, pd2=B1, pd3=deltal, pret=ret; i<nfail; i++, pd1++, pd2++, pd3++, pret++) {
    *pret = *pd3/(*pd1 + *pd2);
  }  

} /* END: getLambdaNP */

static double C_fn(Zt, eventStatus, tGEtstar, n)
double *Zt;
int *eventStatus, *tGEtstar, n;
{
  /*sum(Zt*(t>=tstar)*event_status)*/
  int i, *pi1, *pi2;
  double sum, *pd1;

  pd1 = Zt;
  pi1 = eventStatus;
  pi2 = tGEtstar;
  sum = 0.0;
  for(i=0; i<n; i++, pd1++, pi1++, pi2++) sum += *pd1 * *pi1 * *pi2;

  return(sum);

} /* END: C_fn */

static double D_fn(Zt, M2fn, tGEtstar, n)
double *Zt, *M2fn;
int *tGEtstar, n;
{
  /* ret <- sum(Zt*(t>=tstar)*M2fn) */

  int i, *pi1;
  double sum, *pd1, *pd2;

  pd1 = Zt;
  pd2 = M2fn;
  pi1 = tGEtstar;
  sum = 0.0;
  for(i=0; i<n; i++, pd1++, pi1++, pd2++) sum += *pd1 * *pi1 * *pd2;

  return(sum);

} /* END: D_fn */

static double loglik_fn(Mstar1fn, M1fn, trtEq1, tGEtstar, Zt, lambda, Cfn, Dfn, 
                        beta, logEffect, log1mEffect, n, nfail, tmpVec1, ret_loglike_marg)
double Mstar1fn, *M1fn, *Zt, *lambda, Cfn, Dfn, beta, logEffect, log1mEffect, *tmpVec1, *ret_loglike_marg;
int *trtEq1, *tGEtstar, n, nfail;
{
  /* 
  trteq1   <- trt == 1
  Xge      <- X >= tstar
  tmp      <- Zt*Xge

  r1     = sum(log(lambda))
  r2     = -sum(tmp)*Mstar1fn
  r3     = -t(M1fn)%*%(1-tmp)
  r4     = beta*C.fn(Zt, t=X, tstar, event_status = event_status)
  r5     = -(exp(beta))*D.fn(Zt, t=X, tstar, lambda, t.all, N, n.failure, M2fn=NULL)
  r6     = sum(Zt*trteq1*log(effect_p))
  r7     = sum((1-Zt)*trteq1*log(1-effect_p))
  loglik = r1 + r2 + r3 + r4 + r5 + r6 + r7  
  */

  int i, *pi1, trti;
  double ret, *pd1, *pd2, *pd3, r1, r2, r3, r4, r5, r6, r7, d1, r8, r9, dtmp;

  /* store Zt*Xge in tmpVec1 */
  for (i=0, pd1=Zt, pd2=tmpVec1, pi1=tGEtstar; i<n; i++, pd1++, pd2++, pi1++) {
    *pd2 = *pd1 * *pi1;
  }

  r1 = sumLogdVec(lambda, nfail);
  r2 = -Mstar1fn*sum_dVec(tmpVec1, n);
  r4 = beta*Cfn;
  r5 = -exp(beta)*Dfn;
  r3 = 0.0;
  r6 = 0.0;
  r7 = 0.0;
  r8 = 0.0;
  r9 = 0.0;

  /*
  Could you please fix the problem by modifying the calculation of r8 and r9 in 
  the log likelihood function as follows, that should fix the problem,
  r8       = sum(Zt*trteq1*log(effect_p/(Zt+1e-100) ))
  r9       = sum((1-Zt)*trteq1*log((1-effect_p)/(1-Zt +1e-100))))
  */


  for (i=0, pi1=trtEq1, pd1=Zt, pd2=M1fn, pd3=tmpVec1; i<n; i++, pi1++, pd1++, pd2++, pd3++) {
    trti = *pi1;
    d1   = *pd1;
    r3  += *pd2 * (1.0 - *pd3);
    r6  += d1*trti;
    r7  += (1.0 - d1)*trti; 
    if (ret_loglike_marg && trti) {
      dtmp = 1.0 - d1;
      r8  += d1*(logEffect - log(d1 + NUMERIC_ZERO));
      r9  += dtmp*(log1mEffect - log(dtmp + NUMERIC_ZERO));
    }
  }
  r3   = -r3;
  r6   = r6*logEffect;
  r7   = r7*log1mEffect;

  dtmp = r1 + r2 + r3 + r4 + r5;
  ret  = dtmp + r6 + r7;
  if (ret_loglike_marg) *ret_loglike_marg = dtmp + r8 + r9;
  
  return(ret);

} /* END: loglik_fn */

static void pdf_r_fn(Mstar1fn, M1fn, M2fn, lambdaDeltafn, tGEtstar, beta, eventStatus, ntrt, ret)
double Mstar1fn, *M1fn, *M2fn, *lambdaDeltafn, beta, *ret;
int *eventStatus, ntrt, *tGEtstar;
{
  /*
  r1    <- apply(lambda^Deltafn, 2, prod)[vec]
  tmp1  <- (t >= tstar)[vec]
  r2    <- exp(beta*event_status[vec]*tmp1)
  m1fn  <- M1.fn(lambda,t,t.all,N, n.failure, Mfn=NULL)
  r3    <- exp(-(1-tmp1)*m1fn)
  r4    <- exp((-Mstar1fn-exp(beta)*M2fn)*tmp1)
  pdf.r <- (r1*r2*r3*r4)
  */

  int i, *pi1, *pi2, i2;
  double r1, r2, r3, r4, ebeta, *pret, *pd1, *pd2, *pd3;

  ebeta = exp(beta);
  pi1   = eventStatus;
  pi2   = tGEtstar;
  pd1   = M1fn;
  pd2   = M2fn;
  pd3   = lambdaDeltafn;
  pret  = ret;
  for (i=0; i<ntrt; i++, pi1++, pi2++, pd1++, pd2++, pd3++, pret++) {
    i2    = *pi2;
    r1    = *pd3;
    r2    = exp(beta * *pi1 * i2);
    r3    = exp(-(1.0 - i2)* *pd1);
    r4    = exp((-Mstar1fn - ebeta * *pd2)*i2);
    *pret = r1*r2*r3*r4;
  }

} /* END: pdf_r_fn */

static void pdf_nr_fn(lambdaDeltafn, M1fn, ntrt, ret)
double *lambdaDeltafn, *M1fn, *ret;
int ntrt;
{
  /*
  r1     <- apply(lambda^Deltafn, 2, prod)[1:N]
  r5     <- exp(-M1.fn(lambda,t,t.all,N, n.failure, Mfn=NULL))
  pdf.nr <- r1*r5
  */

  int i;
  double *pd1, *pd2, *pret;

  for (i=0, pd1=lambdaDeltafn, pd2=M1fn, pret=ret; i<ntrt; i++, pd1++, pd2++, pret++) {
    *pret = *pd1 * exp(- *pd2);
  }

} /* END: pdf_nr_fn */

static void getZtNP(effect, pdfr, pdfnr, ntrt, ret)
double effect, *pdfr, *pdfnr, *ret;
int ntrt;
{
  /*
  Zt.trt  <- (effect_p*pdf.r/(effect_p*pdf.r + (1-effect_p)*pdf.nr))[1:n_trt]
  Zt      <- c(Zt.trt, Zt[trtEq0])
  */

  int i;
  double oneMinusEffect, tmp, *pd1, *pd2, *pret;

  oneMinusEffect = 1.0 - effect;

  /* NOTE: Zt does not change for non-treated subjects. Data is ordered */
  for (i=0, pret=ret, pd1=pdfr, pd2=pdfnr; i<ntrt; i++, pret++, pd1++, pd2++) {
    tmp = effect * *pd1;
    *pret = tmp/(tmp + oneMinusEffect * *pd2);
  }
  
} /* END: getZtNP */

static void getPosMat(csDeltafn, n, posVec, np, ret)
int *csDeltafn, n, *posVec, np, **ret;
{
  /* ret should be n by maxPosVec */

  int i, j, k, m;

  k = 0;
  for (i=0; i<n; i++) {
    m = csDeltafn[i];
    if (m) {
      if (k >= np) error("INTERNAL CODING ERROR in getLambdaDelta\n");
      for (j=0; j<m; j++) {
        ret[i][j] = posVec[k];
        k++;
      }
    }
  }

} /* END: getPosMat */

static void getLambdaDelta(lambda, csDeltafn, n, posMat, ret)
double *lambda, *ret;
int *csDeltafn, n, **posMat;
{
  /* csDeltafn are colSums fom Deltafn from R code.
     posMat[i] are c-positions indices of lambda for sub i.
     ret should have length n.
  */

  int i, j, *pi1, m, *posVec;
  double *pret, prod;

  for (i=0, pi1=csDeltafn, pret=ret; i<n; i++, pi1++, pret++) {
    m      = *pi1;
    posVec = posMat[i];
    if (!m) {
      prod = 1.0;
    } else if (m == 1) {
      prod = lambda[posVec[0]];
    } else if (m == 2) {
      prod = lambda[posVec[0]]*lambda[posVec[1]];
    } else if (m == 3) {
      prod = lambda[posVec[0]]*lambda[posVec[1]]*lambda[posVec[2]];
    } else {
      prod = 1;
      for (j=0; j<m; j++) prod = prod*lambda[posVec[j]];
    }
    *pret = prod;
  }

} /* END: getLambdaDelta */

static int EM_mainNP(n, Zt, eventStatus, tGEtstar, nfail, tall, tstar,
                     tallLEtstar, X, trtEq1, ntrt, csDeltafn, posMat, 
                     effect, maxiter, stopTol, print, lambda, deltal,
                     ret_beta, ret_lambda, pdfr, pdfnr, ret_Zt, ret_niter, 
                     ret_loglike, Mfn, Mstarfn, ret_loglike_marg)
double *Zt, *lambda, *ret_lambda, *tall, tstar, *X, stopTol, *ret_beta, 
       *pdfr, *pdfnr, effect, *ret_Zt, *deltal, *ret_loglike, **Mfn, *Mstarfn,
       *ret_loglike_marg;
int n, *eventStatus, *tGEtstar, nfail, *tallLEtstar, ntrt, *csDeltafn, 
    **posMat, maxiter, print, *trtEq1, *ret_niter;
{
  int iter=0, conv=0, Mflag=0;
  double Cfn, Dfn, *M1fn, *M2fn, Mstar1fn, beta, *Afn, *Bfn, *tmpVec1;
  double *lambdaDeltafn, diff, loglik, logEffect, log1mEffect;
  double mybeta=0.0, mybeta0=0.0, loglike=0.0, loglike0=0.0;
 
  /* Error check */
  if ((!Mfn && Mstarfn) || (Mfn && !Mstarfn)) error("INTERNAL CODING ERROR in EM_mainNP"); 
  if (Mfn && Mstarfn) Mflag = 1;
  
  if (!Mflag) {
    Mfn         = dMat_alloc(nfail, n, 0, -9999.0);
    Mstarfn     = dVec_alloc(nfail, 0, -9999.0);
  }
  M1fn          = dVec_alloc(n, 0, -9999.0);
  M2fn          = dVec_alloc(n, 0, -9999.0);
  Afn           = dVec_alloc(nfail, 0, -9999.0);
  Bfn           = dVec_alloc(nfail, 0, -9999.0);
  lambdaDeltafn = dVec_alloc(n, 0, -9999.0);
  tmpVec1       = dVec_alloc(n, 0, -9999.0);

  copy_dVec(ret_lambda, lambda, nfail);
  copy_dVec(ret_Zt, Zt, n);
  M_fn(X, tall, n, nfail, Mfn);
  Mstar_fn(tall, tstar, tallLEtstar, nfail, Mstarfn); 

  logEffect   = log(effect);
  log1mEffect = log(1.0 - effect);

  while(iter < maxiter) {
    iter++;
    Cfn = C_fn(ret_Zt, eventStatus, tGEtstar, n);  
    Mstar1fn = Mstar1_fn(Mstarfn, ret_lambda, nfail);
    M1_fn(Mfn, ret_lambda, n, nfail, M1fn);
    M2_fn(M1fn, n, Mstar1fn, M2fn);

    Dfn    = D_fn(ret_Zt, M2fn, tGEtstar, n);
    mybeta = log(Cfn/Dfn);

    if (!R_FINITE(mybeta)) break;
 
    A_fn(Mstarfn, Mfn, ret_Zt, tGEtstar, n, nfail, tmpVec1, Afn);
    B_fn(Mstarfn, Mfn, ret_Zt, tGEtstar, n, nfail, mybeta, tmpVec1, Bfn);
    getLambdaNP(Afn, Bfn, deltal, nfail, ret_lambda);

    /* Update objects depending on lambda */
    Mstar1fn = Mstar1_fn(Mstarfn, ret_lambda, nfail); 
    M1_fn(Mfn, ret_lambda, n, nfail, M1fn);
    M2_fn(M1fn, n, Mstar1fn, M2fn);
    getLambdaDelta(ret_lambda, csDeltafn, n, posMat, lambdaDeltafn);
    pdf_r_fn(Mstar1fn, M1fn, M2fn, lambdaDeltafn, tGEtstar, mybeta, eventStatus, ntrt, pdfr);
    pdf_nr_fn(lambdaDeltafn, M1fn, ntrt, pdfnr);
    getZtNP(effect, pdfr, pdfnr, ntrt, ret_Zt);

    /* Update obects depending on Zt */
#if EM_STOP_CRITERIA == 0
    Cfn = C_fn(ret_Zt, eventStatus, tGEtstar, n);
    Dfn = D_fn(ret_Zt, M2fn, tGEtstar, n);

    /* Note: old beta gets passed in */
    loglike = loglik_fn(Mstar1fn, M1fn, trtEq1, tGEtstar, ret_Zt, ret_lambda, Cfn, Dfn, 
                        mybeta, logEffect, log1mEffect, n, nfail, tmpVec1, NULL);
    if (iter > 1) {
      diff = fabs(loglike - loglike0);
      if (diff <= stopTol) {  
        conv = 1;
        break;
      }
    }
    if (print > 1) {
      if (iter > 1) {
        Rprintf("Iter=%d, loglike=%g, diff=%g\n", iter, loglike, diff);
      } else {
        Rprintf("Iter=%d, loglike=%g\n", iter, loglike);
      }     
    }  
    loglike0 = loglike;
#else
    if (iter > 1) {
      diff = fabs(mybeta - mybeta0);
      if (diff <= stopTol) {  
        conv   = 1;
        break;
      }
    }
    if (print > 1) {
      if (iter > 1) {
        Rprintf("Iter=%d, beta=%g, diff=%g\n", iter, mybeta, diff);
      } else {
        Rprintf("Iter=%d, beta=%g\n", iter, mybeta);
      }     
    }  
    mybeta0 = mybeta;
#endif  
  
  } /* END: while */

  if (print) {
    if (conv) {
      Rprintf("EM algoritm converged in %d iteration(s)\n", iter);
    } else {
      Rprintf("EM algoritm did not converge\n");
    }
  }

  if (conv) {
    /* Compute marginal log-likelihood */
    *ret_loglike = loglik_fn(Mstar1fn, M1fn, trtEq1, tGEtstar, ret_Zt, ret_lambda, Cfn, Dfn, 
                    mybeta, logEffect, log1mEffect, n, nfail, tmpVec1, ret_loglike_marg);
  }

  if (!Mflag) {
    matrix_free((void **) Mfn, nfail);
    free(Mstarfn);
  }
  free(M1fn);
  free(M2fn);
  free(Afn);
  free(Bfn);
  free(lambdaDeltafn);
  free(tmpVec1);

  *ret_beta    = mybeta;
  *ret_niter   = iter;

  return(conv);

} /* END: EM_mainNP */

void C_EM_mainNP(pn, Zt, eventStatus, tGEtstar, pnfail, tall, ptstar,
                     tallLEtstar, X, trt, pntrt, csDeltafn, posVec, pnposVec, 
                     pMaxPosVec, peffect, pmaxiter, pstopTol, pprint, lambda, deltal,
                     ret_conv, ret_beta, ret_lambda, pdfr, pdfnr, ret_Zt, 
                     ret_loglike, ret_loglike_marg)
double *Zt, *lambda, *ret_lambda, *tall, *ptstar, *X, *pstopTol, *ret_beta, 
       *pdfr, *pdfnr, *peffect, *ret_Zt, *deltal, *ret_loglike, *ret_loglike_marg;
int *pn, *eventStatus, *tGEtstar, *pnfail, *tallLEtstar, *pntrt, *csDeltafn, 
    *posVec, *pnposVec, *pmaxiter, *pprint, *trt, *ret_conv, *pMaxPosVec;
{
  int n, np, maxPosVec, **posMat, niter;

  n         = *pn;
  np        = *pnposVec;
  maxPosVec = *pMaxPosVec;

  /* Matrix for the indices of lambda to multiply together. It will use
     a little more memory than what is actually needed, but makes 
     permuting much simpler. */
  posMat = iMat_alloc(n, maxPosVec, 0, 0);
  getPosMat(csDeltafn, n, posVec, np, posMat);

  *ret_conv = EM_mainNP(n, Zt, eventStatus, tGEtstar, *pnfail, tall, *ptstar,
                     tallLEtstar, X, trt, *pntrt, csDeltafn, posMat,
                     *peffect, *pmaxiter, *pstopTol, *pprint, lambda, deltal,
                     ret_beta, ret_lambda, pdfr, pdfnr, ret_Zt, &niter,
                     ret_loglike, NULL, NULL, ret_loglike_marg);

  matrix_free((void **) posMat, n);

  return;

} /* END:  C_EM_mainNP */

static void orderDataNP(trtPerm, X0, eventStatus0, tGEtstar0, csDeltafn0, posMat0, Zt0,
       n, ntrt, nMaxPos, trt, X, eventStatus, tGEtstar, csDeltafn, posMat, Zt)
double *X, *X0, *Zt0, *Zt;
int *trtPerm, *trt, *eventStatus, *tGEtstar, *eventStatus0, *tGEtstar0, n, ntrt,
     **posMat0, **posMat, *csDeltafn0, *csDeltafn, nMaxPos;
{
  int i, j, trtIndex, nonTrtIndex, flag;

  trtIndex    = 0;
  nonTrtIndex = ntrt; /* Begin non-trt subs at position ntrt */
  flag        = (nMaxPos == 1) ? 1 : 0;

  for (i=0; i<n; i++) {
    if (trtPerm[i]) {
      trt[trtIndex]         = 1;
      X[trtIndex]           = X0[i];
      eventStatus[trtIndex] = eventStatus0[i];
      tGEtstar[trtIndex]    = tGEtstar0[i];
      csDeltafn[trtIndex]   = csDeltafn0[i];
      /*Zt[trtIndex]          = Zt0[i];*/
      if (flag) {
        posMat[trtIndex][0] = posMat0[i][0];
      } else {
        for (j=0; j<nMaxPos; j++) posMat[trtIndex][j] = posMat0[i][j];
      }
      trtIndex++;
    } else {
      trt[nonTrtIndex]         = 0;
      X[nonTrtIndex]           = X0[i];
      eventStatus[nonTrtIndex] = eventStatus0[i];
      tGEtstar[nonTrtIndex]    = tGEtstar0[i];
      csDeltafn[nonTrtIndex]   = csDeltafn0[i];
      /*Zt[nonTrtIndex]          = Zt0[i];*/
      if (flag) {
        posMat[nonTrtIndex][0] = posMat0[i][0];
      } else {
        for (j=0; j<nMaxPos; j++) posMat[nonTrtIndex][j] = posMat0[i][j];
      }
      nonTrtIndex++;
    } 
  }

} /* END: orderDataNP */

static void getInitProbResp(trt, n, ret)
int *trt, n;
double *ret;
{
  int i, *pi1;
  double *pd; 

  for (i=0, pi1=trt, pd=ret; i<n; i++, pi1++, pd++) {
    if (*pi1) {
      *pd = 0.5;
    } else {
      *pd = 0.0;
    }  
  }

} /* END: getInitProbResp */

static int ReRandNP(num_rand, n, eventStatus0, tGEtstar0, nfail, tall, tstar, lambda,
                     tallLEtstar, X0, trt0, ntrt, csDeltafn0, posMat0, beta_obs,
                     effect, maxiter, stopTol, print, nMaxPos, Zt0, deltal, ret_p)
int num_rand, n, nfail, ntrt,  *trt0, *eventStatus0, maxiter, print, *tGEtstar0,
     *tallLEtstar, *csDeltafn0, **posMat0, nMaxPos;
double *X0, effect, tstar, stopTol, *tall, *lambda, beta_obs, *Zt0, *ret_p, *deltal;
{
  double *ret_lambda, ret_beta, *ret_pdfr, *ret_pdfnr, *X, *Zt, *ret_Zt, 
         ret_loglike, ret_loglike2, ret_loglike3;
  int iter, conv, sumNrand=0, sumGTobs=0, ret_niter, sum_niter=0;
  int *trtPerm, *trt,  **posMat, *csDeltafn, *eventStatus, *tGEtstar;

  *ret_p     = -1.0;
  
  trt         = iVec_alloc(n, 0, -9999);
  trtPerm     = iVec_alloc(n, 0, -9999);
  eventStatus = iVec_alloc(n, 0, -9999);
  tGEtstar    = iVec_alloc(n, 0, -9999);
  ret_lambda  = dVec_alloc(nfail, 0, -9999.0);
  ret_pdfr    = dVec_alloc(ntrt, 0, -9999.0);
  ret_pdfnr   = dVec_alloc(ntrt, 0, -9999.0);
  ret_Zt      = dVec_alloc(n, 0, -9999.0);
  X           = dVec_alloc(n, 0, -9999.0);
  csDeltafn   = iVec_alloc(n, 0, -9999);
  posMat      = iMat_alloc(n, nMaxPos, 1, -9999);

  for (iter=0; iter<num_rand; iter++) {
    permute_iVec(trt0, n, trtPerm);

    orderDataNP(trtPerm, X0, eventStatus0, tGEtstar0, csDeltafn0, posMat0, Zt0, 
                n, ntrt, nMaxPos, trt, X, eventStatus, tGEtstar, csDeltafn, posMat, Zt);

    conv = EM_mainNP(n, Zt0, eventStatus, tGEtstar, nfail, tall, tstar,
                     tallLEtstar, X, trt, ntrt, csDeltafn, posMat, 
                     effect, maxiter, stopTol, print, lambda, deltal, 
                     &ret_beta, ret_lambda, ret_pdfr, ret_pdfnr, ret_Zt, 
                     &ret_niter, &ret_loglike, NULL, NULL, NULL);

    if (conv) {
      sumNrand++;
      sum_niter += ret_niter;
      if (ret_beta > beta_obs) sumGTobs++;
    }
  }

  if (sumNrand) *ret_p = 1.0 - ((double) sumGTobs)/((double) sumNrand);

  free(trt);
  free(trtPerm);
  free(eventStatus);
  free(tGEtstar);
  free(ret_lambda);
  free(ret_pdfr);
  free(ret_pdfnr);
  free(ret_Zt);
  free(X);
  free(csDeltafn);
  matrix_free((void **) posMat, n);
 
  return(sumNrand);

} /* END: ReRandNP */

void C_ReRandNP(pnum_rand, pn, X, trt, eventStatus, peffect_p, ptstar, 
               pnfail, tall, lambda, tGEtstar, tallLEtstar, pntrt, pbeta_obs, 
               csDeltafn, posVec, pnposVec, pnMaxPos, Zt,
               pstopTol, pmaxiter, pprint, deltal, ret_nrand, ret_p)
int *pnum_rand, *pn, *trt, *eventStatus, *pmaxiter, *pprint, *ret_nrand, *tGEtstar,
    *pnfail, *tallLEtstar, *pntrt, *csDeltafn, *pnMaxPos, *posVec, *pnposVec;
double *X, *peffect_p, *ptstar, *pstopTol, *ret_p, *tall, *lambda, *Zt,
       *pbeta_obs, *deltal;
{
  int n, np, nMaxPos, **posMat;

  GetRNGstate();

  n       = *pn;
  np      = *pnposVec;
  nMaxPos = *pnMaxPos;

  /* Matrix for the indices of lambda to multiply together. It will use
     a little more memory than what is actually needed, but makes 
     permuting much simpler. */
  posMat = iMat_alloc(n, nMaxPos, 1, -9999);
  getPosMat(csDeltafn, n, posVec, np, posMat);

  *ret_nrand = ReRandNP(*pnum_rand, n, eventStatus, tGEtstar, *pnfail, tall, 
                        *ptstar, lambda, tallLEtstar, X, trt, *pntrt, csDeltafn, 
                        posMat, *pbeta_obs, *peffect_p, *pmaxiter, *pstopTol, 
                        *pprint, nMaxPos, Zt, deltal, ret_p);

  matrix_free((void **) posMat, n);

  PutRNGstate();  

  return;

} /* END: C_ReRandP */


/* Likelihood at Zt = 0, beta = 0 */
static double loglik_fn_0(M1fn, ntrtEq1, lambda, log1mEffect, n, nfail)
double *M1fn, *lambda, log1mEffect;
int ntrtEq1, n, nfail;
{
  /* 
  trteq1   <- trt == 1
  Xge      <- X >= tstar
  tmp      <- Zt*Xge

  r1     = sum(log(lambda))
  r2     = -sum(tmp)*Mstar1fn
  r3     = -t(M1fn)%*%(1-tmp)
  r4     = beta*C.fn(Zt, t=X, tstar, event_status = event_status)
  r5     = -(exp(beta))*D.fn(Zt, t=X, tstar, lambda, t.all, N, n.failure, M2fn=NULL)
  r6     = sum(Zt*trteq1*log(effect_p))
  r7     = sum((1-Zt)*trteq1*log(1-effect_p))
  loglik = r1 + r2 + r3 + r4 + r5 + r6 + r7  
  */

  int i, *pi1, trti;
  double ret, *pd2, r1, r3;

  r1 = sumLogdVec(lambda, nfail);
  r3 = 0.0;
  /*r9 = log1mEffect*ntrtEq1;*/

  for (i=0, pd2=M1fn; i<n; i++, pd2++) r3 += *pd2;
 
  r3  = -r3;
  ret = r1 + r3;
  
  return(ret);

} /* END: loglik_fn */


void C_LRT_NP(pn, Zt, eventStatus, tGEtstar, pnfail, tall, ptstar,
                     tallLEtstar, X, trt, pntrt, csDeltafn, posVec, pnposVec, 
                     pMaxPosVec, peffect, pmaxiter, pstopTol, pprint, lambda, deltal,
                     ret_conv, ret_beta, ret_lambda, pdfr, pdfnr, ret_Zt, 
                     ret_loglike, ret_loglike0, ret_p)
double *Zt, *lambda, *ret_lambda, *tall, *ptstar, *X, *pstopTol, *ret_beta, 
       *pdfr, *pdfnr, *peffect, *ret_Zt, *deltal, *ret_loglike, *ret_loglike0, *ret_p;
int *pn, *eventStatus, *tGEtstar, *pnfail, *tallLEtstar, *pntrt, *csDeltafn, 
    *posVec, *pnposVec, *pmaxiter, *pprint, *trt, *ret_conv, *pMaxPosVec;
{
  int n, np, maxPosVec, **posMat, niter, conv, nfail, i;
  double **Mfn, *Mstarfn, *Afn, *Bfn, *tmpVec1, *M1fn, tstar, mybeta, log1mEffect, test;
  double LL_marg;

  n           = *pn;
  np          = *pnposVec;
  maxPosVec   = *pMaxPosVec;
  nfail       = *pnfail;
  tstar       = *ptstar;
  log1mEffect = log(1.0 - *peffect); 

  posMat  = iMat_alloc(n, maxPosVec, 0, 0);
  Mfn     = dMat_alloc(nfail, n, 0, -9999.0);
  Mstarfn = dVec_alloc(nfail, 0, -9999.0);

  /* Matrix for the indices of lambda to multiply together. It will use
     a little more memory than what is actually needed, but makes 
     permuting much simpler. */
  getPosMat(csDeltafn, n, posVec, np, posMat);

  /* Use marginal likelihoods for LRT test */
  conv = EM_mainNP(n, Zt, eventStatus, tGEtstar, *pnfail, tall, *ptstar,
                     tallLEtstar, X, trt, *pntrt, csDeltafn, posMat,
                     *peffect, *pmaxiter, *pstopTol, *pprint, lambda, deltal,
                     ret_beta, ret_lambda, pdfr, pdfnr, ret_Zt, &niter,
                     ret_loglike, Mfn, Mstarfn, &LL_marg);
  if (!R_FINITE(*ret_loglike)) conv = 0;

  if (conv) {
    tmpVec1 = dVec_alloc(n, 0, -9999.0);
    Afn     = dVec_alloc(nfail, 0, -9999.0);
    Bfn     = dVec_alloc(nfail, 1, 0.0); /* Is a zero vector at Zt = 0 */
    M1fn    = dVec_alloc(n, 0, -9999.0);
   
    /* Set beta and Zt to 0 */ 
    mybeta = 0.0;
    for (i=0; i<n; i++) ret_Zt[i] = 0.0;
    
    A_fn(Mstarfn, Mfn, ret_Zt, tGEtstar, n, nfail, tmpVec1, Afn);
    getLambdaNP(Afn, Bfn, deltal, nfail, ret_lambda);

    /* Update needed objects that depend on lambda */
    M1_fn(Mfn, ret_lambda, n, nfail, M1fn);

    /* Loglike at Zt = 0 and beta = 0 */
    *ret_loglike0 = loglik_fn_0(M1fn, *pntrt, ret_lambda, log1mEffect, n, nfail);

    *ret_loglike = LL_marg;

    test         = 2.0*(*ret_loglike - *ret_loglike0);
    *ret_p       = pchisq(test, 1.0, 0, 0);


/*
double Cfn=0.0, Dfn=0.0;
double Mstar1fn = Mstar1_fn(Mstarfn, ret_lambda, nfail); 
double dtmp = loglik_fn(Mstar1fn, M1fn, trt, tGEtstar, ret_Zt, ret_lambda, Cfn, Dfn, 
                        mybeta, log(*peffect), log1mEffect, n, nfail, tmpVec1, &LL_marg);
Rprintf("%g %g\n", *ret_loglike0, LL_marg);
*/

    
    
    free(tmpVec1);
    free(Afn);
    free(Bfn);
    free(M1fn);
  }

  matrix_free((void **) posMat, n);
  matrix_free((void **) Mfn, nfail);
  free(Mstarfn);

  *ret_conv = conv;

  return;

} /* END:  C_LRT_NP */









