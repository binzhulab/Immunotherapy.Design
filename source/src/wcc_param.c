#include "./source.h"

static void EM_setupObjP(n, X, eventStatus, t1, xt1, xt12, exst1, Xminust1, nevents)
int n, *eventStatus, *exst1, *nevents;
double *X, t1, *xt1, *xt12, *Xminust1;
{
  int i, tmp, evi, sum=0;
  double Xi;

  for (i=0; i<n; i++) { 
     Xi                  = X[i];
     tmp                 = 0;
     evi                 = eventStatus[i];
     if (Xi > t1) tmp    = 1;
     Xminust1[i]         = Xi - t1;
     exst1[i]            = tmp*evi;
     sum                += evi;
     if (Xi < t1) {
       xt1[i] = Xi;
     } else {
       xt1[i] = t1;
     }
     xt12[i]             = Xminust1[i]*tmp;
  }
  *nevents = sum;

  return;

} /* END: EM_setupObjP */

static void EM_setupTrtP(n, trt, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect)
int n, *trt;
double logEff, log1mEff, *trt1LogEffect, *trt1Log1mEffect;
{
  int i, trti;

  for (i=0; i<n; i++) { 
    trti                = trt[i];
    trt1LogEffect[i]    = trti*logEff;
    trt1Log1mEffect[i]  = trti*log1mEff;
  }

  return;

} /* END: EM_setupTrtP */

static void setupZt(n, trt, ret)
int n, *trt;
double *ret;
{
  int i;

  for (i=0; i<n; i++) { 
    if (trt[i]) {
      ret[i] = 0.5;
    } else {
      ret[i] = 0.0;
    }
  }

  return;

} /* END: EM_setupTrtP */

static double loglik_p_new(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
             trt1Log1mEffect, lambda, h_nr, trt, log_Wt)
int n, nevents, *exst1, *trt;
double *X, *Zt, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, lambda, h_nr;
double *log_Wt;
{
  int i;
  double ret, r1=0.0, r2, r3=0.0, r4=0.0, r5=0.0, r6=0.0, r7=0.0, Zti, oneMinusZti;
  double r5_trt0=0.0;

  for (i=0; i<n; i++) { 
    if (trt[i]) {
      r1 += log_Wt[i];
    } else {
      r5 += X[i];
    }
  }

  /* The complete-data likelihood is only correct for control. */
  // ret = r1*log(lambda) + r2*log(h_nr) - h_nr*(r3 + lambda*r4 + r5) + r6 + r7;

  /* This is the full logL because status of responder is unknown. */
  ret = r1 - h_nr * r5;

  return(ret);

} /* END: loglik_p_new */


static void EM_setupTrtP_new(n, trt, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect,
                             trtEq1Rows)
int n, *trt, *trtEq1Rows;
double logEff, log1mEff, *trt1LogEffect, *trt1Log1mEffect;
{
  int i, row=0;

  for (i=0; i<n; i++) { 
    if (trt[i]) {
      trt1LogEffect[i]    = logEff;
      trt1Log1mEffect[i]  = log1mEff;
      trtEq1Rows[row++]   = i;
    } else {
      trt1LogEffect[i]    = 0.0;
      trt1Log1mEffect[i]  = 0.0;
    }
  }

  return;

} /* END: EM_setupTrtP_new */

static void updateZt(nTrtEq1, trtEq1Rows, X, eventStatus, lambda, h_nr, t1, 
                     Xminust1, effect_p, oneMinusEffect, Zt, log_Wt)
int nTrtEq1, *trtEq1Rows, *eventStatus;
double *X, lambda, h_nr, t1, *Xminust1, effect_p, oneMinusEffect, *Zt;
double *log_Wt;
{
  int i, row, evi;
  double Xi, dtmp, vec, vec2, pdfr, pdfnr;

  dtmp = lambda*h_nr;
  for (i=0; i<nTrtEq1; i++) {
     row  = trtEq1Rows[i];
     Xi   = X[row];
     evi  = eventStatus[row];
     if (!evi) {
       vec = 1.0;
     } else {
       vec  = dtmp;
     }
     if (Xi < t1) {
       vec  = h_nr;
       vec2 = -h_nr*Xi;
     } else {
       vec2 = -h_nr*t1 - dtmp*Xminust1[row];
     }  
     pdfr  = vec*exp(vec2);
     vec   = exp(-h_nr*Xi);
     if (evi) vec = vec*h_nr;
     pdfnr   = vec;
     vec     = effect_p*pdfr;
     Zt[row] = vec/(vec + oneMinusEffect*pdfnr);
     log_Wt[row] = log(vec + oneMinusEffect*pdfnr); // See APECM or APECMa
  }
}

static double updateLambda(n, Zt, trt, exst1, xt1, xt12, X, r2, ret_h_nr)
int n, *trt, *exst1;
double *Zt, *xt1, *xt12, *X, r2, *ret_h_nr;
{
   int i;
   double lambda, h_nr, r1, r3, r4, r5, Zti;

   r1 = 0.0; 
   r3 = 0.0;
   r4 = 0.0;
   r5 = 0.0;
   for (i=0; i<n; i++) {
     Zti = Zt[i];
     if (trt[i]) {
       r1 += Zti*exst1[i];
       r3 += Zti*xt1[i];
       r4 += Zti*xt12[i];
     }
     r5 += (1.0 - Zti)*X[i];   
   }
   h_nr   = (r2 - r1)/(r3 + r5);
   lambda = r1/(r4*h_nr);
   if (!R_FINITE(lambda)) error("ERROR: non-finite value for lambda");
   *ret_h_nr = h_nr;

   return(lambda);
}

static int EM_loopP_new(n, nevents, maxiter, eventStatus, exst1, X, xt1, xt12, 
                    effect_p, trt1LogEffect, trt1Log1mEffect, stopTol, t1, print, method,
                    trtEq1Rows, nTrtEq1, trt, 
                    Xminust1, ret_h_nr, ret_lambda, Zt, LAMBDA_OBS, ret_iter, ret_loglike)
int n, nevents, maxiter, *eventStatus, *exst1, print, method, *trtEq1Rows, nTrtEq1, *trt, *ret_iter;
double *X, *xt1, *xt12, effect_p, *trt1LogEffect, *trt1Log1mEffect, stopTol, *Zt,
       t1, *Xminust1, *ret_h_nr, *ret_lambda, LAMBDA_OBS, *ret_loglike;
{
  int iter=0, conv=0, lambdaFlag, lambdaInc=0, lambdaDec=0,
       lambdaIncIter=0, lambdaDecIter=0;
  double dtmp, r1, r3, r2, h_nr, oneMinusEffect;
  double loglik, loglik0=-1.0, lambda=0.0, lambda0, diff;

  double *log_Wt;
  log_Wt = dVec_alloc(n, 0, 0.0);

  r2             = (double) nevents;
  oneMinusEffect = 1.0 - effect_p;
  lambdaFlag     = (LAMBDA_OBS  > DBL_NOTMISS_GT) && (method == 1); 

  /* Get good initial estimates for permutation tests */
  if (LAMBDA_OBS  > DBL_NOTMISS_GT) {  
    updateZt(nTrtEq1, trtEq1Rows, X, eventStatus, LAMBDA_OBS, *ret_h_nr, t1, 
                     Xminust1, effect_p, oneMinusEffect, Zt, log_Wt);
  }

  /* Zt must be 0 for trt = 0 */
  while(iter < maxiter) {
    iter++;
    lambda = updateLambda(n, Zt, trt, exst1, xt1, xt12, X, r2, &h_nr);
    updateZt(nTrtEq1, trtEq1Rows, X, eventStatus, lambda, h_nr, t1, 
                     Xminust1, effect_p, oneMinusEffect, Zt, log_Wt);

#if EM_STOP_CRITERIA == 0
    loglik = loglik_p_new(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
             trt1Log1mEffect, lambda, h_nr, trt, log_Wt);
    if (!R_FINITE(loglik)) error("ERROR: non-finite value for log-likelihood");
#endif

    if (iter > 1) {
#if EM_STOP_CRITERIA == 0
      diff = fabs(loglik - loglik0);
#else
      diff = fabs(lambda - lambda0);
#endif
      if (diff <= stopTol) {  
        conv = 1;
#if EM_STOP_CRITERIA == 0
        *ret_loglike = loglik_p_new(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
             trt1Log1mEffect, lambda, h_nr, trt, log_Wt);
#endif
        break;
      }

      if (lambdaFlag) {
        /* Test for stopping early for re-randomizations tests */
        if (lambda > lambda0) {
          lambdaInc     = 1;
          lambdaDec     = 0;
          lambdaDecIter = 0; 
          if (lambda0 > LAMBDA_OBS) {
            lambdaIncIter++;
          } else {
            lambdaIncIter = 0;
          }
        } else {
          lambdaInc     = 0;
          lambdaDec     = 1;
          lambdaIncIter = 0;
          if (lambda0 < LAMBDA_OBS) {
            lambdaDecIter++;
          } else {
            lambdaDecIter = 0;
          }
        }
        if ((lambdaIncIter >= LAMBDA_MONOITER) || (lambdaDecIter >= LAMBDA_MONOITER)) {
          conv = 1;
          break;
        } 
      }
      if (conv) break;
    }

    if (print > 1) {
      if (iter > 1) {
#if EMSTOP_CRITERIA == 0
        Rprintf("Iter=%d, lambda=%g, loglike=%11.4f, diff=%g\n", iter, lambda, loglik, diff);
#else
        Rprintf("Iter=%d, lambda=%g, diff=%g\n", iter, lambda, diff);
#endif
      } else {
#if EMSTOP_CRITERIA == 0
        Rprintf("Iter=%d, lambda=%g, loglike=%11.4f\n", iter, lambda, loglik);
#else
        Rprintf("Iter=%d, lambda=%g\n", iter, lambda);
#endif
      }     
    } 

    loglik0 = loglik;
    lambda0 = lambda;

  } /* END: while */

  if (print) {
    if (conv) {
      Rprintf("EM algoritm converged in %d iteration(s)\n", iter);
    } else {
      Rprintf("EM algoritm did not converge\n");
    }
  }

  *ret_lambda = lambda;
  *ret_h_nr   = h_nr;
  *ret_iter   = iter;

  free(log_Wt);
  return(conv);

} /* END: EM_loopP_new */


/* Zt is input and output */
static int EM_mainP(n, X, trt, eventStatus, effect_p, t1, stopTol, maxiter, print, method,
                    ret_lambda, ret_h_nr, Zt, ret_loglike)
int n, print, method, maxiter, *trt, *eventStatus;
double *X, *Zt, effect_p, t1, stopTol, *ret_lambda, *ret_h_nr, *ret_loglike;
{
  double *Xminust1, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, logEff, log1mEff;
  int nevents=0, conv=0, *exst1, nTrtEq1, *trtEq1Rows, retIter;

  /* Number of treated subjects */
  nTrtEq1         = sum_iVec(trt, n); 

  Xminust1        = dVec_alloc(n, 0, 0.0);
  exst1           = iVec_alloc(n, 0, 0.0);
  xt1             = dVec_alloc(n, 0, t1);
  xt12            = dVec_alloc(n, 0, 0.0);
  trt1LogEffect   = dVec_alloc(n, 0, 0.0);
  trt1Log1mEffect = dVec_alloc(n, 0, 0.0);
  trtEq1Rows      = iVec_alloc(nTrtEq1, 0, 0);
  logEff          = log(effect_p);
  log1mEff        = log(1.0 - effect_p);
  
  EM_setupObjP(n, X, eventStatus, t1, xt1, xt12, exst1, Xminust1, &nevents);

  EM_setupTrtP_new(n, trt, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect,
                   trtEq1Rows);
  conv = EM_loopP_new(n, nevents, maxiter, eventStatus, exst1, X, xt1, xt12, 
                    effect_p, trt1LogEffect, trt1Log1mEffect, stopTol, t1, print, 0,
                    trtEq1Rows, nTrtEq1, trt,
                    Xminust1, ret_h_nr, ret_lambda, Zt, DBL_MISS, &retIter, ret_loglike);


  free(Xminust1);
  free(exst1);
  free(xt1);
  free(xt12);
  free(trt1LogEffect);
  free(trt1Log1mEffect);
  free(trtEq1Rows);

  return(conv);

} /* END: EM_mainP */

void wcc_C_EM_mainP(pn, X, trt, eventStatus, peffect_p, pt1, pstopTol, pmaxiter, pprint, 
                    ret_conv, ret_lambda, ret_h_nr, Zt, ret_loglike)
int *pn, *pprint, *pmaxiter, *trt, *eventStatus, *ret_conv;
double *X, *Zt, *peffect_p, *pt1, *pstopTol, *ret_lambda, *ret_h_nr, *ret_loglike;
{

  *ret_conv = EM_mainP(*pn, X, trt, eventStatus, *peffect_p, *pt1, *pstopTol, 
                *pmaxiter, *pprint, 0, ret_lambda, ret_h_nr, Zt, ret_loglike);

  return;

} /* END: wcc_C_EM_mainP */

