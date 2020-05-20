#include "./source.h"

static double loglik_p(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
             trt1Log1mEffect, lambda, h_nr)
int n, nevents, *exst1;
double *X, *Zt, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, lambda, h_nr;
{
  int i;
  double ret, r1=0.0, r2, r3=0.0, r4=0.0, r5=0.0, r6=0.0, r7=0.0, Zti, oneMinusZti;

  r2 = (double) nevents;
  
  for (i=0; i<n; i++) {
    Zti         = Zt[i];
    oneMinusZti = 1.0 - Zti;
    r1 += Zti*exst1[i];
    r3 += Zti*xt1[i];
    r4 += Zti*xt12[i];
    r5 += oneMinusZti*X[i];
    r6 += Zti*trt1LogEffect[i];
    r7 += oneMinusZti*trt1Log1mEffect[i];
  }

  ret = r1*log(lambda) + r2*log(h_nr) - h_nr*(r3 + lambda*r4 + r5) + r6 + r7;

  return(ret);

} /* END: loglik_p */

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
             trt1Log1mEffect, lambda, h_nr, trt)
int n, nevents, *exst1, *trt;
double *X, *Zt, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, lambda, h_nr;
{
  int i;
  double ret, r1=0.0, r2, r3=0.0, r4=0.0, r5=0.0, r6=0.0, r7=0.0, Zti, oneMinusZti;

  r2 = (double) nevents;

  for (i=0; i<n; i++) { 
    Zti         = Zt[i];
    oneMinusZti = 1.0 - Zti;
    if (trt[i]) {
      r1 += Zti*exst1[i];
      r3 += Zti*xt1[i];
      r4 += Zti*xt12[i];
      r6 += Zti*trt1LogEffect[i];
      r7 += oneMinusZti*trt1Log1mEffect[i];
    }
    r5 += oneMinusZti*X[i];
    
  }

  ret = r1*log(lambda) + r2*log(h_nr) - h_nr*(r3 + lambda*r4 + r5) + r6 + r7;

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
                     Xminust1, effect_p, oneMinusEffect, Zt)
int nTrtEq1, *trtEq1Rows, *eventStatus;
double *X, lambda, h_nr, t1, *Xminust1, effect_p, oneMinusEffect, *Zt;
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

  r2             = (double) nevents;
  oneMinusEffect = 1.0 - effect_p;
  lambdaFlag     = (LAMBDA_OBS  > DBL_NOTMISS_GT) && (method == 1); 

  /* Get good initial estimates for permutation tests */
  if (LAMBDA_OBS  > DBL_NOTMISS_GT) {  
    updateZt(nTrtEq1, trtEq1Rows, X, eventStatus, LAMBDA_OBS, *ret_h_nr, t1, 
                     Xminust1, effect_p, oneMinusEffect, Zt);
  }

  /* Zt must be 0 for trt = 0 */
  while(iter < maxiter) {
    iter++;
    lambda = updateLambda(n, Zt, trt, exst1, xt1, xt12, X, r2, &h_nr);
    updateZt(nTrtEq1, trtEq1Rows, X, eventStatus, lambda, h_nr, t1, 
                     Xminust1, effect_p, oneMinusEffect, Zt);

#if EM_STOP_CRITERIA == 0
    loglik = loglik_p_new(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
             trt1Log1mEffect, lambda, h_nr, trt);
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
#if EM_STOP_CRITERIA == 1
        *ret_loglike = loglik_p_new(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
             trt1Log1mEffect, lambda, h_nr, trt);
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
        Rprintf("Iter=%d, loglike=%11.4f, diff=%g\n", iter, loglik, diff);
#else
        Rprintf("Iter=%d, lambda=%g, diff=%g\n", iter, lambda, diff);
#endif
      } else {
#if EMSTOP_CRITERIA == 0
        Rprintf("Iter=%d, loglike=%11.4f\n", iter, loglik);
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

void C_EM_mainP(pn, X, trt, eventStatus, peffect_p, pt1, pstopTol, pmaxiter, pprint, 
                    ret_conv, ret_lambda, ret_h_nr, Zt, ret_loglike)
int *pn, *pprint, *pmaxiter, *trt, *eventStatus, *ret_conv;
double *X, *Zt, *peffect_p, *pt1, *pstopTol, *ret_lambda, *ret_h_nr, *ret_loglike;
{

  *ret_conv = EM_mainP(*pn, X, trt, eventStatus, *peffect_p, *pt1, *pstopTol, 
                *pmaxiter, *pprint, 0, ret_lambda, ret_h_nr, Zt, ret_loglike);

  return;

} /* END: C_EM_mainP */

static int ReRandP_new(num_rand, n, X, trt, eventStatus, effect_p, t1, stopTol, 
                       maxiter, print, alpha, method, ret_p, ret_lambda, ret_hnr)
int num_rand, n, *trt, *eventStatus, maxiter, print, method;
double *X, effect_p, t1, stopTol, *ret_p, alpha, *ret_lambda, *ret_hnr;
{
  double *Xminust1, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, logEff, log1mEff;
  double *Zt, lambda, h_nr, minP, maxP, M, lambda_obs, h_nr_obs, alphaBndL, alphaBndU;
  double LAMBDA_OBS, ret_loglike;
  int iter, conv, nevents=0, *exst1, *trtPermute, sumNrand=0, sumGTobs=0;
  int nTrtEq1, *trtEq1Rows, alphaFlag, retIter, totalIter=0;

  *ret_p   = -1.0;
  logEff   = log(effect_p);
  log1mEff = log(1.0 - effect_p);

  /* Number of treated subjects */
  nTrtEq1   = sum_iVec(trt, n);
  alphaFlag = (alpha > 0.0);
  alphaBndL = alpha + 1.0/num_rand;
  alphaBndU = alpha - 1.0/num_rand;

  Xminust1        = dVec_alloc(n, 0, 0.0);
  exst1           = iVec_alloc(n, 0, 0.0);
  xt1             = dVec_alloc(n, 0, t1);
  xt12            = dVec_alloc(n, 0, 0.0);
  trt1LogEffect   = dVec_alloc(n, 0, 0.0);
  trt1Log1mEffect = dVec_alloc(n, 0, 0.0);
  trtPermute      = iVec_alloc(n, 0, 0.0);
  Zt              = dVec_alloc(n, 1, 0.0); 
  trtEq1Rows      = iVec_alloc(nTrtEq1, 0, 0);

  /* Set up objects that do not depend on permuted trt */
  EM_setupObjP(n, X, eventStatus, t1, xt1, xt12, exst1, Xminust1, &nevents);

  /* Set up objects that depend on trt */
  EM_setupTrtP_new(n, trt, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect, trtEq1Rows);

  /* Get observed lambda.
     Initial Zt is used, h_nr and lambda get updated immediately.
  */
  setupZt(n, trt, Zt);
  conv = EM_loopP_new(n, nevents, maxiter, eventStatus, exst1, X, xt1, xt12, 
                    effect_p, trt1LogEffect, trt1Log1mEffect, stopTol, t1, print, 0, 
                    trtEq1Rows, nTrtEq1, trt,
                    Xminust1, &h_nr_obs, &lambda_obs, Zt, DBL_MISS, &retIter, &ret_loglike);
  if (print) Rprintf("Observed lambda = %g\n", lambda_obs);
  if (!conv) error("ERROR: EM algorithm did not converge for the observed data. Try increasing maxiter.");
  *ret_lambda = lambda_obs;
  *ret_hnr    = h_nr_obs;
  LAMBDA_OBS  = lambda_obs;

  print = 0;
  for (iter=0; iter<num_rand; iter++) {
    permute_iVec(trt, n, trtPermute);
    setupZt(n, trtPermute, Zt);

    /* Set up objects that depend on trt */
    EM_setupTrtP_new(n, trtPermute, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect, trtEq1Rows);

    /* Set for each perm */
    h_nr = h_nr_obs;
    conv =  EM_loopP_new(n, nevents, maxiter, eventStatus, exst1, X, xt1, xt12, 
                    effect_p, trt1LogEffect, trt1Log1mEffect, stopTol, t1, print, method,
                    trtEq1Rows, nTrtEq1, trtPermute,
                    Xminust1, &h_nr, &lambda, Zt, LAMBDA_OBS, &retIter, &ret_loglike);    

    if (conv) {

      /*totalIter += retIter;*/
      sumNrand++;
      if (lambda > lambda_obs) sumGTobs++;

      /* Check to stop early. If alpha > 0, then we really only need to consider
         pvalues <= alpha for power calculations. In this case, compute the 
         smallest possible p-value based on the current number of randomizations. */
      if (alphaFlag) {
        M    = (double) (num_rand - iter - 1);
        minP = 1.0 - ((M + sumGTobs)/(M + sumNrand));
        maxP = 1.0 - (sumGTobs/(M + sumNrand));
        if ((minP > alphaBndL) || (maxP < alphaBndU)) break;
      }
    }
  }

  if (sumNrand) *ret_p = 1.0 - ((double) sumGTobs)/((double) sumNrand);

  free(Xminust1);
  free(exst1);
  free(xt1);
  free(xt12);
  free(trt1LogEffect);
  free(trt1Log1mEffect);
  free(trtPermute);
  free(Zt);
  free(trtEq1Rows);

  return(sumNrand);

} /* END: ReRandP_new */


void C_ReRandP_new(pnum_rand, pn, X, trt, eventStatus, peffect_p, pt1, pstopTol, pmaxiter,
               pprint, palpha, pmethod, ret_nrand, ret_p, ret_lambda, ret_hnr)
int *pnum_rand, *pn, *trt, *eventStatus, *pmaxiter, *pprint, *ret_nrand, *pmethod;
double *X, *peffect_p, *pt1, *pstopTol, *ret_p, *palpha, *ret_lambda, *ret_hnr;
{

  /* For random number generation */
  GetRNGstate();

  *ret_nrand = ReRandP_new(*pnum_rand, *pn, X, trt, eventStatus, *peffect_p, *pt1, 
                       *pstopTol, *pmaxiter, *pprint, *palpha, *pmethod,
                       ret_p, ret_lambda, ret_hnr);

  PutRNGstate();  

  return;

} /* END: C_ReRandP_new */

