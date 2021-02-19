
wcc.Pembedded.EM.P <- function(data, time.var="X", trt.var="trt", status.var="event_status",
               effect_p=0.6, t1=1, probResponder=NULL, 
               stopTol=1e-5, maxiter=100000, print=0) {

  tmp <- setUpAndCheckData(data, time.var, trt.var, status.var, probResponder)
  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=NULL, alpha=NULL)

  ret <- wcc.EM.P.main(tmp$time, tmp$treatment, tmp$event_status, tmp$probResponder, effect_p, t1,  
                      stopTol=stopTol, maxiter=maxiter, print=print)

  ret

} # END: wcc.Pembedded.EM.P

################################################################
# EM Algorithm
# num_rep     = number of replications in the bias assessment
# nmax        = maximum sample size
# rand_ratio  = probability of assignment to experimental arm
# effect_P    = Proportion of responders in the treatment arm
# HR          = hazard ratio
# enroll_rate = enrollment rate in subjects per month
# tau         = total study duration
################################################################

wcc.EM.P.main <- function(X, trt, event_status, Zt, effect_p, t1, 
                      stopTol=0.00001, maxiter=1000, print=0) {

  N <- length(trt)
  if (is.null(Zt)) {
    Zt           <- rep(0, N)
    Zt[trt == 1] <- 0.5
  }

  ret_conv    <- as.integer(0)
  ret_lambda  <- as.numeric(-9999)
  ret_h_nr    <- as.numeric(-9999)
  ret_ll      <- as.numeric(-9999)

  tmp <- .C("wcc_C_EM_mainP", as.integer(N), as.numeric(X), as.integer(trt), 
             as.integer(event_status), as.numeric(effect_p), as.numeric(t1), 
             as.numeric(stopTol), as.integer(maxiter), as.integer(print), 
             ret_conv=ret_conv, ret_lambda=ret_lambda, ret_h_nr=ret_h_nr, 
             Zt=as.numeric(Zt), ret_ll=ret_ll, PACKAGE="Immunotherapy.Design")
  
  ret_conv <- as.logical(tmp$ret_conv)
  ret_ll   <- tmp$ret_ll
  if (!ret_conv) ret_ll <- NA
  result <- list(converged=ret_conv, lambda=tmp$ret_lambda, 
                 baseline=tmp$ret_h_nr, probResponder=tmp$Zt, loglike=ret_ll)

  return(result)

} # END: wcc.EM.P.main
