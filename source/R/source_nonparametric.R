
Delta.fn <- function(t, t.fail.o, N, n.failure)
{
   tmp <- matrix(rep(t.fail.o, each=N), nrow=n.failure, ncol=N, byrow=TRUE)
   ret <- t(t(tmp) == t)

   ret
}

M.fn <- function(t, t.all, N, n.failure)
{
  rc  <- -1
  M   <- as.numeric(rep(0, n.failure*N))
  tmp <- .C("C_M_fn", as.numeric(t), as.numeric(t.all), as.integer(N), 
           as.integer(n.failure), retCode=as.integer(rc), 
           retVec=M, PACKAGE="Immunotherapy.Design")
  M   <- matrix(tmp$retVec, nrow=n.failure, ncol=N, byrow=TRUE)

  M
}

M1.fn <- function (lambda, t, t.all, N, n.failure, Mfn=NULL)
{ 
  if (is.null(Mfn)) Mfn = M.fn(t, t.all, N, n.failure)
  as.vector(t(lambda)%*%Mfn)
}

Mstar.fn <- function(tstar, t.all, n.failure)
{
  tvec   <- t.all[1:n.failure]
  tmp1   <- tvec <= tstar
  tmp2   <- pmin(rep(tstar, n.failure), t.all[-1]) - tvec
  M      <- tmp1*tmp2
  dim(M) <- c(n.failure, 1)
  
  M
}

Mstar1.fn <- function(lambda, tstar, t.all, n.failure, Mstarfn=NULL)
{
  if (is.null(Mstarfn)) Mstarfn = Mstar.fn(tstar, t.all, n.failure)
  t(lambda)%*%Mstarfn
}

M2.fn <- function(lambda, t, tstar, t.all, N, n.failure, M1fn=NULL, Mstar1fn=NULL)
{
  if (is.null(M1fn)) M1fn = M1.fn(lambda, t, t.all, N, n.failure)
  if (is.null(Mstar1fn)) Mstar1fn = Mstar1.fn(lambda, tstar, t.all, n.failure)
  ret <- M1fn - rep(Mstar1fn, N)

  ret
}

A.fn <- function(Zt, tstar, t.all, t, N, n.failure, Mstarfn=NULL, Mfn=NULL)
{
  if (is.null(Mstarfn)) Mstarfn = Mstar.fn(tstar, t.all, n.failure)
  if (is.null(Mfn)) Mfn = M.fn(t, t.all, N, n.failure)

  tmp0     <- Zt*(tstar <= t)
  tmp      <- tmp0
  dim(tmp) <- c(N, 1)
  mat      <- Mfn%*%(1-tmp)
  tmp2     <- tmp0*matrix(rep(Mstarfn, each=N), nrow=N, ncol=n.failure, byrow=FALSE)
  M        <- colSums(tmp2) + as.vector(mat)

  M
}

B.fn <- function(beta, t, t.all, tstar, Zt, N, n.failure, Mstarfn=NULL, Mfn=NULL)
{

  if (is.null(Mstarfn)) Mstarfn = Mstar.fn(tstar, t.all, n.failure)
  if (is.null(Mfn)) Mfn = M.fn(t, t.all, N, n.failure)
 
  tmp      <- Zt*(tstar <= t)
  dim(tmp) <- c(N, 1)
  mat      <- Mfn - matrix(rep(Mstarfn, each=N), nrow=n.failure, ncol=N, byrow=TRUE)
  M        <- exp(beta)*as.vector((mat %*% tmp))

  M
}

C.fn <- function(Zt, t, tstar, event_status)
{
  sum(Zt*(t>=tstar)*event_status)
}

D.fn <- function(Zt, t, tstar, lambda, t.all, N, n.failure, M2fn=NULL)
{
  if (is.null(M2fn)) M2fn = M2.fn(lambda, t, tstar, t.all, N, n.failure)
  ret <- sum(Zt*(t>=tstar)*M2fn)

  ret
}

loglik.fn = function(X, trt, event_status,
                     tstar,
                     t.all,
                     beta,
                     lambda, 
                     Zt,
                     effect_p, 
                     N, n.failure)
{

  Mstar1fn = as.vector(Mstar1.fn(lambda, tstar, t.all, n.failure, Mstarfn=NULL))
  M1fn     = M1.fn(lambda, X, t.all, N, n.failure, Mfn=NULL)
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
  loglik
}

pdf.r.fn <- function(lambda, t, Deltafn, event_status, beta, tstar, N, n.failure, t.all)
{

  vec      <- 1:N
  Mstar1fn = as.vector(Mstar1.fn(lambda, tstar, t.all, n.failure, Mstarfn=NULL))
  M2fn     = M2.fn(lambda,t,tstar,t.all,N, n.failure, M1fn=NULL, Mstar1fn=NULL)[vec]

  r1    <- apply(lambda^Deltafn, 2, prod)[vec]

  tmp1  <- (t >= tstar)[vec]
  r2    <- exp(beta*event_status[vec]*tmp1)
  m1fn  <- M1.fn(lambda,t,t.all,N, n.failure, Mfn=NULL)
  r3    <- exp(-(1-tmp1)*m1fn)
  r4    <- exp((-Mstar1fn-exp(beta)*M2fn)*tmp1)
  pdf.r <- (r1*r2*r3*r4)

  pdf.r
}

pdf.nr.fn <- function(lambda, t, Deltafn, N, n.failure, t.all)
{

  r1     <- apply(lambda^Deltafn, 2, prod)[1:N]
  r5     <- exp(-M1.fn(lambda,t,t.all,N, n.failure, Mfn=NULL))
  pdf.nr <- r1*r5

  pdf.nr

}

get_tfailo <- function(X, event_status) {

  xvec          <- X[event_status == 1]
  t.failure     <- unique(xvec)
  tmp           <- order(t.failure)
  t.fail.o      <- t.failure[tmp]
  delta.l       <- as.numeric(table(xvec))

  list(t.fail.o=t.fail.o, delta.l=delta.l)

} # END: get_tfailo

#####################################################
# !!! Data must be ordered before calling EM.main !!!
#####################################################
EM.main.NP <- function(X, trt, event_status, Zt, effect_p, t1,
               lambda, stopTol=1e-4, maxiter=10000, print=FALSE) {

  # X is time (t)

  n             <- length(trt)
  tmp           <- trt == 1
  ntrt          <- sum(tmp)
  if (is.null(Zt)) Zt <- ifelse(tmp, 0.5, 0)
  tmp0          <- get_tfailo(X, event_status)
  tmp           <- tmp0$t.fail.o
  delta.l       <- tmp0$delta.l
  nfail         <- length(tmp)
  tall          <- c(0, tmp)
  if (is.null(lambda)) lambda <- getHazard(X, trt, event_status)
  if (nfail != length(lambda)) stop("ERROR: lambda has the wrong length")
  tGEtstar      <- X >= t1
  tallLEtstar   <- tall <= t1
  tmp0          <- Delta.fn(X, tall[-1], n, nfail)
  csDeltafn     <- colSums(tmp0)
  tmp           <- rep(1:nfail, times=n)
  tmp           <- matrix(tmp, nrow=nfail, ncol=n, byrow=FALSE)
  posVec        <- tmp[tmp0] - 1
  nposVec       <- length(posVec)
  maxnp         <- max(csDeltafn)
  conv          <- as.integer(-1)
  beta          <- as.numeric(-9999)
  ret_lambda    <- as.numeric(rep(-9999, nfail))
  pdfr          <- as.numeric(rep(-9999, ntrt))
  pdfnr         <- as.numeric(rep(-9999, ntrt))
  ret_Zt        <- as.numeric(rep(-9999, n))
  tmp0          <- NULL
  ret_ll        <- as.numeric(-1)
  ret_ll_marg   <- as.numeric(-1)


  tmp <- .C("C_EM_mainNP", as.integer(n), as.numeric(Zt), as.integer(event_status), 
             as.integer(tGEtstar), as.integer(nfail), as.numeric(tall), as.numeric(t1),
             as.integer(tallLEtstar), as.numeric(X), as.integer(trt), as.integer(ntrt), 
             as.integer(csDeltafn), as.integer(posVec), as.integer(nposVec), as.integer(maxnp), 
             as.numeric(effect_p), as.integer(maxiter), as.numeric(stopTol), as.integer(print), 
             as.numeric(lambda), as.numeric(delta.l), 
             ret_conv=conv, ret_beta=beta, ret_lambda=ret_lambda, 
             pdfr=pdfr, pdfnr=pdfnr, ret_Zt=ret_Zt, 
             ret_ll=ret_ll, ret_ll_marg=ret_ll_marg, PACKAGE="Immunotherapy.Design")

  conv          <- tmp$ret_conv == 1
  beta          <- tmp$ret_beta
  mat           <- cbind(tall[-1], tmp$ret_lambda)
  colnames(mat) <- c("EventTime", "Lambda")
  pdfr          <- tmp$pdfr
  pdfnr         <- tmp$pdfn
  pr            <- effect_p*pdfr/(effect_p*pdfr + (1-effect_p)*pdfnr)
  tmp2          <- trt == 0
  pr[tmp2]      <- (tmp$ret_Zt)[tmp2]
  ret_ll        <- tmp$ret_ll
  ret_ll_marg   <- tmp$ret_ll_marg
  if (!conv) {
    ret_ll      <- NA
    ret_ll_marg <- NA
  }  

  result        <- list(converged=conv, logHR=beta, baseline=mat, probResponder=pr,
                        lambda=lambda, csDeltafn=csDeltafn, posVec=posVec, 
                        t.all=tall, tGEtstar=tGEtstar, tallLEtstar=tallLEtstar,
                        ret_lambda=tmp$ret_lambda, ret_Zt=tmp$ret_Zt, delta.l=delta.l,
                        loglike=ret_ll, loglike.marg=ret_ll_marg)

  result

} # END EM.main.NP


getHazard <- function(time, treatment, event_status) {

  control         <- coxph.control()
  control$timefix <- FALSE
  t.fail.o        <- get_tfailo(time, event_status)$t.fail.o
  n.failure       <- length(t.fail.o)
  t.all           <- c(0, t.fail.o)
  t.diff          <- t.all[-1]-t.all[-length(t.all)]
  fit             <- coxph(Surv(time, event_status) ~ treatment, control=control)
  ss              <- survfit(fit)
  cum.hazard      <- unique(-log(ss$surv))

  if (length(cum.hazard) != n.failure) {
    tmp <- cum.hazard != 0
    cum.hazard <- cum.hazard[tmp]
    if (length(cum.hazard) != n.failure) stop("ERROR computing hazard")
  }

  hazard.dis   = c(cum.hazard[1], cum.hazard[-1] - cum.hazard[-length(cum.hazard)])
  hazard       = hazard.dis / t.diff

  hazard

} # END: getHazard

PRIME.EM <- function(data, time.var="X", trt.var="trt", status.var="event_status",
               effect_p=0.6, t1=1, lambda0=NULL, probResponder=NULL, 
               stopTol=1e-4, maxiter=100000, print=0) {

  tmp <- setUpAndCheckData(data, time.var, trt.var, status.var, probResponder,
                           which="NP", lambda0=lambda0)
  checkOptions(effect_p, t1, stopTol, maxiter)

  ret <- EM.main.NP(tmp$time, tmp$treatment, tmp$event_status, tmp$probResponder, effect_p, 
                    t1, tmp$lambda, stopTol=stopTol, maxiter=maxiter, print=print)

  ret[c("converged", "logHR", "baseline", "probResponder", "loglike", "loglike.marg")]

} # END: Pembedded.EM.NP


PRIME.ReRandomizationTest <- function(data, time.var="X", trt.var="trt", 
              status.var="event_status", effect_p=0.6, t1=1, stopTol=1e-4, 
              maxiter=100000, print=0, num_rand=1000,
              min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  # Check for errors and order data
  tmp <- setUpAndCheckData(data, time.var, trt.var, status.var, NULL,
                           which="NP", lambda0=NULL)
  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=num_rand)

  time   <- tmp$time
  trt    <- tmp$treatment
  event  <- tmp$event_status
  lambda <- tmp$lambda
  Zt     <- tmp$probResponder
  checkTrtStatus(trt, event, min.sample.size=min.sample.size, min.n.event=min.n.event,  
                               min.per.trt=min.per.trt)

  # For observed data
  tmp <- EM.main.NP(time, trt, event, tmp$probResponder, effect_p, 
                    t1, lambda, stopTol=stopTol, maxiter=maxiter, print=print)
  logHR.obs <- tmp$logHR
  if (!is.finite(logHR.obs)) stop("ERROR estimating the observed log(HR), check parameter settings.")
  if (!tmp$converged) {
    stop("EM algorithm did not converge for the observed data. Try increasing maxiter and/or stopTol")
  }
  baseline  <- tmp$baseline
  if (print) cat(paste("Observed log(HR) = ", logHR.obs, "\n", sep=""))
  lambda    <- tmp$ret_lambda
  Zt        <- tmp$ret_Zt

  # shuffle the treatment labels in the observed data and obtain the 
  # re-randomization distribution of beta.hat on the shuffled data
  n         <- length(time)
  nfail     <- length(tmp$t.all) - 1
  ntrt      <- sum(trt == 1)
  np        <- length(tmp$posVec)
  maxNPos   <- max(tmp$csDeltafn)
  ret_nrand <- -1
  ret_p     <- -1

  tmp <- .C("C_ReRandNP", as.integer(num_rand), as.integer(n), as.numeric(time), 
            as.integer(trt), as.integer(event), as.numeric(effect_p), as.numeric(t1), 
               as.integer(nfail), as.numeric(tmp$t.all), as.numeric(lambda), 
               as.integer(tmp$tGEtstar), as.integer(tmp$tallLEtstar), as.integer(ntrt), 
               as.numeric(logHR.obs), as.integer(tmp$csDeltafn), as.integer(tmp$posVec),
               as.integer(np), as.integer(maxNPos), as.numeric(Zt),
               as.numeric(stopTol), as.integer(maxiter), as.integer(print), as.numeric(tmp$delta.l),
               ret_nrand=as.integer(ret_nrand), 
               ret_p=as.numeric(ret_p), PACKAGE="Immunotherapy.Design")
  m  <- tmp$ret_nrand
  p  <- tmp$ret_p
  ll <- tmp$ret_ll
  if (m < num_rand) warning("EM algorithm did not converge for all randomizations\n", sep="")
  if (print) {  
    if (m < num_rand) cat("NOTE: EM algorithm did not converge for all randomizations\n", sep="")   
    cat(paste("P-value is based on ", m, " randomizations\n", sep=""))
  }
  if (!m) stop("ERROR: p-value could not be estimated")
  
  result <- list(p.val.rerand=p, baseline=baseline, logHR=logHR.obs)

  return(result)

} # END: PRIME.ReRandomizationTest

Pow.PRIME.ReRandomizationTest <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=10000, print=0,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=num_rand, alpha=alpha, nsim=nsim)

  p.val.all <- rep(NA, nsim) 
  for (i in 1:nsim)
  {  
    if (print) cat(paste("Simulation ", i, "\n", sep=""))
    data.o <- generate_data(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, 
                 enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1)

    ret <- try(PRIME.ReRandomizationTest(data.o, effect_p=effect_p, t1=t1, 
               stopTol=stopTol, maxiter=maxiter, print=0, num_rand=num_rand,
               min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt),
               silent=FALSE)
    if (!("try-error" %in% class(ret))) p.val.all[i] <- ret$p.val.rerand
  }
  m     <- sum(is.na(p.val.all))
  n     <- nsim - m
  if (m) warning(paste("power based only on ", n, " simulated datasets", sep=""))
  if (!n) stop("ERROR: power could not be estimated")
  power <- mean(as.numeric(p.val.all <= alpha), na.rm=TRUE) 
  
  list(power=power, n=n)

} # END: Pow.PRIME.ReRandomizationTest


#####################################################
# !!! Data must be ordered before calling LRT.main !!!
#####################################################
LRT.main.NP <- function(X, trt, event_status, Zt, effect_p, t1,
               lambda, stopTol=1e-4, maxiter=10000, print=FALSE) {

  # X is time (t)

  n             <- length(trt)
  tmp           <- trt == 1
  ntrt          <- sum(tmp)
  if (is.null(Zt)) Zt <- ifelse(tmp, 0.5, 0)
  tmp0          <- get_tfailo(X, event_status)
  tmp           <- tmp0$t.fail.o
  delta.l       <- tmp0$delta.l
  nfail         <- length(tmp)
  tall          <- c(0, tmp)
  if (is.null(lambda)) lambda <- getHazard(X, trt, event_status)
  if (nfail != length(lambda)) stop("ERROR: lambda has the wrong length")
  tGEtstar      <- X >= t1
  tallLEtstar   <- tall <= t1
  tmp0          <- Delta.fn(X, tall[-1], n, nfail)
  csDeltafn     <- colSums(tmp0)
  tmp           <- rep(1:nfail, times=n)
  tmp           <- matrix(tmp, nrow=nfail, ncol=n, byrow=FALSE)
  posVec        <- tmp[tmp0] - 1
  nposVec       <- length(posVec)
  maxnp         <- max(csDeltafn)
  conv          <- as.integer(-1)
  beta          <- as.numeric(-9999)
  ret_lambda    <- as.numeric(rep(-9999, nfail))
  pdfr          <- as.numeric(rep(-9999, ntrt))
  pdfnr         <- as.numeric(rep(-9999, ntrt))
  ret_Zt        <- as.numeric(rep(-9999, n))
  tmp0          <- NULL
  ret_ll        <- as.numeric(-1)
  ret_ll0       <- as.numeric(-1) # Loglike at Zt = 0, beta = 0
  ret_p         <- as.numeric(-1)

  tmp <- .C("C_LRT_NP", as.integer(n), as.numeric(Zt), as.integer(event_status), 
             as.integer(tGEtstar), as.integer(nfail), as.numeric(tall), as.numeric(t1),
             as.integer(tallLEtstar), as.numeric(X), as.integer(trt), as.integer(ntrt), 
             as.integer(csDeltafn), as.integer(posVec), as.integer(nposVec), as.integer(maxnp), 
             as.numeric(effect_p), as.integer(maxiter), as.numeric(stopTol), as.integer(print), 
             as.numeric(lambda), as.numeric(delta.l), 
             ret_conv=conv, ret_beta=beta, ret_lambda=ret_lambda, 
             pdfr=pdfr, pdfnr=pdfnr, ret_Zt=ret_Zt, 
             ret_ll=ret_ll, ret_ll0=ret_ll0, ret_p=ret_p, PACKAGE="Immunotherapy.Design")

  conv          <- tmp$ret_conv == 1
  beta          <- tmp$ret_beta
  mat           <- cbind(tall[-1], tmp$ret_lambda)
  colnames(mat) <- c("EventTime", "Lambda")
  pdfr          <- tmp$pdfr
  pdfnr         <- tmp$pdfn
  pr            <- effect_p*pdfr/(effect_p*pdfr + (1-effect_p)*pdfnr)
  tmp2          <- trt == 0
  pr[tmp2]      <- (tmp$ret_Zt)[tmp2]
  ret_ll        <- tmp$ret_ll
  ret_ll0       <- tmp$ret_ll0
  ret_p         <- tmp$ret_p
  if (!conv) {
    ret_ll  <- NA
    ret_ll0 <- NA
    ret_p   <- NA
  }  

  result <- list(converged=conv, logHR=beta, baseline=mat, probResponder=pr,
                        lambda=lambda, csDeltafn=csDeltafn, posVec=posVec, 
                        t.all=tall, tGEtstar=tGEtstar, tallLEtstar=tallLEtstar,
                        ret_lambda=tmp$ret_lambda, ret_Zt=tmp$ret_Zt, delta.l=delta.l,
                        loglike.max=ret_ll, loglike.0=ret_ll0, p.value=ret_p)

  result

} # END LRT.main.NP

PRIME.LRT <- function(data, time.var="X", trt.var="trt", status.var="event_status",
               effect_p=0.6, t1=1, lambda0=NULL, probResponder=NULL, 
               stopTol=1e-4, maxiter=100000, print=0) {

  tmp <- setUpAndCheckData(data, time.var, trt.var, status.var, probResponder,
                           which="NP", lambda0=lambda0)
  checkOptions(effect_p, t1, stopTol, maxiter)

  ret <- LRT.main.NP(tmp$time, tmp$treatment, tmp$event_status, tmp$probResponder, effect_p, 
                    t1, tmp$lambda, stopTol=stopTol, maxiter=maxiter, print=print)

  ret <- ret[c("p.value", "loglike.max", "loglike.0")]

  ret

} # END: PRIME.LRT

Pow.Pembedded.LRT.NP_main <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, nsim=1000, print=0, powerL=-1, powerU=-1,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  p.val.all <- rep(NA, nsim) 
  nok       <- 0
  pflag     <- (powerL > 0) && (powerU > 0)
  oper      <- "="

  for (i in 1:nsim)
  {  

    if (print) cat(paste("Simulation ", i, "\n", sep=""))
    dat <- generate_data(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, 
                 enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1)

    ret <- try(PRIME.LRT(dat, time.var="X", trt.var="trt", status.var="event_status",
               effect_p=effect_p, t1=t1, lambda0=NULL, probResponder=NULL, 
               stopTol=stopTol, maxiter=maxiter, print=print), silent=TRUE)
    if (!("try-error" %in% class(ret))) {
      p.val.all[i] <- ret$p.value
      nok          <- nok + 1
    } 


    # compute min and max power
    if (pflag) {
      M        <- nsim - i
      tmp      <- sum(p.val.all <= alpha, na.rm=TRUE) 
      denom    <- M + nok
      minPower <- tmp/denom
      maxPower <- (tmp + M)/denom
      if (minPower > powerL) {
        if (M) oper <- ">"
        break
      } else if (maxPower < powerU) {
        if (M) oper <- "<"
        break
      } 
    }  
  }
  m     <- sum(is.na(p.val.all))
  n     <- nsim - m
  if (m && !pflag) warning(paste("power based only on ", n, " simulated datasets", sep=""))
  if (!n) stop("ERROR: power could not be estimated")
  if (!pflag) minPower <- mean(as.numeric(p.val.all <= alpha), na.rm=TRUE) 
  
  list(power=minPower, n=n, oper=oper)

} # END: Pow.Pembedded.LRT.NP_main

Pow.PRIME.LRT <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
                       alpha=0.05, nsim=10000, print=0,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=NULL, alpha=alpha, nsim=nsim)

  ret <- Pow.Pembedded.LRT.NP_main(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, enroll_rate=enroll_rate, 
                       lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter, stopTol=stopTol,
                       alpha=alpha, nsim=nsim, print=print,
                       min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)

  ret[c("power", "n")]

}

N.PRIME.LRT <- function(power=0.8, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
                       alpha=0.05, nsim=10000, min.N=100, max.N=700, 
                       tol.power=0.01, tol.N=1, print=1,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=NULL, alpha=alpha, nsim=nsim)
  if (min.N < 10) stop("ERROR: min.N must be at least 10")
  if (max.N < min.N) stop("ERROR: max.N cannot be smaller than min.N")
  if ((power < 0) || (power > 1)) stop("ERROR: power must be between 0 and 1")
  if (tol.N < 1) stop("ERROR: tol.N must be at least 1")
  PWRL <- power + tol.power + 1e-6
  PWRU <- power - tol.power - 1e-6

  # Test left endpoint
  tmp  <- Pow.Pembedded.LRT.NP_main(nmax=min.N, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, nsim=nsim, 
             powerL=PWRL, powerU=PWRU, 
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
  pwr1 <- tmp$power
  if (print) cat(paste("N = ", min.N, " power ", tmp$oper, " ", pwr1, "\n", sep=""))
  if (abs(power - pwr1) <= tol.power)  return(list(sampleSize=min.N, power=pwr1))
  if (pwr1 > power) {
    warning("The desired power is less than min.N]")
    return(list(sampleSize=min.N, power=pwr1)) 
  }

  # Test midpoint
  Nmid <- floor((min.N + max.N)/2)
  tmp  <- Pow.Pembedded.LRT.NP_main(nmax=Nmid, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter, 
             stopTol=stopTol, alpha=alpha, nsim=nsim, 
             powerL=PWRL, powerU=PWRU,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
  pwrm <- tmp$power
  if (print) cat(paste("N = ", Nmid, " power ", tmp$oper, " ",  pwrm, "\n", sep=""))
  if (abs(power - pwrm) <= tol.power)  return(list(sampleSize=Nmid, power=pwrm))

  step <- 100
  if (pwrm > power) {
    N0    <- min.N
    N1    <- Nmid
    diff1 <- abs(pwr1 - power)
    diff2 <- abs(pwrm - power)
    pwr2  <- pwrm
  } else {
    # Test right endpoint
    tmp  <- Pow.Pembedded.LRT.NP_main(nmax=max.N, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, nsim=nsim, 
             powerL=PWRL, powerU=PWRU,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
    pwr2 <- tmp$power
    if (print) cat(paste("N = ", max.N, " power ", tmp$oper, " ", pwr2, "\n", sep=""))
    if (abs(power - pwr2) <= tol.power)  return(list(sampleSize=max.N, power=pwr2))
    if (pwr2 < power) {
      warning("The desired power is greater than max.N")
      return(list(sampleSize=max.N, power=pwr2))
    }
    N0    <- Nmid
    N1    <- max.N
    diff1 <- abs(pwrm - power)
    diff2 <- abs(pwr2 - power)
    pwr1  <- pwrm
  }
 
  if ((pwr1 == 0) || (pwr2 == 1)) {
    N      <- floor((N0 + N1)/2)
  } else if (diff1 < diff2) {
    N      <- min(N0 + step, N1)
  } else {
    N      <- floor((N0 + N1)/2)
  }
  iter    <- 0
  minDiff <- 1e100
  minN    <- 0
  minPwr  <- -1
  while (1) {
    tmp  <- Pow.Pembedded.LRT.NP_main(nmax=N, rand_ratio=rand_ratio, effect_p=effect_p, 
              enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
              maxiter=maxiter, stopTol=stopTol, alpha=alpha, 
              nsim=nsim, powerL=PWRL, powerU=PWRU,
              min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
    pwr  <- tmp$power
    oper <- tmp$oper
    if (print) cat(paste("N = ", N, " power ", oper, " ",  pwr, "\n", sep=""))
    diff <- abs(power - pwr)
    if (diff < minDiff) {
      minDiff <- diff
      minN    <- N
      minPwr  <- pwr
    }
    if ((diff <= tol.power) || (abs(N1 - N0) <= tol.N)) break
    if (pwr < power) {
      N0 <- N
    } else {
      N1 <- N
    }
    N  <- floor((N0 + N1)/2) 
  }
  if (oper != "=") {
    tmp  <- Pow.Pembedded.LRT.NP_main(nmax=minN, rand_ratio=rand_ratio, effect_p=effect_p, 
              enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
              maxiter=maxiter, stopTol=stopTol, alpha=alpha,
              nsim=nsim, powerL=-1, powerU=-1,
              min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
    minPwr <- tmp$power
    if (print) cat(paste("N = ", N, " power ", tmp$oper, " ",  minPwr, "\n", sep=""))
  }
  
  list(sampleSize=minN, power=minPwr)

} # END: N.PRIME.LRT


