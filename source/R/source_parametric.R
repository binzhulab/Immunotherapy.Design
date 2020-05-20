
Pembedded.EM.P <- function(data, time.var="X", trt.var="trt", status.var="event_status",
               effect_p=0.6, t1=1, probResponder=NULL, 
               stopTol=1e-5, maxiter=100000, print=0) {

  tmp <- setUpAndCheckData(data, time.var, trt.var, status.var, probResponder)
  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=NULL, alpha=NULL)

  ret <- EM.P.main(tmp$time, tmp$treatment, tmp$event_status, tmp$probResponder, effect_p, t1,  
                      stopTol=stopTol, maxiter=maxiter, print=print)

  ret

} # END: Pembedded.EM.P

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

EM.P.main <- function(X, trt, event_status, Zt, effect_p, t1, 
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

  tmp <- .C("C_EM_mainP", as.integer(N), as.numeric(X), as.integer(trt), 
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

} # END: EM.P.main

Pembedded.ReRandomizationTest.P_old <- function(data, time.var="X", trt.var="trt", 
               status.var="event_status", effect_p=0.6, t1=1, stopTol=1e-5, 
               maxiter=10000, print=0, num_rand=10000) {

  # EM estimate parameters of interest and Zt
  EM.est <- Pembedded.EM.P(data, time.var=time.var, trt.var=trt.var, status.var=status.var,
               effect_p=effect_p, t1=t1, probResponder=NULL, 
               stopTol=stopTol, maxiter=maxiter, print=print)
  lambda.obs <- EM.est$lambda

  if (print) cat(paste("Observed lambda = ", lambda.obs, "\n", sep=""))
   
  # shuffle the treatment labels in the observed data and obtain the 
  # re-randomization distribution of lambda.hat on the shuffled data
  nr        <- nrow(data)
  ret_p     <- as.numeric(-1)
  ret_nrand <- as.integer(0)
  prt       <- 0

  tmp <- .C("C_ReRandP", as.integer(num_rand), as.integer(nr), as.numeric(data[, time.var]), 
            as.integer(data[, trt.var]), as.integer(data[, status.var]), as.numeric(effect_p), 
            as.numeric(t1), as.numeric(stopTol), as.integer(maxiter),
            as.integer(prt), as.numeric(lambda.obs), 
            ret_nrand=ret_nrand, ret_p=ret_p, PACKAGE="Immunotherapy.Design")
  m   <- tmp$ret_nrand
  if (print) {  
    if (m < num_rand) cat("NOTE: EM algorithm did not converge for all randomizations\n", sep="")
    cat(paste("P-value is based on ", m, " randomizations\n", sep=""))
  }
  if (!m) stop("ERROR: p-value could not be estimated")
 
  result <- list(p.val.rerand=tmp$ret_p, baseline=EM.est$baseline, lambda=lambda.obs)

  return(result)


} # END: Pembedded.ReRandomizationTest.P_old

Pembedded.ReRandomizationTest.P_main <- function(time, trt, status, effect_p=0.6, t1=1, stopTol=1e-5, 
               maxiter=100000, print=0, num_rand=10000, alpha=-1, method=1,
               min.sample.size=50, min.n.event=5, min.per.trt=0.25, stopOnBadData=0) {
 
  # Check the data
  flag       <- checkTrtStatus(trt, status, min.sample.size=min.sample.size, min.n.event=min.n.event,  
                               min.per.trt=min.per.trt)
  if (flag && stopOnBadData) stop("Bad data generated")

  nr         <- length(time)
  ret_p      <- as.numeric(-1)
  ret_nrand  <- as.integer(0)
  ret_lambda <- as.numeric(-1)
  ret_hnr    <- as.numeric(-1)

  tmp <- .C("C_ReRandP_new", as.integer(num_rand), as.integer(nr), as.numeric(time), 
            as.integer(trt), as.integer(status), as.numeric(effect_p), 
            as.numeric(t1), as.numeric(stopTol), as.integer(maxiter),
            as.integer(print), as.numeric(alpha), as.integer(method),
            ret_nrand=ret_nrand, ret_p=ret_p,
            ret_lambda=ret_lambda, ret_hnr=ret_hnr, PACKAGE="Immunotherapy.Design")

  tmp

} # END: Pembedded.ReRandomizationTest.P_main

Pembedded.ReRandomizationTest.P <- function(data, time.var="X", trt.var="trt", 
               status.var="event_status", effect_p=0.6, t1=1, stopTol=1e-5, 
               maxiter=100000, print=0, num_rand=10000,
               min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  tmp <- setUpAndCheckData(data, time.var, trt.var, status.var, NULL)
  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=num_rand)

  tmp <- Pembedded.ReRandomizationTest.P_main(tmp$time, tmp$treatment, tmp$event_status, 
               effect_p=effect_p, t1=t1, stopTol=stopTol, 
               maxiter=maxiter, print=print, num_rand=num_rand, alpha=-1, method=0,
               min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)

  m   <- tmp$ret_nrand
  if (m < num_rand) warning("EM algorithm did not converge for all randomizations\n", sep="")
  if (print) {  
    if (m < num_rand) cat("NOTE: EM algorithm did not converge for all randomizations\n", sep="")
    cat(paste("P-value is based on ", m, " randomizations\n", sep=""))
  }
  if (!m) stop("ERROR: p-value could not be estimated")
 
  result <- list(p.val.rerand=tmp$ret_p, baseline=tmp$ret_hnr, lambda=tmp$ret_lambda)

  return(result)

} # END: Pembedded.ReRandomizationTest.P

Pow.Pembedded.P_old <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=1000, print=0) {

  p.val.all <- rep(NA, nsim) 
  for (i in 1:nsim)
  {  
    if (print) cat(paste("Simulation ", i, "\n", sep=""))
    data.o <- generate_data(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, 
                 enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1)

    ret <- try(Pembedded.ReRandomizationTest.P(data.o, effect_p=effect_p, t1=t1, 
              stopTol=stopTol, maxiter=maxiter, print=0, num_rand=num_rand), silent=FALSE)
    if (!("try-error" %in% class(ret))) p.val.all[i] <- ret$p.val.rerand
  }
  m     <- sum(is.na(p.val.all))
  n     <- nsim - m
  if (m) warning(paste("power based only on ", n, " simulated datasets", sep=""))
  if (!n) stop("ERROR: power could not be estimated")
  power <- mean(as.numeric(p.val.all <= alpha), na.rm=TRUE) 
  
  list(power=power, n=n)

} # END: Pow.Pembedded.P_old

Pow.Pembedded.P_main <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=1000, print=0, powerL=-1, powerU=-1,
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

    ret <- try(Pembedded.ReRandomizationTest.P_main(dat[, "X"], dat[, "trt"], dat[, "event_status"], 
             effect_p=effect_p, t1=t1, stopTol=stopTol, maxiter=maxiter, print=print, 
             num_rand=num_rand, alpha=alpha, method=1,
             min.sample.size=min.sample.size, min.n.event=min.n.event, 
             min.per.trt=min.per.trt), silent=FALSE)
    if (!("try-error" %in% class(ret))) {
      p.val.all[i] <- ret$ret_p
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

} # END: Pow.Pembedded.P_main


Pow.Pembedded.P <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=10000, print=0,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=num_rand, alpha=alpha, nsim=nsim)
  ret <- Pow.Pembedded.P_main(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, print=print,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
  
  ret

} # END: Pow.Pembedded.P


N.Pembedded.P <- function(power=0.8, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=10000, min.N=100, max.N=700, 
                       tol.power=0.01, tol.N=1, print=1,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  checkOptions(effect_p, t1, stopTol, maxiter, num_rand=num_rand, alpha=alpha, nsim=nsim)
  if (min.N < 10) stop("ERROR: min.N must be at least 10")
  if (max.N < min.N) stop("ERROR: max.N cannot be smaller than min.N")
  if ((power < 0) || (power > 1)) stop("ERROR: power must be between 0 and 1")
  if (tol.N < 1) stop("ERROR: tol.N must be at least 1")
  PWRL <- power + tol.power + 1e-6
  PWRU <- power - tol.power - 1e-6

  # Test left endpoint
  tmp  <- Pow.Pembedded.P_main(nmax=min.N, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, 
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
  tmp  <- Pow.Pembedded.P_main(nmax=Nmid, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter, 
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, 
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
    tmp  <- Pow.Pembedded.P_main(nmax=max.N, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, 
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
    tmp  <- Pow.Pembedded.P_main(nmax=N, rand_ratio=rand_ratio, effect_p=effect_p, 
              enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
              maxiter=maxiter, stopTol=stopTol, alpha=alpha, num_rand=num_rand, 
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
    tmp  <- Pow.Pembedded.P_main(nmax=minN, rand_ratio=rand_ratio, effect_p=effect_p, 
              enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
              maxiter=maxiter, stopTol=stopTol, alpha=alpha, num_rand=num_rand, 
              nsim=nsim, powerL=-1, powerU=-1,
              min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
    minPwr <- tmp$power
    if (print) cat(paste("N = ", N, " power ", tmp$oper, " ",  minPwr, "\n", sep=""))
  }
  
  list(sampleSize=minN, power=minPwr)

} # END: N.Pembedded.P


