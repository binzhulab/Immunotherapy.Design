generate_data <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                          lambda1=0.117, HR=0.5, tau=12*5, t1=1) {
  
  # output data
  dat <- as.data.frame(matrix(0, nrow=nmax, ncol=9))
  names(dat) <- c("id", "trt", "Z", "tau", "enroll_time", "time_to_event", "event_status", "X", "t1")
  
  dat$id = 1:nmax
  
  # treatment allocation: 1 = experimental arm
  # sort by treatment assignment status to make the responders only arise among treated subjects 
  # set.seed(seed)
  dat$trt <- sort(rbinom(n = nmax, size = 1, prob=rand_ratio), decreasing = TRUE)
  n_trt = sum(dat$trt == 1)
  n_cnt = sum(dat$trt == 0)
  
  # simulate the response statue among the experimental arm with proportion effect_P: constraint: proportion of responders = effect_p
  # set.seed(seed)
  Z.trt = rbinom(n = n_trt, size = 1, prob=effect_p)
  Z.cnt = rep(0, n_cnt)
  dat$Z = c(Z.trt, Z.cnt)
  
  # total study duration 
  dat$tau = tau
  
  # enrollment follows expected trajectory from Poisson process with memoryless property
  ## waiting time between two consective subjects follows an Exponential distribution with enrollment rate
  ## arrival time for each subject is converted from the waitig by taking the cumulative sum of the waiting times
  waiting_times <- rexp(n = nmax, rate=enroll_rate)
  #dat$enroll_time <- cumsum(waiting_times)
  dat$enroll_time <- sample(ifelse(cumsum(waiting_times) <= dat$tau, cumsum(waiting_times), 0))

  tmp <- dat$enroll_time > 0
  if (!all(tmp)) {
    dat <- dat[tmp, , drop=FALSE]
    msg <- paste("Invalid enrollment times generated. Data will only have ", nrow(dat),
                 " observations.", sep="")
    warning(msg)
  }

  # t*
  dat$t1 <- t1
  
  # time to event
  # set.seed(seed)
  n_resp  <- sum(dat$Z == 1) # number of responders among the experimental and control groups
  n_nresp <- sum(dat$Z == 0) # number of non-responders among the experimental and control groups
  # dat$time_to_event[dat$Z == 1] <- rexp(n_resp, rate=lambda1*HR)

  dat$time_to_event[dat$Z == 1] <- rpexp(n_resp, rate=c(lambda1, HR*lambda1), t=c(0, dat$t1[1]))
  dat$time_to_event[dat$Z == 0] <- rexp(n_nresp, rate=lambda1)
  
  # build in the administrative censoring 
  dat$event_status <- ifelse(dat$time_to_event <= tau - dat$enroll_time, 1, 0)
  
  # observational time
  dat$X <- ifelse(dat$time_to_event <= tau - dat$enroll_time, dat$time_to_event, tau - dat$enroll_time)
  if (any(dat$X <= 0)) warning("Some observational times are <= 0. Try increasing tau or changing other parameters.")  

  row.names(dat) <- NULL
  
  dat
}

checkTrtStatus <- function(trt, status, min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  ret <- 0

  N <- length(trt)
  if (N < min.sample.size) {
    msg <- paste("WARNING: Sample size < ", min.sample.size, ". This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }
  M     <- ceiling(N*min.per.trt)
  ntrt0 <- sum(trt == 0)
  ntrt1 <- N - ntrt0
  if (!ntrt0) stop("ERROR: all subjects have treatment=1")
  if (!ntrt1) stop("ERROR: all subjects have treatment=0")
  if (ntrt0 < M) {
    msg <- paste("WARNING: Too few subjects with treatment=0. This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }  
  if (ntrt1 < M) {
    msg <- paste("WARNING: Too few subjects with treatment=1. This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }  
  Nev <- sum(status == 1)
  if (!Nev) stop("ERROR: all subjects have event_status=0")
  if (Nev < min.n.event) {
    msg <- paste("WARNING: Too few subjects with event_status=1. This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }

  ret

}

setUpAndCheckData <- function(data, time.var, trt.var, status.var, probResponder,
                              which="P", lambda0=NULL) {

  if ( (!is.data.frame(data)) && (!is.matrix(data)) ) stop("ERROR: data must be a data frame or matrix")
  # data must be ordered by treatment
  data         <- data[order(data[, trt.var], decreasing = TRUE), , drop=FALSE]
  treatment    <- data[, trt.var]
  event_status <- data[, status.var]
  time         <- data[, time.var]
  n            <- length(time)
  if (length(treatment) != n) stop("ERROR: length(treatment) != length(time)")
  if (length(event_status) != n) stop("ERROR: length(event_status) != length(time)")
  tmp <- (!is.finite(time)) | (time < 0)
  tmp[is.na(tmp)] <- TRUE
  if (any(tmp)) stop("ERROR: with time vector")
  if (!all(treatment %in% 0:1)) stop("ERROR: with treatment vector")
  if (!all(event_status %in% 0:1)) stop("ERROR: with event_status vector")
  if (is.null(probResponder)) probResponder <- ifelse(treatment == 1, 0.5, 0)
  if (length(probResponder) != n) stop("ERROR: length(probResponder) != length(time)")
  if ((which == "NP") && !length(lambda0)) lambda0 <- getHazard(time, treatment, event_status)

  list(time=time, treatment=treatment, event_status=event_status, 
       probResponser=probResponder, lambda0=lambda0)

} # END: setUpAndCheckData

checkOptions <- function(effect_p, t1, stopTol, maxiter, num_rand=NULL, alpha=NULL, nsim=NULL) {

  if ((effect_p <= 0) || (effect_p > 1)) stop("ERROR: effect_p must be > 0 and <= 1")
  if (length(t1) != 1) stop("ERROR: t1 must have length 1")
  if (t1 < 0) stop("ERROR: t1 must be non-negative")
  if (maxiter < 1) stop("ERROR: maxiter must be greater than 0")
  if (stopTol <= 0) stop("ERROR: stopTol must be positive")
  
  if (length(num_rand)) {
    if (num_rand < 1) stop("ERROR: num_rand must be greater than 0")
  }
  if (length(alpha)) {
    if ((alpha <= 0) || (alpha >= 1)) stop("ERROR: alpha must be > 0 and < 1")
  }
  if (length(nsim)) {
    if (nsim < 1) stop("ERROR: nsim must be greater than 0")
  }

  NULL

} # checkOptions
