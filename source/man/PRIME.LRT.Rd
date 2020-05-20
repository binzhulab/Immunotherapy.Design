
\name{PRIME.LRT}
\alias{PRIME.LRT}
\title{Likelihood ratio test}
\description{
   PRIME likelihood ratio test                        
}
\usage{ PRIME.LRT(data, time.var="X", trt.var="trt", status.var="event_status",
           effect_p=0.6, t1=1, lambda0=NULL, probResponder=NULL, 
           stopTol=1e-4, maxiter=100000, print=0)
}
\arguments{
  \item{data}{Data frame or matrix containing a time-to-event variable (\code{time.var}),
              a treatment variable (\code{trt.var}), 
              and a censoring variable (\code{status.var}).}
  \item{time.var}{Time-to-event variable name in \code{data}. The default is "X".}
  \item{trt.var}{Binary treatment variable name in \code{data} coded as 
                 0 for controls and 1 for subjects that received treatment.}
  \item{status.var}{Name of the binary censoring variable in \code{data} coded
                   as 0 for censored subjects and 1 for subjects that
                   experienced an event.}
  \item{effect_p}{Proportion of responders in the treatment arm at baseline. 
                  The default is 0.6.}
  \item{t1}{Delayed duration. The default is 1.}
  \item{lambda0}{NULL or vector of initial estimates for the baseline hazards
                 corrsponding to the ordered event times. 
                The default is NULL and will be computed from \code{\link{getHazard}}.}
  \item{probResponder}{NULL or vector of initial probabilities of a subject being a responder.
                       The default is NULL so that the initial probability is 0.5 for
                       treated subjects and 0 for controls.}
  \item{stopTol}{Stopping tolerance. The default is 1e-4.}
  \item{maxiter}{Maximum number of iterations. The default is 100000.}
  \item{print}{0-2 to print information. Larger values will print more information.
               The default is 0.}
}
\details{
 
}
\value{ A list containing the objects:
  \tabular{ll}{
    Name          \tab Description \cr
    p.value       \tab 1 df LRT p-value \cr 
    loglike.max   \tab log-likelihood value at the final estimates \cr
    loglike.0     \tab log-likelihood value under the null setting
                       where HR=1 and probability of being a \cr
                  \tab responder is \code{effect_p}
  }
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="Immunotherapy.Design")
  PRIME.LRT(data)
}
\keyword{ EM, survival } 
\seealso{ \code{\link{Pembedded.ReRandomizationTest.P}} }
