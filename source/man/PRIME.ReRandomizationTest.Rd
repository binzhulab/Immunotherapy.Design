
\name{PRIME.ReRandomizationTest}
\alias{PRIME.ReRandomizationTest}
\title{Randomization test}
\description{
   Compute a randomization test p-value where test statistic is calculated based on PRIME strategy.                     
}
\usage{ PRIME.ReRandomizationTest(data, time.var="X", trt.var="trt", status.var="event_status",
           effect_p=0.6, t1=1, stopTol=1e-4, maxiter=100000, print=0, num_rand=1000,
           min.sample.size=50, min.n.event=5, min.per.trt=0.25)
}
\arguments{
  \item{data}{Data frame or matrix containing a time-to-event variable (\code{time.var}),
              a treatment variable (\code{trt.var}), 
              and a censoring variable (\code{status.var}).}
  \item{time.var}{Observational time variable name in \code{data}.}
  \item{trt.var}{Name of treatment assignment indicator in \code{data} coded as 
                 0 for control subjects and 1 for treated subjects.}
  \item{status.var}{Name of the binary censoring variable in \code{data} coded
                   as 0 for censored subjects and 1 for subjects that
                   experienced an event.}
  \item{effect_p}{Proportion of responders in the treatment arm at baseline.
                  The default is 0.6.}
  \item{t1}{Delayed duration. The default is 1 (month).}
  \item{stopTol}{Stopping tolerance in the EM algorithm. The default is 1e-4.}
  \item{maxiter}{Maximum number of iterations in the EM algorithm. The default is 100000.}
  \item{print}{0-2 to print information. Larger values will print more information.
               The default is 0.}
  \item{num_rand}{The number of replications in the re-randomization test. The default is 1000.}
  \item{min.sample.size}{Minimum sample size. The default is 50.}
  \item{min.n.event}{Minimum number of events. The default is 5.}
  \item{min.per.trt}{Minimum proportion of controls and treated subjects. The default is 0.25.}
}
\details{
  In each randomization, the treatment label is resampled and then the EM algorithm is called.
  The final p-value is based on all re-randomizations in which the EM algorithm converged.
}
\value{ A list containing the objects:
  \tabular{ll}{
    Name          \tab Description \cr
    p.val.rerand  \tab re-randomization test p-value \cr
    baseline      \tab matrix of event times and baseline hazards from observed data \cr
    logHR         \tab estimated log(hazard ratio) of responders versus controls
  }
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="Immunotherapy.Design")
  set.seed(1)
  # Will take much time to complete
  #ret <- PRIME.ReRandomizationTest(data)
}
\keyword{ EM, survival } 
\seealso{ \code{\link{Pembedded.ReRandomizationTest.P}} }
