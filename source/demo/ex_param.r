library(Immunotherapy.Design)
data(data, package="Immunotherapy.Design")

ret <- Pembedded.EM.P(data, print = 2)
ret$lambda
ret$loglike

ret.wcc <- wcc.Pembedded.EM.P(data, print = 2)
ret.wcc$lambda
ret.wcc$loglike

