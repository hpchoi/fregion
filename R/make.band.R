#' @export
make.band.BEc <- function(eigen, conf.level, fd.eval.grid.size=200){
  alpha.level <- 1-conf.level
  pc.to.use <- sum(eigen$values > .Machine$double.eps)
  c.square <- sqrt(eigen$values[1:pc.to.use])
  weights <- eigen$values[1:pc.to.use]/c.square
  xi <- fregion::get.schisq.q.gamma(weights,conf.level) ## Approximate Quantile of Weighted Sum of Chi-square by Gamma
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    evalgrid <- fregion::make.grid(p=fd.eval.grid.size, rangevals=eigen$harmonics$basis$rangeval)
    eigen$vectors <- eval.fd(evalgrid,eigen$harmonics)
  }
  band.eval <- sqrt(apply(t(eigen$vectors[,1:pc.to.use]^2) * c.square * xi,2,sum))
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    return(Data2fd(evalgrid,band.eval,basisobj=eigen$harmonics$basis)) # return as fd object
  } else return(band.eval)                                               # return as vector
}
#' @export
make.band.BEPC <- function(eigen, conf.level, J, fd.eval.grid.size=200){ ## Finite dim. version of FPC based ellipse -> band.
  alpha.level <- 1-conf.level
  pc.to.use <- J
  c.square <- eigen$values[1:pc.to.use]
  xi <- qchisq(conf.level,J) ## Quantile of Chi-square
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    evalgrid <- fregion::make.grid(p=fd.eval.grid.size, rangevals=eigen$harmonics$basis$rangeval)
    eigen$vectors <- eval.fd(evalgrid,eigen$harmonics)
  }
  band.eval <- sqrt(apply(t(eigen$vectors[,1:pc.to.use]^2) * c.square * xi,2,sum))
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    return(Data2fd(evalgrid,band.eval,basisobj=eigen$harmonics$basis)) # return as fd object
  } else return(band.eval)                                               # return as vector
}
#' @export
make.band.Bs <- function(cov, conf.level, sim.size=10000, fd.eval.grid.size=200){
  if (inherits(cov,"bifd")) {
    evalgrid <- fregion::make.grid(p=fd.eval.grid.size, rangevals=cov$sbasis$rangeval)
    cov.m <- eval.bifd(evalgrid,evalbrid,cov) } else {cov.m <- cov}
  crit.Bs <- fregion::get.crit.supnorm.simple(cov.m,n.sim=sim.size,p=conf.level)
  band.eval <- sqrt(diag(cov.m)) * crit.Bs
  if (inherits(cov,"bifd")) {
    return(Data2fd(evalgrid,band.eval,basisobj=cov$sbasis))
  } else return(band.eval)
}
#' @export
make.band.naive.t <- function(cov, conf.level, df, fd.eval.grid.size=200){
  if (inherits(cov,"bifd")) {
    evalgrid <- fregion::make.grid(p=fd.eval.grid.size, rangevals=cov$sbasis$rangeval)
    cov.m <- eval.bifd(evalgrid,evalbrid,cov) } else {cov.m <- cov}
  band.eval <- fregion::qt2(conf.level,df) * sqrt(diag(cov.m))
  if (inherits(cov,"bifd")) {
    return(Data2fd(evalgrid,band.eval,basisobj=cov$sbasis))
  } else return(band.eval)
}
