#' fregion
#'
#' This package contains example codes for inferential tools suggested in \url{http://arxiv.org/abs/1607.07771}, titled 'Geometric Approach to Confidence Regions and Bands for Functional Parameters'.
#' The top level functions below take functional estimate and the covariance of it to perform hypothesis testings, to construct confidence bands, and to construct/visuallize hyper-rectangular regions.
#' The inputs can be either vector/matrix or \link{fd}/\link{bifd} object from \link{fda} package.
#' \itemize{
#'   \item \link{fregion.test}: Perform hypothesis testings based on confidence regions.
#'   \item \link{fregion.band}: Construct (simultaneous) confidence bands using hyper-ellipsoid confidence regions.
#'   \item \link{fregion.rect}: Construct a hyper-rectangular confidence region in Hilbert space. It can be visuzlized using \link{plot.fregion.rect}.
#' }
#' The example below is given in a mean function estimation context but works for other functional estimates as long as the estimators are Gaussian.
#'
#' @examples
#' # 1. Vector/matrix version
#'
#' # Generate a sample
#' p = 200 ; N = 80 ; rangeval = c(0,1)
#' grid = make.grid(p, rangevals=rangeval)
#' mu0 = meanf.poly(grid,c(0,1)) ; names(mu0) = grid
#' mu = meanf.poly(grid,c(0,1.1)) ; names(mu) = grid
#' cov.m = make.cov.m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' x = make.sample(mu,cov.m,N)
#'
#' # Find the estimate and covariance
#' hat.mu = rowMeans(x)
#' hat.cov.m = crossprod(t(x - hat.mu)) / (N-1)
#' e.hat.cov.m = eigen(hat.cov.m)   # <- This is optional and can be provide into the functions instead of hat.cov.m below.
#'
#' # Compare different methods for Hypothesis testings.
#' (a1 <- fregion.test(hat.mu,mu0,hat.cov.m,N,type=c("ALL"),pc.cut=c(1,3,4,5,0.99,0.999)))
#'
#' # Make and visualize/compare confidence bands
#' b <- fregion.band(hat.mu,hat.cov.m,N=N,type=c("Bs","BEc","BEPC.10","naive.t"),conf.level=c(0.95))
#' plot(b)
#'
#' # Make rectangular region and visulize
#' c <- fregion.rect(hat.mu-mu0,hat.cov.m,N=N)
#' plot(c)
#'
#'
#' # 2. fd/bifd version
#'
#' # create basis, convert vector/matrix into fd/bifd objects.
#' require(fda)
#' nbasis <- round(p*0.9)
#' fd.basis <- create.bspline.basis(rangeval,nbasis)
#' mu0.fd <- Data2fd(names(mu0),mu0,fd.basis)
#' mu.fd <- Data2fd(names(mu),mu,fd.basis)
#' x.fd <- Data2fd(rownames(x),x,fd.basis)
#' hat.mu.fd <- mean.fd(x.fd)
#' hat.cov.fd <- var.fd(x.fd)
#' e.hat.cov.fd <- eigen.fd(hat.cov.fd)   # <- This is optional and can be provide into the functions instead of hat.cov.fd below.
#'
#' # Compare different methods for Hypothesis testings.
#' (a1.fd <- fregion.test(hat.mu.fd, mu0.fd, hat.cov.fd, N, type=c("ALL"), pc.cut=c(1,3,4,5,0.99,0.999)))
#'
#' # Make and visualize/compare confidence bands
#' b.fd <- fregion.band(hat.mu.fd,hat.cov.fd,N=N,type=c("Bs","BEc","BEPC.10","naive.t"),conf.level=c(0.95))
#' plot(b.fd)
#'
#' # Make rectangular region and visulize
#' c.fd <- fregion.rect(hat.mu.fd-mu0.fd,hat.cov.fd,N=N)
#' plot(c.fd)
#' @docType package
#' @name fregion
NULL
