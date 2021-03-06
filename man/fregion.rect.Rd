% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fregion.rect.R
\name{fregion.rect}
\alias{fregion.rect}
\title{Builds a Rectangular Confidence Region}
\usage{
fregion.rect(x, cov, N = 1, type = "Rz", conf.level = 0.95,
  pc.cut = 0.999, df = NULL)
}
\arguments{
\item{x}{Functional parameter estimate. It can be either a vector or \link{fd} object from \link{fda}.}

\item{cov}{N * Cov(X), in which X is the functional estimator. It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.}

\item{N}{It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.}

\item{type}{This specifies which rectangular region to be constructed. It should be either one of "Rz", "Rz1", "Rzs", or "Rz1s".}

\item{conf.level}{Confidence level of the region.}

\item{pc.cut}{A numeric value. For integer values, fPC up to those values will be used. If it's a value from 0 to 1, this specifies the proportion of (estimated) variance that should be explained by the fPCs. If it is 0, all the available fPCs will be used as long as the size of eigenvalues are greater than .Machine$double.eps.}

\item{df}{Degree of freedom to use in small sample versions.}
}
\value{
fregion.rect Use \link{plot.fregion.rect} to visualize.
}
\description{
Builds a Rectangular Confidence Region
}
\examples{
# 1. Vector/matrix version

# Generate a sample
p = 200 ; N = 80 ; rangeval = c(0,1)
grid = make.grid(p, rangevals=rangeval)
mu0 = meanf.poly(grid,c(0,1)) ; names(mu0) = grid
mu = meanf.poly(grid,c(0,1.1)) ; names(mu) = grid
cov.m = make.cov.m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
x = make.sample(mu,cov.m,N)

# Find the estimate and covariance
hat.mu = rowMeans(x)
hat.cov.m = crossprod(t(x - hat.mu)) / (N-1)
e.hat.cov.m = eigen(hat.cov.m)   # <- This is optional and can be provide into the functions instead of hat.cov.m below.

# Make rectangular region and visulize
c <- fregion.rect(hat.mu-mu0,hat.cov.m,N=N)
plot(c)


# 2. fd/bifd version

# create basis, convert vector/matrix into fd/bifd objects.
require(fda)
nbasis <- round(p*0.9)
fd.basis <- create.bspline.basis(rangeval,nbasis)
mu0.fd <- Data2fd(names(mu0),mu0,fd.basis)
mu.fd <- Data2fd(names(mu),mu,fd.basis)
x.fd <- Data2fd(rownames(x),x,fd.basis)
hat.mu.fd <- mean.fd(x.fd)
hat.cov.fd <- var.fd(x.fd)
e.hat.cov.fd <- eigen.fd(hat.cov.fd)   # <- This is optional and can be provide into the functions instead of hat.cov.fd below.

# Make rectangular region and visulize
c.fd <- fregion.rect(hat.mu.fd-mu0.fd,hat.cov.fd,N=N)
plot(c.fd)
}
