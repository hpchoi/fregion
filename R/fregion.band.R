#' Makes confidence bands
#'
#' @param x Functional parameter estimate. It can be either a vector or \link{fd} object from \link{fda}.
#' @param cov N * Cov(X), in which X is the functional estimator. It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.
#' @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
#' @param type The band(s) to be constructed.
#' \itemize{
#'   \item BEc : The suggested modified Scheffe style band from hyper-ellipsoie Ec, which uses up to the very last dimension.
#'   \item BEc1 : Another modified Scheffe style band from hyper-ellipsoie Ec1. (Wider than BEc although it does not require smoothness assumption. For comparision purpose only, not recommend for use)
#'   \item Bs : Parametric bootstrap simultaneous confidence band, similar to the one appeard in Degras(2011) (for comparision purpose)
#'   \item naive-t : A collection of point-wise t-intervals. (for comparision purpose)
#'   \item BEPC.x : Scheffe band from 'x' dimensional chi-square ellipse. (for comparision purpose)
#' }
#' @param conf.level A vector of confidence levels for the bands to achieve.
#' @param grid.size This determines on how fine grid the bands will be constructed before converted as an `fd' object. This parameter is used only when 'x' is fd object and 'cov' is bifd object.
#' @param inv.method (Currently Not Used) This determines how the inverse of the sum of chi squares distribution will be achieved. It can be either "gamma.approx", or "inv.imhof". Currently, only "gamma.approx" works.
#' @param prec (Currently Not Used) This determines the accuracy of imhof. It's used only when inv.method is inv.imhof.
#' @param sim.size This determines bootstrap sample size for Bs
#' @return fregion.band Either a collection of vector valued bands or `fd' object whose objectname is changed to fregion.band.
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
#' # Make and visualize/compare confidence bands
#' b <- fregion.band(hat.mu,hat.cov.m,N=N,type=c("Bs","BEc","BEPC.10","naive.t"),conf.level=c(0.95))
#' plot(b)
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
#' # Make and visualize/compare confidence bands.
#' b.fd <- fregion.band(hat.mu.fd,hat.cov.fd,N=N,type=c("Bs","BEc"), conf.level=c(0.95, 0.9))
#' plot(b.fd)
#'
#' @export

fregion.band <- function(x, cov, N=1, type=c("Bs", "BEc"), conf.level=c(0.95), grid.size=200,
                        inv.method="gamma.approx", prec=NULL, sim.size=10000){

  ### 1. Check the data type ###
  if (inherits(x,"fd") & (inherits(cov,"bifd") | inherits(cov,"pca.fd") | inherits(cov,"eigen.fd"))) datatype="fd" else if
     (inherits(x,"numeric") & (inherits(cov,"matrix") | inherits(cov,"list"))) datatype="vector" else stop ("The format of data is unknown")

  ### 2. Evaluate x and cov ###
  # evaluate x and cov if datatype is "fd".
  # Since all functions for generating bands evaluate fd object inside, we do this here and just use vector/matrix version
  if (datatype=="fd") {
    if (!inherits(cov,"bifd")) {
      J <- min(sum(cov$values > 0),dim(cov$harmonics$coefs)[2])
      coef <- cov$harmonics$coefs[,c(1:J)] %*% diag(cov$values[c(1:J)]) %*% t(cov$harmonics$coefs[,c(1:J)])
      cov <- bifd(coef,cov$harmonics$basis,cov$harmonics$basis)
    }
    evalgrid <- make.grid(p=grid.size, rangevals=x$basis$rangeval)
    cov.m <- eval.bifd(evalgrid, evalgrid, cov)
    x.v <- eval.fd(evalgrid, x)
  } else {
    if (inherits(cov,"list")) {
      J <- sum(cov$values > 0)
      cov.m <- cov$vectors[,c(1:J)] %*% diag(cov$values[c(1:J)]) %*% t(cov$vectors[,c(1:J)])
    } else cov.m <- cov;
    x.v <- x
  }
  p <- dim(cov.m)[1]
  if (!isSymmetric(cov.m)) cov.m <- (cov.m + t(cov.m))/2  # force cov.m to be symmetric

  ## Take loop for conf.level
  result <- as.matrix(x.v,ncol=1) ;
  colnames(result) <- c("x")
  ## Take eigen decomposition if BEc or BEc1 is used.
  if (  (sum(c("BEc","BEc1") %in% type) > 0) | (length(grep("BEPC.", type)) > 0) ) {
    eigen.cov.m <- eigen(cov.m) ; eigen.cov.m$values[ eigen.cov.m$values < 0 ] <- 0 # trim negative eigenvalues.
  }

  for (i in c(1:length(conf.level))){
    level <- conf.level[i]
    # 5. Find number of fpc to use.
    if (!(level > 0 & level < 1)) stop("conf.level should have values between 0 and 1")

    # 6. Make bands
    # If 'ALL' is included then run all the tests

    if ("Bs" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("Bs.u.",level), paste0("Bs.l.",level))
      Bs <- make.band.Bs(cov=cov.m,conf.level=level,sim.size=sim.size) / sqrt(N)
      result <- cbind(result, x.v + Bs, x.v - Bs);
      colnames(result) <- tmp.colnames
    }

    if ("BEc" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("BEc.u.",level), paste0("BEc.l.",level))
      BEc <- make.band.BEc(eigen=eigen.cov.m, conf.level=level) / sqrt(N)
      result <- cbind(result, x.v + BEc, x.v - BEc);
      colnames(result) <- tmp.colnames
    }

    if ("BEc1" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("BEc1.u.",level), paste0("BEc1.l.",level))
      BEc <- make.band.BEc1(eigen=eigen.cov.m, conf.level=level) / sqrt(N)
      result <- cbind(result, x.v + BEc1, x.v - BEc1);
      colnames(result) <- tmp.colnames
    }

    if ("naive.t" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("naive.t.u.",level), paste0("naive.t.l.",level))
      naive.t <- make.band.naive.t(cov.m, conf.level=level, df=N-1) / sqrt(N)
      result <- cbind(result, x.v + naive.t, x.v - naive.t);
      colnames(result) <- tmp.colnames
    }

    if ("BRz" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("BRz.u.",level), paste0("BRz.l.",level))
      BRz <- make.band.BRz(eigen=eigen.cov.m, conf.level=level) / sqrt(N)
      result <- cbind(result, x.v + BRz, x.v - BRz);
      colnames(result) <- tmp.colnames
    }

    if ("BRz0" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("BRz0.u.",level), paste0("BRz0.l.",level))
      BRz0 <- make.band.BRz0(eigen=eigen.cov.m, conf.level=level) / sqrt(N)
      result <- cbind(result, x.v + BRz0, x.v - BRz0);
      colnames(result) <- tmp.colnames
    }


    tmpJs.location <- grep("BEPC.",type)
    if (length(tmpJs.location) > 0) {
      tmpJs <- as.integer(sub("BEPC.","",type[tmpJs.location]))
      for (j in tmpJs){
        tmp.colnames <- c(colnames(result), paste0("BEPC.",j,".u.",level), paste0("BEPC.",j,".l.",level))
        BEPC <- make.band.BEPC(eigen=eigen.cov.m, conf.level=level, J=j) / sqrt(N)
        result <- cbind(result, x.v + BEPC, x.v - BEPC);
        colnames(result) <- tmp.colnames
      }
    }

  }
  if (datatype=="fd") {
    result.fd <- Data2fd(evalgrid, result, basisobj=x$basis)
    class(result.fd) <- "fregion.band"
    return(result.fd)
  } else {
    class(result) <- "fregion.band"
    return(result)
  }
}
