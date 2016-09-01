#' Make confidence bands
#'
#' @param x Observed function. It can be either a vector, of fd object in 'fda' package.
#' @param cov The covaraicne operator (as matrix) for the random function x.
#' @param N Sample size. It should be '1' if 'cov' is the covariance operator for 'x' itself. If 'x' is the mean of an i.i.d. sample of size 'N' and the 'cov' is the covariance operator of each unit in the sample, then this 'N' can be used here.
#' @param type This specifies which bands to be constructed. It should be a collection of the followings : "Bs", "BEc", "BEc1".
#' @param conf.level A vector of confidence levels for the bands to achieve.
#' @param grid.size This determines on how fine grid the bands will be constructed. Note that this parameter is used only when 'x' is fd object and 'cov' is bifd object.
#' @param inv.method (Currently Not Used) This determines how the inverse of the sum of chi squares distribution will be achieved. It can be either "gamma.approx", or "inv.imhof".
#' @param prec (Currently Not Used) This determines the accuracy of imhof. It's used only when inv.method is inv.imhof.
#' @param sim.size This determines the simulation size for Bs
#' @return bands returns bands
#' @examples
#' p <- 100 ; N <- 50
#' grid <- make.grid(p)
#' mu <- mean.f.poly(grid,c(0,1)) ; names(mu) <- grid
#' cov.m <- make.cov.m(cov.f=cov.f.st.matern, grid=grid, cov.f.params=c(2/3,1,1/4))
#' X <- make.sample(mu,cov.m,N)
#' hat.mu <- rowMeans(X)
#' hat.cov.m <- crossprod(X - hat.mu) / (N-1)
#' fregion.band(hat.mu, hat.cov.m, N, type=c("Bs","BEc"), conf.level=c(0.90,0.95))
#' @export

fregion.band <- function(x, cov, N=1, type=c("Bs", "BEc"), conf.level=c(0.95), grid.size=200,
                        inv.method="gamma.approx", prec=NULL, sim.size=10000){


  ### 1. Check the data type ###
  if (inherits(x,"fd") & inherits(cov,"bifd")) datatype="fd" else if
     (inherits(x,"numeric") & inherits(cov,"matrix")) datatype="vector" else stop("The format of data is unknown")

  ### 2. Evaluate x and cov ###
  # evaluate x and cov if datatype is "fd".
  # Since all functions for generating bands evaluate fd object inside, we do this here and just use vector/matrix version
  if (datatype=="fd") {
    evalgrid <- make.grid(p=grid.size, rangevals=x$basis$rangeval)
    cov.m <- eval.bifd(evalgrid, evalgrid, cov)
    x.v <- eval.fd(evalgrid, x)
  } else {
    cov.m <- cov ; x.v <- x
  }
  p <- dim(cov.m)[1]


  ## Take loop for conf.level
  result <- as.matrix(x.v,ncol=1) ;
  colnames(result) <- c("x")
  ## Take eigen decomposition if BEc or BEc1 is used.
  if (sum(c("BEc","BEc1") %in% type) > 0) eigen.cov.m <- eigen(cov.m)
  eigen.cov.m$values[eigen.cov.m$values<0] <- 0 # trim negative eigenvalues.

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
    #class(result.fd) <- "fregion.band.fd"
    return(result.fd)
  } else {
    class(result) <- "fregion.band"
    return(result)
  }
}
