#' Make Rectangular Confidence Regoion
#'
#' @param x Observed function. It can be either a vector or fd object in 'fda' package.
#' @param cov The covaraicne operator (as matrix, or 'bifd' object) for the random function x.
#' @param N Sample size. It should be '1' if 'cov' is the covariance operator for 'x' itself. If 'x' is the mean of an i.i.d. sample of size 'N' and the 'cov' is the covariance operator of each unit in the sample, then this 'N' can be used here.
#' @param type This specifies which rectangular region to be constructed. It should be either "Rz", "Rzs", "Rz1", "Rz1s".
#' @param conf.level This specifies confidence level.
#' @param pc.cut This specifies number of fpc to use. For integer values, fpc up to those values will be used. If it's a value from 0 to 1, this specifies the proportion of variance that should be explained by the fpcs. If it is 0, all the available fpcs will be used as long as the size of eigenvalues are greater than .Machine$double.eps.
#' @param df Degree of freedom used for "Rzs" and "Rz1s". If null, N-1 will be used.
#' @examples
#' p <- 100 ; N <- 50
#' grid <- make.grid(p)
#' mu <- mean.f.poly(grid,c(0,1)) ; names(mu) <- grid
#' cov.m <- make.cov.m(cov.f <- cov.f.st.matern, grid=grid, cov.f.params=c(2/3,1,1/4))
#' X <- make.sample(mu,cov.m,N)
#' hat.mu <- rowMeans(X)
#' hat.cov.m <- crossprod(X - hat.mu) / (N-1)
#' fregion.rect(hat.mu, hat.cov.m, N, type="Rz", conf.level= 0.95)
#' @export


fregion.rect <- function(x, cov, N=1, type="Rz", conf.level=0.95, pc.cut=0.999, df=NULL){

  ### 1. Check the data type ###
  if (inherits(x,"fd") & (inherits(cov,"bifd") | inherits(cov,"pca.fd") | inherits(cov,"eigen.fd"))) datatype="fd" else if
  (inherits(x,"numeric") & (inherits(cov,"matrix") | inherits(cov,"list"))) datatype="vector" else stop ("The format of data is unknown")

  ### 2. If covariance is given as it is (not eigen decomposition), take eigen decomposition.
  e.cov <- cov
  if (inherits(cov,"bifd")) e.cov <- eigen.fd(cov)
  if (inherits(cov,"matrix")) e.cov <- eigen(cov)

  # Trim cov to have all non-negative eigenvalues (precautionary)
  e.cov$values[e.cov$values<0] <- 0

  if (is.null(df)) df <- N - 1
  cut=pc.cut
  if (cut < 0) stop("fpc.cut should be some positive fraction from 0 to 1, or integer greater than or equal to 1, or just 0 (to use all available pcs)")
  if (cut == 0 ) cut <- sum(e.cov$values > .Machine$double.eps) else { # if fpc.cut is 0, use all PCs
    if (cut < 1) cut <- get.req.n.pc(cut,e.cov$values) else {
      if (cut != round(cut)) {cut=round(cut) ; print("fpc.cut[",i,"] was rounded to the closest integer")}
    }
  }
  pc.range <- c(1:cut)


  eigenvalues <- e.cov$values[pc.range]

  harmonics <- eigenfunctions <- NULL
  if (datatype=="fd"){
    harmonics <- e.cov$harmonics[pc.range]
  } else {
    eigenfunctions <- e.cov$vectors[,pc.range]
    rownames(eigenfunctions) <- names(x)
    colnames(eigenfunctions) <- paste0("v",pc.range)
  }

  if (datatype=="fd") {coef <- as.vector(inprod(x, harmonics))} else
                     {coef <- as.vector(crossprod(x, eigenfunctions))}
  score <- coef / sqrt(eigenvalues)
  if ("Rz" %in% type) {mult <- get.crit.Rz(eigenvalues,conf.level) ; marginal.coverage <- pnorm2(mult) }
  if ("Rz1" %in% type) {mult <- get.crit.Rz1(eigenvalues,conf.level) ; marginal.coverage <- pnorm2(mult) }
  if ("Rzs" %in% type) {mult <- get.crit.Rzs(eigenvalues,conf.level,df) ; marginal.coverage <- pt2(mult,df) }
  if ("Rz1s" %in% type) {mult <- get.crit.Rz1s(eigenvalues,conf.level,df) ; marginal.coverage <- pt2(mult,df) }
  #boundarybands <- crossprod( t(cov$vectors[,pc.range]), diag(boundaryscores * sqrt(cov$values[pc.range])) )
  #actualcurves <- crossprod( t(cov$vectors[,pc.range], diag(actualscores * sqrt(cov$values[pc.range]))) )
  ME <- mult * sqrt(eigenvalues)
  cum.var.explained <- cumsum(eigenvalues) / sum(e.cov$values[e.cov$values>0])
  var.explained <- eigenvalues / sum(e.cov$values[e.cov$values>0])
  table <- data.frame(PC=factor(pc.range), coefficient=coef, MoE=ME, score=score, abs.score=abs(score),
                     score.bound=mult, var.explained=var.explained, cum.var.explained=cum.var.explained, marginal.coverage=marginal.coverage)
  result <- list(table=table, eigenfunctions=eigenfunctions, harmonics=harmonics, type=type)
  class(result) <- "fregion.rect"
  return(result)
}
