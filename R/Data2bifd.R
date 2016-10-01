#' Creates a bifd object from a matrix.
#'
#' @param sargvals Argument values for `s' in the matrix input y(s,t).
#' @param targvals Argument values for `t' in the matrix input y(s,t).
#' @param y Data (only one) as a matrix.
#' @param sbasisobj basisfd object for `s'
#' @param tbasisobj basisfd object for `t'
#' @param Lfdobjs Number of derivatives used for smoothing along `s'. No smoohting will be done for NULL
#' @param Lfdobjt Number of derivatives used for smoothing along `t'. No smoohting will be done for NULL
#' @param lambdas Weight on the smoothing along `s'
#' @param lambdat Weight on the smoothing along `t'
#' @export

Data2bifd <- function(sargvals=NULL, targvals=NULL, y=NULL,
                      sbasisobj=NULL, tbasisobj=NULL,
                      Lfdobjs=NULL, Lfdobjt=NULL, lambdas=NULL, lambdat=NULL){
  if (is.null(lambdas)) lambdas <- 3e-8/diff(as.numeric(range(sargvals)))
  if (is.null(lambdat)) lambdat <- 3e-8/diff(as.numeric(range(sargvals)))


  Bs <- eval.basis(sargvals,sbasisobj)
  Bt <- eval.basis(targvals,tbasisobj)

  if (!is.null(Lfdobjs)) Ps <- getbasispenalty(sbasisobj,Lfdobj=Lfdobjs) else Ps <- 0
  if (!is.null(Lfdobjt)) Pt <- getbasispenalty(sbasisobj,Lfdobj=Lfdobjt) else Pt <- 0

  coef <- t(chol2inv(chol(crossprod(Bt) + lambdat * Pt)) %*%
            t(Bt) %*% t(y) %*% Bs %*%
            chol2inv(chol(crossprod(Bs) + lambdas * Ps)))
  return(bifd(coef,sbasisobj=sbasisobj,tbasisobj=tbasisobj))
}
