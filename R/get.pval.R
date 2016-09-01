#' @export
get.pval.Enorm <- function(x, sample.size, eigen, fpc.cut=NULL, prec=NULL){ ## prec : precision option for CompQuadForm::imhof
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut)) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (inherits(x,"fd")) {Enorm <- sample.size * sum( inprod(x,eigen$harmonics[1:fpc.cut])^2) } else  ## fd object
                        {Enorm <- sample.size * sum( crossprod(x,eigen$vectors[,1:fpc.cut])^2 ) }            ## vector
  CompQuadForm::imhof(Enorm,eigen$values[1:fpc.cut],epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]] ;
}

#' @export
get.pval.Epc <- function(x, sample.size, eigen, fpc.cut=NULL){
  if (is.null(fpc.cut)) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (inherits(x,"fd")) {Epc <- sample.size * sum( inprod(x, eigen$harmonics[1:fpc.cut])^2 / eigen$values[1:fpc.cut] )} else
                        {Epc <- sample.size * sum( crossprod( x, eigen$vectors[,1:fpc.cut])^2 / eigen$values[1:fpc.cut] )}
  1-pchisq(Epc, fpc.cut)
}

#' @export
get.pval.Ec <- function(x, sample.size, eigen, fpc.cut=NULL, prec=NULL){
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut)) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  c.square <- sqrt(eigen$values[1:fpc.cut])
  if (inherits(x,"fd")) {coef.square <- inprod(x, eigen$harmonics[1:fpc.cut])^2} else
                        {coef.square <- as.vector(crossprod(x, eigen$vectors[,1:fpc.cut])^2)}
  weights <- eigen$values[1:fpc.cut]/c.square/sample.size
  stat <- sum(coef.square / c.square)
  CompQuadForm::imhof(stat, weights,epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]]
}

#' @export
get.pval.Ec1 <- function(x, sample.size, eigen, fpc.cut=NULL, prec=NULL){
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut)) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  c.square <- sqrt(rev(cumsum(rev(eigen$values[1:fpc.cut])))) ## The only difference from Ec
  if (inherits(x,"fd")) {coef.square <- inprod(x, eigen$harmonics[1:fpc.cut])^2} else
                        {coef.square <- as.vector(crossprod(x, eigen$vectors[,1:fpc.cut])^2)}
  weights <- eigen$values[1:fpc.cut]/c.square/sample.size
  stat <- sum(coef.square / c.square)
  CompQuadForm::imhof(stat, weights,epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]]
}

#' @export
get.pval.Rz <- function(x, sample.size, eigen, fpc.cut=NULL){
  if (is.null(fpc.cut)) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (inherits(x,"fd")) {Zs.actual <- abs( sqrt(sample.size) * inprod(x, eigen$harmonics[1:fpc.cut])   / sqrt(eigen$values[1:fpc.cut])) } else
                        {Zs.actual <- abs( sqrt(sample.size) * crossprod(x, eigen$vectors[,1:fpc.cut]) / sqrt(eigen$values[1:fpc.cut])) }
  min(1-exp(sum(eigen$values[1:fpc.cut]) / eigen$values[1:fpc.cut] * log(pnorm2(Zs.actual))))
}

#' @export
get.pval.Rz1 <- function(x, sample.size, eigen, fpc.cut=NULL){
  if (is.null(fpc.cut)) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (inherits(x,"fd")) {Zs.actual <- abs( sqrt(sample.size) * inprod(x, eigen$harmonics[1:fpc.cut])   / sqrt(eigen$values[1:fpc.cut])) } else
                        {Zs.actual <- abs( sqrt(sample.size) * crossprod(x, eigen$vectors[,1:fpc.cut]) / sqrt(eigen$values[1:fpc.cut])) }
  M.star <- max(sqrt(2*pi) * eigen$values[1:fpc.cut] * sf.f1(Zs.actual))
  Zs.star <- sapply(M.star / (sqrt(2*pi) * eigen$values[1:fpc.cut]), sf.f1.inv)
  1 - prod(pnorm2(Zs.star))
}

#' @export
get.pval.Rzs <- function(x, sample.size, eigen, fpc.cut=NULL, hat.cov=NULL, df=NULL){
  if (is.null(hat.cov)) tilde.lambda <- eigen$values else {
    if (inherits(x,"fd")) { tilde.lambda <- diag(inprod(eigen$harmonics$basis, hat.cov$sbasis) %*% hat.cov$coefs %*%
                                                inprod(hat.cov$tbasis, eigen$harmonics$basis))} else
                          { tilde.lambda <- diag(t(eigen$vectors) %*% hat.cov %*% eigen$vectors) } # for vector x
  }
  if (is.null(fpc.cut)) fpc.cut <- sum(tilde.lambda > .Machine$double.eps)
  if (inherits(x,"fd")) {Ts.actual <- abs( sqrt(sample.size) * inprod(x, eigen$harmonics[1:fpc.cut] )   / sqrt(tilde.lambda[1:fpc.cut]))} else
                        {Ts.actual <- abs( sqrt(sample.size) * crossprod(x, eigen$vectors[,1:fpc.cut] ) / sqrt(tilde.lambda[1:fpc.cut]))}
  if (is.null(df)) {df=sample.size-1}
  min(1-exp(sum(eigen$values[1:fpc.cut]) / eigen$values[1:fpc.cut] * log(pt2(Ts.actual,df=df))))
}

#' @export
get.pval.Rz1s <- function(x, sample.size, eigen, fpc.cut=NULL, hat.cov=NULL, df=NULL){
  if (is.null(hat.cov)) tilde.lambda <- eigen$values else {
    if (inherits(x,"fd")) { tilde.lambda <- diag(inprod(eigen$harmonics$basis, hat.cov$sbasis) %*% hat.cov$coefs %*%
                                                  inprod(hat.cov$tbasis, eigen$harmonics$basis))} else
                          { tilde.lambda <- diag(t(eigen$vectors) %*% hat.cov %*% eigen$vectors) } # for vector x
  }
  if (is.null(fpc.cut)) fpc.cut <- sum(tilde.lambda > .Machine$double.eps)
  if (inherits(x,"fd")) {Ts.actual <- abs( sqrt(sample.size) * inprod(x, eigen$harmonics[1:fpc.cut] )   / sqrt(tilde.lambda[1:fpc.cut]))} else
                        {Ts.actual <- abs( sqrt(sample.size) * crossprod(x, eigen$vectors[,1:fpc.cut] ) / sqrt(tilde.lambda[1:fpc.cut]))}
  if (is.null(df)) {df=sample.size-1}
  Zs.actual <- qnorm2(pt2(Ts.actual,df=df))
  M.star <- max(sqrt(2*pi) * eigen$values[1:fpc.cut] * sf.f1(Zs.actual))
  Zs.star <- sapply(M.star / (sqrt(2*pi) * eigen$values[1:fpc.cut]), sf.f1.inv)
  1 - prod(pnorm2(Zs.star))
}
