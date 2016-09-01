#' @export
pnorm2 <- function(z){1-2*(1-pnorm(z))} #### make two-sided version of pnorm with only positive z
#' @export
qnorm2 <- function(p){-qnorm((1-p)/2)}
#' @export
pt2 <- function(t,df){1-2*(1-pt(t,df))}
#' @export
qt2 <- function(p,df){-qt((1-p)/2,df)}
#' @export
make.grid <- function(p=100,rangevals=c(0,1),type="open"){
  if (type=="open") seq(0.5/p, (p-0.5)/p, by=1/p) * (rangevals[2] - rangevals[1]) + rangevals[1] else
    seq(0, 1, by=1/p) * (rangevals[2] - rangevals[1]) + rangevals[1]
}
#' @export
make.cov.m <- function(cov.f=cov.f.st.matern, grid=100, cov.f.params=NULL){  ### Make cov. matrix from cov. function.
  if (length(grid)==1) {grid=make.grid(p=grid)} ## input grid as a single number (as grid size), or as vector (actual grid)
  grid.size <- length(grid)
  cov.m <- matrix(0,nrow=grid.size,ncol=grid.size)
  if (is.null(cov.f.params)) {
    for (i in (1:grid.size)){
      cov.m[i,]=sapply(grid, cov.f, x1=grid[i])
    }
  }
  else{
    for (i in (1:grid.size)){
      cov.m[i,]=sapply(grid, cov.f, x1=grid[i], params=cov.f.params)
    }
  }
  return(cov.m)
}

#' @export
make.sample <- function(mean.v,cov.m,N,dist="rnorm",...){
  p <- length(mean.v)
  if (p != dim(cov.m)[1] | p != dim(cov.m)[2]) stop("Dimensions of mean vector and cov. matrix do not match")
  dist <- get(dist, mode <- "function", envir <- parent.frame())
  Z <- matrix(dist(N*p,...),nrow=p,ncol=N)
  eigen.cov.m <- eigen(cov.m);
  eigen.cov.m$values[eigen.cov.m$values<0] <- 0
  X <- crossprod( t(eigen.cov.m$vectors), crossprod(diag(sqrt(eigen.cov.m$values)), Z)) + mean.v
  rownames(X) <- names(mean.v)
  colnames(X) <- paste0("x",c(1:N))
  return(X)
}

#' @export
get.req.n.pc <- function(proportion, lambdas){ # Finds how many pcs are reuqired to achieve desired variance proportion(prop. as 1-10^(-x))
  satisfied <- (proportion <= cumsum(lambdas) / sum(lambdas))
  min(c(1:length(lambdas))[satisfied])
}

#' @export
get.schisq.q.gamma <- function(weights,prob){  # prob as 1-alpha
  meanvalue <- sum(weights)
  varvalue <- sum(2*weights^2)
  gammascale <- varvalue / meanvalue   # 's'
  gammashape <- meanvalue / gammascale # 'a'
  qgamma(p=prob,shape=gammashape,scale=gammascale,lower.tail=TRUE)
}

#' @export
eigen.fd <- function(cov.fd){
  BtB <- inprod(cov.fd$sbasis,cov.fd$tbasis) ## sbisis and tbasis should be the same.
  B <- chol(BtB)
  coefs <- cov.fd$coefs
  coefsU <- B%*%coefs%*%t(B)
  e.coefsU <- eigen(coefsU)
  vcoefs <-  solve(B) %*% e.coefsU$vectors
  harmonics <- fd(vcoefs,cov.fd$sbasis)
  rtnobj <- list(values=e.coefsU$values, harmonics=harmonics)
  class(rtnobj) <- "eigen.fd"
  return(rtnobj)
}

#' @export
eigen.fd2 <- function(cov.fd){
  W2inv <- inprod(cov.fd$sbasis,cov.fd$tbasis) ## sbisis and tbasis should be the same.
  coefs <- cov.fd$coefs
  A <- coefs %*% W2inv
  e.A <- eigen(A)
  e.A$values <- Re(e.A$values)
  e.A$vectors <- Re(e.A$vectors)
  cut <- sum(e.A$values > .Machine$double.eps)
  vcoefs <- e.A$vectors[,1:cut]
  normv <- diag(t(vcoefs) %*% W2inv %*% vcoefs)
  vcoefs <- vcoefs %*% diag(1/sqrt(normv))
  harmonics <- fd(vcoefs,cov.fd$sbasis)
  rtnobj <- list(values=e.A$values,harmonics=harmonics)
  class(rtnobj) <- "eigen.fd"
  return(rtnobj)
}

#' @export
cov.m.to.cov.fd <- function(cov.m, eval.points, basisobj){
  e.cov.m <- eigen(cov.m)
  basis.evaluated <- eval.basis(eval.points, basisobj)
  M <- t(as.matrix(basis.evaluated)) %*% e.cov.m$vectors
  coefs <- M %*% diag(e.cov.m$values) %*% t(M)
  bifd(coefs,basisobj,basisobj)
}
