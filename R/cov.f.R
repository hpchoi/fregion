#' @export
cov.f.st.matern <- function(x1,x2,params=c(1/2,1,1)){
  nu <- params[1] ; l <- params[2] ; sigma <- params[3]  #params=c(nu,l,sigma)
  d  <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}

#' @export
cov.f.st.matern.warp.power <- function(x1,x2,params=c(1/2,1,1,10)){ #params <- c(nu,l,sigma,power)
  cov.f.st.matern(warp.f.power(x1,params[4]),warp.f.power(x2,params[4]),params=params[1:3])
}

#' @export
cov.f.BM <- function(x1,x2,params=c(1)){ params[1]^2 * (abs(x1)+abs(x2)-abs(x2-x1))/2 } # min(x2,x1)
