#' @export
covf.st.matern <- function(x1,x2,params=c(1/2,1,1)){
  nu <- params[1] ; l <- params[2] ; sigma <- params[3]  #params=c(nu,l,sigma)
  d  <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}

#' @export
covf.st.matern.warp.power <- function(x1,x2,params=c(1/2,1,1,10)){ #params <- c(nu,l,sigma,power)
  covf.st.matern(warpf.power(x1,params[4]),warpf.power(x2,params[4]),params=params[1:3])
}
