#' @export
mean.f.poly <- function(x,params=c(0,1)){  params[1] + params[2]*(10*x^3 - 15*x^4 + 6*x^5)} #params=c(shift,scale)

#' @export
mean.f.peak <- function(x,params=c(0,1,2,16,1)){
  peak <- params[3] - params[4]*abs(x-0.5)^(params[5])
  peak[peak<0] <- 0   #params=c(shift,scale,peak-top,peak-sharpness,peak-power )
  params[2] * (params[1] + peak)
}
