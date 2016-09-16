#' Visualizes confidence bands constructed from \link{fregion.band}.
#'
#' @param band A 'fregion.band' object, the output from \link{fregion.band} funciton.
#' @param col a vector of colors to be used for the bands.
#' @param lty a vector of colors to be used for the bands.
#' @param center Whether to include the functional estimate or not.
#' @export

plot.fregion.band <- function(band,col=NULL,lty=NULL,center=TRUE){
  if (is.null(dim(band))) {type <- "fd" ; class(band) <- "fd"} else type <- "vector"

  if (type=="vector") bandnames <- colnames(band) else bandnames <- colnames(band$coefs)
  bandnames <- sub(".u","",bandnames)
  bandnames <- sub(".l","",bandnames)
  bandnames <- bandnames[c(T,F)]

  if (is.null(col)) col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628')
  if (is.null(lty)) lty <- c(5,3,4,2,6)

  if (type=="vector") {
    matplot(rownames(band),band[,-1],lty=rep(lty,each=2),col=rep(col,each=2),type="l",xlab="T",ylab="bands")
  } else {
    fda::plot.fd(band[-1],lty=rep(lty,each=2),col=rep(col,each=2),xlab="T",ylab="bands")
  }

  if (center) {
    if (type=="vector") {
      lines(rownames(band),band[,1],col="black",lty=1)
    } else {
      fda::plot.fd(band[1],col="black",lty=1,add=TRUE)
    }
    legend(x="topleft",legend=c("Estimate",bandnames[-1]),col=c(1,col),lty=c(1,lty))
  } else {
    legend(x="topleft",legend=bandnames[-1],col=c(col),lty=c(lty))
  }
}

