#' Perform a hypothesis test based on a functional confidence region
#'
#' @param x Observed function. It can be either a vector, of fd object in 'fda' package.
#' @param x0 The function under the null hypothesis. zero function is assumed when not given.
#' @param cov The covaraicne operator (as matrix) for the random function x. The eigen decomposition of covariance operator can be used instead.
#' @param N Sample size. It should be '1' if 'cov' is the covariance operator for 'x' itself. If 'x' is the mean of an i.i.d. sample of size 'N' and the 'cov' is the covariance operator of each unit in the sample, then this 'N' can be used here.
#' @param type This specifies which regions to be used for the tests. It should be a collection of the followings : "Enorm", "Epc", "Ec", "Ec1", "Rz", "Rz1", "Rzs", or "Rz1s".
#' @param pc.cut It takes a vector of number of fpc to use in each testing. For integer values, fpc up to those values will be used. If it's a value from 0 to 1, this specifies the proportion of variance that should be explained by the fpcs. If it is 0, all the available fpcs will be used as long as the size of eigenvalues are greater than .Machine$double.eps.
#' @param prec This determines the accuracy of imhof. One may try to modify this if p-value achieved in Ellipsoid form other than "Epc" gives negative value. It should the the form of c(epsabs, epsrel, limit).
#' @param hat.cov This is estimated covariance operator, which will be used when 'cov' is given as a parametric one and small sample version is used. Usually not needed.
#' @param df Degree of freedom to use in small sample versions with suffix 's'
#' @return p-value achieved from it.
#' @examples
#' ## Generate Sample
#' p <- 100 ; N <- 50
#' grid <- make.grid(p)
#' mu0 <- mean.f.poly(grid,c(0,1)) ; names(mu0) <- grid
#' mu <- mean.f.poly(grid,c(0,1.1)) ; names(mu) <- grid
#' cov.m <- make.cov.m(cov.f <- cov.f.st.matern, grid=grid, cov.f.params=c(1/2,1,1/4))
#' X <- make.sample(mu,cov.m,N)
#' ## Find Estimates
#' hat.mu <- rowMeans(X)
#' hat.cov.m <- crossprod(t(X - hat.mu)) / (N-1)
#' ## Perform Test
#' fregion.test(hat.mu,mu0,hat.cov.m,N,type=c("ALL"),pc.cut=c(3,5,0.99,0.999,0))
#' @export

fregion.test <- function(x, x0=NULL, cov, N=1, type=c("Ec"), pc.cut=c(0.99), prec=NULL, hat.cov=NULL, df=NULL){

  ### 1. Check the data type ###
  if (inherits(x,"fd") & (inherits(cov,"bifd") | inherits(cov,"pca.fd") | inherits(cov,"eigen.fd"))) datatype="fd" else if
     (inherits(x,"numeric") & (inherits(cov,"matrix") | inherits(cov,"list"))) datatype="vector" else stop ("The format of data is unknown")

  ### 2. If covariance is given as it is (not eigen decomposition), take eigen decomposition.
  e.cov <- cov
  if (inherits(cov,"bifd")) e.cov <- eigen.fd(cov)
  if (inherits(cov,"matrix")) e.cov <- eigen(cov)

  #if (class(cov) == "matrix") {
  #  if (dim(cov)[1] != dim(cov)[2]) stop("Covariance matrix is not square")
  #  cov <- eigen(cov) # take eigen decomposition
  #} else if (class(cov) == "list") {
  #  if  ((sum(names(cov) == c("values", "vectors"))!=2)) stop("The format of 'cov' is not known -- a list, but not an eigen decomposition of a matrix")
  #} else if (class(cov) == "bifd") {
  #  cov.fd <- cov
  #} else if (class(cov) == "bifd")
  # 3. For all vector valued data, check whether dimensions match.
  # for matrix 'x', or 'x0', change them to vectors
  # if (class(x)=="matrix") { if (min(dim(x))!=1) stop("x should be a vector") else x <- as.vector(x) }
  # if (class(x0)=="matrix") { if (min(dim(x0))!=1) stop("x0 should be a vector") else x0 <- as.vector(x0) }
  # if x0 is null, make zero function. If not, check type and dimension
  # if (is.null(x0)) x0 <- numeric(p) else if (!(class(x0)=="numeric" & length(x0)==p)) stop("x0 doesn't seem to be a vector of the same length with x")

  # Trim cov to have all non-negative eigenvalues
  e.cov$values[e.cov$values < .Machine$double.eps] <- 0

  # 4. Center the function x
  x <- x - x0  # this works both for "fd" and vector.

  ################ Take loop for fpc.cut ##############
  result <- NULL
  for (i in c(1:length(pc.cut))){
    cut <- pc.cut[i]
    # 5. Find number of fpc to use.
    if (cut < 0) stop("fpc.cut should be some positive fraction from 0 to 1, or integer greater than or equal to 1, or just 0 (to use all available pcs)")
    if (cut == 0 ) cut <- sum(e.cov$values > .Machine$double.eps) else { # if fpc.cut is 0, use all PCs
       if (cut < 1) cut <- get.req.n.pc(cut,e.cov$values) else {
          if (cut != round(cut)) {cut=round(cut) ; print("fpc.cut[",i,"] was rounded to the closest integer")}
       }
    }
  #####################################################

    # 6. Run tests ; If 'ALL' is included then run all the tests
    if ( sum(c("all", "All", "ALL") %in% type) > 0 ) type <- c("Enorm","Epc","Ec","Ec1","Rz","Rz1","Rzs","Rz1s")
    if ("Enorm" %in% type)  pval.Enorm <- get.pval.Enorm(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut,prec=prec)
    if ("Epc" %in% type)    pval.Epc   <- get.pval.Epc(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut)
    if ("Ec" %in% type)     pval.Ec    <- get.pval.Ec(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut,prec=prec)
    if ("Ec1" %in% type)    pval.Ec1   <- get.pval.Ec1(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut,prec=prec)
    if ("Rz" %in% type)     pval.Rz    <- get.pval.Rz(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut)
    if ("Rz1" %in% type)    pval.Rz1   <- get.pval.Rz1(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut)
    if ("Rzs" %in% type)    pval.Rzs   <- get.pval.Rzs(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut,hat.cov=hat.cov,df=df)
    if ("Rz1s" %in% type)   pval.Rz1s  <- get.pval.Rz1s(x=x,sample.size=N,eigen=e.cov,fpc.cut=cut,hat.cov=hat.cov,df=df)

    names.pval <- ls(pattern=glob2rx("pval.*"))
    pvalues <- sapply(names.pval,get,inherits=FALSE,envir=1)
    names(pvalues) <- sub("pval.","",names(pvalues))
    pvalues <- c(pc.cut[i],cut,pvalues)
    result <- rbind(result,pvalues)
    rownames(result)[i] <- i
  }
  colnames(result)[c(1,2)] <- c("pc.cut","pc.used")
  return(result)
}
