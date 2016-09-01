
## What this package is about

This package contains examplary codes for methods suggested in http://arxiv.org/abs/1607.07771 (Geometric Approach to Confidence Regions and Bands for Functional Parameters)


## Some top-level functions (with help available)

- fregion.test() : hypothesis testings
- fregion.band() : make confidence bands
- fregion.rect() : construct rectangular confidence region in Hilbert space.

- plot.fregion.band() : Make a plot for confidence bands. This works only for the vector version since for 'fd' object, the fda package will handle itslef (although it misses some features this function provides.)
- plot.fregion.rect() : Provides visualization tool for rectangular regions.

## More base-line functions

- get.pval.*
- get.crit.*
- make.band.*

These functions work behind the top-level functions but may be found to be more useful.

## Input values

Most functions need 1) an estimate (as vector, or 'fd' object of fda package) and 2) covariance (or eigen decomposition of it) of the estimator. The covariance can be either matrix or 'bifd' object of fda package. The eigen decomposition can take the typical 'list' type (for matrix), 'pca.fd' object of fda package, of `eigen.fd' object from this package, which is similar to 'pca.fd'.

## Other functions
There're more functions in this package (but without documentation) to support those main functions shown above.

## Example Codes

### vector/matrix version

#### Generate Sample
p = 200 ; N = 80 ; rangeval = c(0,1) 

grid = make.grid(p, rangevals=rangeval)

mu0 = mean.f.poly(grid,c(0,1)) ; names(mu0) = grid

mu = mean.f.poly(grid,c(0,1.1)) ; names(mu) = grid

cov.m = make.cov.m(cov.f = cov.f.st.matern, grid=grid, cov.f.params=c(2/2,1,1))

e.cov.m = eigen(cov.m)

X = make.sample(mu,cov.m,N)

#### Find Estimates
hat.mu = rowMeans(X)

hat.cov.m = crossprod(t(X - hat.mu)) / (N-1)

e.hat.cov.m = eigen(hat.cov.m)

#### Plot sample and the sample mean 
matplot(grid,X,type="l")

lines(grid,hat.mu,lwd=2)

#### Hypothesis Tests
(a1=fregion.test(hat.mu,mu0,hat.cov.m,N,type=c("ALL"),pc.cut=c(1,3,4,5,0.99,0.999)))

#### Make confidence bands
b = fregion.band(hat.mu,hat.cov.m,N=N,type=c("Bs","BEc","BEPC.10","naive.t"),conf.level=c(0.95))

plot(b)

#### Make ellipsoie region and visulize
c = fregion.rect(hat.mu,hat.cov.m,N=N)

plot(c)


### 'fda'version

#### create basis 
nbasis = p

fd.basis = create.bspline.basis(rangeval,nbasis)

#### convert vectors in to fd objects
mu0.fd = Data2fd(names(mu0),mu0,fd.basis)

mu.fd = Data2fd(names(mu),mu,fd.basis)

X.fd = Data2fd(rownames(X),X,fd.basis)

cov.fd =  cov.m.to.cov.fd(cov.m, grid, fd.basis)


hat.mu.fd = mean.fd(X.fd)

hat.cov.fd = var.fd(X.fd)

hat.cov.fd.pca = pca.fd(X.fd,nharm=nbasis)

e.hat.cov.fd = eigen.fd(hat.cov.fd)


#### Plot sample and the sample mean ###
plot(X.fd)

plot(hat.mu.fd,lwd=2,add=TRUE)

#### Hypothesis Tests
(a1.fd=fregion.test(hat.mu.fd,mu0.fd,hat.cov.fd,N,type=c("ALL"),pc.cut=c(1,3,4,5,0.99,0.999)))
(a1.fd.v1=fregion.test(hat.mu.fd,mu0.fd,hat.cov.fd.pca,N,type=c("ALL"),pc.cut=c(1,3,4,5,0.99,0.999)))
(a1.fd.true=fregion.test(hat.mu.fd,mu0.fd,cov.fd,N,type=c("ALL"),pc.cut=c(1,3,4,5,0.99,0.999,0.99999)))


#### Make confidence bands
b.fd = fregion.band(hat.mu.fd,hat.cov.fd,N=N,type=c("Bs","BEc","BEPC.10","naive.t"),conf.level=c(0.95))

plot(b.fd)

##### compare b.fd with b
par(mfrow=c(1,2)) ;plot(b) ;plot(b.fd); par(mfrow=c(1,1))


#### Make ellipsoie region and visulize
c.fd = fregion.rect(hat.mu.fd,hat.cov.fd,N=N)
plot(c.fd)

c.fd.true = fregion.rect(hat.mu.fd, cov.fd, N=N, pc.cut=0)
plot(c.fd.true)


