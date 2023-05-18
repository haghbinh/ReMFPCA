require(fda)
require(ReMFPCA)
bs1 <- create.fourier.basis(c(0,2*pi),5)
bs2 <- create.bspline.basis(c(0,1),7)
bs3 <- create.exponential.basis(c(0,2),3)

####### 1-D Basis ######## (similar to the fd features) 
mdbs1 <- Basismfd(bs1)

mdbs1$basis
mdbs1$dimSupp
mdbs1$nbasis
mdbs1$supp
mdbs1$gram
mdbs1$eval(1:7/10)
image(as.matrix(mdbs1$gram))

####### 2-D Basis ######## (fd cannot handle this)
mdbs2 <- Basismfd(bs1,bs2)
mdbs2$basis
mdbs2$dimSupp
mdbs2$nbasis
mdbs2$supp
mdbs2$gram
arg_mdbs <- list(1:10,1:9/10)
mdbs2$eval(arg_mdbs)
image(as.matrix(mdbs2$gram))
