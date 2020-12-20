#' @title Generate sample from the special distribution
#' @description The pdf of distribution if fe=3/4*(1-x^2),|x|<=1
#' Devroye and Gyorfi give the following algorithm for simulation from this distribution.
#' @param n the sample size
#' @return a random sample of size \code{n} from fe
#' @importFrom stats runif
#' @useDynLib StatComp20056
#' @examples
#' \dontrun{
#' rangfe(10)
#' }
#' @export
ranfe<-function(n){## the function to generate random variates from fe
  u1<-runif(n,-1,1);u2<-runif(n,-1,1);u3<-runif(n,-1,1)
  x1<-rep(0,n)
  for(i in 1:n){
    if(abs(u3[i])>=abs(u2[i])&&abs(u3[i])>=abs(u1[i])){
      x1[i]<-u2[i]
    }else{
      x1[i]<-u3[i]
    }
  }
  return(x1)
}

#' @title  Compute the sample skewness statistic.
#' @description The calculated statistic is used for skewness test of normality 
#' and has an asymptotic distribution of N(0,n/6).
#' @param x the sample vector
#' @return sample skewness statistic
#' @examples
#' \dontrun{
#' x<-rnorm(100)
#' sktest<-(abs(sk(x))>= qnorm(.975, 0, sqrt(6/n)))
#' }
#' @export
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

#' @title Power of the skewness test of symmetric Beta distributions.
#' @description Estimate the power of the skewness test of normality against 
#' symmetric Beta(alpha,alpha) distributions by using the sample skewness 
#' statistic function sk(). The value of power is the proportion of null hypothesis 
#' (skewness of distribution is 0) rejected in m replicates.
#' @param para parameter vector of symmetric Beta distributions
#' @param level significance level 
#' @param n size of sample drawn from Beta distribution for each experiment
#' @param m number of repeated experiments
#' @return the power of skewness test
#' @importFrom stats qnorm rbeta
#' @examples
#' \dontrun{
#' alpha<-c(seq(0.5,10,0.5),11:30)
#' pwr<-power(alpha)
#' plot(para1, pwr, type = "b",xlab = bquote(alpha), ylim = c(0,0.1))
#' }
#' @export
power<-function(para,level=0.05,n=50,m=2500){
  N <- length(para)
  pwr <- numeric(N)
  #critical value for the skewness test
  cv <- qnorm(1-level/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
  for (j in 1:N){
    e <- para[j]
    sktests <- numeric(m)
    for (i in 1:m){
      x <- rbeta(n, para[j],para[j])
      sktests[i] <- as.integer(abs(sk(x)) >= cv)
    }
    pwr[j] <- mean(sktests)
  }
  return(pwr)
}

#' @title Power of the skewness test of student t distributions.
#' @description Estimate the power of the skewness test of normality against 
#' student t(nu) distributions by using the sample skewness 
#' statistic function sk(). The value of power is the proportion of null hypothesis 
#' (skewness of distribution is 0) rejected in m replicates.
#' @param para parameter vector of t distributions
#' @param level significance level 
#' @param n size of sample drawn from t distribution for each experiment
#' @param m number of repeated experiments
#' @return the power of skewness test
#' @importFrom stats qnorm rt
#' @examples
#' \dontrun{
#' alpha<-c(seq(0.5,10,0.5),11:30)
#' pwr<-power_stu(alpha)
#' plot(para1, pwr, type = "b",xlab = bquote(alpha), ylim = c(0,0.1))
#' }
#' @export
power_stu<-function(para,level=0.05,n=50,m=2500){
  N <- length(para)
  pwr <- numeric(N)
  #critical value for the skewness test
  cv <- qnorm(1-level/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
  for (j in 1:N){
    e <- para[j]
    sktests <- numeric(m)
    for (i in 1:m){
      x <- rt(n, para[j])
      sktests[i] <- as.integer(abs(sk(x)) >= cv)
    }
    pwr[j] <- mean(sktests)
  }
  return(pwr)
}

#' @title  Equal variance test for two samples(equal size)
#' @description It is applied for independent random samples when the random variables are similarly 
#' distributed and sample sizes are equal. For two independent random samples of 
#' equal size, the maximum number of extreme observation values is calculated and compared with 5.  
#' The null hypothesis is equal variance.
#' @param x the one sample vector
#' @param y the other sample vector
#' @return  the value 1(reject H0) or 0 (do not reject H0)
#' @examples
#' \dontrun{
#' x<-rnorm(20,0,2)
#' y<-rnorm(20,1,4)
#' count5test(x,y)
#' }
#' @export
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

#' @title  F test of equal variance.
#' @description Under the null hypothesis (H0), the variance ratio follows the F-distribution.
#' Large values will reject the null hypothesis.
#' @param x the one sample vector
#' @param y the other sample vector
#' @param level Significant level
#' @param n Sample size
#' @importFrom stats qf var
#' @return  the value 1(reject H0) or 0 (do not reject H0)
#' @examples
#' \dontrun{
#' x<-rnorm(20,0,2)
#' y<-rnorm(20,1,4)
#' ftest(x,y)
#' }
#' @export
ftest<-function(x,y,level=0.055,n=20){
  X <- x - mean(x)
  Y <- y - mean(y)
  t<-1
  # return 1 (reject) or 0 (do not reject H0)
  if(var(X)/var(Y)<qf(1-level/2,n-1,n-1)&&var(X)/var(Y)>qf(level/2,n-1,n-1)){t<-0}
  return(t)
}

#' @title  Find the maximum number of extreme observations
#' @description Count the number of extreme observations in each sample that 
#' are not in the other sample. Calculate the number of maximum extreme points.
#' @param x the one sample vector
#' @param y the other sample vector
#' @return  the maximum number of extreme observations from two samples
#' @examples
#' \dontrun{
#' x<-1:10
#' y<-3:12
#' maxout(x,y)
#' }
#' @export
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}

#' @title  Permutation test for equal variance for two samples.
#' @description A permutation test for equal variance based on the maximum number
#'  of extreme points that applies when sample sizes are not necessarily equal. 
#'  Calculate original test statistic and permutation test statistics, and obtain the 
#'  permutation test p-value defined as proportion of permutation test statistics 
#'  larger than the original statistic.
#' @param x the one sample vector
#' @param y the other sample vector
#' @return  permutation test p-value
#' @examples
#' \dontrun{
#' set.seed(11)
#' b<-numeric(1000)
#' for(i in 1:1000){
#' x<--rnorm(15,0,1)
#' y<-rnorm(20,1,2)
#' b[i]<-eqtest.permu(x,y)<=0.05
#' }
#' print(mean(b))
#' }
#' @export
eqtest.permu<-function(x,y){
  R<-999
  z<-c(x,y)
  n<-length(x);N<-length(z)
  reps<-numeric(R)
  m0<-maxout(x,y)
  
  for(i in 1:R){
    #generate indices k for the first sample
    k <- sample(N, size = n, replace = FALSE)
    x1 <- z[k]
    y1 <- z[-k] #complement of x1
    reps[i] <- maxout(x1,y1)
  }
  phat<-mean(c(m0,reps)>=m0)
  return(phat)
}

#' @title  Nearest Neighbor Statistic for Equal Distribution Test
#' @description Find the nearest k neighbors by function "nn2".According to 
#' whether each sample point in the sample 1 and its corresponding k neighbors
#' are from the same sample, assign 1 or 0 to points of sample 1 and take the 
#' average as a statistic.
#' @param z A data set composed of two multidimensional samples
#' @param ix index used during bootstrap
#' @param sizes A vector of two sample sizes
#' @param k the number of chosen neighbors
#' @return  nearest neighbor test statistic
#' @examples
#' \dontrun{
#'   n1<-10;n2<-20;n<-n1+n2
#'   sizes<-c(n1,n2);k=3
#'   x<-rnorm(10,0,1)
#'   y<-rnorm(20,0,2)
#'   z <- c(x, y)
#'   ix<-1:n
#'   tn<-Tn(z,ix,sizes,k)
#' }
#' @export
Tn <- function(z, ix, sizes,k) {# k denote the number of chosen neighbors
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);# we consider the case that z is 1-dimension
  z <- z[ix, ]
  NN <- nn2(data=z, k=k+1) # the first column is itself!
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < (n1 +0.5)); i2 <- sum(block2 > (n1+.5))
  (i1 + i2) / (k * n)
}

#' @title  Equal Distribution Test by Nearest Neighbor Method
#' @description Use Bootstrap permutation method to obtain the sample of statistics
#'  and the p value for equal distribution test
#' @param z A data set composed of two multidimensional samples
#' @param sizes A vector of two sample sizes
#' @param k the number of chosen neighbors
#' @param R number of bootstrap
#' @return  original nearest neighbor statistic and equal distribution test p-value
#' @importFrom RANN nn2
#' @importFrom boot boot
#' @examples
#' \dontrun{
#' n1<-10;n2<-20;n<-n1+n2
#' sizes<-c(n1,n2);k=3
#' m <- 500;
#' p.values<-numeric[m]
#' p.values <- matrix(NA,m,3)
#' for(i in 1:m){
#' x<-rnorm(10,0,1)
#' y<-rnorm(20,1,2)
#' z <- c(x, y)
#' p.values[i] <- eqdist.nn(z,sizes,k)$p.value
#' alpha <- 0.1;
#' pow <- mean(p.values<alpha)
#' }
#' }
#' @export
eqdist.nn <- function(z,sizes,k,R){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

#' @title Standard Laplace Distribution Density
#' @description This function is defined to make it easier to write
#' @param x a real number
#' @return  standard Laplace distribution density at x
#' @examples
#' \dontrun{
#' lap(0)
#' lap(0.5)
#' }
#' @export
lap<-function(x){
  return(exp(-abs(x))/2)
}

#' @title  A random walk Metropolis Sampler for Generating the standard Laplace distribution
#' @description Generate  a sample that obeys standard Laplace distribution by 
#' Random walk Metropolis, which uses proposal distribution normal (X_t,sigma^2). 
#' The data generation process records the acceptance rate.
#' @param sigma the standard deviation chosen for a normal distribution
#' @param x0 the initial value of random walk
#' @param N size of generated sample
#' @return  A sample that obeys standard Laplace distribution and acceptance rates.
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' rw1 <- rw.Metropolis(1, 10, 5000)
#' plot(rw1,type="l",xlab="sigma=1",ylab="X", ylim=range(rw1))
#' refline <- c(log(0.05),-log(0.05))
#' abline(h=refline)
#' }
#' @export
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  u <- runif(N)
  k <- 0
  
  y <- rnorm(1, x0, sigma)
  if (u[1] <= (lap(y) / lap(x0)))
    x[1] <- y else {
      x[1] <- x0
      k <- k + 1
    }
  
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (lap(y) / lap(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=(1-k/N)))
}

#' @title Calculate the Gelman-Rubin statistic
#' @description The scalar summary statistic psi_ij is the mean of the ith chain up to time j.
#' The convergence rate is monitored by calculating Gelman-Rubin statistics of 
#' multiple data chains at different times.
#' @param psi A matrix of scalar summary statistic with n rows and k columns
#' @return  A sample that obeys standard Laplace distribution and acceptance rates.
#' @examples
#' \dontrun{
#' x0<-c(1,5,10)
#' for(i in 1: 3){
#' rw[i,] <- rw.Metropolis(1, x0[i], 1000)#generate random walk Metropolis sample
#' }
#' psi <- t(apply(X, 1, cumsum))
#' for (i in 1:nrow(psi))
#' psi[i,] <- psi[i,] / (1:ncol(psi))
#' print(Gelman.Rubin(psi))
#' }
#' @export
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

#' @title Observed data likelihood
#' @description As for the A-B-O blood type problem, observed data include nA, nB and nOO. 
#' Under these given conditions, the likelihood function of frequency p and q (and r=1-p-q) 
#' is shown.
#' @param theta a Gene frequency vector (p,q)
#' @return  Log of Likelihood of frequencies under observed data
#' @examples
#' \dontrun{
#' loglike_obs(c(0.3,0.3))
#' }
#' @export
loglike_obs<-function(theta){
  nA<-444;nB<-132;nOO<-361;nAB<-63
  p<-theta[1];q<-theta[2]
  r<-1-p-q
  t<-nA*log(p^2+2*p*r)+nB*log(q^2+2*q*r)+nAB*log(2*p*q)+nOO*log(r^2)
  return(t)
}

#' @title Epanechnikov kernel function
#' @description As shown in the expression. It is also a density function.
#' @param x real number
#' @return  Epanechnikov density(value) at x
#' @examples
#' \dontrun{
#' epan(0)
#' }
epan<-function(x){
  (3/4)*(1-x^2)*(abs(x)<=1)
}

#' @title Kernel estimation of binary joint density
#' @description Multiplicative kernel and spherical-symmetric kernel can be selected 
#' to estimate the joint pdf of two variables according to their binary sample matrix.
#' @param x the x-coordinate of point at which the joint density is to be obtained
#' @param y the y-coordinate of point at which the joint density is to be obtained
#' @param hx the bandwidth on the x-component
#' @param hy the bandwidth on the y-component
#' @param mat sample data
#' @param method the kernel we choose to estimate the joint pdf
#' @return  the joint density at point (x,y)
#' @examples
#' \dontrun{
#' attach(faithful)
#' esti_multi(x=2,y=56,mat=mat1,method = "multi")
#' xlab<-seq(0,6,0.1)
#' ylab<-seq(40,100,1)
#' z_density<-matrix(0,length(xlab),length(ylab))
#' for(i in 1:length(xlab)){
#' for(j in 1:length(ylab)){
#'  z_density[i,j]<-esti_multi(x=xlab[i],y=ylab[j],mat=mat1,method="multi")
#'  }}
#'  persp(xlab, ylab,z_density,phi = 30, theta = 40, col ="blue", border = 0, 
#'   zlim = c(0, 0.08),zlab = "Multiplicative KDE")
#' }
#' @export
esti_multi<-function(x,y,hx=0.5,hy=0.5,mat,method=c("multi","sphe")){
  n<-nrow(mat)
  temp<-numeric(n)
  if(method=="multi"){
    for(i in 1:n){
      temp[i]<-epan((x-mat[i,1])/hx)*epan((y-mat[i,2])/hy)/(hx*hy)
    }
    return(mean(temp))
  }
  if(method=="sphe"){
    for(i in 1:n){
      u<-(x-mat[i,1])/hx
      v<-(y-mat[i,2])/hy
      temp[i]<-(1-(u^2+v^2))*((u^2+v^2)<=1)/(pi/2)/(hx*hy)
    }
    return(mean(temp))
  }
}

#' @title Gibbs sampler for Bivariate distribution
#' @description In each iteration, samples at specific component locations are 
#' taken one by one according to the conditional distribution under other component samples.
#' In order to obtain the sanple from the goal distribution, discard the previous 
#' part of the unstable sample
#' @param n The desired sample size
#' @param x1 initial value of one variable
#' @param x2 initial value of one variable
#' @param plot whether to generate scatter diagram of sample
#' @param paras a vector which includes means, standard deviations and the correlation coefficient
#' @return  A list that includes sample chain and mean, covariance matrix and correlation coefficient
#' @importFrom stats cor cov
#' @examples
#' \dontrun{
#' out1<-Gibbs_samp(3000,0,0)
#' out1
#' out2<-Gibbs_samp(1000,0,0,plot="F")
#' out2
#' }
#' @export
Gibbs_samp<-function(n=3000,x1,x2,plot=c("T","F"),paras=c(1,2,1,2,0.5)){
  mu1<-paras[1];mu2<-paras[2];sigma1<-paras[3];sigma2<-paras[4];ro<-paras[5]
  m<-1000
  mat<-matrix(0,n+m+1,2)
  mat[1,]<-c(x1,x2)
  for(i in 2:n+1){
    x2<-mat[i-1,2]
    x1<-rnorm(1,mu1+ro*sigma1/sigma2*(x2-mu2),sqrt(1-ro^2)*sigma1)
    mat[i,1]<-x1
    x2<-rnorm(1,mu2+ro*sigma2/sigma1*(x1-mu1),sqrt(1-ro^2)*sigma2)
    mat[i,2]<-x2
  }
  mat1<-mat[-(1:(m+1)),]
  mu<-apply(mat1, 2,mean)
  sigma<-cov(mat1)
  corr<-cor(mat1)[1,2]
  if(plot=="T"){
    plot(mat1[,1],mat1[,2])
  }
  return(list(chain=mat1,mu=mu,sigma=sigma,corr=corr))
}