## -----------------------------------------------------------------------------
set.seed(1)
x<-rnorm(20)
y<-rt(20,5)
layout(matrix(1:2,1,2))
qqnorm(x,xlim=c(-2,2),ylim=c(-2,2))
mtext("random samples from N(0,1)",side=3,line=0.5,at=1.3,cex=0.9)
qqnorm(y,xlim=c(-2,2),ylim=c(-2,2))
mtext("random samples from t(5)",side=3,line=0.5,at=1.3,cex=0.9)

## -----------------------------------------------------------------------------
library(knitr)
set.seed(2)
z1<-round(runif(10),3);z2<-round(rnorm(10),3);z3<-round(rexp(10),3);
m<-matrix(c(1:10,z1,z2,z3),10,4)
m<-as.data.frame(m)
names(m)<-c("No","uniform","normal","exponential")
##knitr::kable(m, format="html")
##print(xtable(m), type="html")
kable(m, format = "html",caption = "Different Distributions")

## -----------------------------------------------------------------------------
set.seed(11)
a<-2;b<-2
u1<-runif(500)
x1<-b*(1-u1)^(-1/a)## generate samples from Preto(2,2)
hist(x1,prob=T,main="Density histogram of the sample  with the Pareto(2, 2) density")
x2<-seq(b,max(x1),0.1)
lines(x2,a*b^a*x2^(-a-1))## theoretical density curve

## -----------------------------------------------------------------------------
set.seed(22)
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
x1<-ranfe(n=10000)
hist(x1,prob=T)
##x2<-seq(-1,1,0.01)
##lines(x2,3/4*(1-x2^2))## theoretical density curve

## -----------------------------------------------------------------------------
set.seed(33)
n<-10000;r<-4;beta<-2
x1<-rgamma(n,r,beta)
x2<-rep(0,n)
for(i in 1: n){
  x2[i]<-rexp(1,x1[i])
}
hist(x2,prob=T,main="the empirical and theoretical (Pareto) distributions ")
x3<-seq(0,15,0.01)
lines(x3,r*beta^r/(beta+x3)^(1+r))

## -----------------------------------------------------------------------------
set.seed(11)
x<-runif(1e4,min=0,max=pi/3)
theta.hat<-mean(sin(x)*pi/3)
theta<-cos(0)-cos(pi/3)
print(theta.hat);print(theta)

## -----------------------------------------------------------------------------
set.seed(22)
m<-1e4
u1<-runif(m/2,min=0,max=1);u2<-runif(m/2,min=0,max=1);u3<-1-u1
##simple MC
x1<-c(u1,u2)
theta.hat1<-mean(exp(x1))
##antithetic variate approach
x2<-c(u1,u3)
theta.hat2<-mean(exp(x2))
print(theta.hat1);print(theta.hat2)

## -----------------------------------------------------------------------------
var1<-var(exp(x1))/m
var2<-var((exp(u1)+exp(1-u1))/2)/(m/2)
redu<-1-var2/var1

## -----------------------------------------------------------------------------
x1<-seq(1,6,0.01)
g<-function(t){
  return(t^2/sqrt(2*pi)*exp(-t^2/2))
}
plot(x1,g(x1),type="l",ylim=c(0,0.6),xlab="x",ylab="g(x)",main="Plot of g, f1 and f2")
f1<-function(t){
  return(exp(1-t))
}
lines(x1,f1(x1),col="blue")
f2<-function(t){
  return(2/sqrt(2*pi)*exp(-(t-1)^2/2))
}
lines(x1,f2(x1),col="red")
legend("topright", legend = c("g","f1","f2"),col=c("black","blue","red"),lty=1)

m <- 10000
set.seed(111)
theta.hat <- se <- numeric(2)
## f1=Exp distribution, inverse transform method
u<-runif(m)
x<-1-log(u)
fg<-g(x)/f1(x)
theta.hat[1]<-mean(fg)
se[1]<-round(var(fg),4)
## f2 is half of normal distribution
v<-rnorm(m,1,1)
y<-numeric(m)
for(i in 1:m){
  y[i]<-v[i]*(v[i]>=1)+(2-v[i])*(v[i]<1)
}
fg<-g(y)/f2(y)
theta.hat[2]<-mean(fg)
se[2]<-round(var(fg),4)
print(theta.hat)
print(se)

## -----------------------------------------------------------------------------
g<-function(x){
  exp(-x)/(1+x^2)
}
theta_I<-0.5257801
se_I<- 0.0970314
set.seed(222)
m<-10000
theta.hat<-se.hat<-numeric(5)
for(j in 1:5){
  u<-runif(m/5,(j-1)/5,j/5)
  temp<- -log(exp(-(j-1)/5)-(u*(exp(-(j-1)/5)-exp(-j/5))))
  theta.hat[j]<-mean(g(temp)/(exp(-temp)/(exp(-(j-1)/5)-exp(-j/5))))
  se.hat[j]<-var(g(temp)/(exp(-temp)/(exp(-(j-1)/5)-exp(-j/5))))*(m/5-1)
}
##print(theta.hat)
theta_hat<-sum(theta.hat)
se_hat<-sum(se.hat)
print(c(theta_hat,sqrt(se_hat)))
print(c(theta_I,se_I))

## -----------------------------------------------------------------------------
mu<-0;sigma2<-1^2
m<-1e4;n<-20
set.seed(333)
ybar<-sigmahat<-numeric(m)
for(i in 1:m){
  y<-rnorm(n,mu,sigma2)
  ybar[i]<-mean(y)
  sigmahat[i]<-sd(y)
}
y_bar<-mean(ybar)
sigma_hat<-mean(sigmahat)
##print(c(y_bar,sigma_hat))

## -----------------------------------------------------------------------------
cp_1<-0
for(i in 1:m){
  if(mu>ybar[i]-qt(0.975,n-1)*sigmahat[i]/sqrt(n)&&mu<ybar[i]+qt(0.975,n-1)*sigmahat[i]/sqrt(n)) cp_1<-cp_1+1
}
cp_1<-cp_1/m
cp_1

## -----------------------------------------------------------------------------
m<-1e4;n<-20
mean<-2##mean is equal to df
set.seed(444)
ybar<-sigmahat<-numeric(m)
for(i in 1:m){
  y<-rchisq(n,2)
  ybar[i]<-mean(y)
  sigmahat[i]<-sd(y)
}
y_bar<-mean(ybar)
sigma_hat<-mean(sigmahat)
##print(c(y_bar,sigma_hat))

## -----------------------------------------------------------------------------
cp_2<-0
for(i in 1:m){
  if(mean>ybar[i]-qt(0.975,n-1)*sigmahat[i]/sqrt(n)&&mean<ybar[i]+qt(0.975,n-1)*sigmahat[i]/sqrt(n)) cp_2<-cp_2+1
}
cp_2<-cp_2/m
cp_2

## -----------------------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
set.seed(11)
level<-0.05
n <- 50
m <- 2500
para1<- seq(0.5,10,0.5)
para2<- 1:40
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
#plot power vs epsilon
pwr1<-power(para1)
pwr2<-power(para2)
plot(para1, pwr1, type = "b",xlab = bquote(alpha), ylim = c(0,0.1))
abline(h = .05, lty = 3)
plot(para2, pwr2, type = "b",xlab = bquote(alpha), ylim = c(0,0.1))
abline(h = .05, lty = 3)


## -----------------------------------------------------------------------------
set.seed(12)
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
nu<-seq(0.5,20,0.5)
pwr<-power_stu(nu)
plot(nu, pwr, type = "b",xlab = bquote(nu), ylim = c(0,1))
abline(h = .05, lty = 3)
#print(power_stu(50:60))

## -----------------------------------------------------------------------------
set.seed(21)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
ftest<-function(x,y,level=0.055,n=20){
  X <- x - mean(x)
  Y <- y - mean(y)
  t<-1
  # return 1 (reject) or 0 (do not reject H0)
  if(var(X)/var(Y)<qf(1-level/2,n-1,n-1)&&var(X)/var(Y)>qf(level/2,n-1,n-1)){t<-0}
  return(t)
}
n1 <- n2 <- 20
mu1 <- mu2 <- 0
sigma1 <-1
sigma2 <-1.5
m <- 10000
counttest<-numeric(m)
Ftest<-numeric(m)
for(i in 1:m){
  x <- rnorm(20, 0, sigma1)
  y <- rnorm(20, 0, sigma2)
  counttest[i]<-count5test(x, y)
  Ftest[i]<-ftest(x,y)
}
power_count <- mean(counttest)
power_F<-mean(Ftest)
print(power_count)
print(power_F)

## -----------------------------------------------------------------------------
set.seed(21)
#change the sample size
n<-c(10,20,30,50,100,500)
N<-length(n)
power_count<-numeric(N)
power_F<-numeric(N)
for(j in 1: N){
  counttest<-numeric(m)
  Ftest<-numeric(m)
  for(i in 1:m){
  x <- rnorm(n[j], 0, sigma1)
  y <- rnorm(n[j], 0, sigma2)
  counttest[i]<-count5test(x, y)
  Ftest[i]<-ftest(x,y,level=0.055,n=n[j])
  }
  power_count[j] <- mean(counttest)
  power_F[j]<-mean(Ftest)
}
print(power_count)
print(power_F)
plot(n,power_count,type = "b",xlab = "size of sample", ylim = c(0,1))
lines(n,power_F,lty=3,type = "b")
legend("bottomright", legend = c("CountFive test","F test"),lty=c(1,3),lwd=2)

## -----------------------------------------------------------------------------
##Example 6.8
set.seed(33)
d<-1
level<-0.05
n <- c(10, 20, 30,50, 100, 500) #sample sizes
cv <- qchisq(1-level,d*(d+1)*(d+2)/6)

p.reject <- numeric(length(n)) #to store sim. results
m <- 10000 #num. repl. each sim.
#sigma<-matrix(c(1,0.5,0.5,1),2,2)
#mean<-numeric(d)
for (i in 1:length(n)) {
  sktests <- numeric(m) #test decisions
  for (j in 1:m) {
  x <- rnorm(n[i])
  #x<-mvrnorm(n[i],mean,sigma)
  #test decision is 1 (reject) or 0
  sigma_hat<-var(x)
  sigma_inv<-1/var(x)
  mean_hat<-mean(x)
  b1d<-(t(x-mean_hat))^3%*%(x-mean_hat)^3*((sigma_inv)^3)/(n[i]^2)
  sktests[j] <- as.integer(b1d*n[i]/6>=cv)
 }
p.reject[i] <- mean(sktests) #proportion rejected
}
print(p.reject)

## ----eval=FALSE---------------------------------------------------------------
#  #Example 6.10
#  set.seed(33)
#  d<-1
#  alpha <- .1
#  n <- 30
#  m <- 2500
#  epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
#  N <- length(epsilon)
#  pwr <- numeric(N)
#  #critical value for the skewness test
#  cv <- qchisq(1-alpha,d*(d+1)*(d+2)/6)
#  for (j in 1:N) { #for each epsilon
#  e <- epsilon[j]
#  sktests <- numeric(m)
#  for (i in 1:m) { #for each replicate
#    sigma <- sample(c(1, 10), replace = TRUE,
#    size = n, prob = c(1-e, e))
#    x <- rnorm(n, 0, sigma)
#  
#    sigma_hat<-var(x)
#    sigma_inv<-1/var(x)
#    mean_hat<-mean(x)
#    b1d<-0
#    b1d<-(t(x-mean_hat))^3%*%(x-mean_hat)^3*((sigma_inv)^3)/(n^2)
#    sktests[i] <- as.integer(b1d*n/6>= cv)
#  }
#   pwr[j] <- mean(sktests)
#  }
#  print(pwr)
#  #plot power vs epsilon
#  plot(epsilon, pwr, type = "b",
#  xlab = bquote(epsilon), ylim = c(0,1))
#  abline(h = .1, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap) #for the law data

## -----------------------------------------------------------------------------
n<-nrow(law)
LSAT<-law$LSAT
GPA<-law$GPA
theta.hat<-cor(LSAT,GPA)
#compute the jackknife replicates, leave-one-out estimates
theta.jack <- numeric(n)
for (i in 1:n)
{theta.jack[i] <- cor(LSAT[-i],GPA[-i])}
#estimate bias
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
print(bias)
#estimate stabdard error
se <- sqrt((n-1) *mean((theta.jack - mean(theta.jack))^2))
print(se)

## -----------------------------------------------------------------------------
library(boot)
#mode(aircondit)
n<-nrow(aircondit)
time<-aircondit[1:n,]

## -----------------------------------------------------------------------------
set.seed(22)
boot.obj<-boot(time,statistic = function(t, i){mean(t[i])},1000)
print(boot.obj)
#boot.obj$t0 is original mean
ci<-boot.ci(boot.obj, type = c("basic", "norm", "perc"))
print(ci)

## -----------------------------------------------------------------------------
hist(boot.obj$t,xlab = "Bootstrap sample",main="Histogram of Bootstrap sample",breaks = seq(0,270,10))

## -----------------------------------------------------------------------------
library(bootstrap)
lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat)
n <- nrow(scor) # number of rows (data size)
# Jackknife
theta_j <- rep(0, n)
for (i in 1:n) {
x <- scor [-i,]
lambda <- eigen(cov(x))$values
theta_j[i] <- lambda[1] / sum(lambda)
# the i-th entry of theta_j is the i-th "leave-one-out" estimation of theta
}
bias_jack <- (n - 1) * (mean(theta_j) - theta_hat)
# the estimated bias of theta_hat, using jackknife
se_jack <- (n - 1) * sqrt(var(theta_j) / n)
# the estimated se of theta_hat, using jackknife
print(c(bias_jack,se_jack))

## ----eval=FALSE---------------------------------------------------------------
#  data<-read.table("./data.txt",quote = "")
#  data<-t(as.matrix(data))
#  magnetic<-data[,1]
#  chemical<-data[,2]
#  r<-nrow(data)
#  c<-ncol(data)
#  #library(DAAG); attach(ironslag)
#  #a <- seq(10, 40, .1) #sequence for plotting fits
#  e1 <- e2 <- e3 <- e4 <- matrix(0,r,r)
#  k<-r*(r-1)/2
#  #index i<j
#  for(i in 1:(r-1)){
#    for(j in (i+1):r){
#      y <- magnetic[-c(i,j)]
#      x <- chemical[-c(i,j)]
#  
#      J1 <- lm(y ~ x)
#      yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i,j)]
#      e1[i,j] <- mean((magnetic[c(i,j)] - yhat1)^2)
#  
#      J2 <- lm(y ~ x + I(x^2))
#      yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(i,j)] + J2$coef[3] * chemical[c(i,j)]^2
#      e2[i,j] <- mean((magnetic[c(i,j)] - yhat2)^2)
#  
#      J3 <- lm(log(y) ~ x)
#      logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
#      yhat3 <- exp(logyhat3)
#      e3[i,j] <- mean((magnetic[c(i,j)] - yhat3)^2)
#  
#      J4 <- lm(log(y) ~ log(x))
#      logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
#      yhat4 <- exp(logyhat4)
#      e4[i,j] <- mean((magnetic[c(i,j)] - yhat4)^2)
#    }
#  }
#  sum_e1<-sum(e1)/k
#  sum_e2<-sum(e2)/k
#  sum_e3<-sum(e3)/k
#  sum_e4<-sum(e4)/k

## ----eval=FALSE---------------------------------------------------------------
#  print(c(sum_e1,sum_e2,sum_e3,sum_e4))

## ----eval=FALSE---------------------------------------------------------------
#  J2

## ----eval=FALSE---------------------------------------------------------------
#  #par(mfrow = c(2, 2)) #layout for graphs
#  plot(J2$fit, J2$res) #residuals vs fitted values
#  abline(0, 0) #reference line
#  qqnorm(J2$res) #normal probability plot
#  qqline(J2$res) #reference line
#  #par(mfrow = c(1, 1)) #restore display

## -----------------------------------------------------------------------------
#count the maximum of extreme points
set.seed(11)
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}
#define function of replicate
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

## ----eval=FALSE---------------------------------------------------------------
#  # test the performance of eqtest.permu
#  set.seed(11)
#  a<-numeric(1000)
#  for(i in 1:1000){
#    x<--rnorm(15,0,4)
#    y<-rnorm(20,1,4)
#    a[i]<-eqtest.permu(x,y)<=0.05
#  }
#  #reject rate of 1000 replicates
#  print(mean(a))

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(11)
#  b<-numeric(1000)
#  for(i in 1:1000){
#    x<--rnorm(15,0,1)
#    y<-rnorm(20,1,2)
#    b[i]<-eqtest.permu(x,y)<=0.05
#  }
#  #reject rate of 1000 replicates
#  print(mean(b))

## -----------------------------------------------------------------------------
# Test whether distributions are same
# Different methods include NN, energy, and ball methods
# Generate experiment data z

set.seed(111)
#NN
library(RANN)
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
#"boot" function to carry out the permutation
library(boot)
R<-499
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
#use function eqdist.nn

#Energy method
library(energy)
#use function eqdist.etest

#Ball method
library(Ball)
# use function bd.test

## ----eval=FALSE---------------------------------------------------------------
#  # Now we conduct the power comparison
#  n1<-10;n2<-20;n<-n1+n2
#  sizes<-c(n1,n2);k=3
#  m <- 500;
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x<-rnorm(10,0,1)
#    y<-rnorm(20,0,2)
#    z <- c(x, y)
#    p.values[i,1] <- eqdist.nn(z,sizes,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=sizes,R=499)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,R=499,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(111)
#  n1<-10;n2<-20;n<-n1+n2
#  sizes<-c(n1,n2);k=3
#  m <- 500;
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x<-rnorm(10,0,1)
#    y<-rnorm(20,1,2)
#    z <- c(x, y)
#    p.values[i,1] <- eqdist.nn(z,sizes,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=sizes,R=499)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,R=499,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(111)
#  n1<-20;n2<-20;n<-n1+n2
#  sizes<-c(n1,n2);k=3
#  m <- 500;
#  e<-0.5#epsilon
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x<-rt(20,1)
#    sigma <- sample(c(1, 2), replace = TRUE, size = 20, prob = c(1-e, e))
#    y <- rnorm(20, 0, sigma)
#    z <- c(x, y)
#    p.values[i,1] <- eqdist.nn(z,sizes,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=sizes,R=499)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,R=499,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## -----------------------------------------------------------------------------
set.seed(11)
#define standard Laplace density
lap<-function(x){
  return(exp(-abs(x))/2)
}
#define the random walk sampler
#The function returns the sample x and acceptance rates.
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
N <- 5000
sigma <- c(.05, 1, 4, 16)
x0 <- 10
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
print(c(rw1$k, rw2$k, rw3$k, rw4$k))

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2)) #display 4 graphs together
refline <- c(log(0.05),-log(0.05))
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
for (j in 1:4) {
  plot(rw[,j],type="l",xlab=bquote(sigma == .(round(sigma[j],3))),ylab="X", ylim=range(rw[,j]))
  abline(h=refline)
}
#par(mfrow=c(1,1)) #reset to default

## -----------------------------------------------------------------------------
set.seed(22)
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

## -----------------------------------------------------------------------------
normal.chain <- function(sigma, N, X1) {
  #generates a Metropolis chain for Normal(0,1)
  #with Normal(X[t], sigma) proposal distribution
  #and starting value X1
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  for (i in 2:N) {
  xt <- x[i-1]
  y <- rnorm(1, xt, sigma) #candidate point
  r1 <- lap(y) * dnorm(xt, y, sigma)
  r2 <- lap(xt) * dnorm(y, xt, sigma)
  r <- r1 / r2
  if (u[i] <= r) x[i] <- y else
    x[i] <- xt
  }
  return(x)
}

## ----eval=FALSE---------------------------------------------------------------
#  sigma <-1 #parameter of proposal distribution
#  k <- 4 #number of chains to generate
#  n <- 15000 #length of chains
#  b <- 1000 #burn-in length
#  x0 <- c(-10, -5, 5, 10)#different initial values
#  
#  #generate the MC chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k){
#    X[i, ] <- normal.chain(sigma, n, x0[i])
#  }
#  
#  #compute diagnostic statistics, that is mean of the first j elements in the k chain
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#  psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  #plot psi for the four chains
#  #par(mfrow=c(2,2))
#  for (i in 1:k){
#    plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
#  }
#  #par(mfrow=c(1,1)) #restore default
#  
#  #plot the sequence of R-hat statistics
#  rhat <- rep(0, n)
#  for (j in (b+1):n)
#  rhat[j] <- Gelman.Rubin(psi[,1:j])
#  plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)
#  

## -----------------------------------------------------------------------------
K<-c(4:25,100,500,1000)
n<-length(K)
#a>0 and a^2<k
pre<-numeric(n)
for(i in 1:n){
  k<-K[i]
  f<-function(a){
    pt(sqrt(a^2*k/(k+1-a^2)),k)-pt(sqrt(a^2*(k-1)/(k-a^2)),k-1)
  }
  pre[i]<-f(0.1)*f(sqrt(k)-0.1)
  if(pre[i]>=0) print(K[i])
}

## -----------------------------------------------------------------------------
Ak<-rep(NA,n)
for(i in 1:22){
  Ak[i]<-uniroot(f,c(0.1,sqrt(K[i])-0.1))$root
}
print(Ak)

## -----------------------------------------------------------------------------
for(i in 23:25){
  Ak[i]<-uniroot(f,c(1,2))$root
}
print(Ak)

## ----warning=FALSE------------------------------------------------------------
#赋初值
nA<-444;nB<-132;nOO<-361;nAB<-63
loglike_obs<-function(theta){
  p<-theta[1];q<-theta[2]
  r<-1-p-q
  t<-nA*log(p^2+2*p*r)+nB*log(q^2+2*q*r)+nAB*log(2*p*q)+nOO*log(r^2)
  return(t)
}

theta0<-c(0.33,0.33)
m<-10
record<-matrix(m*3,m,3)
for(i in 1:m){
#E-step
#这里应当求出nAO和nBO这两个缺失值的估计，它们是上述函数最大值时的p0和q0
p0<-theta0[1];q0<-theta0[2];r0<-1-sum(theta0)
nAO<-nA*(2*p0*r0)/(p0^2+2*p0*r0)
nBO<-nB*(2*q0*r0)/(q0^2+2*q0*r0)
loglike_all<-function(theta){
  p<-theta[1];q<-theta[2]
  r<-1-p-q
  t<-nA*log(p^2)+nAO*log(2*r/p)+nB*log(q^2)+nBO*log(2*r/q)*is.na(2*r/q)+nOO*log(r^2)+nAB*log(2*p*q)*is.na(2*p/q)
  return(-t)
}
#M-step
res1<-optim(c(0.33,0.33),loglike_all)#这里的初值是极小值函数optim的初值设定
theta1<-res1$par
p1<-theta1[1];q1<-theta1[2];r1<-1-p1-q1
loglike_obs<-nA*log( (p1)^2+2*p1*r1) + nB*log( (q1)^2+2*q1*r1) + nOO*log( (r1)^2) + nAB*log(2*p1*q1)
record[i,]<-c(theta1,loglike_obs)
theta0<-theta1
}
colnames(record)=c("p","q","logML_obs")
knitr::kable(record)

## -----------------------------------------------------------------------------
nA<-444;nB<-132;nOO<-361;nAB<-63
m<-10

record<-matrix(m*3,m,3)
#根据理论推导的（全数据）条件似然的极大值点表达式，写出下列递推函数
#并同时计算log max likelihood values(for observed data)
#自变量p0,q0，输出结果是p1,q1
recur<-function(p0,q0){
  r0<-1-p0-q0
  nAO<-nA*(2*p0*r0)/(p0^2+2*p0*r0)
  nBO<-nB*(2*q0*r0)/(q0^2+2*q0*r0)
  c1<-2*nA+nAB-nAO
  c2<-2*nB+nAB-nBO
  c3<-2*nOO+nAO+nBO
  p1<-c1/(c1+c2+c3)
  q1<-c2/(c1+c2+c3)
  r1<-1-p1-q1
  loglike_obs<-nA*log( (p1)^2+2*p1*r1) + nB*log( (q1)^2+2*q1*r1) + nOO*log( (r1)^2) + nAB*log(2*p1*q1)
  return(c(p1,q1,loglike_obs))
}
#下面写出对数最大似然的值（对于观测数据）

p0<-0.3;q0<-0.3
for(i in 1:m){
  record[i,]<-c(recur(p0,q0))
  p0<-record[i,1]
  q0<-record[i,2]
}
colnames(record)=c("p","q","logML_obs")
knitr::kable(record)

## ----warning=FALSE------------------------------------------------------------
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
#lapply function
lap_models<-lapply(formulas,lm)
lap_models

## -----------------------------------------------------------------------------
#loops
loop_models<-vector("list",length(formulas))
for(i in seq_along(formulas)){
  loop_models[[i]]<-lm(formulas[[i]])
}
loop_models

## -----------------------------------------------------------------------------
set.seed(1)
trials <- replicate(100, t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
pvalue<-sapply(trials,function(x) x[]$p.value)
print(pvalue)

## -----------------------------------------------------------------------------
f4=function(x,fun4,fun4_value){
    row4=length(x[,1]);col4=length(x[1,])##行列数
    i4=vapply(x,fun4,fun4_value)##对矩阵每个单独值，类似lapply函数单独求值
    print(matrix(i4,ncol=col4))##返回矩阵形式的值
  }
matrix(1:6,ncol=2)
f4(matrix(1:6,ncol=2),function(x){x^2},c(0))

## -----------------------------------------------------------------------------
library(Rcpp)
#metro function returns rate and chain
library(StatComp20056)
#sourceCpp("./metro.cpp")
sigma<-c(0.05,1,4,16)
x0<-5;N<-5000
a1<-as.numeric(metro(N,sigma[1],x0))
a2<-as.numeric(metro(N,sigma[2],x0))
a3<-as.numeric(metro(N,sigma[3],x0))
a4<-as.numeric(metro(N,sigma[4],x0))
#acceptance rates with different variances
print(c(a1[1],a2[1],a3[1],a4[1]))
#compare different variances
#use R plot
#par(mfrow=c(2,2)) #display 4 graphs together
refline <- c(log(0.05),-log(0.05))
a<- cbind(a1,a2,a3,a4)
for (j in 1:4) {
  plot(a[-1,j],type="l",xlab=bquote(sigma == .(round(sigma[j],3))),ylab="X", ylim=range(a[-1,j]))
  abline(h=refline)
}

## -----------------------------------------------------------------------------
quan<-(1:N)/N
theo1<-log(2*quan[1:(N/2)])
theo2<--rev(theo1)
theo<-c(theo1,theo2)
#par(mfrow=c(2,2))
for(j in 1:4){
  qqplot(a[-1,j],theo,xlab=bquote(sigma == .(round(sigma[j],3))),ylab = "Sample quantile")
  abline(coef = c(0,1))
}

## -----------------------------------------------------------------------------
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
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
# acceptance rates
#print(c(rw1$k, rw2$k, rw3$k, rw4$k))
rw<-cbind(rw1$x,rw2$x,rw3$x,rw4$x)
#par(mfrow=c(2,2))
#qqplot
for(j in 1:4){
  qqplot(rw[,j],theo,xlab=bquote(sigma == .(round(sigma[j],3))),ylab = "Sample quantile")
  abline(coef = c(0,1))
}

## -----------------------------------------------------------------------------
#par(mfrow=c(1,2))
qqplot(a[-1,2],rw[,2],xlab = "Sample quantile (Cpp function) with sigma=1",ylab = "Sample quantile (R function)")
abline(coef = c(0,1))
qqplot(a[-1,2],rw[,2],xlab = "Sample quantile (Cpp function) with sigma=4",ylab = "Sample quantile (R function)")
abline(coef = c(0,1))

## -----------------------------------------------------------------------------
library(microbenchmark)
ts1<-microbenchmark(rw_R=rw.Metropolis(sigma[1],x0,N),metro_C=metro(N,sigma[1],x0))
ts2<-microbenchmark(rw_R=rw.Metropolis(sigma[2],x0,N),metro_C=metro(N,sigma[2],x0))
ts3<-microbenchmark(rw_R=rw.Metropolis(sigma[3],x0,N),metro_C=metro(N,sigma[3],x0))
ts4<-microbenchmark(rw_R=rw.Metropolis(sigma[4],x0,N),metro_C=metro(N,sigma[4],x0))
summary(ts1)[,c(1,3,5,6)]
summary(ts2)[,c(1,3,5,6)]
summary(ts3)[,c(1,3,5,6)]
summary(ts4)[,c(1,3,5,6)]

