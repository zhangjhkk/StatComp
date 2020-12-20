## -----------------------------------------------------------------------------
library(kedd)
epan<-function(x){
  (3/4)*(1-x^2)*(abs(x)<=1)
}
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
      u<-(x-data[i,1])/hx
      v<-(y-data[i,2])/hy
      temp[i]<-(1-(u^2+v^2))*((u^2+v^2)<=1)/(pi/2)/(hx*hy)
    }
    return(mean(temp))
  }
}
attach(faithful)
mat1<-cbind(eruptions,waiting)
n<-nrow(data)
esti_multi(x=2,y=56,mat=mat1,method = "multi")
xlab<-seq(0,6,0.1)
ylab<-seq(40,100,1)
z_density<-matrix(0,length(xlab),length(ylab))
for(i in 1:length(xlab)){
for(j in 1:length(ylab)){
  z_density[i,j]<-esti_multi(xlab[i],ylab[j])
  }}
persp(xlab, ylab,z_density,phi = 30, theta = 40, col ="blue", border = 0, 
zlim = c(0, 0.08),zlab = "Multiplicative KDE")


## -----------------------------------------------------------------------------
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
out1<-Gibbs_samp(3000,0,0)
out1
out2<-Gibbs_samp(1000,0,0,plot="F")
out2

