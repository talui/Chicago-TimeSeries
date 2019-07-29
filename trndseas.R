trndseas=function(y,seas,lam,degtrnd){
  
# requires the R-package 'pracma'
  
# fits  a trend plus seasonal for the "best" Box-Cox 
# transformation.
  
# input: y, observed series; seas, seasons
  
# input: lam, the grid of Box-Cox transformations (lambda values)
  
# input: degtrnd, degree of the polynomial trend, if
# degtrnd=0, then the fitted trend is constant.
  
# output:  coef, regression coefficients - the
# first degtrnd+1 values for the trend part and the
# rest associated with the seasonals
  
# output: fit, fitted y-values; res, residuals,
  
# output: trend, fitted trend; season, fitted seasonals
  
# output: rsq, adjusted r-square values for different lambda in the
  
# output: lamopt, the value of lambda (among those supplied 
# in the vector lam) at which r-square is maximum.

m=length(lam)
n=length(y)

# Part of design matrix for estimating trend
if(degtrnd>0) {
   tm=seq(1/n,1,by=1/n)
   x1=poly(tm,degree=degtrnd,raw=TRUE)
   x1=cbind(rep(1,n),x1)
   } else {
    x1=as.matrix(rep(1,n),ncol=1)
   }

# Part of design matrix for estimating seasonality
x2=NULL
if(seas>1){
sn=rep(1:seas,length.out=n)
x2=factor(sn,levels=unique(sn),ordered=TRUE)
x2=model.matrix(~x2-1)
m2=ncol(x2)
m21=m2-1
x2=x2[,1:m21]-matrix(rep(x2[,m2],m21),ncol=m21,nrow=nrow(x2),byrow=F)
}

x=cbind(x1,x2)  # design matrix

xx=t(x)%*%x
rsq=rep(1,m)
m1=ncol(x1)     #degtrnd+1
m11=m1+1
mx=ncol(x)      # degtrnd+1+seas-1

for(i in 1:m) {
  if (lam[i]==0) {
    yt=log(y)
  } else {
    yt=y^lam[i]
   }
  xy=t(x)%*%yt
  coef=solve(xx,xy)
  fit=x%*%coef
  res=yt-fit
  ssto=(n-1)*var(yt)
  sse=t(res)%*%res
  rsq[i]=1-((n-1)/(n-mx))*sse/ssto
  }

  ii=which.max(rsq)
  lamopt=lam[ii]
  if (lamopt==0) {
    yt=log(y)
  } else {
    yt=y^lamopt
   }
  xy=t(x)%*%yt
  coef=solve(xx,xy)
  fit=x%*%coef
  trnd=x1%*%coef[1:m1]
  season=NULL
  if(seas>1){
  season=c(coef[m11:mx],-sum(coef[m11:mx]))
  }
  res=yt-fit

  result=list(coef=coef,fitted=fit,trend=trnd,residual=res,season=season,rsq=rsq,lamopt=lamopt)
  return(result)
}
  