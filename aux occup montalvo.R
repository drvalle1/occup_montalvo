tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#------------------------------------------------
sample.z=function(xmat.occ,betas,nloc,nspp,nrep,y,xmat.det,gammas,alpha.s){
  alpha.s.mat=matrix(alpha.s,nloc,nspp,byrow=T)
  media.occ=alpha.s.mat+xmat.occ%*%betas
  prob.occ=pnorm(media.occ)
  prob.det=matrix(1,nloc,nspp)
  for (i in 1:nrep){
    media.det=xmat.det[,,i]%*%gammas
    tmp=1-pnorm(media.det)
    prob.det=prob.det*tmp
  }
  tmp=prob.det*prob.occ
  prob.fim=tmp/(tmp+(1-prob.occ))
  
  #sample z
  z=matrix(1,nloc,nspp)
  max.y=apply(y,c(1,2),max)
  cond=max.y==0
  z[cond]=rbinom(sum(cond),size=1,prob=prob.fim[cond])
  z
}
#----------------------------------------------
sample.zstar=function(z,xmat.occ,betas,nloc,nspp,alpha.s){
  alpha.s.mat=matrix(alpha.s,nloc,nspp,byrow=T)
  media=alpha.s.mat+xmat.occ%*%betas
  lo.mat=matrix(ifelse(z==1,0,-1000),nloc,nspp)
  hi.mat=matrix(ifelse(z==0,0, 1000),nloc,nspp)
  zstar=tnorm(n=nloc*nspp,lo=lo.mat,hi=hi.mat,mu=media,sig=1)
  matrix(zstar,nloc,nspp)
}
#--------------------------------------------
sample.ystar=function(nrep,xmat.det,gammas,y,nnloc,nspp){
  ystar=y
  for (j in 1:nrep){
    media=xmat.det[,,j]%*%gammas
    lo.mat=matrix(ifelse(y[,,j]==1,0,-1000),nloc,nspp)
    hi.mat=matrix(ifelse(y[,,j]==0,0, 1000),nloc,nspp)
    ystar[,,j]=tnorm(n=nloc*nspp,lo=lo.mat,hi=hi.mat,mu=media,sig=1)
  }
  ystar
}
#--------------------------------------------
sample.alpha.s=function(m.alpha,tau2.alpha,xmat.occ,zstar,nspp,
                        w,betas,nloc){
  prec=nloc+(1/tau2.alpha)
  var1=1/prec
  
  err=colSums(zstar-xmat.occ%*%betas)
  pmedia=(1/tau2.alpha)*m.alpha+err
  rnorm(nspp,mean=var1*pmedia,sd=sqrt(var1))
}
#----------------------------------------
sample.m.alpha=function(nspp,alpha.s,tau2.alpha){
  prec=(nspp/tau2.alpha)+(1/100)
  var1=1/prec
  pmedia=sum(alpha.s)/tau2.alpha
  rnorm(1,mean=var1*pmedia,sd=sqrt(var1))
}
#----------------------------------------
sample.tau2.alpha=function(nspp,alpha.s,m.alpha,tau2.a,tau2.b){
  a1=(nspp/2)+tau2.a
  err2=sum((alpha.s-m.alpha)^2)
  b1=tau2.b+(err2/2)
  1/rgamma(1,a1,b1)
}
#----------------------------------------
sample.gammas=function(ystar,xmat.det,z,m.gamma,tau2.gamma,nparam.det,nspp,nrep){
  invTau=diag(1/tau2.gamma)
  gammas=matrix(NA,nparam.det,nspp)
  for (i in 1:nspp){
    cond=z[,i]==1
    for (j in 1:nrep){
      if (j==1) {
        xmat.det1=xmat.det[cond,,j]
        ystar1=ystar[cond,i,j]
      }
      if (j> 1) {
        xmat.det1=rbind(xmat.det1,xmat.det[cond,,j])
        ystar1=c(ystar1,ystar[cond,i,j])
      }
    }
    xtx=t(xmat.det1)%*%xmat.det1    
    var1=solve(xtx+invTau)
    pmedia1=t(xmat.det1)%*%ystar1+invTau%*%m.gamma
    gammas[,i]=rmvnorm(1,var1%*%pmedia1,var1)
  }
  gammas
}
#--------------------------------------------
sample.m.gamma=function(gammas,tau2.gamma,nspp,nparam.det){
  soma.gammas=rowSums(gammas)
  invTau=diag(1/tau2.gamma)
  prec=nspp*invTau + diag(rep(1/100,nparam.det))
  var1=solve(prec)
  pmedia=invTau%*%soma.gammas
  t(rmvnorm(1,mean=var1%*%pmedia,var1))
}
#--------------------------------------------
sample.tau2.gamma=function(gammas,m.gamma,nspp,tau2.a,tau2.b,nparam.det){
  a1=(nspp+2*tau2.a)/2
  m.gamma.mat=matrix(m.gamma,nparam.det,nspp)
  err2=rowSums((gammas-m.gamma.mat)^2)
  1/rgamma(nparam.det,a1,(err2/2)+tau2.b)
}
#----------------------------------------
sample.betas=function(m.betas,tau2.betas,xmat.occ,zstar,nspp,nparam.occ,
                      w,xtx.occ,alpha.s){
  betas=matrix(NA,nparam.occ,nspp)
  prec=1/tau2.betas
  
  #to speed things up
  prec1=xtx.occ+diag(prec,nparam.occ)
  var1=solve(prec1)
  
  for (i in 1:nspp){
    w1=w[i]
    pmedia=t(xmat.occ)%*%(zstar[,i]-alpha.s[i])+(1/tau2.betas)*m.betas[,w1]
    betas[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
  }
  betas
}
#----------------------------------------
sample.m.betas=function(w,betas,tau2.betas,nparam.occ,ngr){
  m.betas=matrix(NA,nparam.occ,ngr)
  for (i in 1:ngr){
    cond=w==i
    nw=sum(cond)
    const=(nw/tau2.betas)+(1/100)
    var1=diag(rep(1/const,nparam.occ))
    if (nw==0) soma.betas=matrix(0,nparam.occ,1)
    if (nw==1) soma.betas=betas[,cond]
    if (nw >1) soma.betas=matrix(rowSums(betas[,cond]),nparam.occ,1)
    pmedia=(1/tau2.betas)*soma.betas
    m.betas[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
  }
  m.betas
}
#-----------------------------------
sample.w=function(tau2.betas,betas,ltheta,w,ngr,m.betas,nparam.occ){
  p1=(-nparam.occ/2)*log(tau2.betas)
  for (i in 1:nspp){
    err2=(betas[,i]-m.betas)^2
    var.part=(-1/(2*tau2.betas))
    p2=colSums(var.part*err2)
    lprob=p1+p2+ltheta

    #do we need a new group?    
    max.w=max(w)
    if (max.w<ngr){
      lprob=lprob[1:max.w]

      #probability of new group
      ind=max.w+1
      p1=(-nparam.occ/2)*log(100+tau2.betas)
      err2=betas[,ind]^2
      var.part=-1/(2*(100+tau2.betas))
      p2=sum(var.part*err2)
      lprob=c(lprob,c(p1+p2+ltheta[ind]))
    }
    tmp=lprob-max(lprob)
    tmp=exp(tmp)
    prob=tmp/sum(tmp)
    ind=rmultinom(1,size=1,prob=prob)
    w[i]=which(ind==1)
  }
  w
}
#----------------------------------------------
sample.theta=function(gamma1,w,ngr){
  nk=rep(0,ngr)
  tmp=table(w)
  nk[as.numeric(names(tmp))]=tmp
  theta=v=rep(NA,ngr)
  aux=1
  for(i in 1:(ngr-1)){
    nk.maior=nk[(i+1):ngr]
    v[i]=rbeta(1,nk[i]+1,sum(nk.maior)+gamma1)
    theta[i]=v[i]*aux
    aux=aux*(1-v[i])
  }
  theta[ngr]=aux
  theta
}
#--------------------------------------------
