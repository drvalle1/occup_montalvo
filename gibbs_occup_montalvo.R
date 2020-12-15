gibbs_occup=function(y,xmat.occ,xmat.det,ngr,tau2.a,tau2.b,
                     gamma1,ngibbs,nburn,tau2.betas){

  #initial values
  nspp=ncol(y[,,1])
  nloc=nrow(xmat.occ)
  nparam.occ=ncol(xmat.occ)
  nparam.det=ncol(xmat.det[,,1])
  ystar=y
  ystar=ifelse(y==1,1,-1)
  m.betas=matrix(0,nparam.occ,ngr)
  w=sample(1:ngr,size=nspp,replace=T)
  alpha.s=rep(0,nspp)
  betas=matrix(0,nparam.occ,nspp)
  theta=rep(1/ngr,ngr)
  z=apply(y,c(1,2),max)
  zstar=matrix(ifelse(z==1,1,-1),nloc,nspp)
  m.gamma=rep(0,nparam.det)
  tau2.gamma=rep(1,nparam.det)
  m.alpha=0
  tau2.alpha=1
  gammas=matrix(0,nparam.det,nspp)

  #MCMC settings
  store.alpha.s=matrix(NA,ngibbs,nspp)
  store.m.alpha=matrix(NA,ngibbs,1)
  store.tau2.alpha=matrix(NA,ngibbs,1)
  store.gammas=matrix(NA,ngibbs,nparam.det*nspp)  
  store.m.gamma=matrix(NA,ngibbs,nparam.det)
  store.tau2.gamma=matrix(NA,ngibbs,nparam.det)
  store.betas=matrix(NA,ngibbs,nparam.occ*nspp)
  store.m.betas=matrix(NA,ngibbs,nparam.occ*ngr)
  store.w=matrix(NA,ngibbs,nspp)
  store.theta=matrix(NA,ngibbs,ngr)
  store.llk=matrix(NA,ngibbs,1)
  
  options(warn=2)

  #start gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    print(table(w))
    
    #update latent variables
    z=sample.z(xmat.occ=xmat.occ,betas=betas,nloc=nloc,alpha.s=alpha.s,
               nspp=nspp,nrep=nrep,y=y,xmat.det=xmat.det,gammas=gammas)
    # z=z.true
    zstar=sample.zstar(z=z,xmat.occ=xmat.occ,betas=betas,nloc=nloc,nspp=nspp,
                       alpha.s=alpha.s)
    # zstar=zstar.true
    ystar=sample.ystar(nrep=nrep,xmat.det=xmat.det,gammas=gammas,
                       y=y,nspp=nspp)
    # ystar=ystar.true

    #sample betas and associated prior parameters
    betas=sample.betas(m.betas=m.betas,tau2.betas=tau2.betas,xmat.occ=xmat.occ,
                       zstar=zstar,nspp=nspp,nparam.occ=nparam.occ,
                       w=w,xtx.occ=xtx.occ,alpha.s=alpha.s)
    # betas=betas.true
    m.betas=sample.m.betas(w=w,betas=betas,tau2.betas=tau2.betas,nparam.occ=nparam.occ,ngr=ngr)
    
    #update intercept and associated prior parameters
    alpha.s=sample.alpha.s(m.alpha=m.alpha,tau2.alpha=tau2.alpha,
                           xmat.occ=xmat.occ,zstar=zstar,nspp=nspp,
                           w=w,betas=betas,nloc=nloc)
    m.alpha=sample.m.alpha(nspp=nspp,alpha.s=alpha.s,tau2.alpha=tau2.alpha)
    tau2.alpha=sample.tau2.alpha(nspp=nspp,alpha.s=alpha.s,
                                 m.alpha=m.alpha,tau2.a=tau2.a,tau2.b=tau2.b)
    
    #update gammas and associated prior parameters
    gammas=sample.gammas(ystar=ystar,xmat.det=xmat.det,z=z,m.gamma=m.gamma,tau2.gamma=tau2.gamma,
                         nparam.det=nparam.det,nspp=nspp,nrep=nrep)
    m.gamma=sample.m.gamma(gammas=gammas,tau2.gamma=tau2.gamma,nspp=nspp,nparam.det=nparam.det)
    tau2.gamma=sample.tau2.gamma(gammas=gammas,m=m.gamma,nspp=nspp,tau2.a=tau2.a,
                                 tau2.b=tau2.b,nparam.det=nparam.det)
    
    #sample other parameters
    w=sample.w(tau2.betas=tau2.betas,betas=betas,ltheta=log(theta),w=w,
               ngr=ngr,m.betas=m.betas,nparam.occ=nparam.occ)
    theta=sample.theta(gamma1=gamma1,w=w,ngr=ngr)
    # theta=rep(1/ngr,ngr)

    llk=get.llk(alpha.s=alpha.s,nloc=nloc,nspp=nspp,betas=betas,
                xmat.occ=xmat.occ,xmat.det=xmat.det,y=y,gammas=gammas)  
      
    #re-order w from time to time
    if (i<nburn & i%%50==0){
      k=table(w)
      k1=rep(0,ngr)
      k1[as.numeric(names(k))]=k
      ind=order(k1,decreasing=T)
      theta=theta[ind]
      m.betas=m.betas[,ind]
      wnew=w
      for (j in 1:ngr){
        cond=w==ind[j]
        wnew[cond]=j
      }
      w=wnew
    }
  
    #store results
    store.betas[i,]=betas
    store.m.betas[i,]=m.betas
    store.alpha.s[i,]=alpha.s
    store.m.alpha[i,]=m.alpha
    store.tau2.alpha[i,]=tau2.alpha
    store.gammas[i,]=gammas    
    store.m.gamma[i,]=m.gamma
    store.tau2.gamma[i,]=tau2.gamma
    store.w[i,]=w
    store.theta[i,]=theta
    store.llk[i]=llk
  }
  
  list(betas=store.betas,m.betas=store.m.betas,
       alpha.s=store.alpha.s,m.alpha=store.m.alpha,tau2.alpha=store.tau2.alpha,
       gammas=store.gammas,m.gamma=store.m.gamma,tau2.gamma=store.tau2.gamma,
       w=store.w,theta=store.theta,
       llk=store.llk)
}