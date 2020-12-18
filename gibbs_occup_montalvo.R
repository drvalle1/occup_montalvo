gibbs_occup=function(y,xmat.occ,xmat.det,ngr,tau2.a,tau2.b,
                     gamma1,ngibbs,nburn,tau2.betas2,
                     ind.reff,ind.gr.reff){

  #initial values
  nspp=ncol(y[,,1])
  nloc=nrow(xmat.occ)
  nparam.occ=ncol(xmat.occ)
  nparam.occ1=length(ind.reff)
  nparam.occ2=length(ind.gr.reff)
  nparam.det=ncol(xmat.det[,,1])
  ystar=y
  ystar=ifelse(y==1,1,-1)
  m.betas1=rep(0,nparam.occ1)
  m.betas2=matrix(0,nparam.occ2,ngr)
  tau2.betas1=rep(1,nparam.occ1)
  w=sample(1:ngr,size=nspp,replace=T)
  betas=matrix(0,nparam.occ,nspp)
  theta=rep(1/ngr,ngr)
  z=apply(y,c(1,2),max)
  zstar=matrix(ifelse(z==1,1,-1),nloc,nspp)
  m.gamma=rep(0,nparam.det)
  tau2.gamma=rep(1,nparam.det)
  gammas=matrix(0,nparam.det,nspp)

  #MCMC settings
  store.m.betas2=matrix(NA,ngibbs,nparam.occ2*ngr)
  store.m.betas1=matrix(NA,ngibbs,nparam.occ1)
  store.tau2.betas1=matrix(NA,ngibbs,nparam.occ1)
  store.gammas=matrix(NA,ngibbs,nparam.det*nspp)  
  store.m.gamma=matrix(NA,ngibbs,nparam.det)
  store.tau2.gamma=matrix(NA,ngibbs,nparam.det)
  store.betas=matrix(NA,ngibbs,nparam.occ*nspp)
  store.w=matrix(NA,ngibbs,nspp)
  store.theta=matrix(NA,ngibbs,ngr)
  store.llk=matrix(NA,ngibbs,1)
  
  options(warn=2)

  #start gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    print(table(w))
    
    #update latent variables
    z=sample.z(xmat.occ=xmat.occ,betas=betas,nloc=nloc,
               nspp=nspp,nrep=nrep,y=y,xmat.det=xmat.det,gammas=gammas)
    # z=z.true
    zstar=sample.zstar(z=z,xmat.occ=xmat.occ,betas=betas,nloc=nloc,nspp=nspp)
    # zstar=zstar.true
    
    ystar=sample.ystar(nrep=nrep,xmat.det=xmat.det,gammas=gammas,
                       y=y,nspp=nspp)
    # ystar=ystar.true

    #sample betas and associated prior parameters
    betas=sample.betas(m.betas1=m.betas1,tau2.betas1=tau2.betas1,
                       m.betas2=m.betas2,tau2.betas2=tau2.betas2,
                       xmat.occ=xmat.occ,
                       zstar=zstar,nspp=nspp,nparam.occ=nparam.occ,
                       w=w,xtx.occ=xtx.occ,nparam.occ2=nparam.occ2)

    # betas=betas.true
    m.betas2=sample.m.betas2(w=w,betas=betas,tau2.betas2=tau2.betas2,nparam.occ2=nparam.occ2,
                             ngr=ngr,ind.gr.reff=ind.gr.reff)
    # tmp=m.betas2
    # tmp[,1:3]=m.betas2.true
    # m.betas2=tmp                             
    
    m.betas1=sample.m.betas1(nspp=nspp,betas=betas,tau2.betas1=tau2.betas1,
                             ind.reff=ind.reff,nparam.occ1=nparam.occ1)
    # m.betas1=m.betas1.true
    
    tau2.betas1=sample.tau2.betas1(nspp=nspp,betas=betas,ind.reff=ind.reff,
                                   nparam.occ1=nparam.occ1,
                                   m.betas1=m.betas1,tau2.a=tau2.a,tau2.b=tau2.b)

    #update gammas and associated prior parameters
    gammas=sample.gammas(ystar=ystar,xmat.det=xmat.det,z=z,m.gamma=m.gamma,tau2.gamma=tau2.gamma,
                         nparam.det=nparam.det,nspp=nspp,nrep=nrep)
    m.gamma=sample.m.gamma(gammas=gammas,tau2.gamma=tau2.gamma,nspp=nspp,nparam.det=nparam.det)
    tau2.gamma=sample.tau2.gamma(gammas=gammas,m=m.gamma,nspp=nspp,tau2.a=tau2.a,
                                 tau2.b=tau2.b,nparam.det=nparam.det)
    
    #sample other parameters
    w=sample.w(tau2.betas2=tau2.betas2,betas=betas,ltheta=log(theta),w=w,
               ngr=ngr,m.betas2=m.betas2,nparam.occ2=nparam.occ2,
               ind.gr.reff=ind.gr.reff)
    # w=w.true
    theta=sample.theta(gamma1=gamma1,w=w,ngr=ngr)
    # theta=rep(1/ngr,ngr)

    llk=get.llk(nloc=nloc,nspp=nspp,betas=betas,
                xmat.occ=xmat.occ,xmat.det=xmat.det,y=y,gammas=gammas)  
      
    #re-order w from time to time
    if (i<nburn & i%%50==0){
      k=table(w)
      k1=rep(0,ngr)
      k1[as.numeric(names(k))]=k
      ind=order(k1,decreasing=T)
      theta=theta[ind]
      m.betas2=m.betas2[,ind]
      wnew=w
      for (j in 1:ngr){
        cond=w==ind[j]
        wnew[cond]=j
      }
      w=wnew
    }
  
    #store results
    store.betas[i,]=betas
    store.m.betas2[i,]=m.betas2
    store.m.betas1[i,]=m.betas1
    store.tau2.betas1[i,]=tau2.betas1
    store.gammas[i,]=gammas    
    store.m.gamma[i,]=m.gamma
    store.tau2.gamma[i,]=tau2.gamma
    store.w[i,]=w
    store.theta[i,]=theta
    store.llk[i]=llk
  }
  
  list(betas=store.betas,m.betas1=store.m.betas1,m.betas2=store.m.betas2,
       tau2.betas1=store.tau2.betas1,
       gammas=store.gammas,m.gamma=store.m.gamma,tau2.gamma=store.tau2.gamma,
       w=store.w,theta=store.theta,
       llk=store.llk)
}