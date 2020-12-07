rm(list=ls())
set.seed(32)

#basic settings
nloc=1000
nrep=5
nspp=150
ngr=3
nparam.occ=4
nparam.det=3

#latent group membership
w.true=w=sample(1:ngr,size=nspp,replace=T)

#OCCUPANCY

#parameters
seq1=seq(from=-2,to=2,by=1)
m.betas.true=m.betas=matrix(sample(seq1,size=nparam.occ*ngr,replace=T),
                            nparam.occ,ngr)
tau2.betas=0.1 #set by user

#get betas
betas=matrix(NA,nparam.occ,nspp)
for (i in 1:nspp){
  mu1=m.betas[,w[i]]
  betas[,i]=rnorm(nparam.occ,mean=mu1,sd=sqrt(tau2.betas))
}
betas.true=betas

#visualize betas
rango=range(betas)
par(mfrow=c(3,1),mar=rep(1,4))
for (i in 1:ngr) image(betas[,w==i],zlim=rango)

#get intercepts
alpha.s.true=alpha.s=runif(nspp,min=-0.2,max=0.2)
alpha.mat=matrix(alpha.s,nloc,nspp,byrow=T)

#get covariates for occupancy
xmat.occ=matrix(rnorm(nloc*nparam.occ),nloc,nparam.occ)

#generate occupancy status
media=alpha.mat+xmat.occ%*%betas
zstar.true=zstar=matrix(rnorm(nloc*nspp,mean=media,sd=1),nloc,nspp)
z.true=z=matrix(ifelse(zstar>0,1,0),nloc,nspp)

#---------------------------------
#---------------------------------
#DETECTION
tau2.gamma.true=tau2.gamma=runif(nparam.det,min=0,max=0.3)
m.gamma.true=m.gamma=runif(nparam.det,min=-1,max=1)
gammas=matrix(NA,nparam.det,nspp)
for (i in 1:nparam.det){
  gammas[i,]=rnorm(nspp,mean=m.gamma[i],sd=sqrt(tau2.gamma[i]))
}
gammas.true=gammas

#get covariates for detection
xmat.det=array(NA,dim=c(nloc,nparam.det,nrep))
for (i in 1:nrep){
  tmp=rnorm(nloc*(nparam.det-1))
  tmp1=cbind(1,matrix(tmp,nloc,nparam.det-1))
  xmat.det[,,i]=tmp1
}

#generate observations
media=ystar=y=array(NA,dim=c(nloc,nspp,nrep))
for (i in 1:nrep){
  media[,,i]=xmat.det[,,i]%*%gammas
  ystar[,,i]=rnorm(nloc*nspp,mean=media[,,i],sd=1)
  #site has to be occupied and species has to have been detected
  y[,,i]=ifelse(ystar[,,i]>0 & z.true==1,1,0) 
}
ystar.true=ystar

#export data
setwd('U:\\GIT_models\\occup_montalvo\\fake data')
y1=matrix(y,nloc*nspp*nrep,1)
write.csv(y1,'fake data y.csv',row.names=F)
write.csv(xmat.occ,'fake data xmat occ.csv',row.names=F)
xmat.det1=matrix(xmat.det,nloc*nparam.det*nrep,1)
write.csv(xmat.det1,'fake data xmat det.csv',row.names=F)