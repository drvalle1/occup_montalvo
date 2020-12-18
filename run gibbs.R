# rm(list=ls())
library('mvtnorm')
set.seed(35)

#basic settings
nrep=5

#get functions
setwd('U:\\GIT_models\\occup_montalvo')
source('aux occup montalvo.R')
source('gibbs_occup_montalvo.R')

#get design matrix for occupancy
setwd('U:\\GIT_models\\occup_montalvo\\fake data')
xmat.occ=data.matrix(read.csv('fake data xmat occ.csv',as.is=T))
xtx.occ=t(xmat.occ)%*%xmat.occ
nparam.occ=ncol(xmat.occ)
nloc=nrow(xmat.occ)

#get data
tmp=read.csv('fake data y.csv',as.is=T)
nspp=nrow(tmp)/(nrep*nloc); nspp
y=array(tmp$V1,dim=c(nloc,nspp,nrep))

#design matrix for detection
tmp=read.csv('fake data xmat det.csv',as.is=T)
nparam.det=nrow(tmp)/(nloc*nrep); nparam.det
xmat.det=array(tmp$V1,dim=c(nloc,nparam.det,nrep))

#basic settings
ngr=10
ngibbs=1000
nburn=ngibbs/2
tau2.betas2=0.5

#priors
tau2.a=0.1; tau2.b=0.1
gamma1=0.01
ind.reff=1:2
ind.gr.reff=3:ncol(xmat.occ)

#run gibbs
mod1=gibbs_occup(y=y,xmat.occ=xmat.occ,xmat.det=xmat.det,ngr=ngr,
                 tau2.a=tau2.a,tau2.b=tau2.b,gamma1=gamma1,
                 ind.reff=ind.reff,ind.gr.reff=ind.gr.reff,
                 ngibbs=ngibbs,nburn=nburn,tau2.betas2=tau2.betas2)