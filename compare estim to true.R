compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango,col='red',lwd=2)  
}

#how many groups
plot(mod1$theta[ngibbs,],type='h')

#compare group membership
aux=data.frame(w.true=w.true,w.estim=mod1$w[ngibbs,])
k=table(aux); k
seq1=c(3,1,2)
k[,seq1]

#compare betas
betas.estim=matrix(mod1$betas[ngibbs,],nparam.occ,nspp)
compare1(estim=betas.estim,true=betas.true)

#compare m.betas
ngr=10
m.betas.estim=matrix(mod1$m.betas[ngibbs,],nparam.occ,ngr)
ind=seq1
compare1(m.betas.estim[,ind],m.betas.true)

#compare alpha.s
alpha.s.estim=mod1$alpha.s[ngibbs,]
compare1(estim=alpha.s.estim,true=alpha.s.true)

#compare gammas
gammas.estim=matrix(mod1$gammas[ngibbs,],nparam.det,nspp)
compare1(estim=gammas.estim,true=gammas.true)

#compare m.gammas
m.gammas.estim=mod1$m.gamma[ngibbs,]
compare1(estim=m.gammas.estim,true=m.gamma.true)

#compare tau2.gammas
tau2.gammas.estim=mod1$tau2.gamma[ngibbs,]
compare1(estim=tau2.gammas.estim,true=tau2.gamma.true)