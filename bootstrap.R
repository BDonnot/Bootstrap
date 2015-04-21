# install.packages("actuar")
library(actuar)

##shape
shape = 1

##basic functions for bootstraping
naif = function(data)
{
  sample(data,replace = T)
}

smooth = function(data,sd = 0.1)
{
  naif(data)+rnorm(length(data),0,sd=sd)
}

param = function(data,hat_scale)
{
  rpareto(length(data),shape = shape, scale = hat_scale )
}

empiricQuantile = function(data,probs)
{
  quantile(data,probs = probs)
}

invFun = function(hat_scale,probs)
{
  qpareto(probs,shape = shape, scale = hast_scale)
}

getParamHill = function(data)
{
  #return the parametric estimation of the coefficients shape and scale
  res = list()
  betaHat = 1/mean(log(data))
  res$hast_scale = betaHat**2
  return(res)
}
##more evoluate functions

basicParam = function(data,NsBootstrap,probs,funBoot)
{
  #on simule les echantillons bootstrap (naivement)
  matBoot = sapply(1:NsBootstrap,function(el) funBoot(data)) #chaque colonne est une simulation
  
  #on applique la fonction de quantile empirique :
  quantiles = apply(matBoot,2,function(col) empiricQuantile(col,probs = probs)) 
  #meme ligne : meme quantile
  
  resQuant = apply(quantiles,1,mean)
  
  confidenceInter = apply(quantiles,1,function(line) quantile(line,probs = c(0.05,0.95)))
  
  res = list()
  res$resQuant = resQuant
  res$confidenceInter = confidenceInter
  return(res)
}

naifParam = function(data,NsBootstrap,probs)
{
  res = basicParam(data,NsBootstrap,probs,naif)
  return(res)
}

smoothParam = function(data,NsBootstrap,probs,sdSmooth = 0.1)
{
  smoothSpecifiq = function(x) 
  {
    smooth(x,sd = sdSmooth)
  }
  res = basicParam(data,NsBootstrap,probs,smoothSpecifiq)
  return(res)
}


#les donnees
x = rpareto(30,shape=shape,scale = 1)

naifParam(x,1000,c(0.01,0.99))
smoothParam(x,1000,c(0.01,0.99),0.1)


x = rpareto(10000000,shape=1,scale = 2)
betaHat = (1/mean(log(x)))**2
betaHat
