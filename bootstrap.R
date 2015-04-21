# install.packages("actuar")
library(actuar)
library(e1071) 


##shape
scale = 1

##basic functions for bootstraping
naif = function(data)
{
  sample(data,replace = T)
}

smooth = function(data,sd = 0.1)
{
  naif(data)+rnorm(length(data),0,sd=sd)
}

param = function(data,hat_shape)
{
  rpareto(length(data),shape = hat_shape, scale = scale )+1
}

empiricQuantile = function(data,probs)
{
  quantile(data,probs = probs,type = 1)
}

invFun = function(hat_shape,probs)
{
  qpareto(probs,shape = hat_shape, scale = scale)+1
}

getParamEMV = function(data)
{
  #return the parametric estimation of the coefficients shape and scale
  res = list()
  betaHat = 1/mean(log(data))
  res$hat_shape = betaHat
  return(res)
}

getParamHill = function(data,k = 10)
{
  dataS = log(sort(data,decreasing = T)[1:k])
  res = list()
  res$hat_shape = 1/(mean(dataS[1:(k-1)]) - dataS[k])
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

paramBootstrap = function(data,NsBootstrap,probs)
{
  shape_hat = getParamHill(data)$hat_shape
  paramSpecifiq = function(x) 
  {
    param(x,shape_hat)
  }
  res = basicParam(data,NsBootstrap,probs,paramSpecifiq)
  return(res)
}

asymptotiqueEst = function(data,probs)
{
  dataSorted = sort(data)
  n = length(data)
  nalpha = qnorm(c(.975),0,1)
  i1 = round(n*probs[1]-sqrt(n)*nalpha*sqrt(probs[1]*(1-probs[1])),0 )
  i2 = round(n*probs[1]+sqrt(n)*nalpha*sqrt(probs[1]*(1-probs[1])),0 )
  res = list()
  res$min = dataSorted[i1]
  res$max = dataSorted[i2]
  return(res)
}

#les donnees
shape = 2
x = rpareto(30,shape=shape,scale = scale)+1

probs = c(.5,.6,.7,.75,.90,.95,.99)
qpareto(probs,shape = shape, scale = 1)+1
naifParam(x,10000,probs)
smoothParam(x,10000,probs,log(length(x))/length(x) )
smoothParam(x,10000,probs,1/sqrt(length(x)) )
paramBootstrap(x,10000,probs,1/sqrt(length(x)) )
qpareto(probs,shape = shape, scale = 1)+1


p = .75
asymptotiqueEst(x,p)

c(quantile(x,p)-(p*(1-p))/dpareto(qpareto(p,shape = shape, scale = scale),shape = shape, scale = scale)**2*1.96/sqrt(length(x))
  ,quantile(x,p)+(p*(1-p))/dpareto(qpareto(p,shape = shape, scale = scale),shape = shape, scale = scale)**2*1.96/sqrt(length(x))
)



x = rpareto(30000,shape=2,scale =1)+1
getParamEMV(x)
getParamHill(x,sqrt(length(x)))

# x = rpareto(10000000,shape=12,scale =1)+1
# betaHat = 1/mean(log(x))
# betaHat
# hist(x)
