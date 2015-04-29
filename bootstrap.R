# install.packages("actuar")
library(actuar)
library(e1071) 


##assign the scale (unchanged)
scale = 1
probsInter = c(2.5,97.5)
names(probsInter) = paste0(as.character(probsInter),"%")
probsInter = probsInter/100
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

getParamEMV = function(data,kHill)
{
  #return the parametric estimation of the coefficients shape and scale
  #kHill unused
  res = list()
  betaHat = 1/mean(log(data))
  res$hat_shape = betaHat
  return(res)
}

getParamHill = function(data,k = 10)
{
  k = round(k)
  dataS = log(sort(data,decreasing = T)[1:k])
  res = list()
  res$hat_shape = 1/(mean(dataS[1:(k-1)]) - dataS[k])
  return(res)
}

##more evoluate functions
basicParam = function(data,NsBootstrap,probs,funBoot,funparam=empiricQuantile)
{
  #on simule les echantillons bootstrap (naivement)
  matBoot = sapply(1:NsBootstrap,function(el) funBoot(data)) #chaque colonne est une simulation
  
  #on applique la fonction de quantile empirique :
  quantiles = apply(matBoot,2,function(col) funparam(col,probs)) 
  row.names(quantiles) = probs
  #meme ligne : meme quantile
  
  resQuant = apply(quantiles,1,mean)
  
  confidenceInter = apply(quantiles,1,function(line) quantile(line,probs = probsInter))
  
  plotFun =function(line,main)
  {
    hist(line,100,main=main,xlab = "")
    abline(v = quantile(line,probs = probsInter),col = "blue",lty = 2)
  }
  nline = 1
  for(p in probs)
  {
    plotFun(quantiles[nline,],main = names(probs)[nline])
    nline = nline +1
    
  }
  
  res = list()
  res$resQuant = resQuant
  res$confidenceInter = confidenceInter
  return(res)
}

naifBootstrap = function(data,NsBootstrap,probs)
{
  #bootstrap naif
  #1a
  res = basicParam(data,NsBootstrap,probs,naif)
  return(res)
}

smoothBootstrap = function(data,NsBootstrap,probs,sdSmooth = 0.1)
{
  #smooth bootstrap
  #1b
  smoothSpecifiq = function(x) 
  {
    smooth(x,sd = sdSmooth)
  }
  res = basicParam(data,NsBootstrap,probs,smoothSpecifiq)
  return(res)
}

paramBootstrap = function(data,NsBootstrap,probs,funParam = getParamHill)
{
  #bootstrap parametrique en utilisant Hill pour determiner le parametre
  #parametrique = 1c
  shape_hat = funParam(data)$hat_shape
  paramSpecifiq = function(x) 
  {
    param(x,shape_hat)
  }
  res = basicParam(data,NsBootstrap,probs,paramSpecifiq)
  return(res)
}


bootstrapModifa = function(data,
                           NsBootstrap,
                           probs,
                           funParam = getParamHill,
                           sdSmooth = log(length(data))/length(data),
                           kHill = sqrt(length(data)))
{
  #2a
  funBootstrap = function(x,probs)
  {
    estShape = funParam(x,kHill)$hat_shape
    return(invFun(estShape,probs))
  }
  res = basicParam(data,NsBootstrap,probs,naif,funBootstrap)
  return(res)
}

bootstrapModifb = function(data,
                           NsBootstrap,
                           probs,
                           funParam = getParamHill,
                           sdSmooth = log(length(data))/length(data),
                           kHill = sqrt(length(data)))
{
  #2b
  funBootstrap = function(x,probs)
  {
    estShape = funParam(x,kHill)$hat_shape
    return(invFun(estShape,probs))
  }
  smoothSpecifiq = function(x) 
  {
    smooth(x,sd = sdSmooth)
  }
  res = basicParam(data,NsBootstrap,probs,smoothSpecifiq,funBootstrap)
  return(res)
}

asymptotiqueEst = function(data,probs)
{
  #estimateur d'erreur asymptotique, cf formule de Clara
  #intervalle de confiance 'Clara'
  res = c()
  for(p in probs)
  {
    dataSorted = sort(data)
    n = length(data)
    nalpha = qnorm(c(.975),0,1)
    i1 = round(n*p-sqrt(n)*nalpha*sqrt(p*(1-p)),0 )
    i2 = round(n*p+sqrt(n)*nalpha*sqrt(p*(1-p)),0 )
    res = cbind(res,c(data[i1],data[i2]))
  }
  colnames(res) = probs
  rownames(res) = names(probsInters)
  return(res)
}

getTheoricBounds = function(data,probs,shape)
{
  #intervalles de confiance asymptotique, en connaissant vraiment shape et scale (oracle)
  #intervalle de confiance 'Matthieu'
  res = c()
  for(p in probs)
  {
    inter = matrix(c(quantile(data,p)-(p*(1-p))/(dpareto(qpareto(p,shape = shape, scale = scale),shape = shape, scale = scale)*1.96*sqrt(length(data)))
                     ,quantile(data,p)+(p*(1-p))/(dpareto(qpareto(p,shape = shape, scale = scale),shape = shape, scale = scale)*1.96*sqrt(length(data)))
    ),nrow = 2)
    
    res = cbind(res,inter)
  }
  rownames(res) = c("5%","95%")
  colnames(res) = probs
  return(res)
}



#les donnees
shape = 2
set.seed(1)
x = rpareto(30,shape=shape,scale = scale)+1
max=1-(1/length(x))
probs = c(.75,.90,max)
names(probs) = c(".75",".90","max")
qpareto(probs,shape = shape, scale = 1)+1
naifBootstrap(x,10000,probs)

smoothBoostrap(x,10000,probs,log(length(x))/length(x) )
smoothBootstrap(x,10000,probs,1/sqrt(length(x)) )
paramBootstrap(x,10000,probs)

bootstrapModifa(x,1000,probs)
bootstrapModifb(x,1000,probs,log(length(x))/length(x))

#les donnees
x = rpareto(100,shape=shape,scale = scale)+1
max=1-(1/length(x))
probs = c(.75,.90,max)
names(probs) = c(".75",".90","max")
qpareto(probs,shape = shape, scale = 1)+1
naifBootstrap(x,10000,probs)

smoothBoostrap(x,10000,probs,log(length(x))/length(x) )
smoothBootstrap(x,10000,probs,1/sqrt(length(x)) )
paramBootstrap(x,10000,probs)

bootstrapModifa(x,1000,probs,kHill=10)
bootstrapModifa(x,1000,probs,kHill=20)
bootstrapModifb(x,1000,probs,log(length(x))/length(x))

#il y a probablement un probleme avec Hill ...
p = .75
asymptotiqueEst(x,probs)

getTheoricBounds(x,probs,scale)


x = rpareto(30,shape=2,scale =1)+1
getParamEMV(x)
getParamHill(x,5)

# x = rpareto(10000000,shape=12,scale =1)+1
# betaHat = 1/mean(log(x))
# betaHat
# hist(x)
