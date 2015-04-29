# install.packages("actuar")
library(actuar)
library(e1071) 

##Plot of the pareto distribution
xplot = seq(0.001,5,length = 1e5)
plot(xplot,dpareto(xplot,scale = 1, shape = 1),
     xlab = "",
     ylab = "",
     main = "Densite Pareto pour different paremetre de forme",
     type = 'l')
lines(xplot,dpareto(xplot,scale = 1, shape = 2),col = 'blue')
lines(xplot,dpareto(xplot,scale = 1, shape = 3),col = 'red')
lines(xplot,dpareto(xplot,scale = 1, shape = 5),col = 'darkgreen')
legend("topright",
       legend = c("beta = 1","beta = 2","beta = 3","beta = 5"),
       col = c("black","red","blue","darkgreen"),
       lty = 1)
##assign the scale (unchanged)
scale = 1
probsInter = c(2.5,97.5)
names(probsInter) = paste0(as.character(probsInter),"\\%")
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

getParamHill = function(data,kHill = 10)
{
  k = round(kHill,0)
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
  row.names(quantiles) = names(probs)
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

paramBootstrap = function(data,NsBootstrap,probs,funParam = getParamHill,...)
{
  #bootstrap parametrique en utilisant Hill pour determiner le parametre
  #parametrique = 1c
  shape_hat = funParam(data,...)$hat_shape
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
  result = list()
  result$resQuant = apply(res,2,mean)
  names(result$resQuant) = names(probs)
  colnames(res) = names(probs)
  rownames(res) = names(probsInter)
  result$confidenceInter = res
  
  return(result)
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
  result = list()
  result$resQuant = apply(res,2,mean)
  names(result$resQuant) = names(probs)
  colnames(res) = names(probs)
  rownames(res) = names(probsInter)
  result$confidenceInter = res
  return(result)
}

#library(xtable)
getResults = function(listOfRes)
{
  #fonction pour exporter directemment les resultats en tableau latex
  #normalement ca se fait avec xtable, mais la ca ne marchait pas ...
  result = matrix("",ncol = length(listOfRes)+1,nrow = length(listOfRes[[1]][[1]])+1)
  result[1,] = c("",names(listOfRes))
  result[,1] = c("",names(listOfRes[[1]][[1]]))
  n = length(listOfRes[[1]][[1]])+1
  j = 2
  for(res in listOfRes)
  {
    tempVal = res[[1]]
    tempIC = res[[2]]
    result[2:n,j] = paste0("\\begin{tabular}{c}$",
                     round(tempVal,2),"$\\\\",
                     apply(tempIC,2,function(IC) paste("{\\small $[",
                                                       round(IC[1],2),",",
                                                       round(IC[2],2),"]$ }","\\end{tabular}" ) ) )
    j = j+1
  }
  result1 = apply(result,1,function(x) paste(x,sep = "&",collapse = "&") )
  return(cat(result1,sep = "\\\\\n"))
}

#les donnees
shape = 2
set.seed(1)
x = rpareto(30,shape=shape,scale = scale)+1
max=1-(1/length(x))
probs = c(.75,.90,max)
names(probs) = c("75\\%","90\\%","max")
qpareto(probs,shape = shape, scale = 1)+1

#partie 1
resTh = getTheoricBounds(x,probs,scale)
resAs = asymptotiqueEst(x,probs)
resNaif = naifBootstrap(x,10000,probs)
resSmooth = smoothBootstrap(x,10000,probs,1/sqrt(length(x)) )
res = list(resTh,resAs,resNaif,resSmooth)
names(res) = c("IC Oracle","IC Asympt.","IC Naif","IC Lisse")
getResults(res)

#smoothBootstrap(x,10000,probs,1/sqrt(length(x)) )

#partie 2
set.seed(1)
x = rpareto(100,shape=shape,scale = scale)+1
max=1-(1/length(x))
probs = c(.75,.90,max)
names(probs) = c("75\\%","90\\%","max")
qpareto(probs,shape = shape, scale = 1)+1

#1c
resParamEMV = paramBootstrap(data=x,NsBootstrap=10000,probs=probs,funParam=getParamEMV,kHill = NULL)
resParamHill = paramBootstrap(data=x,NsBootstrap=10000,probs=probs,funParam=getParamHill,kHill = length(x)/3)

#2a
resAEMV = bootstrapModifa(x,1000,probs,funParam=getParamEMV,kHill = NULL)
resAHill = paramBootstrap(data=x,NsBootstrap=10000,probs=probs,funParam=getParamHill,kHill = length(x)/3)










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
