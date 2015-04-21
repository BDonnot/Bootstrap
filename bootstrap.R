install.packages("actuar")
library(actuar)


##basic functions for bootstraping
naif = function(data)
{
  sample(data,replace = T)
}

smooth = function(data)
{
  naif(data)+rnorm(length(data))
}

param = function(data,hat_shape,hat_scale)
{
  rpareto(length(data),shape = hat_shape, scale = hat_scale )
}

empiricQuantile = function(data,probs)
{
  quantile(data,probs = probs)
}

invFun = function(hat_shape,hat_scale,probs)
{
  qpareto(probs,shape = hat_shape, scale = hast_scale)
}

##more evoluate functions

naifParam = function(data,NsBootstrap,probs)
{
  #on simule les echantillons bootstrap (naivement)
  matBoot = sapply(1:NsBootstrap,function(el) naif(data)) #chaque colonne est une simulation
  
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



#les donnees
x = rpareto(30,shape=1,scale = 1)

qpareto(c(0.5,0.6),1,1)
