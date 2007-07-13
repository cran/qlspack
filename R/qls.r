#QLS

#family="gaussian","binomial","poisson"
#correlation="ar1","exchangeable","markov","tridiagonal","ex.fam"
#std.err="san.se","jack","j1s","fij"

qls<- function(formula,data,id,family=gaussian,time="NA",
             correlation="ar1",std.err="san.se"){
  
  if(is.character(family)){
    family<- switch(family,
                  "gaussian"="gaussian",
                  "Gaussian"="gaussian",
                  "binomial"="binomial",
                  "Binomial"="binomial",
                  "Poisson"="poisson")
  }else{ 
    if(is.function(family)){
      family<- family()[[1]]
    }else{
       print(family)
       stop("'family' not recognized")
    }    
  }
     
  epsilon<-1e-6
  dataset<- proc(data=data,time=time,id=id)
  data<- dataset$data
  m<- data$m
  n<- data$n
  id<- data$id
  time<- data$time
  autotime<- data$autotime
  
  alpha.pre<- 0 
  fit.pre<- geeglm(formula=formula,data=data,id=id,family=family,
                   corstr="independence",std.err=std.err)
  beta.pre<- fit.pre$coef
  
  yname<- formula[[2]]
  y<- data[,which(colnames(data)==yname)]
  xname<- strsplit(as.character(formula)," ")[[3]]
  xlen<- length(xname)
  xname<- xname[seq(1,xlen,2)]
  xlen<- length(xname)
  x<- matrix(0,nrow=length(y),ncol=xlen)
  for(j in 1:xlen){
    x[,j]<- data[,which(colnames(data)==xname[j])]
  } 
  x<-cbind(rep(1,length(y)),x)

  resid<- residual(x,y,beta.pre,family=family)
 
  alpha<- switch(correlation,
                 ar1= ar1.one(resid=resid,time=time,id=id),
                 exchangeable= exch.one(resid=resid, time=time, id=id),
                 markov=markov.one(resid=resid,time=time,id=id),
                 tridiagonal= tri.one(resid=resid,time=time,id=id),
                 ex.fam= exfam.one(resid=resid,id=id)
                 )

  iteration<- 0
  while(sum((alpha-alpha.pre)^2)>epsilon){
    iteration<- iteration+1
    geefit<- gee.fixed(formula=formula,data=data,id=id,alpha=alpha,std.err=std.err,
                       family=family,time=time,correlation=correlation)
    beta<- geefit[[1]]
    resid<- residual(x,y,beta,family=family)
    alpha.pre<- alpha
    alpha<- switch(correlation,
                   ar1= ar1.one(resid=resid,time=time,id=id),
                   exchangeable= exch.one(resid=resid, time=time, id=id),
                   markov=markov.one(resid=resid,time=time,id=id),
                   tridiagonal= tri.one(resid=resid,time=time,id=id),
                   ex.fam= exfam.one(resid=resid,id=id)
                   )
   
  }    
  ru<- switch(correlation,
              ar1= ar1.two(alpha,resid,time,id),
              exchangeable= exch.two(alpha,resid,time,id),
              markov= markov.two(alpha,resid,time,id),
              tridiagonal=tri.two(alpha,resid,time,id),
              ex.fam= exfam.two(alpha,id=id)      
              )
  
  qlsfit<- gee.fixed(formula=formula,data,id=id,alpha=ru,
                  family=family,time=time,correlation=correlation,
                  std.err=std.err)
  qlsfit$geese$alpha<- ru
  parn<- length(ru)
  qlsfit$geese$zcor.names<- rep(0,parn)
  for (i in 1:parn){
    qlsfit$geese$zcor.names[i]<- paste("alpha",i,sep=":")
  }
  names(qlsfit$geese$alpha)<- qlsfit$geese$zcor.names
  qlsfit$geese$model$corstr<- correlation
  qlsfit$call<-match.call()
  qlsfit
  return(qlsfit)
}
