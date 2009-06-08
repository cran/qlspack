#QLS

#family="gaussian","binomial","poisson"
#correlation="ar1","exchangeable","markov","tridiagonal","fam","ex.fam"
#std.err="san.se","jack","j1s","fij"

qls<- function(formula,data,id,family=gaussian,time=NULL,
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
  dt.fm<- data.frame(data)
  dataset<- data.proc(data=dt.fm,formula,time=time,id=id,del.n=0)
  #print(formula) 
  data<- dataset$data
  #print(data)
  m<- dataset$m
  n<- dataset$n
  id<- dataset$id
  time<- dataset$time
  autotime<- dataset$autotime
  
  alpha.pre<- 0 
  fit.pre<- geeglm(formula=formula,data=dt.fm,id=id,family=family,
                   corstr="independence",std.err=std.err)
  #cat("fit.pre:","\n")
  #print(fit.pre)
  beta<- as.matrix(coef(fit.pre))
  #print(beta)
  #print(data)
  data.end<- ncol(data)
  xx<- data[,-data.end]
  yy<- data[,data.end]
  #print(xx)
  #print(yy)
  resid<- residual(xx,yy,beta,family=family)
  #print(resid)
  #print(id)
  
  #data.proc(data=resid,formula=resid~resid,time=NULL,id=id,del.n=2)
  
  alpha<- switch(correlation,
                 ar1= ar1.one(resid=resid,time=time,id=id),
                 exchangeable= exch.one(resid=resid,time=time, id=id),
                 markov=markov.one(resid=resid,time=time,id=id),
                 tridiagonal= tri.one(resid=resid,time=time,id=id),
                 fam = fam.one(resid=resid,time=NULL,id=id),
                 ex.fam= exfam.one(resid=resid,time=NULL,id=id)
                )
  #cat("alpha0",alpha,"\n")
  
  iteration<- 0
  while(sum((alpha-alpha.pre)^2)>epsilon){
    iteration<- iteration+1
    #cat("iteration",iteration,"\n")
    geefit<- gee.fixed(formula=formula,data=dt.fm,id=id,alpha=alpha,std.err=std.err,
                       family=family,time=time,correlation=correlation)
    beta<- as.matrix(coef(geefit))
    #cat("beta",beta,"\n")
    resid<- residual(xx,yy,beta,family=family)
    alpha.pre<- alpha
    alpha<- switch(correlation,
                   ar1= ar1.one(resid=resid,time=NULL,id=id),
                   exchangeable= exch.one(resid=resid,time=NULL, id=id),
                   markov=markov.one(resid=resid,time=time,id=id),
                   tridiagonal= tri.one(resid=resid,time=NULL,id=id),
                   fam = fam.one(resid=resid,id=id),
                   ex.fam= exfam.one(resid=resid,id=id)
                   )
   #cat("alpha",alpha,"\n\n")
  }    
  ru<- switch(correlation,
              ar1= ar1.two(alpha,resid,time,id),
              exchangeable= exch.two(alpha,resid,time,id),
              markov= markov.two(alpha,resid,time,id),
              tridiagonal=tri.two(alpha,resid,time,id),
              fam = fam.two(alpha,id=id),
              ex.fam= exfam.two(alpha,id=id)      
              )
  #cat("ru",ru,"\n\n")
  
  qlsfit<- gee.fixed(formula=formula,data=dt.fm,id=id,alpha=ru,
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
