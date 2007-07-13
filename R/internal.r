#########################Getting Info from the Data#################

#This function will remain all the subjects.
#This function will not help order the subjects.

cluster.size<- function(id){
  clid<- unique(id)
  m<- length(unique(id))
  n<- rep(0,m)
  autotime<- rep(0,0)
  for(i in 1:m){
    n[i]<- length(which(id==clid[i]))
    autotime<- c(autotime,1:n[i])
  }
  id<- rep(1:m,n)
  return(list(m=m,n=n,id=id,autotime=autotime))
}

####################################################################

#################################Data Process#######################

#This function will not delete any subjects.

proc<- function(data="NA",time="NA",id){
  
  data<- as.data.frame(data)
  data.end<- ncol(data)
  data<- data.frame(data,time,id)
  colnames(data)[(data.end+1):(data.end+2)]<- c("time","id")
  index<- order(data$id,data$time)
  data<- data[index,]

  cluster<- cluster.size(data$id)
  m<- cluster$m
  n<- cluster$n
  data$id<- cluster$id
  if(data$time[1]=="NA"){
    time<- cluster$autotime
  }else{
    time<- data$time
  }
  autotime<- cluster$autotime
  
  return(list(data=data,time=time,autotime=autotime,id=id,m=m,n=n))
}

####################################################################

#This function will delete subjects with only one observations.

data.process<- function(data="NA",time="NA",id){
  
  data<- as.data.frame(data)
  data.end<- ncol(data)
  data<- data.frame(data,time,id)
  colnames(data)[(data.end+1):(data.end+2)]<- c("time","id")
  index<- order(data$id,data$time)
  data<- data[index,]

  cluster<- cluster.size(data$id)
  m<- cluster$m
  n<- cluster$n
  data$id<- cluster$id

  index<- which(n>1)
  index2<- which(data$id==index[1])
  for(i in 2:m){
    index2<-c(index2,which(data$id==index[i]))  
  }
  data<- data[index2,]


  cluster.new<- cluster.size(data$id)
  id.new<- cluster.new$id
  m.new<- cluster.new$m
  n.new<- cluster.new$n
  data.new<- data[,1:data.end]
  if(data$time[1]=="NA"){
    time.new<- cluster.new$autotime
  }else{
    time.new<- data$time
  }
  autotime<- cluster.new$autotime
  
  return(list(data=data.new,time=time.new,autotime=autotime,id=id.new,m=m.new,n=n.new))
}

##############################################################################

#This function will delete subjects with fewer than 2 observations.

data.proc.exfam<- function(data="NA",time="NA",id){
  
  data<- as.data.frame(data)
  colnames(data)<- "data"
  data.end<- ncol(data)
  data<- data$data
  data<- data.frame(data,time,id)
  colnames(data)[(data.end+1):(data.end+2)]<- c("time","id")
  index<- order(data$id,data$time)
  data<- data[index,]

  cluster<- cluster.size(data$id)
  m<- cluster$m
  n<- cluster$n
  data$id<- cluster$id

  index<- which(n>2)
  index2<- which(data$id==index[1])
  for(i in 2:m){
    index2<-c(index2,which(data$id==index[i]))  
  }
  data<- data[index2,]


  cluster.new<- cluster.size(data$id)
  id.new<- cluster.new$id
  m.new<- cluster.new$m
  n.new<- cluster.new$n
  data.new<- data[,1:data.end]
  if(data$time[1]=="NA"){
    time.new<- cluster.new$autotime
  }else{
    time.new<- data$time
  }
  autotime<- cluster.new$autotime
  
  return(list(data=data.new,time=time.new,autotime=autotime,id=id.new,m=m.new,n=n.new))
}

#################################################################################

##########################Calculating Residuals##################################

residual<- function(x,y,beta,family="gaussian"){
  x<- as.matrix(x)
  y<- as.matrix(y)
  beta<- as.matrix(beta)
  u<- switch(family, 
        gaussian=x%*%beta,
        binomial=exp(x%*%beta)/(1+exp(x%*%beta)),
        poisson=exp(x%*%beta)
      )
  h<- switch(family,
        gaussian=1,
        binomial=u*(1-u),
        poisson=u     
      ) 
  return((y-u)/sqrt(h)) 
}

###########################Correlation Matrix############################

#AR1 Structure

cormax.ar1<- function(alpha,id,time="NA"){
  cluster<- cluster.size(id)
  n<- cluster$n
  
  n.max<- max(n)
  cor.max<- diag(1,n.max)
  lowertri<- rep(0,0)
  for(j in (n.max-1):1){
    lowertri<- c(lowertri,1:j)
  }
  cor.max[lower.tri(cor.max)]<- alpha^lowertri
  cor.max[upper.tri(cor.max)]<- alpha^lowertri[length(lowertri):1]
  return(cor.max)
}

#####################################################################

#Exchangeable Structure

cormax.exch<- function(alpha, id, time="NA"){
  cluster<- cluster.size(id)
  n<- cluster$n
  time<- cluster$autotime
  
  n.max<- max(n)
  cor.max<- diag(1,n.max)
  cor.max[lower.tri(cor.max)]<- rep(alpha,n.max*(n.max-1)/2)
  cor.max[upper.tri(cor.max)]<- rep(alpha,n.max*(n.max-1)/2) 
  return(cor.max)
}

####################################################################

#Markov Structure

cormax.markov<- function(alpha,time){
  t.max<- max(time)
  cor.max<- diag(1,t.max)
  lowertri<- rep(0,0)
  for(j in (t.max-1):1){
    lowertri<- c(lowertri,1:j)
  }
  cor.max[lower.tri(cor.max)]<- alpha^lowertri
  cor.max[upper.tri(cor.max)]<- alpha^lowertri[length(lowertri):1]
  return(cor.max)
} 

####################################################################

#Tridiagonal Structure

cormax.tri<- function(alpha,id,time="NA"){
  cluster<- cluster.size(id)
  n<- cluster$n
  
  n.max<- max(n)
  cor.max<- diag(1,n.max)
  lowertri<- rep(0,0)
  for(j in (n.max-1):1){
    lowertri<- c(lowertri,alpha,rep(0,(j-1)))
  }
  cor.max[lower.tri(cor.max)]<- lowertri
  cor.max[upper.tri(cor.max)]<- lowertri[length(lowertri):1]
  return(cor.max)

}

#####################################################################

#Extended Familiar Structure

cormax.exfam<- function(alpha,id,time="NA"){
  cluster<- cluster.size(id)
  n<- cluster$n

  n.max<- max(n)

  ru1<- alpha[1]
  ru2<- alpha[2]
  gam<- alpha[3]
  alp<- alpha[4]
  
  m1<- cormax.ar1(gam,rep(1,2))
  m2<- cbind(rep(ru1,(n.max-2)),rep(ru2,(n.max-2)))
  cor.max<- rbind(m1,m2)
  m3<- cormax.exch(alp,rep(1,(n.max-2)))
  m4<- rbind(t(m2),m3)
  cor.max<- cbind(cor.max,m4)
  return(cor.max)
}

#################################################################################

#############################Creating zcor#######################################

gen.zcor<- function(cor.max,id,time="NA",markov=FALSE){
  
  cluster<- cluster.size(id)
  id<- cluster$id
  m<- cluster$m
  n<- cluster$n
  if((time[1]=="NA")||(markov==FALSE)){time=cluster$autotime}
  
  covariate.max<- diag(1,max(time)) 
  zcor<- rep(0,0)
  
  for(i in 1:m){
    #cat("i=",i,"\n")
    i.index<- which(id==i)
    t.i<- time[i.index]
    if(n[i]==1){
       zcor<- c(zcor,1[lower.tri(1)])
    }else{
      covariate.i<- covariate.max[t.i,]
      #cat("dim cov.i", dim(covariate.i),"\n")
      cor.i<- covariate.i%*%cor.max%*%t(covariate.i)
      #print(cor.i)
      zcor<- c(zcor,cor.i[lower.tri(cor.i)])      
    }
  }
  zcor
}

#################################################################################

##########################GEE for Different structures###########################

#AR1 Structure

gee.ar1.fixed<- 
function(formula,data,id,alpha,family="gaussian",time="NA",std.err=std.err){
  cor.max<- cormax.ar1(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                 family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

#################################################################################

#Exchangeable Structure

gee.exch.fixed<- 
function(formula,data,id,alpha,family="gaussian",time="NA",std.err=std.err){
  cor.max<- cormax.exch(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                 family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

#################################################################################

#Markov Structure
  
gee.markov.fixed<- 
function(formula,data,id,alpha,family="gaussian",time="NA",std.err=std.err){
  cor.max<- cormax.markov(alpha=alpha,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=TRUE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                 family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

#################################################################################

#Tridiagonal Structure

gee.tri.fixed<- 
function(formula, data,id, alpha, family="gaussian",time="NA",std.err=std.err){
  cor.max<- cormax.tri(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula, data=data,id=id,std.err=std.err,
                 family=family, corstr="fixed", zcor=zcor)
  return(geefit)
}

##################################################################################

#Extended Familiar Structure

gee.exfam.fixed<- function(formula,data,id,alpha,family="gaussian",time="NA",std.err=std.err){
  cor.max<- cormax.exfam(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                  family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

##################################################################################

#GEE for all Structures.

gee.fixed<- 
function(formula,data,id,alpha,family="gaussian",time="NA",correlation,std.err=std.err){
  switch(correlation,
         ar1=gee.ar1.fixed(formula=formula,data=data,std.err=std.err,
                           id=id,alpha=alpha,family=family,time=time),
         exchangeable= gee.exch.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                      alpha=alpha,family=family,time=time),
         markov= gee.markov.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                  alpha=alpha,family=family,time=time),
         tridiagonal= gee.tri.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                    alpha=alpha,family=family,time=time),
         ex.fam= gee.exfam.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                    alpha=alpha,family=family,time=time)
         )
}
