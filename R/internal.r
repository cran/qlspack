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


#This function will delete subjects with less or equal to #=del.n observations.

data.proc<- function(data,formula,time=NULL,id,del.n){
  
  dat<- data.frame(data)
  col.name<- names(dat)
  #cat("1\n")
  #print(dat)
  
  cluster<- cluster.size(id)
  m<- cluster$m
  n<- cluster$n
  id<- cluster$id
  if(length(time)==0){
    time<- cluster$autotime
  }
  autotime<- cluster$autotime
  index<- order(id,time)  
  #cat("index",index,"\n")
  #print(dat)
  #cat("ncol.dat",ncol(dat),"\n")
  if(ncol(dat)==1){
    dat<- dat[index,]
  }else{
    dat<- dat[index,]
  }
  dat<- data.frame(dat)
  names(dat)<- col.name
  
  
  del<- which(n<=del.n)
  if(length(del)>0){
    n<- n[-del]
    m<- length(n)
    mtch<- match(id,del)
    del.id<- which(mtch!="NA")
    #cat("ncol(dat)",ncol(dat),"\n")
    dat<- dat[-del.id,]
    dat<- data.frame(dat)
    names(dat)<- col.name
    row.names(dat)<- 1:nrow(dat)
    time<- time[-del.id]
    autotime<- autotime[-del.id]
    id<- rep(1:m,n)
  }
  
  formula<- as.formula(formula)
  fml<- as.formula(paste("~",formula[3],"+",formula[2],sep="")) 
  #print(fml)
  #print(dat)
  dat<- model.matrix(fml,data=dat)
  
  return(list(data=dat,time=time,autotime=autotime,id=id,m=m,n=n))
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

cormax.ar1<- function(alpha,id,time=NULL){
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

cormax.exch<- function(alpha, id, time=NULL){
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

cormax.tri<- function(alpha,id,time=NULL){
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

#Familiar Structure

cormax.fam<- function(alpha,id,time=NULL){
  cluster<- cluster.size(id)
  n<- cluster$n
  
  n.max<- max(n)
  r1<- alpha[1]
  r2<- alpha[2]
  cor.max<- matrix(r2,n.max,n.max)
  cor.max[,1]<- r1
  cor.max[1,]<- r1
  diag(cor.max)<- 1
  return(cor.max)
}



#####################################################################

#Extended Familiar Structure

cormax.exfam<- function(alpha,id,time=NULL){
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

gen.zcor<- function(cor.max,id,time=NULL,markov=FALSE){
  
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
function(formula,data,id,alpha,family="gaussian",time=NULL,std.err=std.err){
  cor.max<- cormax.ar1(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                 family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

#################################################################################

#Exchangeable Structure

gee.exch.fixed<- 
function(formula,data,id,alpha,family="gaussian",time=NULL,std.err=std.err){
  cor.max<- cormax.exch(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                 family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

#################################################################################

#Markov Structure
  
gee.markov.fixed<- 
function(formula,data,id,alpha,family="gaussian",time=NULL,std.err=std.err){
  cor.max<- cormax.markov(alpha=alpha,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=TRUE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                 family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

#################################################################################

#Tridiagonal Structure

gee.tri.fixed<- 
function(formula, data,id, alpha, family="gaussian",time=NULL,std.err=std.err){
  cor.max<- cormax.tri(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula, data=data,id=id,std.err=std.err,
                 family=family, corstr="fixed", zcor=zcor)
  return(geefit)
}

##################################################################################

#Familiar Structure

gee.fam.fixed<- 
function(formula, data,id, alpha, family="gaussian",time=NULL,std.err=std.err){
  cor.max<- cormax.fam(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula, data=data,id=id,std.err=std.err,
                 family=family, corstr="fixed", zcor=zcor)
  return(geefit)
}

###################################################################################

#Extended Familiar Structure

gee.exfam.fixed<- function(formula,data,id,alpha,family="gaussian",time=NULL,std.err=std.err){
  cor.max<- cormax.exfam(alpha=alpha,id=id,time=time)
  zcor<- gen.zcor(cor.max,id=id,time=time,markov=FALSE)
  geefit<- geeglm(formula=formula,data=data,id=id,std.err=std.err,
                  family=family,corstr="fixed",zcor=zcor)
  return(geefit)
}

##################################################################################

#GEE for all Structures.

gee.fixed<- 
function(formula,data,id,alpha,family="gaussian",time=NULL,correlation,std.err=std.err){
  switch(correlation,
         ar1=gee.ar1.fixed(formula=formula,data=data,std.err=std.err,
                           id=id,alpha=alpha,family=family,time=time),
         exchangeable= gee.exch.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                      alpha=alpha,family=family,time=time),
         markov= gee.markov.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                  alpha=alpha,family=family,time=time),
         tridiagonal= gee.tri.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                    alpha=alpha,family=family,time=time),
         fam = gee.fam.fixed(formula=formula,data=data,id=id,std.err=std.err,
                              alpha=alpha,family=family,time=time),
         ex.fam= gee.exfam.fixed(formula=formula,data=data,id=id,std.err=std.err,
                                    alpha=alpha,family=family,time=time)
         )
}
