#####################Stage One and Stage Two Estimates########################

###################################AR1 Stucture#############################

#Stage One

ar1.one<- function(resid,time=NULL,id){
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(resid)),resid),formula=resid~resid,time=time,id=id,del.n=1)
  resid<- data1$data[,2]
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n

  a<- rep(0,m)
  b<- rep(0,m)
  c<- rep(0,m)
  d<- rep(0,m)
  z<- matrix(0,nrow=max(n),ncol=m)
  for(i in 1:m){
    z[1:n[i],i]<- resid[which(id==i)]
    a[i]<- sum(z[1:(n[i]-1),i]^2+z[2:n[i],i]^2)
    b[i]<- sum(z[1:(n[i]-1),i]*z[2:n[i],i])
    c[i]<- sum((z[1:(n[i]-1),i]-z[2:n[i],i])^2)
    d[i]<- sum((z[1:(n[i]-1),i]+z[2:n[i],i])^2)
  }
  a<- sum(a)
  b<- sum(b)
  c<- sum(c)
  d<- sum(d)
  alpha<- (a-sqrt(c*d))/(2*b)
  alpha
}

#Stage Two
ar1.two<- function(alpha,resid=NULL,time=NULL,id=NULL){
  2*alpha/(1+alpha^2)
}

#####################################################################

######################Exchageable Stucture###########################

#Stage One

exch.one<- function(resid,time=NULL,id){
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(resid)),resid),formula=resid~resid,time=time,id=id,del.n=1)
  resid<- data1$data[,2]
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n
 
  z<- matrix(0,nrow=max(n),ncol=m)
  a<- rep(0,m)
  b<- rep(0,m)
  for(i in 1:m){
    if(n[i]>1){
      z[1:n[i],i]<- resid[which(id==i)]
      a[i]<- sum(z[(1:n[i]),i]^2)
      b[i]<- sum(z[(1:n[i]),i])^2   
    }
  }
  f<- function(alpha){
    sum(a)-sum(b*(1+alpha^2*(n-1))/((1+alpha*(n-1))^2))
  }
  lower<- -1/(max(n)-1)+0.00001
  upper<- 0.99999
  alpha<- uniroot(f,c(lower,upper),tol=1e-10,maxiter=10000000)$root
  return(alpha)  
}

#Stage Two

exch.two<- function(alpha,resid=NULL,time=NULL,id){
  resid<- rep(1,length(id))
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(id)),resid),formula=resid~resid,time=NULL,id=id,del.n=1)
  resid<- data1$data
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n

  a<- sum(n*(n-1)*alpha*(alpha*(n-2)+2)/(1+alpha*(n-1))^2)
  b<- sum(n*(n-1)*(1+alpha^2*(n-1))/(1+alpha*(n-1))^2)
  a/b
}


########################################################################

############################Markov Structure############################

#Stage One
markov.one<- function(resid,time=NULL,id){
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(resid)),resid),formula=resid~resid,time=time,id=id,del.n=1)
  resid<- data1$data[,2]
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n

  z<- matrix(0,nrow=max(n),ncol=m)
  t<- matrix(0,nrow=max(n),ncol=m)
  for(i in 1:m){
    z[1:n[i],i]<- resid[which(id==i)]
    t[1:n[i],i]<- time[which(id==i)]
  }
  z2<- z[2:max(n),]
  z1<- z[1:(max(n)-1),]
  e<- abs(t[2:max(n),]-t[1:(max(n)-1),])
  if(max(n)==2){
    e<- t(as.matrix(e))
  }
  f<- function(alpha){
    temp<-e*(alpha^e)*(alpha^(2*e)*z2*z1-(alpha^e)*(z2^2+z1^2)+z2*z1)/(1-alpha^(2*e))^2
    x<- rep(0,m)
    for(i in 1:m){
      x[i]<- sum(temp[1:(n[i]-1),i])
    }
    return(sum(x))
  }
  alpha<- uniroot(f,c(0.00000001,0.99999999),tol=1e-10,maxiter=10000000)$root
  return(alpha)$root
}

#Stage Two
markov.two<- function(alpha,resid=NULL,time=NULL,id){
  
  resid<- rep(1,length(id))
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(id)),resid),formula=resid~resid,time=time,id=id,del.n=1)
  resid<- data1$data
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n
 
  t<- matrix(0,nrow=max(n),ncol=m)
  for(i in 1:m){
    t[1:n[i],i]<- time[which(id==i)]
  }
  e<- abs(t[2:max(n),]-t[1:(max(n)-1),])
  if(max(n)==2){
    e<- t(as.matrix(e))
  }
  g<- function(ru){
    temp<- (2*e*alpha^(2*e-1)-(ru^e)*e*(alpha^(e-1)+alpha^(3*e-1)))/(1-alpha^(2*e))^2
    x<- rep(0,m)
    for(i in 1:m){
      x[i]<- sum(temp[1:(n[i]-1),i])
    }
    return(sum(x))
  }
  return(uniroot(g,c(0.00001,0.99999),tol=1e-10,maxiter=10000000)$root)
}

#######################################################################

#########################Tridiagonal Structure#########################

#Stage One
tri.one<- function(resid,time=NULL,id){
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(resid)),resid),formula=resid~resid,time=time,id=id,del.n=1)
  resid<- data1$data[,2]
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n
  
  z<- matrix(0,nrow=max(n),ncol=m)
  funR<- function(a,i){
    R<- diag(1,n[i],n[i])
    dR<- matrix(0,n[i],n[i])
    for(j in 1:(n[i]-1)){
      R[j,j+1]<- a
      R[j+1,j]<- a
      dR[j,j+1]<- 1
      dR[j+1,j]<- 1
    }
    return(solve(R)%*%dR%*%solve(R))
  }
  for(i in 1:m){
    z[1:n[i],i]<- resid[which(id==i)]
  }
  f<- function(a){
    x<- rep(0,m)
    for(i in 1:m){
      x[i]<- t(z[1:n[i],i])%*%funR(a,i)%*%z[1:n[i],i]
    }
    return(sum(x))
  }
  lower<- -1/(2*sin(pi*(max(n)-1)/(2*(max(n)+1))))+0.000001
  upper<- 1/(2*sin(pi*(max(n)-1)/(2*(max(n)+1))))-0.000001
  return(uniroot(f,c(lower,upper),tol=1e-10,maxiter=10000000)$root) 
}

########################################################################

#Stage Two
tri.two<- function(alpha,resid=NULL,time=NULL,id){
  
  resid<- rep(1,length(id))
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(id)),resid),formula=resid~resid,time=time,id=id,del.n=1)
  resid<- data1$data
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n

  funRu<- function(ru,i){
    Ra<- diag(1,n[i],n[i])
    Rru<- diag(1,n[i],n[i])
    dR<- matrix(0,n[i],n[i])
    for(j in 1:(n[i]-1)){
      Ra[j,j+1]<- alpha
      Ra[j+1,j]<- alpha
      Rru[j,j+1]<- ru
      Rru[j+1,j]<- ru
      dR[j,j+1]<- 1
      dR[j+1,j]<- 1
    }
    return(sum(diag(solve(Ra)%*%dR%*%solve(Ra)%*%Rru)))
  }
  g<- function(ru){
    x<- rep(0,m)
    for(i in 1:m){
      x[i]<- funRu(ru,i)
    }
    return(sum(x))
  }
  lower<- -1/(2*sin(pi*(max(n)-1)/(2*(max(n)+1))))+0.00001
  upper<- 1/(2*sin(pi*(max(n)-1)/(2*(max(n)+1))))-0.00001
  return(uniroot(g,c(lower,upper),tol=1e-10,maxiter=10000000)$root) 
}

#########################################################################

####################Familiar Structure################################

fam.one<- function(resid,time=NULL,id){
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(resid)),resid),formula=resid~resid,time=NULL,id=id,del.n=2)
  resid<- data1$data[,2]
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n
  t<- max(n)-1
  
  z<- matrix(0,nrow=max(n),ncol=m)
  for(i in 1:m){
    z[1:n[i],i]<- resid[which(id==i)]
  }
  
  f<- function(ru){
    x<- rep(0,m)
    for(i in 1:m){
      corr<- cormax.fam(ru,id=rep(1,n[i]),time=NULL) 
      x[i]<- t(z[1:n[i],i])%*%solve(corr)%*%z[1:n[i],i]
    }
    return(sum(x))
  }
  g<- function(ru){
     r1<- ru[1]
     r2<- ru[2]
     if((r1^2<(1+(t-1)*r2)/t)&&(r2>-1/(t-1))&&(r2<1)){
       return(f(ru))
     }else{
       return(f(rep(0,2))+100000)}
  }          
  optim(rep(0,2),g,method="BFGS")$par
}

#Stage Two
fam.two<- function(alpha,time=NULL,id){
  
  resid<- rep(1,length(id))
  resid<- data.frame(resid)
  data<- data.proc(data=resid,formula=resid~resid,time=time,id=id,del.n=2)
  resid<- data$resid
  time<- data$autotime
  id<- data$id
  m<- data$m
  n<- data$n
  
  mat<- matrix(0,2,2)
  tr.tmp<- rep(0,m)
  for(j in 1:2){
    for(k in 1:2){
      for(i in 1:m){
        tmp<- cormax.fam(alpha,rep(1,n[i]))%*%
              (cormax.fam(diag(1,2,2)[j,],rep(1,n[i]))-cormax.fam(c(0,0),rep(1,n[i])))%*%
              cormax.fam(alpha,rep(1,n[i]))%*%
              (cormax.fam(diag(1,2)[k,],rep(1,n[i]))-cormax.fam(c(0,0),rep(1,n[i])))
        tr.tmp[i]<- sum(diag(tmp))        
      }
      mat[j,k]<- sum(tr.tmp)
    }
  }
  c<- rep(0,2)
  for(j in 1:2){
    for(i in 1:m){
      tmp<- cormax.fam(alpha,rep(1,n[i]))%*%
            (cormax.fam(diag(1,2,2)[j,],rep(1,n[i]))-cormax.fam(c(0,0),rep(1,n[i])))%*%
            solve(cormax.fam(alpha,rep(1,n[i])))%*%
            cormax.fam(c(0,0),rep(1,n[i]))
      tr.tmp[i]<- sum(diag(tmp))
    }
    c[j]<- -sum(tr.tmp)  
  }
  ru<- solve(mat,c)
  return(ru)
}



####################Extended Familiar Structure#######################

#Stage One Estimate

exfam.one<- function(resid,time=NULL,id){
  resid<- data.frame(resid)
  data1<- data.proc(data=cbind(rep(1,length(resid)),resid),formula=resid~resid,time=time,id=id,del.n=3)
  resid<- data1$data[,2]
  time<- data1$time
  id<- data1$id
  m<- data1$m
  n<- data1$n
  t<- max(n)-2
  
  z<- matrix(0,nrow=max(n),ncol=m)
  for(i in 1:m){
    z[1:n[i],i]<- resid[which(id==i)]
  }
  
  f<- function(ru){
    x<- rep(0,m)
    for(i in 1:m){
      corr<- cormax.exfam(ru,id=rep(1,n[i]),time=NULL) 
      x[i]<- t(z[1:n[i],i])%*%solve(corr)%*%z[1:n[i],i]
    }
    return(sum(x))
  }
  g<- function(ru){
     r1<- ru[1]
     r2<- ru[2]
     gam<- ru[3]
     alp<- ru[4]
     if((alp<1)&&(alp>-3/(t-1))&&(gam<1)&&(gam>-1)&&
        (t*(r1^2+r2^2)<3+2*(t-1)*alp-gam^2)&&
        (t*(r1^2+r2^2-2*gam*r1*r2)<(1-gam^2)*(1+(t-1)*alp))){
       return(f(ru))
     }else{
       return(f(rep(0,4))+100000)}
  }          
  optim(rep(0,4),g,method="BFGS")$par
}


############################################################

#Partial R Partial alp

palp1<- function(t){
  m1<- diag(0,2)
  m2<- rbind(rep(1,t),rep(0,t))
  m3<- cbind(rep(1,t),rep(0,t))
  m4<- diag(0,t)
  rbind(cbind(m1,m2),cbind(m3,m4))
}

palp2<- function(t){
  m1<- diag(0,2)
  m2<- rbind(rep(0,t),rep(1,t))
  m3<- cbind(rep(0,t),rep(1,t))
  m4<- diag(0,t)
  rbind(cbind(m1,m2),cbind(m3,m4))
}

palp3<- function(t){
  m1<- matrix(1,2,2)
  diag(m1)<- rep(0,2)
  m2<- matrix(0,2,t)
  m3<- matrix(0,t,2)
  m4<- diag(0,t)
  rbind(cbind(m1,m2),cbind(m3,m4))
}


palp4<- function(t){
  m1<- diag(0,2)
  m2<- matrix(0,2,t)
  m3<- matrix(0,t,2)
  m4<- matrix(1,t,t)
  diag(m4)<- rep(0,t)
  rbind(cbind(m1,m2),cbind(m3,m4))
}

palp<- function(t,j){
  if(j==1) return(palp1(t))
  if(j==2) return(palp2(t))
  if(j==3) return(palp3(t))
  if(j==4) return(palp4(t)) 
}

########################################################################

#Stage Two Estimate#

exfam.two<- function(alpha,id,time=NULL){
  
  resid<- rep(1,length(id))
  resid<- data.frame(resid)
  data<- data.proc(data=cbind(rep(1,length(id)),resid),formula=resid~resid,time=time,id=id,del.n=3)
  resid<- data$data
  time<- data$time
  id<- data$id
  m<- data$m
  n<- data$n
  t<- max(n)-2
  #cat("t=",t,"\n")
  #cat("m=",m,"\n")
  #cat("n=",n,"\n")
  
  M<- matrix(0,4,4)
  y<- rep(0,4)
  for(j in 1:4){
    x<- rep(0,m)
    for(i in 1:m){
      R<- cormax.exfam(alpha,id=rep(1,n[i]))
      A<- solve(R)%*%palp((n[i]-2),j)%*%solve(R)
      ru<- rep(0,4)
      x[i]<- sum(diag(A%*%cormax.exfam(ru,id=rep(1,n[i]))))
    }
    y[j]<- -sum(x[i])
    #cat("y",j,y[j],"\n")
    for(k in 1:4){
      for(i in 1:m){
        R<- cormax.exfam(alpha,id=rep(1,n[i]))
        A<- solve(R)%*%palp((n[i]-2),j)%*%solve(R)
        ru<- rep(0,4)
        ru[k]<- 1
        x[i]<- sum(diag(A%*%cormax.exfam(ru,id=rep(1,n[i]))))
      }
      M[j,k]<- sum(x[i])+y[j]
      #cat("m",j,k,M[j,k],"\n")
    }
  }
  return(solve(M,y))
}
