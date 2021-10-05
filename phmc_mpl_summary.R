summary.phmc_mpl=function(obj,se="M2HM2",full=FALSE, ORHR = TRUE){
  col.names = c("Estimate","Std. Error", "z-value", "Pr(>|z|)")
  col.names.OR = c("Odds Ratio","Lower CI", "Upper CI", "Pr(>|z|)")
  col.names.HR = c("Hazard Ratio","Lower CI", "Upper CI", "Pr(>|z|)")
  p=obj$dimensions$p
  q=obj$dimensions$q
  m=obj$dimensions$m
  if(se=="M2HM2"){
    seB=obj$se$se_H[1:p]
    seG=obj$se$se_H[(p+1):(p+q)]
    seT=obj$se$se_H[(p+q+1):(p+q+m)]
    
    matxB = cbind(obj$beta,seB,obj$beta/seB,2*(1-pnorm(abs(obj$beta/seB))))
    colnames(matxB)=col.names
    rownames(matxB)=colnames(obj$data$Z)
    matOR = cbind(exp(obj$beta), exp(obj$beta - 1.96*seB),exp(obj$beta + 1.96*seB), 2*(1-pnorm(abs(obj$beta/seB))))
    colnames(matOR)=col.names.OR
    rownames(matOR)=colnames(obj$data$Z)
    
    matxG = cbind(obj$gamma,seG,obj$gamma/seG,2*(1-pnorm(abs(obj$gamma/seG))))
    colnames(matxG)=col.names
    rownames(matxG)=colnames(obj$data$X)
    matHR = cbind(exp(obj$gamma), exp(obj$gamma - 1.96*seG),exp(obj$gamma + 1.96*seG), 2*(1-pnorm(abs(obj$gamma/seG))))
    colnames(matHR)=col.names.HR
    rownames(matHR)=colnames(obj$data$X)
    
    matxT = cbind(obj$theta,seT,obj$theta/seT,2*(1-pnorm(abs(obj$theta/seT))))
    colnames(matxT)=col.names
  }
  if(se=="M2QM2"){
    seB=obj$se$se_Q[1:p]
    seG=obj$se$se_Q[(p+1):(p+q)]
    seT=obj$se$se_Q[(p+q+1):(p+q+m)]
    matxB = cbind(obj$beta,seB,obj$beta/seB,2*(1-pnorm(abs(obj$beta/seB))))
    colnames(matxB)=col.names
    rownames(matxB)=colnames(obj$data$Z)
    matOR = cbind(exp(obj$beta), exp(obj$beta - 1.96*seB),exp(obj$beta + 1.96*seB), 2*(1-pnorm(abs(obj$beta/seB))))
    colnames(matOR)=col.names.OR
    rownames(matOR)=colnames(obj$data$Z)
    
    matxG = cbind(obj$gamma,seG,obj$gamma/seG,2*(1-pnorm(abs(obj$gamma/seG))))
    colnames(matxG)=col.names
    rownames(matxG)=colnames(obj$data$X)
    matHR = cbind(exp(obj$gamma), exp(obj$gamma - 1.96*seG),exp(obj$gamma + 1.96*seG), 2*(1-pnorm(abs(obj$gamma/seG))))
    colnames(matHR)=col.names.HR
    rownames(matHR)=colnames(obj$data$X)
    
    matxT = cbind(obj$theta,seT,obj$theta/seT,2*(1-pnorm(abs(obj$theta/seT))))
    colnames(matxT)=col.names
  }
  
  out=list(Beta=matxB,Gamma=matxG,Theta=matxT, OR=matOR, HR=matHR, ORHR = ORHR, inf=list(call=obj$call,full=full,data=obj$data,ploglik=obj$ploglik,convergence=obj$convergence,smooth=obj$smooth,active=obj$constraint_info))
  class(out) = "summary.phmc_mpl"
  out

}

print.phmc_mpl=function(x){
  cat("\nLog-likelihood : ",x$ploglik[1],"\n",sep="")
  cat("\nSmoothing value: ",x$smooth[1],"\n" ,sep = "")
  cat("\nConvergence: ",x$convergence[4],"\n",sep="")
  cat("\nLogistic regression parameters :\n")
  vect=c(x$Beta)
  print(vect)
  cat("\nProportional hazards regression parameters :\n")
  vectg=c(x$Gamma)
  print(vectg)
  cat("\nBaseline hazard parameters : \n")
  print(x$Theta)
  cat("\n")
}

print.summary.phmc_mpl=function(x,...){
  inf = x$inf
  print(inf$call)
  cat("\n-----\n")

  cat("Mixture Cure Proportional Hazards Model Fitted Using MPL","\n")
  cat("Penalized log-likelihood  :  ",inf$ploglik,"\n",sep="")
  cat(ifelse(inf$convergence[2]==1,
             "Fixed smoothing value     :  ",
             "Estimated smoothing value :  "),
      inf$smooth,"\n",sep="")
  cat(ifelse(inf$convergence[1]==1,
               "Convergence : Yes, ",
               "Convergence : No, "),
      inf$convergence[4]," iterations\n",sep="")
  
  #cat("Data : ",inf$data$name,"\n",sep="")  
  #cat("No. of obs. : ",length(inf$data$time),"\n",sep="")  
  #cat("No. of events : ",sum(inf$data$censoring==1),"\n",sep="")  
  #cat("No. of right cens. : ",sum(inf$data$censoring==0),"\n",sep="") 
  #cat("No. of left cens. : ",sum(inf$data$censoring==2),"\n",sep="") 
  #cat("No. of interval cens. : ",sum(inf$data$censoring==3),"\n",sep="") 
  
  cat("\n-----\n")
  
  if(x$ORHR == TRUE){
    cat("\nIncidence odds ratios : \n",sep="")
    printCoefmat(x$OR, P.values=TRUE, has.Pvalue=TRUE)
    
    cat("\nLatency hazard ratios : \n",sep="")
    printCoefmat(x$HR, P.values=TRUE, has.Pvalue=TRUE) 
    
  }else{
    cat("Logistic regression parameters : \n",sep="")
    printCoefmat(x$Beta, P.values=TRUE, has.Pvalue=TRUE)
    
    cat("\nProportional hazards regression parameters : \n",sep="")
    printCoefmat(x$Gamma, P.values=TRUE, has.Pvalue=TRUE) 
    
  }
  
  
  cat("\n-----\n")
  
  if(inf$full){
    cat("\nBaseline hazard parameter vector : \n",sep="")
    printCoefmat(x$Theta, P.values=TRUE, has.Pvalue=TRUE)
    cat("\nActive constraint in theta : ")
    cat(ifelse(is.null(inf$active$ac),
           "None",
           inf$active$ac),"\n",sep="")
  }
}




#coef functions
coef.summary.phmc_mpl=function(obj,parameters="Gamma"){
  obj[[parameters]]
}

coef.phmc_mpl=function(obj,parameters="gamma"){
  obj[[parameters]]
}

#plot function
#baseline hazard, cumulative baseline hazard, baseline survival
plot.phmc_mpl=function(x,se="M2HM2",ask=TRUE,which=1:3){
  which.plot=rep(TRUE,3)
  if(!is.null(which)){which.plot[-which]=FALSE}
  if(sum(which.plot)==1){ask=FALSE}
  if(ask){oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  }
  
  control=x$control
  knots=x$knots
  pos=x$theta<x$control$min.theta
  n.x=1000
  V_x_X = seq(x$knots$Alpha[1],max(x$knots$Alpha),length=n.x)
  pq=(x$dimensions$p+x$dimensions$q)
  covar=x$covar[[se]]
  t.covar=covar[-c(1:pq),-c(1:pq)]
  
  anyplot=function(j,se,V_x_X,knots=x$knots){
    if(j==1){
      Ppsi=basis_phmc(V_x_X,knots)
    }
    if(j>1){
      Ppsi=basis_phmc(V_x_X,knots,which=2)
    }
    var.Hh0=diag(Ppsi%*%(t.covar)%*%t(Ppsi))
    pos.var=var.Hh0>0
    Hh0=c(Ppsi%*%matrix(x$theta,ncol=1))[pos.var]
    V_x_X=V_x_X[pos.var]
    sd.Hh0=sqrt(var.Hh0[pos.var])
    upper=Hh0+1.96*sd.Hh0
    lower=Hh0-1.96*sd.Hh0
    lower[lower<0]=0
    if(j==3){
      Hh0=exp(-Hh0)
      upper=exp(-upper)
      lower=exp(-lower)
    }
    xlim=range(knots$Alpha)
    ylim=c(ifelse(j<3,0,min(upper[V_x_X<xlim[2]])),ifelse(j<3,max(upper[V_x_X<xlim[2]]),max(lower[V_x_X<xlim[2]])))
    plot(1,1,xlim=xlim,ylim=ylim,main=paste("Estimate of the",c(" baseline hazard"," cumulative baseline hazard"," baseline survival")[j]," function",sep=""),
         xlab="Survival time",ylab=c(expression(h[0]*(t)),expression(H[0]*(t)),expression(S[0]*(t)))[j],type="n")
    lines(V_x_X,Hh0,lwd=1.1)
    lines(V_x_X,lower,col="grey")
    lines(V_x_X,upper,col="grey")
    
    }
  
  if(which.plot[1]){anyplot(1,se=se,V_x_X)}
  if(which.plot[2]){anyplot(2,se=se,V_x_X)}
  if(which.plot[3]){anyplot(3,se=se,V_x_X)}
  
}



#predict
#hazard or survival function
predict.phmc_mpl = function(object,se="M2QM2",type="hazard",cov=NULL,i=NULL,time=NULL,prob=0.95){
  
  beta=object$beta
  gamma=object$gamma
  theta=object$theta
  p=object$dimensions$p
  q=object$dimensions$q
  m=object$dimensions$m
  covar=object$covar[[se]]
  if(length(i)>1){warning("Only the first observation will be considered.\n",call.=FALSE)}
  if(is.null(time)){
    n.x=1000
    V_x_X=seq(object$knots$Alpha[1],max(object$knots$Alpha),length=n.x)
  }else{
    n.x=length(time)
    V_x_X = time
  }
  
  if(is.null(i)&is.null(cov)){
    xT=object$data$X
  }else if(!is.null(i) & is.null(cov)){
    xT=object$data$X[i[1],]
  }else{
    xT=matrix(cov,nrow=1)
  }
  
  
  Mu=c(exp(xT%*%gamma))
  
  out=data.frame(time=V_x_X, mid=NA, se=NA, low=NA, high=NA)
  
  if(type=="hazard"){
    psi=basis_phmc(V_x_X,knots=object$knots,order=object$control$order,which=1)
    h0=psi%*%theta
    out$mid=Mu*h0
    #th_se=sqrt(diag(covar))[(p+q+1):(p+q+m)]
    out$se=sqrt(diag(psi%*%(covar[(p+q+1):(p+q+m),(p+q+1):(p+q+m)])%*%t(psi)))
    out$low=out$mid-1.96*out$se
    out$low[out$low<0]=0
    out$high=out$mid+1.96*out$se
  }else{
    Psi=basis_phmc(V_x_X,knots=object$knots,order=object$control$order,which=2)
    H0=Psi%*%theta
    out$mid=exp(-Mu*H0)
    out$se=sqrt(diag((matrix(rep((exp(-H0)),m),ncol=m)*(-Psi))%*%(covar[(p+q+1):(p+q+m),(p+q+1):(p+q+m)])%*%t(matrix(rep((exp(-H0)),m),ncol=m)*(-Psi))))
    out$low=out$mid-1.96*out$se
    out$low[out$low<0]=0
    out$high=out$mid+1.96*out$se
  }
  
  times=c(object$data$time[,1])
  attributes(out)$inf=list(i=i[1],user.time=!is.null(time),prob=prob,upper.value=quantile(times,prob),max=max(times),m=m,risk=(type=="hazard"))
  colnames(out)[2]=type
  class(out)=c("predict.phmc_mpl","data.frame")
  out
}


phmc_piecewiseCI = function(obj, z_pos, x_pos, lvls = c(0,1), CI = TRUE){
  X_temp_save = NULL
  
  if(is.null(x_pos)){
    cov_obj = apply(obj$data$X,2,mean)
    predict_obj = predict(obj, type = "survival", cov = cov_obj)
    time = predict_obj$time
    s = rep(predict_obj$survival,length(lvls))
    s = matrix(s, nrow = 1000, ncol = length(lvls), byrow = FALSE)
    X_temp_save = matrix(rep(cov_obj,length(lvls)), nrow = length(lvls), byrow = TRUE)
  }else{
    s = matrix(0,nrow = 1000, ncol = length(lvls))
    for(l in 1:length(lvls)){
      cov_obj = apply(obj$data$X,2,mean)
      cov_obj[x_pos] = lvls[l]
      predict_obj = predict(obj, type = "survival", cov = cov_obj)
      time = predict_obj$time
      s[,l] = predict_obj$survival
      X_temp_save = rbind(X_temp_save, cov_obj)
      
    }
  }
  pop.survival = s
  Z_temp_save = NULL
  if(is.null(z_pos)){
    Z_tmp = matrix(apply(obj$data$Z,2,mean),nrow=1)
    ZB = Z_tmp%*%obj$beta
    pi = as.numeric(exp(ZB)/(1+exp(ZB)))
    pi_save = rep(pi, length(lvls))
    pop.survival = pi*s + (1 - pi)
    Z_temp_save = matrix(rep(Z_tmp,length(lvls)), nrow = length(lvls), byrow = TRUE)
    
  }else{
    pi_save = NULL
    for(l in 1:length(lvls)){
      Z_tmp = matrix(apply(obj$data$Z,2,mean),nrow=1)
      Z_tmp[z_pos] = lvls[l]
      ZB = Z_tmp%*%obj$beta
      pi = as.numeric(exp(ZB)/(1+exp(ZB)))
      pi_save = c(pi_save, pi)
      pop.survival[,l] = pi*s[,l] + (1 - pi)
      Z_temp_save = rbind(Z_temp_save, Z_tmp)
    }
    
  }
  
  CI_save = matrix(0, nrow = 1000, ncol = length(lvls)*2)
  
  for(l in 1:length(lvls)){
    X_temp = X_temp_save[l,]
    Z_temp = Z_temp_save[l,]
    pi = pi_save[l]
    mu_X=exp(X_temp%*%obj$gamma)
    cov_eta=obj$covar$M2QM2
    
    UL_population.survival.logit=rep(0,length(time))
    LL_population.survival.logit=rep(0,length(time))
    
    for(p in 2:length(time)){
      t=time[p]
      St=s[p,l]
      logitSpop=log(pop.survival[p,l]/(1-pop.survival[p,l]))
      Psi_t=basis_phmc(t,knots=obj$knots,order=obj$control$order,which=2)
      
      dlogitSt_dbeta=solve(pop.survival[p,l]*(1-pop.survival[p,l]))%*%(St-1)%*%pi%*%(1-pi)%*%Z_temp
      
      dlogitSt_dgamma=-solve(pop.survival[p,l]*(1-pop.survival[p,l]))%*%(pi%*%St%*%(-log(St)))%*%X_temp
      
      dlogitSt_dtheta=-mu_X%*%pi%*%St%*%solve(pop.survival[p,l]*(1-pop.survival[p,l]))%*%Psi_t
      
      dlogitSt_deta=matrix(c(dlogitSt_dbeta,dlogitSt_dgamma,dlogitSt_dtheta),nrow=1)
      
      var_logitSt=dlogitSt_deta%*%cov_eta%*%t(dlogitSt_deta)
      se_logitSt=sqrt(var_logitSt)
      
      
      
      LL_population.survival.logit[p]=logitSpop-1.96*se_logitSt
      UL_population.survival.logit[p]=logitSpop+1.96*se_logitSt
      
    }
    
    LL=exp(LL_population.survival.logit)/(1+exp(LL_population.survival.logit))
    UL=exp(UL_population.survival.logit)/(1+exp(UL_population.survival.logit))
    
    CI_save[,l]=LL
    CI_save[,(l+2)]=UL
    
  }
  
  
  out = list(pop.survival = pop.survival, CI_save = CI_save, time = time)
  return(out)
  
}



