

phmc_mpl.control=function(basis="msplines",maxIter=c(150,75e+02,1e+04),epsilon=c(1e-16, 1e-10),kappa=1/0.6,conv_limit=1e-6,
                          smooth=0,min.theta=1e-10,n.knots=c(8),order=3,range.quant=c(0.075,0.9),
                          ties="epsilon",seed=NULL,new_constraint=0, new_criteria=0,new_knots=0){
  
  basis="msplines"
  basis = basis.name_mpl(basis)
  maxIter = c(ifelse(is.null(smooth),ifelse(maxIter[1]>0,as.integer(maxIter[1]),1.5e+2),1L),
              ifelse(maxIter[2]>0,as.integer(maxIter[2]),7.5e+2),
              ifelse(length(maxIter)==2,1e+6,
                     ifelse(maxIter[3]>ifelse(maxIter[2]>0,as.integer(maxIter[2]),7.5e+2),
                            as.integer(maxIter[3]),1e+4))) 
  
  penalty=order-1
  conv_limit=ifelse(conv_limit>0 & conv_limit<1,conv_limit,1e-6)
  if(is.null(n.knots)|sum(n.knots)<1|length(n.knots)!=1){
    n.knots    = if(basis!='msplines'){stop("Choose msplines basis.")}else{c(8)}
  }
  if(all(range.quant<=1) & all(range.quant>=0) & length(range.quant)==2){
    range.quant = range.quant[order(range.quant)]
  }else{range.quant = c(0.075,.9)}
  min.theta    = ifelse(min.theta>0 & min.theta<1e-3,min.theta,1e-10)
  order        = ifelse(order>0 & order<6,as.integer(order),3L)    
  kappa        = ifelse(kappa>1, kappa, 1/.6)
  if(!is.null(smooth)){
    smooth = ifelse(smooth<0,0,smooth)
  }else{smooth=0}
  
  out = list(basis = basis, smooth = smooth, maxIter = maxIter, conv_limit = conv_limit,
             order = order, penalty = penalty, n.knots = n.knots, range.quant = range.quant,
             min.theta = min.theta, ties = ties, seed = as.integer(seed), kappa = kappa, 
             epsilon = epsilon, new_constraint=new_constraint,new_criteria=new_criteria,new_knots=new_knots)
  
  class(out) = "phmc_mpl.control"
  out
}

basis.name_mpl <- function(k){
  if(k == "m" | k == "msplines" | k == "mspline")
  {"msplines"
  }else{
    stop("Unkown basis choice", call. = FALSE)
  }
}

penalty.order_mpl <- function(p,basis,order){
  p = as.integer(p)
  switch(basis,
         'msplines' = order-1)
}



phmc_mpl = function(ph.formula,pi.formula,data,control,...){
  mc = match.call(expand.dots = FALSE)
  m=match(c("ph.formula","pi.formula","data"), names(mc),0)
  mc = mc[c(1,m)]
  if(m[1]==0){stop("A Cox regression formula argument is required.")}
  data.name = if(m[3]!=0){deparse(match.call()[[4]])}else{"-"}      
  form=lapply(list(mc$ph.formula,mc$pi.formula),as.formula)
  ph_mf=model.frame(form[[1]])
  pi_mf=model.frame(form[[2]])
  
  if(any(is.na(ph_mf))) stop("Missing observations in the proportional hazards regression model variables.")
  if(any(is.na(pi_mf))) stop("Missing observations in the logistic regression model variables.")
  if(nrow(ph_mf)==0) stop("No non-missing observations.")
  
  surv = model.extract(ph_mf,"response")
  type=attr(surv,"class")
  if(!inherits(surv, "Surv")) {stop("Response must be a survival object.")}
  if(attr(surv,which="type")=="right"){
    left=surv[,1]
    right=rep(NA,nrow(surv))
    icase=which(surv[,2]==1)
    right[icase]=surv[icase,1]
    surv=Surv(left,right,type="interval2")
  } #else if(type=="interval"){
  #stop("\nPlease create the survival object using the option type='interval2' in the Surv function.\n")
  #}
  t_1 = surv[,1L]
  t_2 = surv[,2L]
  n = length(t_1)
  
  censor = matrix(NA, nrow=n, ncol=4)
  for(ct in 1:4){
    censor[,ct]=(surv[,3L]==(ct-1))
  }
  n.censor=apply(censor,2,sum)
  censorTF=n.censor>0   
  
  
  # control arguments
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(phmc_mpl.control)) 
    m <- pmatch(names(extraArgs), controlargs, nomatch=0L)
    if (any(m==0L))
      stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m==0L]),
           domain = NA, call. = F)
  }    
  if (missing(control)) control <- phmc_mpl.control(...)
  
  
  # ties 
  t_1.obs  = t_1[censor[,2]]   
  ties     = duplicated(t_1.obs)
  if(any(ties)){
    if(control$ties=="epsilon"){
      if(length(control$seed)>0){
        old <- .Random.seed
        on.exit({.Random.seed <<- old})
        set.seed(control$seed)
      }
      t_1.obs[ties] = t_1.obs[ties]+runif(sum(ties),-1e-11,1e-11)
      t_1[censor[,2]] = t_1.obs
    }else{    
      t_1.obs = t_1.obs[!ties]
      n.obs   = length(t_1.obs)
    }
  }
  
  
  #create X
  ph_mt=attr(ph_mf,"terms")
  X=model.matrix(ph_mt,data)
  X=X[,!apply(X, 2, function(x) all(x==x[1])), drop=FALSE]
  
  #create Z
  pi_mt=attr(pi_mf,"terms")
  if(ncol(pi_mf)==0){
    Z=model.matrix(pi_mt, data = data)
  }else{
    Z=model.matrix(pi_mt,data=pi_mf)
  }
  
  #centre X and Z
  mean_j=apply(X,2,mean)
  X=X-rep(mean_j,each=n)
  mean_t=apply(Z,2,mean)
  Z=cbind(Z[,1],Z[,-1]-rep(mean_t[-1],each=n))
  
  q = ncol(X) #proportional hazard regression number of parameters
  p = ncol(Z) #logistic regression number of parameters (including intercept term)
  
  #get knots sequence
  #knots=knots_phmc(c(t_1[censor[,4]],t_2[censor[,4]],t_1[censor[,2]]-1e-3,t_1[censor[,2]]+1e-3,t_1[censor[,1]],t_1[censor[,3]]),n.knots=control$n.knots)
  #t_t=(t_2[censor[,4]]+t_1[censor[,4]]/2)
  #t_t1=t_1[censor[,3]]/2
  #knots=knots_phmc(c((t_1[censor[,4]]+t_2[censor[,4]])/2,t_1[censor[,2]],t_1[censor[,3]]/2),n.knots=control$n.knots)
  #knots=knots_phmc(c(t_1[censor[,4]],t_2[censor[,4]],t_1[censor[,2]]-1e-3,t_1[censor[,2]]+1e-3,t_1[censor[,3]]),n.knots=control$n.knots)
  
  if(control$new_knots==1){
    knots=knots_phmc(t_1,t_2,censor=censor,n.knots=control$n.knots)
  }else if(control$new_knots==0){
    knots=knots_phmc2(c(t_1[censor[,4]],t_2[censor[,4]],t_1[censor[,2]]-1e-3,t_1[censor[,2]]+1e-3,t_1[censor[,1]],t_1[censor[,3]]),n.knots=control$n.knots)
  }else if(control$new_knots==2){
    knots=knots_phmc3(t_1,t_2,censor=censor,n.knots=control$n.knots)
  }
  
  
  m=knots$m
  lambda=control$smooth
  
  #logistic regression covariates
  Z_o = as.matrix(Z[censor[,2],]) #n1 x p
  Zt_o = t(Z_o) # p x n1
  Z_r = as.matrix(Z[censor[,1],]) #n2 x p
  Zt_r = t(Z_r) #p x n2
  Z_l = as.matrix(Z[censor[,3],]) #n3 x p
  Zt_l = t(Z_l) #p x n3
  Z_i = as.matrix(Z[censor[,4],]) #n4 x p
  Zt_i = t(Z_i) #p x n4
  
  #ph regression covariates
  X_o = as.matrix(X[censor[,2],]) #n1 x q
  Xt_o = t(X_o) # q x n1
  X_r = as.matrix(X[censor[,1],]) #n2 x q
  Xt_r = t(X_r) #q x n2
  X_l = as.matrix(X[censor[,3],]) #n3 x q
  Xt_l = t(X_l) #q x n3
  X_i = as.matrix(X[censor[,4],]) #n4 x q
  Xt_i = t(X_i) #q x n4
  
  #get R for penalty term
  R = penalty_phmc(knots)
  Rstar=rbind(matrix(0,p+q,p+q+m),cbind(matrix(0,m,p+q),R))
  
  #initialise parameter vectors
  beta = matrix(0,nrow=p,ncol=1)
  gamma = matrix(0,nrow=q,ncol=1)
  theta = matrix(1,nrow=m,ncol=1)
  df=-1
  epsilon=control$epsilon
  
  #get elements needed for updating beta when we go around again
  #observed
  pi_Z_o = exp(Z_o%*%beta)/(1+exp(Z_o%*%beta)) 
  #right censored
  pi_Z_r = exp(Z_r%*%beta)/(1+exp(Z_r%*%beta)) 
  Psi_r = basis_phmc(t_1, knots = knots,which=2)[censor[,1],]
  H0_r = Psi_r%*%theta
  mu_r = exp(X_r%*%gamma) 
  H_r = H0_r * mu_r 
  S_r = exp(-H_r)
  S_r[S_r==1]=1-epsilon[1]
  #left censored
  pi_Z_l = exp(Z_l%*%beta)/(1+exp(Z_l%*%beta)) 
  #interval censored
  pi_Z_i = exp(Z_i%*%beta)/(1+exp(Z_i%*%beta)) 
  
  
  #get values needed to calculate the log likelihood
  #observed
  psi_o = basis_phmc(t_1, knots = knots,which=1)[censor[,2],]
  h0_o = psi_o%*%theta
  Psi_o = basis_phmc(t_1, knots = knots,which=2)[censor[,2],]
  H0_o =  Psi_o%*%theta
  mu_o = exp(X_o%*%gamma)
  x_o.gamma=X_o%*%gamma
  H_o = H0_o * mu_o
  S_o = exp(-H_o)
  S_o[S_o==1]=1-epsilon[1]
  #right censored
  #left censored
  Psi_l = basis_phmc(t_1, knots = knots,which=2)[censor[,3],]
  H0_l =  Psi_l%*%theta
  mu_l = exp(X_l%*%gamma)
  H_l = H0_l*mu_l 
  S_l = exp(-H_l)
  S_l[S_l==1]=1-epsilon[1]
  #interval censored
  Psi1_i = basis_phmc(t_1, knots = knots,which=2)[censor[,4],]
  H01_i =  Psi1_i%*%theta
  Psi2_i = basis_phmc(t_2, knots = knots,which=2)[censor[,4],]
  H02_i =  Psi2_i%*%theta 
  mu_i = exp(X_i%*%gamma)
  H1_i = H01_i * mu_i
  S1_i = exp(-H1_i)
  H2_i = H02_i * mu_i
  S2_i = exp(-H2_i)
  S1S2_i = S1_i-S2_i
  S1_i[S1_i==1]=1-epsilon[1]
  S2_i[S2_i==1]=1-epsilon[1]
  S1S2_i = S1_i-S2_i  
  S1S2_i[S1S2_i<epsilon[2]]=epsilon[2]
  #penalty term
  Rtheta = R%*%theta
  thetaRtheta = t(theta)%*%Rtheta
  TwoLRtheta = lambda*2*Rtheta
  
  
  full.iter = 0
  
  
  for(iter in 1:control$maxIter[1]){
    #get initial log likelihood value to ensure that iterations increase it
    log_lik = sum(log(pi_Z_o)) + sum(log(h0_o)) +sum(x_o.gamma) + sum(log(S_o)) +
      sum(log(1-pi_Z_r+pi_Z_r*S_r)) + 
      sum(log(pi_Z_l)) + sum(log(1-S_l)) + 
      sum(log(pi_Z_i)) + sum(log(S1S2_i)) -
      lambda*thetaRtheta
    
    
    for(k in 1:control$maxIter[2]){
      beta_score = 
        Zt_o %*% (1-pi_Z_o) + 
        Zt_r %*% ((pi_Z_r*(1-pi_Z_r)*(S_r-1))/(1-pi_Z_r+pi_Z_r*S_r)) + 
        Zt_l %*% (1-pi_Z_l) + 
        Zt_i %*% (1-pi_Z_i) 
      beta_hessian = 
        -1*Zt_o %*% diag(c(pi_Z_o*(1-pi_Z_o))) %*% Z_o + 
        Zt_r %*% diag(c(((S_r-1)*((pi_Z_r*(1-pi_Z_r)^3 - pi_Z_r^3*(1-pi_Z_r)*S_r)/((1-pi_Z_r + pi_Z_r*S_r)^2))))) %*% Z_r - 
        Zt_l %*% diag(c((pi_Z_l*(1-pi_Z_l)))) %*% Z_l - 
        Zt_i %*% diag(c((pi_Z_i*(1-pi_Z_i)))) %*% Z_i 
      
      beta_OLD = beta 
      log_lik_OLD = log_lik
      omega1=1
      beta_update = solve(-beta_hessian)%*%beta_score
      beta = beta_OLD + omega1*beta_update
      
      #get elements needed to get updated log-likelihood (using new beta)
      #observed
      pi_Z_o = exp(Z_o%*%beta)/(1+exp(Z_o%*%beta))
      #right censored
      pi_Z_r = exp(Z_r%*%beta)/(1+exp(Z_r%*%beta))
      #left censored
      pi_Z_l = exp(Z_l%*%beta)/(1+exp(Z_l%*%beta))
      #interval censored
      pi_Z_i = exp(Z_i%*%beta)/(1+exp(Z_i%*%beta))
      
      log_lik = sum(log(pi_Z_o)) + sum(log(h0_o)) +sum(x_o.gamma) + sum(log(S_o)) +
        sum(log(1-pi_Z_r+pi_Z_r*S_r)) + 
        sum(log(pi_Z_l)) + sum(log(1-S_l)) + 
        sum(log(pi_Z_i)) + sum(log(S1S2_i)) -
        lambda*thetaRtheta
      
      #if likelihood decreases
      if(log_lik < log_lik_OLD){
        i = 0 #want to be able to restrict number of iterations this does
        omega1 = 1/control$kappa
        while(log_lik < log_lik_OLD){
          beta = beta_OLD + omega1*beta_update
          #get elements needed to get updated log-likelihood (using new beta)
          #observed
          pi_Z_o = exp(Z_o%*%beta)/(1+exp(Z_o%*%beta))
          #right censored
          pi_Z_r = exp(Z_r%*%beta)/(1+exp(Z_r%*%beta))
          #left censored
          pi_Z_l = exp(Z_l%*%beta)/(1+exp(Z_l%*%beta))
          #interval censored
          pi_Z_i = exp(Z_i%*%beta)/(1+exp(Z_i%*%beta))
          
          log_lik = sum(log(pi_Z_o)) + sum(log(h0_o)) +sum(x_o.gamma) + sum(log(S_o)) +
            sum(log(1-pi_Z_r+pi_Z_r*S_r)) + 
            sum(log(pi_Z_l)) + sum(log(1-S_l)) + 
            sum(log(pi_Z_i)) + sum(log(S1S2_i)) -
            lambda*thetaRtheta
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/control$kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>500){break}
        }
      }
      
      gamma_score = 
        Xt_o %*% (1-H_o) - #observed term ; q x 1
        Xt_r %*% ((pi_Z_r*S_r*H_r)/(1-pi_Z_r+pi_Z_r*S_r)) + 
        Xt_l %*% ((H_l*S_l)/(1-S_l)) -
        Xt_i %*% ((H1_i*S1_i-H2_i*S2_i)/S1S2_i)
      gamma_hessian = # q x q
        -1*Xt_o %*% diag(c(H_o),n.censor[2],n.censor[2]) %*% X_o -
        
        Xt_r %*% diag(c((pi_Z_r*S_r*H_r)/(1-pi_Z_r+pi_Z_r*S_r)),n.censor[1],n.censor[1]) %*% X_r +
        Xt_r %*% diag(c((pi_Z_r*(1-pi_Z_r)*S_r*(H_r^2))/((1-pi_Z_r+pi_Z_r*S_r)^2)),n.censor[1],n.censor[1]) %*% X_r -
        Xt_l %*% diag(c((S_l*(H_l^2))/((1-S_l)^2)),n.censor[3],n.censor[3]) %*% X_l +
        Xt_l %*% diag(c((H_l*S_l)/(1-S_l)),n.censor[3],n.censor[3]) %*% X_l -
        Xt_i %*% diag(c((S1_i*H1_i-S2_i*H2_i)/S1S2_i),n.censor[4],n.censor[4]) %*% X_i -
        Xt_i %*% diag(c((S1_i*S2_i*(H1_i-H2_i)^2)/(S1S2_i^2)),n.censor[4],n.censor[4]) %*% X_i
      
      gamma_OLD = gamma
      log_lik_OLD = log_lik
      omega2=1
      
      
      gamma_update = solve(-gamma_hessian)%*%gamma_score
      gamma = gamma_OLD + omega2*gamma_update
      
      #get elements needed to get updated log-likelihood (using new gamma)
      #observed
      x_o.gamma=X_o%*%gamma
      mu_o = exp(X_o%*%gamma)
      H_o = H0_o * mu_o
      S_o = exp(-H_o)
      S_o[S_o==1]=1-epsilon[1]
      #right censored
      mu_r = exp(X_r%*%gamma)
      H_r = H0_r * mu_r
      S_r = exp(-H_r)
      S_r[S_r==1]=1-epsilon[1]
      #left censored
      mu_l = exp(X_l%*%gamma)
      H_l = H0_l * mu_l
      S_l = exp(-H_l)
      S_l[S_l==1]=1-epsilon[1]
      #interval censored
      mu_i = exp(X_i%*%gamma)
      H1_i = H01_i * mu_i
      S1_i = exp(-H1_i)
      H2_i = H02_i * mu_i
      S2_i = exp(-H2_i)
      S1S2_i = S1_i-S2_i
      S1_i[S1_i==1]=1-epsilon[1]
      S2_i[S2_i==1]=1-epsilon[1]
      S1S2_i = S1_i-S2_i  
      S1S2_i[S1S2_i<epsilon[2]]=epsilon[2]
      
      log_lik = sum(log(pi_Z_o)) + sum(log(h0_o)) +sum(x_o.gamma) + sum(log(S_o)) +
        sum(log(1-pi_Z_r+pi_Z_r*S_r)) + 
        sum(log(pi_Z_l)) + sum(log(1-S_l)) + 
        sum(log(pi_Z_i)) + sum(log(S1S2_i)) -
        lambda*thetaRtheta
      
      if(log_lik < log_lik_OLD){
        i = 0
        omega2 = 1/control$kappa
        while(log_lik<log_lik_OLD){
          gamma = gamma_OLD + omega2*gamma_update
          #get elements needed to get updated log-likelihood (using new gamma)
          #observed
          x_o.gamma=X_o%*%gamma
          mu_o = exp(X_o%*%gamma)
          H_o = H0_o * mu_o
          S_o = exp(-H_o)
          S_o[S_o==1]=1-epsilon[1]
          #right censored
          mu_r = exp(X_r%*%gamma)
          H_r = H0_r * mu_r
          S_r = exp(-H_r)
          S_r[S_r==1]=1-epsilon[1]
          #left censored
          mu_l = exp(X_l%*%gamma)
          H_l = H0_l * mu_l
          S_l = exp(-H_l)
          S_l[S_l==1]=1-epsilon[1]
          #interval censored
          mu_i = exp(X_i%*%gamma)
          H1_i = H01_i * mu_i
          S1_i = exp(-H1_i)
          H2_i = H02_i * mu_i
          S2_i = exp(-H2_i)
          S1S2_i = S1_i-S2_i
          S1_i[S1_i==1]=1-epsilon[1]
          S2_i[S2_i==1]=1-epsilon[1]
          S1S2_i = S1_i-S2_i  
          S1S2_i[S1S2_i<epsilon[2]]=epsilon[2]
          
          log_lik = sum(log(pi_Z_o)) + sum(log(h0_o)) +sum(x_o.gamma) + sum(log(S_o)) +
            sum(log(1-pi_Z_r+pi_Z_r*S_r)) + 
            sum(log(pi_Z_l)) + sum(log(1-S_l)) + 
            sum(log(pi_Z_i)) + sum(log(S1S2_i)) -
            lambda*thetaRtheta
          
          #update value of omega2
          if(omega2>=1e-2){
            omega2 = omega2/control$kappa
          }else if(omega2<1e-2 & omega2>=1e-5){
            omega2 = omega2*(5e-2)
          }else if(omega2<1e-5){
            omega2 = omega2*(1e-5)
          }
          i = i+1
          if(i>500){break}
          
        }
      }
      
      #get elements needed for updating theta (use updated beta and gamma)
      #observed
      psit_o = t(psi_o) 
      Psit_o = t(Psi_o) 
      #right censored
      Psit_r = t(Psi_r) 
      mu_r = exp(X_r%*%gamma)
      H_r = H0_r * mu_r
      S_r = exp(-H_r) 
      S_r[S_r==1]=1-epsilon[1]
      #left censored
      Psit_l = t(Psi_l) 
      #interval censored
      Psi1t_i = t(Psi1_i)
      Psi2t_i = t(Psi2_i) 
      
      omega3=1
      theta_OLD = theta 
      log_lik_OLD = log_lik
      #positive elements of theta score vector minus negative elements of penalty first d
      gradtheta_A = psit_o%*%(1/h0_o) +
        Psit_l%*%((S_l*mu_l)/(1-S_l)) +
        Psi2t_i%*%((S2_i*mu_i)/(S1S2_i)) - 
        TwoLRtheta*(TwoLRtheta<0)+0.3
      #negative elements of theta score vector plus positive elements of penalty
      gradtheta_B = Psit_o%*%mu_o +
        Psit_r%*%((pi_Z_r*S_r*mu_r)/(1-pi_Z_r+pi_Z_r*S_r)) +
        Psi1t_i%*%((S1_i*mu_i)/(S1S2_i)) +
        TwoLRtheta*(TwoLRtheta>0)+0.3
      gradtheta=gradtheta_A-gradtheta_B
      s=theta/gradtheta_B
      theta_update=s*gradtheta
      theta = theta_OLD+omega3*theta_update
      theta[theta<epsilon[2]]=epsilon[2]
      
      
      #get elements needed to get updated log-likelihood (using new theta)
      #observed
      h0_o = psi_o%*%theta
      H0_o =  Psi_o%*%theta
      H_o = H0_o * mu_o
      S_o = exp(-H_o)
      S_o[S_o==1]=1-epsilon[1]
      #right censored
      H0_r = Psi_r%*%theta
      H_r = H0_r * mu_r 
      S_r = exp(-H_r)
      S_r[S_r==1]=1-epsilon[1]
      #left censored
      H0_l =  Psi_l%*%theta
      H_l = H0_l * mu_l 
      S_l = exp(-H_l)
      S_l[S_l==1]=1-epsilon[1]
      #interval censored
      H01_i =  Psi1_i%*%theta
      H02_i =  Psi2_i%*%theta 
      H1_i = H01_i * mu_i
      S1_i = exp(-H1_i)
      H2_i = H02_i * mu_i
      S2_i = exp(-H2_i)
      S1S2_i = S1_i-S2_i
      S1_i[S1_i==1] = 1-epsilon[1]
      S2_i[S2_i==1] = 1-epsilon[1]
      S1S2_i = S1_i-S2_i  
      S1S2_i[S1S2_i<epsilon[2]] = epsilon[2]
      #penalty term
      Rtheta=R%*%theta
      thetaRtheta = t(theta)%*%Rtheta
      
      log_lik = sum(log(pi_Z_o)) + sum(log(h0_o)) +sum(x_o.gamma) + sum(log(S_o)) +
        sum(log(1-pi_Z_r+pi_Z_r*S_r)) + 
        sum(log(pi_Z_l)) + sum(log(1-S_l)) + 
        sum(log(pi_Z_i)) + sum(log(S1S2_i)) -
        lambda*thetaRtheta
      
      if(log_lik<log_lik_OLD){
        i = 0
        omega3 = omega3/control$kappa
        while(log_lik<log_lik_OLD){
          theta = theta_OLD+omega3*theta_update
          theta[theta<epsilon[2]]=epsilon[2]
          #get elements needed to get updated log-likelihood (using new theta)
          #observed
          h0_o = psi_o%*%theta
          H0_o =  Psi_o%*%theta
          H_o = H0_o * mu_o
          S_o = exp(-H_o)
          S_o[S_o==1]=1-epsilon[1]
          #right censored
          H0_r = Psi_r%*%theta
          H_r = H0_r * mu_r 
          S_r = exp(-H_r)
          S_r[S_r==1]=1-epsilon[1]
          #left censored
          H0_l =  Psi_l%*%theta
          H_l = H0_l * mu_l 
          S_l = exp(-H_l)
          S_l[S_l==1]=1-epsilon[1]
          #interval censored
          H01_i =  Psi1_i%*%theta
          H02_i =  Psi2_i%*%theta 
          H1_i = H01_i * mu_i
          S1_i = exp(-H1_i)
          H2_i = H02_i * mu_i
          S2_i = exp(-H2_i)
          S1S2_i = S1_i-S2_i
          S1_i[S1_i==1] = 1-epsilon[1]
          S2_i[S2_i==1] = 1-epsilon[1]
          S1S2_i = S1_i-S2_i  
          S1S2_i[S1S2_i<epsilon[2]] = epsilon[2]
          #penalty term
          Rtheta=R%*%theta
          thetaRtheta = t(theta)%*%Rtheta
          TwoLRtheta = lambda*2*Rtheta
          
          log_lik = sum(log(pi_Z_o)) + sum(log(h0_o)) +sum(x_o.gamma) + sum(log(S_o)) +
            sum(log(1-pi_Z_r+pi_Z_r*S_r)) + 
            sum(log(pi_Z_l)) + sum(log(1-S_l)) + 
            sum(log(pi_Z_i)) + sum(log(S1S2_i)) -
            lambda*thetaRtheta
          
          #update value of omega3
          if(omega3>=1e-2){
            omega3 = omega3/control$kappa
          }else if(omega3<1e-2 & omega3>=1e-5){
            omega3 = omega3*(5e-2)
          }else if(omega3<1e-5){
            omega3 = omega3*(1e-5)
          }
          i = i+1
          if(i>500){
            break
          }
          
        }
      }
      activeconstraint=NULL
      for(t in 1:m){
        if(control$new_constraint==1){
          if((theta[t]<(1e-03) & gradtheta[t]<=(-1e-02)) | (theta[t]<(1e-05) & gradtheta[t]<0)){
            activeconstraint=c(activeconstraint,t)
          }
        }else{
          if((theta[t]<(1e-03) & gradtheta[t]<=(-1e-02))){
            activeconstraint=c(activeconstraint,t)
          }
        }
        
        
      }
      
      if(!is.null(activeconstraint)){
        if(control$new_criteria==1){
          if(all(c(abs(beta-beta_OLD),abs(gamma-gamma_OLD))<control$conv_limit) & all(((gradtheta_A[-activeconstraint]/gradtheta_B[-activeconstraint])-1)<(1e-3))){
            save=TRUE
            break
          }else{
            save=FALSE
            #print(full.iter+k)
          }
        }else{
          if(all(c(abs(beta-beta_OLD),abs(gamma-gamma_OLD),abs(theta-theta_OLD))<control$conv_limit)){
            save=TRUE
            break
          }else{
            save=FALSE
            #print(full.iter+k)
          }
        }
        
      }else{
        if(control$new_criteria==1){
          if(all(c(abs(beta-beta_OLD),abs(gamma-gamma_OLD))<control$conv_limit) & all(((gradtheta_A/gradtheta_B)-1)<(1e-3))){
            save=TRUE
            break
          }else{
            save=FALSE
            print(full.iter+k)
          }
        }else{
          if(all(c(abs(beta-beta_OLD),abs(gamma-gamma_OLD),abs(theta-theta_OLD))<control$conv_limit)){
            save=TRUE
            break
          }else{
            save=FALSE
            print(full.iter+k)
          }
        }
      }
      
      
      
      
      
      
      
    } #end inner loop
    #H matrix
    H=HRinv=matrix(0,p+q+m,p+q+m)
    H[1:p,1:p] = #beta beta
      Zt_o %*% diag(c(pi_Z_o*(1-pi_Z_o)),n.censor[2],n.censor[2]) %*% Z_o - 
      Zt_r %*% diag(c((S_r-1)*pi_Z_r*(1-pi_Z_r)*(((1-pi_Z_r)^2 - S_r*(pi_Z_r)^2)/((1-pi_Z_r + pi_Z_r*S_r)^2))),n.censor[1],n.censor[1]) %*% Z_r + 
      Zt_l %*% diag(c(pi_Z_l*(1-pi_Z_l)),n.censor[3],n.censor[3]) %*% Z_l + 
      Zt_i %*% diag(c(pi_Z_i*(1-pi_Z_i)),n.censor[4],n.censor[4]) %*% Z_i
    
    H[1:p,(p+1):(p+q)] = #beta gamma
      Zt_r %*% diag(c((H_r*S_r*pi_Z_r*(1-pi_Z_r))/((1-pi_Z_r+pi_Z_r*S_r)^2)),n.censor[1],n.censor[1])%*%X_r
    
    H[1:p,(p+q+1):(p+q+m)] = #beta theta
      Zt_r %*% diag(c((mu_r*S_r*pi_Z_r*(1-pi_Z_r))/((1-pi_Z_r+pi_Z_r*S_r)^2)),n.censor[1],n.censor[1])%*%Psi_r
    
    H[(p+1):(p+q),1:p] = t(H[1:p,(p+1):(p+q)]) #gamma beta
    
    H[(p+1):(p+q),(p+1):(p+q)] = #gamma gamma
      Xt_o %*% diag(c(H_o),n.censor[2],n.censor[2]) %*% X_o +
      Xt_r %*% diag(c((pi_Z_r*S_r*H_r)/(1-pi_Z_r+pi_Z_r*S_r)),n.censor[1],n.censor[1]) %*% X_r -
      Xt_r %*% diag(c((pi_Z_r*(1-pi_Z_r)*S_r*(H_r^2))/((1-pi_Z_r+pi_Z_r*S_r)^2)),n.censor[1],n.censor[1]) %*% X_r +
      Xt_l %*% diag(c(((S_l*(H_l^2))/((1-S_l)^2))),n.censor[3],n.censor[3]) %*% X_l -
      Xt_l %*% diag(c(H_l*S_l/(1-S_l)),n.censor[3],n.censor[3]) %*% X_l +
      Xt_i %*% diag(c((S1_i*H1_i-S2_i*H2_i)/S1S2_i),n.censor[4],n.censor[4]) %*% X_i +
      Xt_i %*% diag(c((S1_i*S2_i*((H1_i-H2_i)^2))/(S1S2_i^2)),n.censor[4],n.censor[4]) %*% X_i
  
     H[(p+1):(p+q),(p+q+1):(p+q+m)] = #gamma theta
      Xt_o%*%diag(c(mu_o),n.censor[2],n.censor[2])%*%Psi_o +
      Xt_r%*%diag(c(S_r*mu_r*(pi_Z_r*S_r + (1-pi_Z_r)*(1-H_r))/((1-pi_Z_r+pi_Z_r*S_r)^2)),n.censor[1],n.censor[1])%*%Psi_r +
      Xt_l%*%diag(c(S_l*mu_l*(H_l+S_l-1)/((1-S_l)^2)),n.censor[3],n.censor[3])%*%Psi_l +
      Xt_i%*%diag(c(mu_i*S1_i/S1S2_i),n.censor[4],n.censor[4])%*%Psi1_i -
      Xt_i%*%diag(c(mu_i*S2_i/S1S2_i),n.censor[4],n.censor[4])%*%Psi2_i +
      Xt_i%*%diag(c(mu_i*S1_i*S2_i*(H1_i-H2_i)/(S1S2_i^2)),n.censor[4],n.censor[4])%*%(Psi1_i-Psi2_i)
    
    H[(p+q+1):(p+q+m),1:p]=t(H[1:p,(p+q+1):(p+q+m)]) #theta beta
    
    H[(p+q+1):(p+q+m),(p+1):(p+q)]=t(H[(p+1):(p+q),(p+q+1):(p+q+m)]) #theta gamma
    
    H[(p+q+1):(p+q+m),(p+q+1):(p+q+m)]= #theta theta
      t(psi_o)%*%diag(c(1/(h0_o^2)),n.censor[2],n.censor[2])%*%psi_o -
      t(Psi_r)%*%diag(c((mu_r^2)*(S_r*(1-pi_Z_r)*pi_Z_r)/((1-pi_Z_r+S_r*pi_Z_r)^2)),n.censor[1],n.censor[1])%*%Psi_r +
      t(Psi_l)%*%diag(c((mu_l^2)*S_l/((1-S_l)^2)),n.censor[3],n.censor[3])%*%Psi_l +
      t(Psi1_i-Psi2_i)%*%diag(c((mu_i^2)*(S1_i*S2_i)/((S1S2_i)^2)),n.censor[4],n.censor[4])%*%(Psi1_i-Psi2_i)
    
    lambda_old   = lambda
    df_old       = df
    sigma2_old   = 1/(2*lambda_old)
    activeconstraint=NULL
    for(t in 1:m){
      if(control$new_constraint==1){
        if((theta[t]<(1e-03) & gradtheta[t]<=(-1e-02)) | (theta[t]<(1e-05) & gradtheta[t]<0)){
          activeconstraint=c(activeconstraint,t)
        }
      }else{
        if((theta[t]<(1e-03) & gradtheta[t]<=(-1e-02))){
          activeconstraint=c(activeconstraint,t)
        }
      }
      
      
    }
    
    pos=c(rep(TRUE,p+q+m))
    pos[p+q+activeconstraint]=FALSE
    
    #temp           = try(chol2inv(chol(H[pos,pos]+(1/sigma2_old)*Rstar[pos,pos])),silent=T)  
    #if(class(temp)!="try-error"&!any(is.infinite(temp))){
    #  HRinv[pos,pos]=temp
    #}else{HRinv[pos,pos]=MASS::ginv(H[pos,pos])}
    HRinv[pos,pos]=solve(H[pos,pos]+(1/sigma2_old)*Rstar[pos,pos])
    df           = m-sum(diag(HRinv%*%Rstar))/sigma2_old
    if(df<0){
      stop("Negative variance in estimation of smoothing parameter. Number of knots may be mis-specified.")
    }
    sigma2       = c(t(theta)%*%R%*%theta/df)    
    lambda       = 1/(2*sigma2)
    TwoLRtheta     = lambda*2*Rtheta 
    
    full.iter      = full.iter+k
    if(full.iter>control$maxIter[3]){break}
    if((k<control$maxIter[2])&
       (abs(df-df_old)<1)
    ){break}
    
    
  } #end outer loop
  
  lambda=lambda_old
  
  M_2=H+2*lambda*Rstar
  
  Q = matrix(NA,n,p+q+m)
  
  if(censorTF[1]){Q[censor[,1],1:p] = rep((pi_Z_r*(1-pi_Z_r)*(S_r-1))/(1-pi_Z_r+pi_Z_r*S_r),p)*Z_r}
  if(censorTF[2]){Q[censor[,2],1:p] = rep((1-pi_Z_o),p)*Z_o}
  if(censorTF[3]){Q[censor[,3],1:p] = rep((1-pi_Z_l),p)*Z_l}
  if(censorTF[4]){Q[censor[,4],1:p] = rep((1-pi_Z_i),p)*Z_i}
  
  if(censorTF[1]){Q[censor[,1],(p+1):(p+q)] = rep(-1*((pi_Z_r*S_r*H_r)/(1-pi_Z_r+pi_Z_r*S_r)),q)*X_r}
  if(censorTF[2]){Q[censor[,2],(p+1):(p+q)] = rep((1-H_o),q)*X_o}
  if(censorTF[3]){Q[censor[,3],(p+1):(p+q)] = rep(((H_l*S_l)/(1-S_l)),q)*X_l}
  if(censorTF[4]){Q[censor[,4],(p+1):(p+q)] = rep(-1*((H1_i*S1_i-H2_i*S2_i)/S1S2_i),q)*X_i}
  
  
  if(censorTF[1]){Q[censor[,1],-c(1:(p+q))] = rep(-((pi_Z_r*S_r*mu_r)/(1-pi_Z_r+pi_Z_r*S_r)),m)*Psi_r}
  if(censorTF[2]){Q[censor[,2],-c(1:(p+q))] = rep((1/h0_o),m)*psi_o-rep(mu_o,m)*Psi_o}
  if(censorTF[3]){Q[censor[,3],-c(1:(p+q))] = rep(((S_l*mu_l)/(1-S_l)),m)*Psi_l}
  if(censorTF[4]){Q[censor[,4],-c(1:(p+q))] = rep(((S2_i*mu_i)/(S1S2_i)),m)*Psi2_i-
    rep(((S1_i*mu_i)/(S1S2_i)),m)*Psi1_i}
  
  Sp = Q-matrix(rep(c(rep(0,p+q),TwoLRtheta),n),n,byrow=T)/n
  Q = t(Sp)%*%Sp
  
  Minv_2=corr=matrix(0,p+q+m,p+q+m)
  diag(corr)=rep(1,m+q+p)
  corr[!pos,]=0
  #temp = try(chol2inv(chol(M_2[pos,pos])),silent=T) 
  #if(class(temp)!="try-error"){
  #  Minv_2[pos,pos] = temp
  #  cov_H=corr%*%(Minv_2%*%H%*%Minv_2)%*%t(corr)
  #  cov_Q=corr%*%(Minv_2%*%Q%*%Minv_2)%*%t(corr)
  #  se_H=sqrt(diag(cov_H))
  #  se_Q=sqrt(diag(cov_Q))
    
  Minv_2[pos,pos]=solve(M_2[pos,pos])
  cov_H=corr%*%(Minv_2%*%H%*%Minv_2)%*%t(corr)
  cov_Q=corr%*%(Minv_2%*%Q%*%Minv_2)%*%t(corr)
  se_H=sqrt(diag(cov_H))
  se_Q=sqrt(diag(cov_Q))
  
  
  convergence=c(save,iter,k,full.iter)
  ZB = Z%*%beta
  pi_estimate = as.numeric(exp(ZB)/(1+exp(ZB)))
  
  #output
  fit=list(beta=beta,gamma=gamma,theta=theta*(theta>control$min.theta))
  fit$data=list(time = surv, censoring=surv[,3L], X = X, Z=Z,name = data.name)
  fit$pi_estimate = pi_estimate
  fit$knots=knots
  fit$call=match.call()
  fit$constraint_info=list(theta=theta,gradtheta=gradtheta,ac=activeconstraint)
  fit$beta_gradient=beta_score
  fit$gamma_gradient=gamma_score
  fit$covar=list(H=H,M2HM2=cov_H,M2QM2=cov_Q,M2=M_2)
  fit$se=list(se_H=se_H,se_Q=se_Q)
  fit$ploglik=log_lik
  fit$loglik=log_lik+lambda*thetaRtheta
  fit$dimensions=list(p=p, q=q, m=m)
  fit$smooth=lambda
  fit$convergence=convergence
  fit$control=control
  class(fit)="phmc_mpl"
  fit
  
}

