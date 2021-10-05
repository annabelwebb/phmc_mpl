basis_phmc = function(events, knots, basis="msplines", order=3, which=1){
  n=length(events)
  Alpha=knots$Alpha
  Delta=knots$Delta
  n.Alpha=length(Alpha)
  m=ifelse(basis=="msplines",n.Alpha+order-2,knots$m)
  M_Psi_nm=M_psi_nm=matrix(0,n,m)
  seq1n = 1:n
  Alpha_star=as.numeric(c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L)))
  M_psi_nm=M_Psi_nm=cbind(M_psi_nm,0) #why add a column?
  if(which==1){
    Alpha_star_x = sapply(events,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
    M_psi_nm[(Alpha_star_x-1L)*n+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
    if(order>1){
      for(ow in 2:order){
        uw_x = Alpha_star_x-ow+1
        for(pw in 0:(ow-1)){
          pos_x = (uw_x+pw-1L)*n+seq1n
          M_psi_nm[pos_x]=(ow/((ow-1)*(Alpha_star[1:m+ow]-Alpha_star[1:m])))[uw_x+pw]*((events-Alpha_star[uw_x+pw])*M_psi_nm[pos_x]+(Alpha_star[uw_x+pw+ow]-events)*M_psi_nm[pos_x+n])
        }
      }
    }
    M_psi_nm = M_psi_nm[,1:m,drop=FALSE]
    M_psi_nm
  }
  
  else{
    rank.events=rank(events)
    events=events[order(events)]
    Alpha_x=sapply(events,function(y,lim=Alpha[-1L])sum(lim<y)+1L)
    up_u     = cumsum(tabulate(Alpha_x,n.Alpha-1))
    for(uw in 1:(m-order+1)){M_Psi_nm[min(n,up_u[uw]+1):n,uw] = 1} 
    Alpha_star2 = c(rep(Alpha[1],order),Alpha,rep(Alpha[n.Alpha],order))    
    factor_v    = c((Alpha_star2[(order+2):length(Alpha_star2)]-Alpha_star2[1:(length(Alpha_star2)-order-1)])/
                      (order+1),rep(0,order-1))
    first=basis_phmc(events,knots,order=order+1,which=1)
    M_psi2_nm   = cbind(basis_phmc(events,knots,basis=basis,order=order+1,which=1),matrix(0,n,order-1))
    pos_xo  = rep((Alpha_x-1L)*n,1)+seq1n
    pos_xo1 = rep(pos_xo,order)+rep(1:order,each=n)*n
    for(ow in 0:(order-1)){
      M_Psi_nm[pos_xo+ow*n] = apply(matrix(M_psi2_nm[pos_xo1+ow*n]*
                                             factor_v[rep(Alpha_x,order)+rep((1:order)+ow,each=n)],ncol=order),1,sum)
    }
    M_Psi_nm = M_Psi_nm[rank.events,1:m,drop=FALSE]
    M_Psi_nm
  }
}


basis_phmc2 = function(events, knots, basis="msplines", order=3, which=1){
  n=length(events)
  Alpha=knots$Alpha
  #Delta=knots$Delta
  n.Alpha=length(Alpha)
  m=ifelse(basis=="msplines",n.Alpha+order-2,knots$m)
  M_Psi_nm=M_psi_nm=matrix(0,n,m)
  seq1n = 1:n
  Alpha_star=as.numeric(c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L)))
  M_psi_nm=M_Psi_nm=cbind(M_psi_nm,0) #why add a column?
  if(which==1){
    Alpha_star_x = sapply(events,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
    M_psi_nm[(Alpha_star_x-1L)*n+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
    if(order>1){
      for(ow in 2:order){
        uw_x = Alpha_star_x-ow+1
        for(pw in 0:(ow-1)){
          pos_x = (uw_x+pw-1L)*n+seq1n
          M_psi_nm[pos_x]=(ow/((ow-1)*(Alpha_star[1:m+ow]-Alpha_star[1:m])))[uw_x+pw]*((events-Alpha_star[uw_x+pw])*M_psi_nm[pos_x]+(Alpha_star[uw_x+pw+ow]-events)*M_psi_nm[pos_x+n])
          M_psi_nm=matrix(M_psi_nm,nrow=4)
        }
      }
    }
    M_psi_nm = matrix(M_psi_nm,nrow=4)
    M_psi_nm = M_psi_nm[,1:m,drop=FALSE]
    M_psi_nm
  }
  
  else{
    rank.events=rank(events)
    events=events[order(events)]
    Alpha_x=sapply(events,function(y,lim=Alpha[-1L])sum(lim<y)+1L)
    up_u     = cumsum(tabulate(Alpha_x,n.Alpha-1))
    for(uw in 1:(m-order+1)){M_Psi_nm[min(n,up_u[uw]+1):n,uw] = 1} 
    Alpha_star2 = c(rep(Alpha[1],order),Alpha,rep(Alpha[n.Alpha],order))    
    factor_v    = c((Alpha_star2[(order+2):length(Alpha_star2)]-Alpha_star2[1:(length(Alpha_star2)-order-1)])/
                      (order+1),rep(0,order-1))
    first=basis_phmc2(events,knots,order=order+1,which=1)
    M_psi2_nm   = cbind(basis_phmc2(events,knots,basis=basis,order=order+1,which=1),matrix(0,n,order-1))
    pos_xo  = rep((Alpha_x-1L)*n,1)+seq1n
    pos_xo1 = rep(pos_xo,order)+rep(1:order,each=n)*n
    for(ow in 0:(order-1)){
      M_Psi_nm[pos_xo+ow*n] = apply(matrix(as.numeric(M_psi2_nm[pos_xo1+ow*n])*
                                             as.numeric(factor_v[rep(Alpha_x,order)+rep((1:order)+ow,each=n)]),ncol=order),1,sum)
    }
    M_Psi_nm = M_Psi_nm[rank.events,1:m,drop=FALSE]
    M_Psi_nm
  }
}

