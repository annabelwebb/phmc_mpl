
penalty_phmc = function(knots, basis="msplines", order=3){
  m = knots$m
  R = matrix(0,m,m)
  Alpha = knots$Alpha
  n.Alpha=length(Alpha)
  Alpha_star = c(rep(Alpha[1],order-1),Alpha,rep(Alpha[n.Alpha],order-1))
  if(basis=="msplines"){
    seq1n=1:(n.Alpha-1)
    n.Alpha_star=length(Alpha_star)
    Alpha_star_x=sapply(Alpha[-1],function(y,lim=Alpha[-1])sum(lim<y)+1)+order-1
    M_d2f_mm = matrix(0,n.Alpha-1,n.Alpha+order-1L)
    M_d2f_mm[(Alpha_star_x-1L)*(n.Alpha-1)+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
    for(ow in 2L:order){
      pw   = 1L:ow 
      uw_x = Alpha_star_x-ow+1L
      for(pw in 0:(ow-1L)){
        M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]=
          (ow/(Alpha_star[1:(n.Alpha+ow)+ow]-Alpha_star[1:(n.Alpha+ow)]))[uw_x+pw]*
          (M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]-M_d2f_mm[(uw_x+pw)*(n.Alpha-1)+seq1n])
      }
    }
    M_d2f_mm = M_d2f_mm[,1:m,drop=FALSE]
    for(uw in 1:m){
      for(vw in uw:m){
        R[uw,vw] = R[vw,uw] = 
          sum((M_d2f_mm[,uw]*M_d2f_mm[,vw])*(Alpha[-1]-Alpha[-n.Alpha]))
      }
    }
  }else{stop("Choose msplines basis.")}
  R
}

