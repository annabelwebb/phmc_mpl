
knots_phmc = function(t1,t2,censor,basis="msplines",n.knots=c(8), range.quant=c(0.025,0.90),order=3){
  
  events=c((t1[censor[,4]]+t2[censor[,4]])/2,t1[censor[,4]],t2[censor[,4]],t2[censor[,3]]/2,t2[censor[,3]],t1[censor[,2]])
  internal=quantile(events,seq(range.quant[1],range.quant[2],length.out = n.knots))
  external=quantile(c(t1[censor[,2]]-1e-3,t1[censor[,2]]+1e-3,t1[censor[,1]],t1[censor[,4]],t2[censor[,4]],t2[censor[,3]]),seq(0,1,length.out = 3))
  Alpha=c(external[1]+1e-4,internal,external[3]+1e-4)
  n.Alpha  = length(Alpha)
  if(basis=="msplines"){
    m = n.Alpha+order-2
    list(m=m, Alpha=Alpha, Delta=rep(1,m))
    #t2[censor[,3]]
  }else{
    stop("Choose msplines basis.")
  }
}


knots_phmc2 = function(events,basis="msplines",n.knots=c(8), range.quant=c(0.075,0.9),order=3){
  #events=c((t1[censor[,4]]+t2[censor[,4]])/2,t1[censor[,3]]/2,t1[censor[,1]])
  #internal=quantile(events,seq(range.quant[1],range.quant[2],length.out = n.knots))
  #external=quantile(c(t1,t2),seq(0,1,length.out = 2))
  Alpha=quantile(events,seq(0,1,length.out=(n.knots[1]+2)))
  #Alpha=c(external[1],internal,external[2])
  n.Alpha  = length(Alpha)
  if(basis=="msplines"){
    m = n.Alpha+order-2
    list(m=m, Alpha=Alpha, Delta=rep(1,m))
  }else{
    stop("Choose msplines basis.")
  }
}

knots_phmc3 = function(t1,t2,censor,basis="msplines",n.knots=c(8),range.quant=c(0.01,0.9),order=3){
  events=c((t1[censor[,4]]+t2[censor[,4]])/2,t2[censor[,3]]/2,t1[censor[,2]])
  main=quantile(events,seq(range.quant[1],range.quant[2],length.out = n.knots))

  
  
  transitionpoint=max(events)
  external=quantile(c(t1[censor[,2]]-1e-3,t1[censor[,2]]+1e-3,t1[censor[,1]],t1[censor[,4]],t2[censor[,4]],t2[censor[,3]]),seq(0,1,length.out = 3))

  
  tail_range=c(max(main),transitionpoint,max(external))
  tail=quantile(tail_range,c(0.2,0.5))
  
  Alpha=c(external[1]+1e-4,main,tail,external[3]+1e-4)
  n.Alpha  = length(Alpha)
  
  #list(transitionpoint,max_interval,max(main),max(external))
  if(basis=="msplines"){
    m = n.Alpha+order-2
    list(m=m, Alpha=Alpha, Delta=rep(1,m))
  }else{
    stop("Choose msplines basis.")
  }
}



