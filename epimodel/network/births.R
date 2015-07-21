births <- function(dat, at) {
  #Variables
  birth.rate <- dat$param$birth.rate
  n <- network.size(dat$nw)
  
  births.mmu <- rpois(1, birth.rate/8)
  births.mmr <- rpois(1, birth.rate/8)
  births.nfu <- rpois(1, birth.rate/4)
  births.nfr <- rpois(1, birth.rate/4)
  births.nmu <- rpois(1, birth.rate/8)
  births.nmr <- rpois(1, birth.rate/8)
  births.total <- births.mmu + births.mmr + births.nfu + births.nfr + births.nmu + births.nmr
  
  #Introduce and activate new vertices
  #probably a faster way to do this but can't
  #figure out syntax
  for (i in 1:births.mmu) {
    dat$nw <- add.vertices(dat$nw, nv=1, vattr = 
                             list(list(sex=1, loc=0, mig_stat=1)))
  }
  for (i in 1:births.mmr) {
    dat$nw <- add.vertices(dat$nw, nv=1, vattr = 
                             list(list(sex=1, loc=1, mig_stat=1)))
  }
  for (i in 1:births.nfu) {
    dat$nw <- add.vertices(dat$nw, nv=1, vattr = 
                             list(list(sex=0, loc=0, mig_stat=0)))
  }
  for (i in 1:births.nfr) {
    dat$nw <- add.vertices(dat$nw, nv=1, vattr = 
                             list(list(sex=0, loc=1, mig_stat=0)))
  }
  for (i in 1:births.nmu) {
    dat$nw <- add.vertices(dat$nw, nv=1, vattr = 
                             list(list(sex=1, loc=0, mig_stat=0)))
  }
  for (i in 1:births.nmr) {
    dat$nw <- add.vertices(dat$nw, nv=1, vattr = 
                             list(list(sex=1, loc=1, mig_stat=0)))
  }
  
  newNodes <- (n+1):(n+births.total)
  dat$nw <- activate.vertices(dat$nw, onset=at, terminus=Inf, v=newNodes)
  
  #Set other vertex attributes
  if (births.total > 0) {
    #dat controlled attributes
    dat$attr$active <- c(dat$attr$active, rep(1, births.total))
    dat$attr$status <- c(dat$attr$status, rep("s", births.total))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, births.total))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, births.total))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, births.total))
    
    
  }
  
  return(dat)
  
}