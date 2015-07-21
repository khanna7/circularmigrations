deaths <- function(dat, at) {
  #Parameters
  death.rate.gen <- dat$param$death.rate.gen
  death.rate.late <- dat$param$death.rate.aids + death.rate.gen
  late.cutoff <- dat$param$late.cutoff 
  n <- network.size(dat$nw)
  
  #Decide who dies
  death.vec <- NULL
  for (i in 1:n) {
    if (dat$attr$infTime[i] < late.cutoff) {
      death.vec[i] <- rbinom(1, 1, death.rate.gen)
    }
    else {
      death.vec[i] <- rbinom(1, 1, death.rate.late)
    }
  }
  
  death.ids <- which(death.vec == 1)
  dat$nw <- delete.vertices(dat$nw, death.ids)
  dat$attr$active <- dat$attr$active[-1*death.ids]
  dat$attr$status <- dat$attr$status[-1*death.ids]
  dat$attr$infTime <- dat$attr$infTime[-1*death.ids]
  dat$attr$entrTime <- dat$attr$entrTime[-1*death.ids]
  dat$attr$exitTime <- dat$attr$exitTime[-1*death.ids]
  
  return(dat)
   
  
}