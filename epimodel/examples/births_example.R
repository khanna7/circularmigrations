births_example <- function(dat, at) {
  birth.rate <- dat$param$birth.rate
  nBirths <- rpois(1, birth.rate)
  n <- network.size(dat$nw)
  if (nBirths > 0) {    
    dat$nw <- add.vertices(dat$nw, nv = nBirths)
    newNodes <- (n + 1):(n + nBirths)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
  }
  
  # Update attributes
  if (nBirths > 0) {
    dat$attr$active <- c(dat$attr$active, rep(1, nBirths))
    dat$attr$status <- c(dat$attr$status, rep("s", nBirths))  
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))
    dat$attr$age <- c(dat$attr$age, rep(18, nBirths))
  }
  
  # Summary statistics
  if (at == 2) {
    dat$epi$b.flow <- c(0, nBirths)
  } else {
    dat$epi$b.flow[at] <- nBirths
  }
  
  return(dat)
}