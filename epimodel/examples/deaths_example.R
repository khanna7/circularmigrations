deaths_example <- function(dat, at) {
  
  death.prob <- dat$param$death.prob
  active.vertices <- which(dat$attr$active == 1)
  num.elig <- length(active.vertices)
  
  if (num.elig > 0) {
    
    deaths <- which(rbinom(num.elig, 1, death.prob) == 1)
    deathIds <- active.vertices[deaths]
    num.deaths <- length(deathIds)
    
    if (num.deaths > 0) {
      dat$attr$active[deathIds] <- 0
      dat$attr$exitTime[deathIds] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = deathIds, deactivate.edges = TRUE)
    }
  }
  if (at == 2) {
    dat$epi$d.flow <- c(0, num.deaths)
  } else {
    dat$epi$d.flow[at] <- num.deaths
  }
  
  return(dat)
}