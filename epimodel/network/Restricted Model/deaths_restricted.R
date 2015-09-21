deaths <- function (dat, at) {
  death.rate.gen <- dat$param$death.rate.gen
  death.rate.late <- dat$param$death.rate.aids + death.rate.gen
  late.cutoff <- dat$param$late.cutoff
  active.vertices <- which(dat$attr$active == 1)
  n <- length(active.vertices)
  late.stage.deaths <- 0
  death.vec <- NULL
  if (n > 0) {
    for (i in 1:n) {
      if (dat$attr$status[active.vertices[i]] == "s") {
        death.vec[i] <- rbinom(1, 1, death.rate.gen)
      } else {
        if ((at - dat$attr$infTime[active.vertices[i]]) > 552) {
          death.vec[i] <- 1
          late.stage.deaths <- late.stage.deaths + 1
        } else {
          if (((at - dat$attr$infTime[active.vertices[i]]) <= 552) && ((at - dat$attr$infTime[active.vertices[i]]) >= late.cutoff)) {
            death.vec[i] <- rbinom(1, 1, death.rate.late)
          } else {
            death.vec[i] <- rbinom(1, 1, death.rate.gen)
          }
        }
      }
    }
    
    #death.ids <- which(death.vec == 1)
    death.ids <- NULL
    for (i in 1:n) {
      if (death.vec[i] == 1) {
        death.ids <- c(death.ids, active.vertices[i])
      }
    }
    m <- length(death.ids)
    
    if (at == 2) {
      dat$epi$tot.deaths <- c(0, m)
      dat$epi$late.stage.deaths <- c(0, late.stage.deaths)
      dat$epi$natural.deaths <- c(0, m-late.stage.deaths)
    } else {
      dat$epi$tot.deaths[at] <- m
      dat$epi$late.stage.deaths[at] <- late.stage.deaths
      dat$epi$natural.deaths[at] <- m-late.stage.deaths
    }
    
    if (m > 0) {
      dat$attr$active[death.ids] <- 0
      dat$attr$exitTime[death.ids] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf, v = death.ids, deactivate.edges = TRUE)
    }
  }
  
  return(dat)
}