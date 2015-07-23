#EpiModel code for infection module with my minor edits to take
#spacial separation into account

infection <- function (dat, at) 
{
  active <- dat$attr$active
  status <- dat$attr$status
  modes <- dat$param$modes
  mode <- idmode(dat$nw)
  inf.prob <- dat$param$inf.prob
  inf.prob.m2 <- dat$param$inf.prob.m2
  act.rate <- dat$param$act.rate
  nw <- dat$nw
  tea.status <- dat$control$tea.status
  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  nInf <- nInfM2 <- totInf <- 0
  if (nElig > 0 && nElig < nActive) {
    
    #del returns a data frame with columns time step, id of susceptible
    #partner, and id of infected partner in discordant partnerships
    
    del <- discord_edgelist(dat, idsInf, idsSus, at)
    
    #delete rows of del corresponding to partnerships in which the 
    #two members are in separate locations.
    #Begin Nathan's edits
    
    n <- length(del$sus)
    if (n > 0) {
    invalid.edges <- NULL 
    for (i in 1:n) {
      if (get.vertex.attribute(nw, "loc")[del$sus[i]] != 
            get.vertex.attribute(nw, "loc")[del$inf[i]]) {
        invalid.edges <- c(invalid.edges, i)
      }
    }
    #invalid.edges <- -invalid.edges
    n <- length(invalid.edges)
    #this is not elegant, but don't know how to 
    #delete multiple rows at once
    for (i in 1:n) {
      del <- del[-invalid.edges[i],]
      invalid.edges <- c(invalid.edges[1:i], invalid.edges[(i+1):n]-1)
    }
    #End Nathan's edits
    }
    if (!(is.null(del))) {
      del$infDur <- at - dat$attr$infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1
      linf.prob <- length(inf.prob)
      if (is.null(inf.prob.m2)) {
        del$transProb <- ifelse(del$infDur <= linf.prob, 
                                inf.prob[del$infDur], inf.prob[linf.prob])
      }
      else {
        del$transProb <- ifelse(del$sus <= nw %n% "bipartite", 
                                ifelse(del$infDur <= linf.prob, inf.prob[del$infDur], 
                                       inf.prob[linf.prob]), ifelse(del$infDur <= 
                                                                      linf.prob, inf.prob.m2[del$infDur], inf.prob.m2[linf.prob]))
      }
      if (!is.null(dat$param$inter.eff) && at >= dat$param$inter.start) {
        del$transProb <- del$transProb * (1 - dat$param$inter.eff)
      }
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate, act.rate[del$infDur], 
                            act.rate[lact.rate])
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      nInf <- sum(mode[idsNewInf] == 1)
      nInfM2 <- sum(mode[idsNewInf] == 2)
      totInf <- nInf + nInfM2
      if (totInf > 0) {
        if (tea.status == TRUE) {
          nw <- activate.vertex.attribute(nw, prefix = "testatus", 
                                          value = "i", onset = at, terminus = Inf, 
                                          v = idsNewInf)
        }
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
        form <- get_nwparam(dat)$formation
        fterms <- get_formula_terms(form)
        if ("status" %in% fterms) {
          nw <- set.vertex.attribute(nw, "status", dat$attr$status)
        }
      }
      if (any(names(nw$gal) %in% "vertex.pid")) {
        del$sus <- get.vertex.pid(nw, del$sus)
        del$inf <- get.vertex.pid(nw, del$inf)
      }
    }
  }
  if (totInf > 0) {
    del <- del[!duplicated(del$sus), ]
    if (at == 2) {
      dat$stats$transmat <- del
    }
    else {
      dat$stats$transmat <- rbind(dat$stats$transmat, del)
    }
  }
  if (at == 2) {
    dat$epi$si.flow <- c(0, nInf)
    if (modes == 2) {
      dat$epi$si.flow.m2 <- c(0, nInfM2)
    }
  }
  else {
    dat$epi$si.flow[at] <- nInf
    if (modes == 2) {
      dat$epi$si.flow.m2[at] <- nInfM2
    }
  }
  dat$nw <- nw
  return(dat)
}