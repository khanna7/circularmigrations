initialize.net.mig <- function (x, param, init, control, s) 
{
  if (control$start == 1) {
    dat <- list()
    dat$param <- param
    dat$init <- init
    dat$control <- control
    dat$attr <- list()
    dat$stats <- list()
    dat$temp <- list()
    if (class(x$fit) == "network") {
      nw <- simulate(x$formation, basis = x$fit, coef = x$coef.form.crude, 
                     constraints = x$constraints)
    }
    else {
      nw <- simulate(x$fit)
    }
    modes <- ifelse(nw %n% "bipartite", 2, 1)
    if (control$depend == TRUE) {
      if (class(x$fit) == "stergm") {
        nw <- network.collapse(nw, at = 1)
      }
      nw <- sim_nets(x, nw, nsteps = 1, control)
    }
    if (control$depend == FALSE) {
      nw <- sim_nets(x, nw, nsteps = control$nsteps, control)
    }
    nw <- activate.vertices(nw, onset = 1, terminus = Inf)
    dat$nw <- nw
    dat$nwparam <- list(x[-which(names(x) == "fit")])
    dat$param$modes <- modes
    dat <- init_status.net(dat)
    if (control$use.pids == TRUE) {
      dat$nw <- init_pids(dat$nw, dat$control$pid.prefix)
    }
    form <- get_nwparam(dat)$formation
    fterms <- get_formula_terms(form)
    dat <- copy_toall_attr(dat, at = 1, fterms)
    dat$temp$t1.tab <- get_attr_prop(dat$nw, fterms)
    dat <- get_prev.net(dat, at = 1)
  }
  else {
    dat <- list()
    dat$nw <- x$network[[s]]
    dat$param <- x$param
    dat$control <- control
    dat$nwparam <- x$nwparam
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$stats <- sapply(x$stats, function(var) var[[s]])
    dat$temp <- list()
  }
  
  #Nathan's modifications for nonzero inital time since infection
  
  if (control$start == 1) {
  
  init.infected.urban.women <- sample(2501:3750, 100, replace = FALSE)
  init.infected.rural.women <- sample(3751:5000, 100, replace = FALSE)
  init.infected.total <- c(init.infected.urban.women, init.infected.rural.women)
  init.infTime <- runif(length(init.infected.total), min = 0, max = 552)
  for (i in 1:200) {
    dat$attr$status[init.infected.total[i]] <- "i"
    dat$attr$infTime[init.infected.total[i]] <- -init.infTime[i]
  }
  
  }
  
  
  return(dat)
}