get_prev_example_4 <- function (dat, at) 
{
  active <- dat$attr$active
  modes <- dat$param$modes
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active == 
                                                              1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL
  status <- l$status
  if (modes == 2) {
    mode <- idmode(dat$nw)[active == 1]
  }
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- dat$control$epi.by
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }
  if (modes == 1) {
    if (at == 1) {
      dat$epi <- list()
      dat$epi$s.num <- sum(status == "s")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == 
                                                       "s" & get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num <- sum(status == "i")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num", ebun[i])]] <- sum(status == 
                                                       "i" & get(ebn) == ebv[i])
        }
      }
      if (dat$control$type == "SIR") {
        dat$epi$r.num <- sum(status == "r")
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("r.num", ebun[i])]] <- sum(status == 
                                                         "r" & get(ebn) == ebv[i])
          }
        }
      }
      dat$epi$num <- length(status)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num", ebun[i])]] <- sum(get(ebn) == 
                                                     ebv[i])
        }
      }
    }
    else {
      dat$epi$s.num[at] <- sum(status == "s")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == 
                                                           "s" & get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num[at] <- sum(status == "i")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num", ebun[i])]][at] <- sum(status == 
                                                           "i" & get(ebn) == ebv[i])
        }
      }
      if (dat$control$type == "SIR") {
        dat$epi$r.num[at] <- sum(status == "r")
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("r.num", ebun[i])]][at] <- sum(status == 
                                                             "r" & get(ebn) == ebv[i])
          }
        }
      }
      dat$epi$num[at] <- length(status)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num", ebun[i])]][at] <- sum(get(ebn) == 
                                                         ebv[i])
        }
      }
    }
  }
  else {
    if (at == 1) {
      dat$epi <- list()
      dat$epi$s.num <- sum(status == "s" & mode == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == 
                                                       "s" & mode == 1 & get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num <- sum(status == "i" & mode == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num", ebun[i])]] <- sum(status == 
                                                       "i" & mode == 1 & get(ebn) == ebv[i])
        }
      }
      if (dat$control$type == "SIR") {
        dat$epi$r.num <- sum(status == "r" & mode == 
                               1)
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == 
                                                         "r" & mode == 1 & get(ebn) == ebv[i])
          }
        }
      }
      dat$epi$num <- sum(mode == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num", ebun[i])]] <- sum(mode == 
                                                     1 & get(ebn) == ebv[i])
        }
      }
      dat$epi$s.num.m2 <- sum(status == "s" & mode == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num.m2", ebun[i])]] <- sum(status == 
                                                          "s" & mode == 2 & get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num.m2 <- sum(status == "i" & mode == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num.m2", ebun[i])]] <- sum(status == 
                                                          "i" & mode == 2 & get(ebn) == ebv[i])
        }
      }
      if (dat$control$type == "SIR") {
        dat$epi$r.num.m2 <- sum(status == "r" & mode == 
                                  2)
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("r.num.m2", ebun[i])]] <- sum(status == 
                                                            "r" & mode == 2 & get(ebn) == ebv[i])
          }
        }
      }
      dat$epi$num.m2 <- sum(mode == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num.m2", ebun[i])]] <- sum(mode == 
                                                        2 & get(ebn) == ebv[i])
        }
      }
    }
    else {
      dat$epi$s.num[at] <- sum(status == "s" & mode == 
                                 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == 
                                                           "s" & mode == 1 & get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num[at] <- sum(status == "i" & mode == 
                                 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num", ebun[i])]][at] <- sum(status == 
                                                           "i" & mode == 1 & get(ebn) == ebv[i])
        }
      }
      if (dat$control$type == "SIR") {
        dat$epi$r.num[at] <- sum(status == "r" & mode == 
                                   1)
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == 
                                                             "r" & mode == 1 & get(ebn) == ebv[i])
          }
        }
      }
      dat$epi$num[at] <- sum(mode == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num", ebun[i])]][at] <- sum(mode == 
                                                         1 & get(ebn) == ebv[i])
        }
      }
      dat$epi$s.num.m2[at] <- sum(status == "s" & mode == 
                                    2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num.m2", ebun[i])]][at] <- sum(status == 
                                                              "s" & mode == 2 & get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num.m2[at] <- sum(status == "i" & mode == 
                                    2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num.m2", ebun[i])]][at] <- sum(status == 
                                                              "i" & mode == 2 & get(ebn) == ebv[i])
        }
      }
      if (dat$control$type == "SIR") {
        dat$epi$r.num.m2[at] <- sum(status == "r" & mode == 
                                      2)
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("r.num.m2", ebun[i])]][at] <- sum(status == 
                                                                "r" & mode == 2 & get(ebn) == ebv[i])
          }
        }
      }
      dat$epi$num.m2[at] <- sum(mode == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num.m2", ebun[i])]][at] <- sum(mode == 
                                                            2 & get(ebn) == ebv[i])
        }
      }
    }
  }
  
  if (at == 1) {
    dat$epi$edge_count <- network.edgecount(network.extract(dat$nw, at = 1))
    dat$epi$size <- network.size(network.extract(dat$nw, at = 1))
    dat$epi$mean_deg <- 2*dat$epi$edge_count/length(which(dat$attr$active == 1))
    dat$epi$prev.as.percent <- dat$epi$i.num[at]/dat$epi$size
    dat$epi$mixing_matrix_urban <- mixingmatrix(network.extract(dat$nw, at = 1), "rest.urb")$matrix[1, 1]
  } else {
    dat$epi$edge_count[at] <- network.edgecount(network.extract(dat$nw, at = at))
    dat$epi$size[at] <- network.size(network.extract(dat$nw, at = at))
    dat$epi$mean_deg[at] <- 2*dat$epi$edge_count[at]/length(which(dat$attr$active == 1))
    dat$epi$prev.as.percent[at] <- dat$epi$i.num[at]/dat$epi$size[at]
    dat$epi$mixing_matrix_urban[at] <- mixingmatrix(network.extract(dat$nw, at=at), "rest.urb")$matrix[1, 1]
  }
  
  active.vertices <- which(active == 1)
  mig.men <- intersect(active.vertices, which(get.vertex.attribute(dat$nw, "type") == "B-MM"))
  non.mig.men.rural <- intersect(active.vertices, which(get.vertex.attribute(dat$nw, "type") == "A-MR"))
  non.mig.men.urban <- intersect(active.vertices, which(get.vertex.attribute(dat$nw, "type") == "C-MU"))
  females <- intersect(active.vertices, which(get.vertex.attribute(dat$nw, "sex") == 0))
  infected <- which(dat$attr$status == "i")
  if (at == 1) {
    dat$epi$mig.men <- length(mig.men)
    dat$epi$non.mig.men <- length(non.mig.men.rural) + length(non.mig.men.urban)
    dat$epi$women <- length(females)
    dat$epi$inf.women <- length(intersect(which(dat$attr$status == "i"), females))
    dat$epi$inf.non.mig.men <- length(intersect(infected, non.mig.men.rural)) + length(intersect(infected, non.mig.men.urban))
    dat$epi$inf.mig.men <- length(intersect(infected, mig.men))
  } else {
    dat$epi$mig.men[at] <- length(mig.men)
    dat$epi$non.mig.men[at] <- length(non.mig.men.rural) + length(non.mig.men.urban)
    dat$epi$women[at] <- length(females)
    dat$epi$inf.women[at] <- length(intersect(which(dat$attr$status == "i"), females))
    dat$epi$inf.non.mig.men[at] <- length(intersect(infected, non.mig.men.rural)) + length(intersect(infected, non.mig.men.urban))
    dat$epi$inf.mig.men[at] <- length(intersect(infected, mig.men))
  }
  
  #output <- data.frame(#time.step = c(1:at),
  #prev = dat$epi$prev.as.percent,
  #pop = dat$epi$size,
  #mean_deg = dat$epi$mean_deg,
  #incidence = dat$epi$si.flow,
  #births = dat$epi$b.total,
  #deaths.tot = dat$epi$tot.deaths,
  #deaths.late = dat$epi$late.stage.deaths,
  #deaths.natural = dat$epi$natural.deaths,
  #mig.men = dat$epi$mig.men,
  #non.mig.men = dat$non.mig.men,
  #women = dat$epi$women,
  #infected = dat$epi$i.num) 
  #inf.women = dat$epi$inf.women,
  #inf.mig.men = dat$epi$inf.mig.men,
  #inf.non.mig.men = dat$epi$inf.non.mig.men)
  
  #write.csv(output, "output.csv")
  #time.step <- data.frame(time.step = c(1:at))
  #write.csv(time.step, "time.step.csv")
  
#   prev <- data.frame(prev = dat$epi$prev.as.percent)
#   write.csv(prev, filename=paste(nsim_, "prev_4.csv", sep=""))
#   
#   pop <- data.frame(pop = dat$epi$size)
#   write.csv(pop, filename=paste(nsim_, "pop_4.csv", sep=""))
#   
#   mean_deg <- data.frame(mean_deg = dat$epi$mean_deg)
#   write.csv(mean_deg, filename=paste(nsim_, "mean_deg_4.csv", sep=""))
#   
#   incidence <- data.frame(incidence = dat$epi$si.flow)
#   write.csv(incidence, filename=paste(nsim_, "incidence_4.csv", sep=""))
#   
#   births <- data.frame(births = dat$epi$b.total)
#   write.csv(births, filename=paste(nsim_, "births_4.csv", sep=""))
#   
#   deaths.tot <- data.frame(deaths.tot = dat$epi$tot.deaths)
#   write.csv(deaths.tot, filename=paste(nsim_, "deaths_tot_4.csv", sep=""))
#   
#   late.stage.deaths <- data.frame(late.stage.deaths = dat$epi$late.stage.deaths)
#   write.csv(late.stage.deaths, filename=paste(nsim_, "late_deaths_4.csv", sep=""))
  
  #natural.deaths <- data.frame(natural.deaths = dat$epi$natural.deaths)
  #write.csv(natural.deaths, "natural_deaths.csv")
  
  #mig.men <- data.frame(mig.men = dat$epi$mig.men)
  #write.csv(mig.men, "mig_men.csv")
  
  #non.mig.men <- data.frame(non.mig.men = dat$epi$non.mig.men)
  #write.csv(non.mig.men, "non_mig_men.csv")
  
  #women <- data.frame(women = dat$epi$women)
  #write.csv(women, "women.csv")
  
#   inf <- data.frame(inf = dat$epi$i.num)
#   write.csv(inf, filename=paste(nsim, "inf_4.csv", sep=""))
  
  #inf.women <- data.frame(inf.women = dat$epi$inf.women)
  #write.csv(inf.women, "inf_women.csv")
  
  
  #inf.mig.men <- data.frame(inf.mig.men = dat$epi$inf.mig.men)
  #write.csv(inf.mig.men, "inf_mig_men.csv")
  
  #inf.non.mig.men <- data.frame(inf.non.mig.men = dat$epi$inf.non.mig.men)
  #write.csv(inf.non.mig.men, "inf_non_mig_men.csv")
  
  
  
  return(dat)
}