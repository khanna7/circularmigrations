#modified from EpiModel source code for get_prev.net

get_prev <- function (dat, at) 
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
  
#Begin Nathan's code. From the old summary_statistics function.
#Moved within get_prev so that this code would be called
#at the end of the simulation
  nw.active <- which(dat$attr$active == 1)
  nw.sex <- intersect(nw.active, get.vertex.attribute(dat$nw, "sex"))
  nw.loc <- intersect(nw.active, get.vertex.attribute(dat$nw, "loc"))
  nw.mig_stat <- intersect(nw.active, get.vertex.attribute(dat$nw, "mig_stat"))
  
  nw.men <- which(nw.sex == 1)
  nw.women <- which(nw.sex == 0)
  nw.migrant <- which(nw.mig_stat == 1)
  nw.nonmigrant <- which(nw.mig_stat == 0)
  nw.urban <- which(nw.loc == 0)
  nw.rural <- which(nw.loc == 1)
  nw.infected <- which(dat$attr$status == "i")
  #total and infected MMU
  
  if (at == 1) {
    dat$epi$n.mmu <- length(intersect(nw.men, intersect(nw.migrant, nw.urban)))
    dat$epi$inf.mmu <- length(intersect(nw.men, intersect(nw.migrant, intersect(nw.urban, nw.infected))))
  } else {
    dat$epi$n.mmu[at] <- length(intersect(nw.men, intersect(nw.migrant, nw.urban)))
    dat$epi$inf.mmu[at] <- length(intersect(nw.men, intersect(nw.migrant, intersect(nw.urban, nw.infected))))
  }
  
  #total and infected MMR
  if (at == 1) {
    dat$epi$n.mmr <- length(intersect(nw.men, intersect(nw.migrant, nw.rural)))
    dat$epi$inf.mmr <- length(intersect(nw.men, intersect(nw.migrant, intersect(nw.rural, nw.infected))))
  } else {
    dat$epi$n.mmr[at] <- length(intersect(nw.men, intersect(nw.migrant, nw.rural)))
    dat$epi$inf.mmr[at] <- length(intersect(nw.men, intersect(nw.migrant, intersect(nw.rural, nw.infected))))
  }
  
  #total and infected NMU
  if (at == 1) {
    dat$epi$n.nmu <- length(intersect(nw.men, intersect(nw.nonmigrant, nw.urban)))
    dat$epi$inf.nmu <- length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected))))
  } else {
    dat$epi$n.nmu[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, nw.urban)))
    dat$epi$inf.nmu[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected))))
  }
  
  #total and infected NMR
  if (at == 1) {
    dat$epi$n.nmr <- length(intersect(nw.men, intersect(nw.nonmigrant, nw.rural)))
    dat$epi$inf.nmr <- length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected))))
  } else {
    dat$epi$n.nmr[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, nw.rural)))
    dat$epi$inf.nmr[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected))))
  }
  
  #total and infected NFU
  if (at == 1) {
    dat$epi$n.nfu <- length(intersect(nw.women, intersect(nw.nonmigrant, nw.urban)))
    dat$epi$inf.nfu <- length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected))))
  } else {
    dat$epi$n.nfu[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, nw.urban)))
    dat$epi$inf.nfu[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected))))
  }
  
  #total and infected NFR
  if (at == 1) {
    dat$epi$n.nfr <- length(intersect(nw.women, intersect(nw.nonmigrant, nw.rural)))
    dat$epi$inf.nfr <- length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected))))
  } else {
    dat$epi$n.nfr[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, nw.rural)))
    dat$epi$inf.nfr[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected))))
  }

  #Prevalence
  if (at == 1) {
    dat$epi$prevalence <- length(nw.infected)/length(nw.active)
  } else {
    dat$epi$prevalence[at] <- length(nw.infected)/length(nw.active)
  }
  
  #dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf, v = nw.active, deactivate.edges = TRUE)
  if (at == 1) {
    dat$epi$mean_deg <- network.edgecount(dat$nw)/network.size(dat$nw)
    dat$epi$size <- network.size(dat$nw)
    dat$epi$edges <- network.edgecount(dat$nw)
  } else {
    dat$epi$mean_deg[at] <- network.edgecount(dat$nw)/network.size(dat$nw)
    dat$epi$size[at] <- network.size(dat$nw)
    dat$epi$edges[at] <- network.edgecount(dat$nw)
  }
  
  return(dat)
}