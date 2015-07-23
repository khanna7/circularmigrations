#call last in the simulation
summary_statistics <- function(dat, at) {
  
  nw.sex <- get.vertex.attribute(dat$nw, "sex")
  nw.loc <- get.vertex.attribute(dat$nw, "loc")
  nw.mig_stat <- get.vertex.attribute(dat$nw, "mig_stat")
  
  nw.men <- which(nw.sex == 1)
  nw.women <- which(nw.sex == 0)
  nw.migrant <- which(nw.mig_stat == 1)
  nw.nonmigrant <- which(nw.mig_stat == 0)
  nw.urban <- which(nw.loc == 0)
  nw.rural <- which(nw.loc == 1)
  nw.infected <- which(dat$attr$status == "i")
  #total and infected MMU
  
  if (at == 2) {
    dat$epi$n.mmu <- c(0, length(intersect(nw.men, intersect(nw.migrant, mw.urban))))
    dat$epi$inf.mmu <- c(0, length(intersect(nw.men, intersect(nw.migrant, intersect(nw.urban, nw.infected)))))
  } else {
    dat$epi$n.mmu[at] <- length(intersect(nw.men, intersect(nw.migrant, mw.urban)))
    dat$epi$inf.mmu[at] <- length(intersect(nw.men, intersect(nw.migrant, intersect(nw.urban, nw.infected))))
  }
  
  #total and infected MMR
  if (at == 2) {
    dat$epi$n.mmr <- c(0, length(intersect(nw.men, intersect(nw.migrant, mw.rural))))
    dat$epi$inf.mmr <- c(0, length(intersect(nw.men, intersect(nw.migrant, intersect(nw.rural, nw.infected)))))
  } else {
    dat$epi$n.mmr[at] <- length(intersect(nw.men, intersect(nw.migrant, mw.rural)))
    dat$epi$inf.mmr[at] <- length(intersect(nw.men, intersect(nw.migrant, intersect(nw.rural, nw.infected))))
  }
  
  #total and infected NMU
  if (at == 2) {
    dat$epi$n.nmu <- c(0, length(intersect(nw.men, intersect(nw.nonmigrant, mw.urban))))
    dat$epi$inf.nmu <- c(0, length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected)))))
  } else {
    dat$epi$n.nmu[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, mw.urban)))
    dat$epi$inf.nmu[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected))))
  }
  
  #total and infected NMR
  if (at == 2) {
    dat$epi$n.nmr <- c(0, length(intersect(nw.men, intersect(nw.nonmigrant, mw.rural))))
    dat$epi$inf.nmr <- c(0, length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected)))))
  } else {
    dat$epi$n.nmr[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, mw.rural)))
    dat$epi$inf.nmr[at] <- length(intersect(nw.men, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected))))
  }
  
  #total and infected NFU
  if (at == 2) {
    dat$epi$n.nfu <- c(0, length(intersect(nw.women, intersect(nw.nonmigrant, mw.urban))))
    dat$epi$inf.nfu <- c(0, length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected)))))
  } else {
    dat$epi$n.nfu[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, mw.urban)))
    dat$epi$inf.nfu[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.urban, nw.infected))))
  }
  
  #total and infected NFR
  if (at == 2) {
    dat$epi$n.nfr <- c(0, length(intersect(nw.women, intersect(nw.nonmigrant, mw.rural))))
    dat$epi$inf.nfr <- c(0, length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected)))))
  } else {
    dat$epi$n.nfr[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, mw.rural)))
    dat$epi$inf.nfr[at] <- length(intersect(nw.women, intersect(nw.nonmigrant, intersect(nw.rural, nw.infected))))
  }
}