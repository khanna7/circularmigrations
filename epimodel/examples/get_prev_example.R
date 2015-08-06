get_prev_example <- function (dat, at) 
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
  } else {
    dat$epi$edge_count[at] <- network.edgecount(network.extract(dat$nw, at = at))
    dat$epi$size[at] <- network.size(network.extract(dat$nw, at = at))
    dat$epi$mean_deg[at] <- 2*dat$epi$edge_count[at]/length(which(dat$attr$active == 1))
  }
  return(dat)
}