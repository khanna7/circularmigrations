migration <- function(dat, at) {
  
  migration.rate <- dat$param$migration.rate
  mig_stat_vec <- get.vertex.attribute(dat$nw, "mig_stat")
  loc_vec <- get.vertex.attribute(dat$nw, "loc")
  migrants <- which(mig_stat_vec == 1) 
  n = length(migrants)
  temp <- rbinom(n, 1, migration.rate)
  
  for (i in 1:n) {
    if (temp[i] == 1) {
      if (loc_vec[migrants[i]] == 0) {
        set.vertex.attribute(dat$nw, "loc", 1, migrants[i])
      }
      else {
        set.vertex.attribute(dat$nw, "loc", 0, migrants[i])
      }
        
    }
  }
  
  return(dat)
}