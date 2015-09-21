verbose_example <- function (x, type, s = 1, at = 2) 
{
  if (type == "startup") {
    if (x$verbose == TRUE) {
      cat("\n* Starting Network Simulation")
    }
  }
  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (x$control$verbose.int == 0 && at == x$control$nsteps) {
        cat("\nSim = ", s, "/", x$control$nsims, sep = "")
      }
      if (x$control$verbose.int > 0 && (at%%x$control$verbose.int == 
                                          0)) {
        cat("\f")
        cat("\nEpidemic Simulation")
        cat("\n----------------------------")
        cat("\nSimulation: ", s, "/", x$control$nsims, 
            sep = "")
        cat("\nTimestep: ", at, "/", x$control$nsteps, 
            sep = "")
        if (x$param$modes == 1) {
          cat("\nIncidence:", x$epi$si.flow[at])
          cat("\nEdge Count: ", x$epi$edge_count[at])
          cat("\nMean Degree: ", x$epi$mean_deg[at])
          cat("\nSize: ", x$epi$size[at])
        }
        if (x$param$modes == 2) {
          cat("\nIncidence:", x$epi$si.flow[at] + x$epi$si.flow.m2[at])
        }
        if (x$control$type == "SIR") {
          if (x$param$modes == 1) {
            cat("\nRecoveries:", x$epi$ir.flow[at])
          }
          if (x$param$modes == 2) {
            cat("\nRecoveries:", x$epi$ir.flow[at] + 
                  x$epi$ir.flow.m2[at])
          }
        }
        if (x$control$type == "SIS") {
          if (x$param$modes == 1) {
            cat("\nRecoveries:", x$epi$is.flow[at])
          }
          if (x$param$modes == 2) {
            cat("\nRecoveries:", x$epi$is.flow[at] + 
                  x$epi$is.flow.m2[at])
          }
        }
        if (x$param$modes == 1) {
          cat("\nPrevalence:", x$epi$i.num[at]/x$epi$size[at])
        }
        if (x$param$modes == 2) {
          cat("\nPrevalence:", x$epi$i.num[at] + x$epi$i.num.m2[at])
        }
        if (x$control$type %in% c("SI", "SIS")) {
          if (x$param$modes == 1) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$i.num[at])
          }
          if (x$param$modes == 2) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$s.num.m2[at] + 
                  x$epi$i.num[at] + x$epi$i.num.m2[at])
          }
        }
        if (x$control$type == "SIR") {
          if (x$param$modes == 1) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$i.num[at] + 
                  x$epi$r.num[at])
          }
          if (x$param$modes == 2) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$i.num[at] + 
                  x$epi$r.num[at] + x$epi$s.num.m2[at] + 
                  x$epi$i.num.m2[at] + x$epi$r.num.m2[at])
          }
        }
        if (x$param$vital == TRUE) {
          if (x$param$modes == 1) {
            cat("\nBirths:", x$epi$b.flow[at])
            cat("\nDeaths, susceptibles:", x$epi$ds.flow[at])
            cat("\nDeaths, infecteds:", x$epi$di.flow[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$epi$dr.flow[at])
            }
          }
          if (x$param$modes == 2) {
            cat("\nBirths:", x$epi$b.flow[at] + x$epi$b.flow.m2[at])
            cat("\nDeaths, susceptible:", x$epi$ds.flow[at] + 
                  x$epi$ds.flow.m2[at])
            cat("\nDeaths, infected:", x$epi$di.flow[at] + 
                  x$epi$di.flow.m2[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$epi$dr.flow[at] + 
                    x$epi$dr.flow.m2[at])
            }
          }
        }
        cat("\n----------------------------")
      }
    }
  }
}