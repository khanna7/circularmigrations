rm(list=ls())

library(EpiModel) 
load("estimation_restricted.RData")

source("migration_restricted.R")
source("infection_restricted.R")
source("births_restricted.R")
source("deaths_restricted.R")
source("get_prev_restricted_example_4.R")
source("initialize_net_restricted.R")
source("verbose_restricted_example.R")
nsim <- sim10

#control <- control.net(type="SI", nsims = 1, nsteps = 200, initialize.FUN = initialize.net.mig, deaths.FUN = deaths, births.FUN = births, edges_correct.FUN = edges_correct, infection.FUN = infection, migration.FUN = migration, get_prev.FUN = get_prev, verbose.FUN = verbose, depend = TRUE)
control <- control.net(type = "SI", nsims = 10, nsteps = 5000, initialize.FUN = initialize.net.mig,
                       deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
                       migration.FUN = migration,
                       get_prev.FUN = get_prev_example_4,
                       verbose.FUN = verbose_example, depend = TRUE)
#status.vector <- rbinom(5000, 1, 0.1)
#status.vector <- replace(status.vector, status.vector == 0, "s")
#status.vector <- replace(status.vector, status.vector == 1, "i")
init <- init.net(i.num = 0)

mod <- netsim(est_migrations, param, init, control)