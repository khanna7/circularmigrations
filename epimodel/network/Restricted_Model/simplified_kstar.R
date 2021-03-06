nw_migrations <- network.initialize(5000, directed = FALSE)

sex.vec <- c(rep(1, 2500), rep(0, 2500))
loc.vec <- c(rep(0, 1250), rep(1, 1250), rep(0, 1250), rep(1, 1250))
mig.vec <- c(rep(1, 625), rep(0, 625), rep(1, 625), rep(0, 625), rep(0, 2500))
type.vec <- c(rep("B-MM", 625), rep("C-MU", 625), rep("B-MM", 625), rep("A-MR", 625), rep("E-FU", 1250), rep("D-FR", 1250))

restriction.vec.urban <- c(rep(0, 625), rep(1, 625), rep(0, 625), rep(2, 625), rep(0, 1250), rep(3, 1250))
restriction.vec.rural <- c(rep(0, 625), rep(1, 625), rep(0, 625), rep(2, 625), rep(3, 1250), rep(0, 1250))


nw_migrations <- set.vertex.attribute(nw_migrations, "sex", sex.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "loc", loc.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "mig_stat", mig.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "type", type.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "rest.urb", restriction.vec.urban)
nw_migrations <- set.vertex.attribute(nw_migrations, "rest.rur", restriction.vec.rural)
formation <- ~edges + kstar(2, "loc") #nodemix("type", base = c(-1, -2, -3, -4, -5, -6, -7, -9, -10, -11, -12, -13, -14, -15)) + kstar(2, "rest.urb") + kstar(2, "rest.rur")
target.stats <- c(2500, 0)#, rep(0, 2), rep(0, 3), 417, rep(0, 2), 0, 833, 417, 0, 0, rep(0, 2))

est_migrations <- netest(nw_migrations, formation = formation, target.stats = target.stats, coef.diss = dissolution_coefs(~offset(edges), 100, 1/(45*52)))

source("migration_restricted.R")
source("infection_restricted.R")
source("births_restricted.R")
source("deaths_restricted.R")
source("get_prev_restricted_example_4.R")
source("initialize_net_restricted.R")
source("verbose_restricted_example.R")

nsim <- 1

#control <- control.net(type="SI", nsims = 1, nsteps = 200, initialize.FUN = initialize.net.mig, deaths.FUN = deaths, births.FUN = births, edges_correct.FUN = edges_correct, infection.FUN = infection, migration.FUN = migration, get_prev.FUN = get_prev, verbose.FUN = verbose, depend = TRUE)
control <- control.net(type = "SI", nsims = 1, nsteps = 5000, initialize.FUN = initialize.net.mig,
                       deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
                       migration.FUN = migration,
                       get_prev.FUN = get_prev_example_4,
                       verbose.FUN = verbose_example, depend = TRUE)
#status.vector <- rbinom(5000, 1, 0.1)
#status.vector <- replace(status.vector, status.vector == 0, "s")
#status.vector <- replace(status.vector, status.vector == 1, "i")
init <- init.net(i.num = 0)

mod <- netsim(est_migrations, param, init, control)
