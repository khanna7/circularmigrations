nw_migrations <- network.initialize(5000, directed = FALSE)

sex.vec <- c(rep(1, 2500), rep(0, 2500))
loc.vec <- c(rep(0, 1250), rep(1, 1250), rep(0, 1250), rep(1, 1250))
mig.vec <- c(rep(1, 625), rep(0, 625), rep(1, 625), rep(0, 625), rep(0, 2500))
type.vec <- c(rep("B-MM", 625), rep("C-MU", 625), rep("B-MM", 625), rep("A-MR", 625), rep("E-FU", 1250), rep("D-FR", 1250))

nw_migrations <- set.vertex.attribute(nw_migrations, "sex", sex.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "loc", loc.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "mig_stat", mig.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "type", type.vec)
formation <- ~edges #+ nodemix("type", base = c(-1, -2, -3, -4, -5, -6, -7, -9, -10, -11, -12, -13, -14, -15))
target.stats <- c(3000) #0, rep(0, 2), rep(0, 3), 500, rep(0, 2), 0, 1000, 500, 0, 0)

est_migrations <- netest(nw_migrations, formation = formation, target.stats = target.stats, coef.diss = dissolution_coefs(~offset(edges), 100))

#dx <- netdx(est_migrations, nsteps = 1500, nsims = 1)

param <- param.net(death.rate.gen = 1/(45*52), death.rate.aids = 1/40, migration.rate = 1/30, birth.rate = 8*1250/(45*52), late.cutoff = 512, inf.prob = c(rep(1-(1-.0007*26)^2.4, 12), rep(1-(1-.0007)^2.4, 500), rep(1-(1-.0007*7)^2.4, 40)))

source('~/circularmigrations/epimodel/network/migration.R')
source('~/circularmigrations/epimodel/network/infection.R')
source('~/circularmigrations/epimodel/network/births.R')
source('~/circularmigrations/epimodel/network/deaths.R')
source('~/circularmigrations/epimodel/network/get_prev.R')
source('~/circularmigrations/epimodel/network/initialize_net.R')

control <- control.net(type="SI", nsims = 1, nsteps = 200, initialize.FUN = initialize.net.mig, deaths.FUN = deaths, births.FUN = births, edges_correct.FUN = edges_correct, infection.FUN = infection, migration.FUN = migration, get_prev.FUN = get_prev, verbose.FUN = verbose, depend = TRUE)
#status.vector <- rbinom(5000, 1, 0.1)
#status.vector <- replace(status.vector, status.vector == 0, "s")
#status.vector <- replace(status.vector, status.vector == 1, "i")
init <- init.net(i.num = 0)

mod <- netsim(est_migrations, param, init, control)
