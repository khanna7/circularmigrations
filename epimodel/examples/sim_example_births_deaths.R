nw <- network.initialize(500, directed = FALSE)

formation <- ~edges

target.stats <- c(300)

est_nw <- netest(nw, formation = formation, target.stats = target.stats,
                 coef.diss = dissolution_coefs(~offset(edges), 25))

param <- param.net(inf.prob = 0, death.prob = 0.001, birth.rate = 0.5)

control <- control.net(type = "SI", nsims = 3, nsteps = 1000, deaths.FUN = deaths_example, births.FUN = births_example)

init <- init.net(i.num = 0)

mod <- netsim(est_nw, param, init, control)