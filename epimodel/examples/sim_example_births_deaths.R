source('~/circularmigrations/epimodel/examples/births_example.R')
source('~/circularmigrations/epimodel/examples/deaths_example.R')

nw <- network.initialize(500, directed = FALSE)

formation <- ~edges

target.stats <- c(300)

est_nw <- netest(nw, formation = formation, target.stats = target.stats,
                 coef.diss = dissolution_coefs(~offset(edges), 25, 0.02))

param <- param.net(inf.prob = 0, death.prob = 0.02, birth.rate = 10)

control <- control.net(type = "SI", nsims = 1, nsteps = 1000, deaths.FUN = deaths_example, 
                       births.FUN = births_example, get_prev.FUN = get_prev_example,
                       verbose.FUN = verbose_example, depend = TRUE)

init <- init.net(i.num = 0)

mod_2_percent <- netsim(est_nw, param, init, control)