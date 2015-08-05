#code to generate a 500 node, 300 edge network, no vital 
#dynamics, and a mean partnership duration of 25 time steps

nw <- network.initialize(500, directed = FALSE)

formation <- ~edges

target.stats <- c(300)

est_nw <- netest(nw, formation = formation, target.stats = target.stats,
             coef.diss = dissolution_coefs(~offset(edges), 25))
dx <- netdx(est_nw, nsteps = 5000, nsims = 1)
dx
