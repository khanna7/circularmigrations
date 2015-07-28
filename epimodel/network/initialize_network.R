nw_migrations <- network.initialize(5000, directed = FALSE)

sex.vec <- c(rep(1, 2500), rep(0, 2500))
loc.vec <- c(rep(0, 1250), rep(1, 1250), rep(0, 1250), rep(1, 1250))
mig.vec <- c(rep(1, 625), rep(0, 625), rep(1, 625), rep(0, 625), rep(0, 2500))

nw_migrations <- set.vertex.attribute(nw_migrations, "sex", sex.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "loc", loc.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "mig_stat", mig.vec)

summary(nw_migrations~nodemix(c("sex", "mig_stat", "loc"))) #to see order of terms for base argument

formation <- ~edges + offset(nodemix(c("sex", "mig_stat", "loc"), base = c(4, 8, 11, 12, 16, 17))) + nodemix(c("sex", "mig_stat", "loc"), base = c(-4, -8, -11, -12, -16, -17))
#formation offsets all but the 6 possible types of ties to 0

target.stats <- c(5000*1.2/2, rep(5000*1.2/12, 6))

est_migrations <- netest(nw_migrations, formation = formation, 
                         target.stats = target.stats,
                         coef.diss = dissolution_coefs(~offset(edges), 100),
                         coef.form = rep(-Inf, 15))

dx <- netdx(est_migrations, nsims = 1, nsteps = 5000) #one long run just to get sense of how it's working

dx

  
  