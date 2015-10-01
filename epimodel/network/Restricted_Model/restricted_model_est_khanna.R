#3 contacts per week
#with death correction
#migration rate = 1/30

rm(list=ls())
library(EpiModel)

N <- 1000*5

N.M.R<- 1/8*N
N.M.M <- 1/4*N # Ditto migrant
N.M.U<- 1/8*N # Ditto urban
N.M<-N.M.R+N.M.M+N.M.U # Total number of males

N.F.R<- 1/4*N # Number of rural females
N.F.M<-0 # Ditto migrant
N.F.U<-1/4*N # Ditto urban
N.F<-N.F.R+N.F.M+N.F.U # Total number of females

nw_migrations <- network.initialize(5000, directed = FALSE)

sex.vec <- c(rep(1, 2500), rep(0, 2500))
loc.vec <- c(rep(0, 1250), rep(1, 1250), rep(0, 1250), rep(1, 1250))
mig.vec <- c(rep(1, 625), rep(0, 625), rep(1, 625), rep(0, 625), rep(0, 2500))
type.vec <- c(rep("A-MR", 625), rep("B-MM", 625), rep("B-MM", 625),
              rep("C-MU", 625), rep("D-FR", 1250), rep("E-FU", 1250))

# restriction.vec.urban <- c(rep(0, 625), rep(1, 625), rep(0, 625), rep(2, 625), rep(0, 1250), rep(3, 1250))
# restriction.vec.rural <- c(rep(0, 625), rep(1, 625), rep(0, 625), rep(2, 625), rep(3, 1250), rep(0, 1250))

nw_migrations <- set.vertex.attribute(nw_migrations, "sex", sex.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "loc", loc.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "mig_stat", mig.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "type", type.vec)
# nw_migrations <- set.vertex.attribute(nw_migrations, "rest.urb", restriction.vec.urban)
# nw_migrations <- set.vertex.attribute(nw_migrations, "rest.rur", restriction.vec.rural)

## migrant-male mixing
nw_migrations %v% "migmalemix.rural" <- rep(c("A-X", "B-MMR.MIX",
                                              "C-Y", "B-MMR.MIX", "D-Z"),
                                            #s for correct ordering
                                            # left out female migrant cell, 5 FM
                                            # because
                                            # it was confusing
                                            # change from masters
                                            c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)

## migrant-male mixing
nw_migrations %v% "migmalemix.urban" <- rep(c("1X", "2MMU.MIX",
                                              "3Y", "4Z", "2MMU.MIX"),
                                            #s for correct ordering
                                            # left out female migrant cell, 5 FM
                                            # because
                                            # it was confusing
                                            # change from masters
                                            c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)


formation <- ~edges + nodemix("type", base = c(-1, -2, -3, -4, -5, -6, -7, -9, -10, -11, -12, -13, -14, -15)) + 
  kstar(2, "migmalemix.rural") + kstar(2, "migmalemix.urban")
target.stats <- c(2500, 
                  0, rep(0, 2), rep(0, 3), 417, rep(0, 2), 0, 833, 417, 0, 0, 
                  rep(0, 2))

est_migrations <- netest(nw_migrations, 
                         formation = formation, 
                         target.stats = target.stats, 
                         coef.diss = dissolution_coefs(~offset(edges), 100, 1/(45*52)))


#dx <- netdx(est_migrations, nsteps = 1500, nsims = 1)

param <- param.net(death.rate.gen = 1/(45*52), death.rate.aids = 1/40, 
                   migration.rate = 1/30, birth.rate = 8*625/(45*52), 
                   late.cutoff = 512, inf.prob = c(rep(1-(1-.0007*26)^3, 12), 
                                                   rep(1-(1-.0007)^3, 500), rep(1-(1-.0007*7)^3, 40)))
save.image(file='estimation_restricted.RData')
# source('Projects/circularmigrations/epimodel/network/Restricted_Model/migration_restricted.R')
# source('Projects/circularmigrations/epimodel/network/Restricted_Model/infection_restricted.R')
# source('Projects/circularmigrations/epimodel/network/Restricted_Model/births_restricted.R')
# source('Projects/circularmigrations/epimodel/network/Restricted_Model/deaths_restricted.R')
# source('Projects/circularmigrations/epimodel/network/Restricted_Model/get_prev__restricted_example_4.R')
# source('Projects/circularmigrations/epimodel/network/Restricted_Model/initialize_net_restricted.R')
# source('Projects/circularmigrations/epimodel/network/Restricted_Model/verbose_restricted_example.R')
# 
# 
# #control <- control.net(type="SI", nsims = 1, nsteps = 200, initialize.FUN = initialize.net.mig, deaths.FUN = deaths, births.FUN = births, edges_correct.FUN = edges_correct, infection.FUN = infection, migration.FUN = migration, get_prev.FUN = get_prev, verbose.FUN = verbose, depend = TRUE)
#control <- control.net(type = "SI", nsims = 1, nsteps = 5000, initialize.FUN = initialize.net.mig,
#deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
#migration.FUN = migration,
#get_prev.FUN = get_prev_example_5,
#verbose.FUN = verbose_example, depend = TRUE)
# #status.vector <- rbinom(5000, 1, 0.1)
# #status.vector <- replace(status.vector, status.vector == 0, "s")
# #status.vector <- replace(status.vector, status.vector == 1, "i")
#init <- init.net(i.num = 0)
# 
#mod <- netsim(est_migrations, param, init, control)