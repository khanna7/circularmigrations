source('~/Projects/circularmigrations/epimodel/network/Restricted_Model/estimation_restricted.R', echo=TRUE)
formation
list.vertex.attributes(nw_migrations)
table(nw_migration %v% "rest.rur")
nw_migrations
head(nw_migrations%v%"rest.rur")
table(nw_migrations%v%"rest.rur")
table(nw_migrations%v%"rest.urb")
N <- 1000*5
N.M.R<- 1/8*N
N.M.M <- 1/4*N # Ditto migrant
N.M.U<- 1/8*N # Ditto urban
N.M<-N.M.R+N.M.M+N.M.U # Total number of males
N.F.R<- 1/4*N # Number of rural females
N.F.M<-0 # Ditto migrant
N.F.U<-1/4*N # Ditto urban
N.F<-N.F.R+N.F.M+N.F.U # Total number of females
### Construct the empty network.
n0<-network.initialize(N,bipartite=FALSE,directed=FALSE)
### Add Atributes
## Type
n0 %v% "type" <- rep(c("A-MR", "B-MM", "C-MU", "D-FR", "F-FU"),
#s for correct ordering
# left out female migrant cell, 5 FM
# because
# it was confusing
# change from masters
c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)
## In order these are
### Number-Male-Rural, Number-Male-Migrant, Number-Male-Urban,
### Number-Female-Rural, Number-Female-Migrant, Number-Female-Urban
## location
n0 %v% "location" <- c(rep(0, N.M.R), #0=rural
rep(0, N.M.M/2), #1=urban
rep(1, N.M.M/2),
rep(1, N.M.U),
rep(0, N.F.R),
rep(1, N.F.U)
)
## infection status
n0 %v% "inf.status" <- rep(0, N)
#set.vertex.attribute(n0, "inf.status", 1, c(2501,3751)) # 1 rural and urban female inf.
set.vertex.attribute(n0, "inf.status", 1, c(2501:2550,
3751:3800)) # for testing
## time of infection
n0 %v% "inf.time" <- rep(NA, N)
#set.vertex.attribute(n0, "inf.time", 0, c(2501,3751)) # 1 rural and urban female inf.
set.vertex.attribute(n0, "inf.time", 0, c(2501:2550,
3751:3800)) # for testing
## time since infection
n0 %v% "time.infected" <- rep(NA, N)
#set.vertex.attribute(n0, "time.infected", 0, c(2501,3751)) # 1 rural and urban female inf.
##lifespan.inf.ind <- 616
lifespan.inf.ind <- 12+500+40 # changed after discussion w. steve
init.time.infected <- runif(length(c(2501:2550, 3751:3800)),
min=0, max=lifespan.inf.ind)
set.vertex.attribute(n0, "time.infected", init.time.infected,
c(2501:2550, 3751:3800)) # for testing
## infector.id
n0 %v% "infector.id" <- rep(NA, N)
#set.vertex.attribute(n0, "infector.id",
#                     0, c(2501,3751)) # 1 rural and urban female inf.
set.vertex.attribute(n0, "infector.id", 0, c(2501:2550,
3751:3800)) # for testing
## migrant-male mixing
n0 %v% "migmalemix.rural" <- rep(c("A-X", "B-MMR.MIX",
"C-Y", "B-MMR.MIX", "D-Z"),
#s for correct ordering
# left out female migrant cell, 5 FM
# because
# it was confusing
# change from masters
c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)
## migrant-male mixing
n0 %v% "migmalemix.urban" <- rep(c("1X", "2MMU.MIX",
"3Y", "4Z", "2MMU.MIX"),
#s for correct ordering
# left out female migrant cell, 5 FM
# because
# it was confusing
# change from masters
c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)
## Behavior
duration <- 100
### Mean-Statistics for Network
mean.deg <- 1
n.edges <- N*mean.deg/2 # change mean degree to 1
#MF.type.mix<-c(450,300,0,450,0,300)
MF.type.mix <- c(0,
rep(0,2),
rep(0,3),
(1/6*n.edges), (1/3*n.edges), rep(0,2),
0, (1/3*n.edges), (1/6*n.edges), 0)
## try with leaving out a non-structural zero
## the following works with term 7 left out from
## the mixing matrix
## MF.type.mix <- c(0,
##                  rep(0,2),
##                  rep(0,3),
##                  (1/3*n.edges), rep(0,2),
##                  0, (1/3*n.edges), (1/6*n.edges), 0, 0)
##  the trick to have the correct 'restricted partnership structure'
## only seems to work when a structural zero is left out. otherwise,
## it fails.
## no male has degree > 2, as it should be
twopath.stat <- (2/3)*n.edges*(1/2)
kstar.loc.stat <- (2/3)*n.edges
kstar.mmr.mix <- 1/3*n.edges
kstar.mmu.mix <- 1/3*n.edges
kstar.mmr.mix.2star <- 0
kstar.mmu.mix.2star <- 0
kstar.mmr.mix.3star <- 0
ls()
n0
list.vertex.attributes(n0)
table(n0%v%"migmalemix.urban")
head(n0%v%"migmalemix.urban")
?ergm.terms
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
type.vec <- c(rep("B-MM", 625), rep("C-MU", 625), rep("B-MM", 625),
rep("A-MR", 625), rep("E-FU", 1250), rep("D-FR", 1250))
# restriction.vec.urban <- c(rep(0, 625), rep(1, 625), rep(0, 625), rep(2, 625), rep(0, 1250), rep(3, 1250))
# restriction.vec.rural <- c(rep(0, 625), rep(1, 625), rep(0, 625), rep(2, 625), rep(3, 1250), rep(0, 1250))
nw_migrations <- set.vertex.attribute(nw_migrations, "sex", sex.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "loc", loc.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "mig_stat", mig.vec)
nw_migrations <- set.vertex.attribute(nw_migrations, "type", type.vec)
# nw_migrations <- set.vertex.attribute(nw_migrations, "rest.urb", restriction.vec.urban)
# nw_migrations <- set.vertex.attribute(nw_migrations, "rest.rur", restriction.vec.rural)
n0 %v% "migmalemix.rural" <- rep(c("A-X", "B-MMR.MIX",
"C-Y", "B-MMR.MIX", "D-Z"),
#s for correct ordering
# left out female migrant cell, 5 FM
# because
# it was confusing
# change from masters
c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)
nw_migrations %v% "migmalemix.rural" <- rep(c("A-X", "B-MMR.MIX",
"C-Y", "B-MMR.MIX", "D-Z"),
#s for correct ordering
# left out female migrant cell, 5 FM
# because
# it was confusing
# change from masters
c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)
nw_migrations %v% "migmalemix.urban" <- rep(c("1X", "2MMU.MIX",
"3Y", "4Z", "2MMU.MIX"),
#s for correct ordering
# left out female migrant cell, 5 FM
# because
# it was confusing
# change from masters
c(N.M.R, N.M.M, N.M.U, N.F.R, N.F.U)
)
#3 contacts per week
#with death correction
#migration rate = 1/30
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
type.vec <- c(rep("B-MM", 625), rep("C-MU", 625), rep("B-MM", 625),
rep("A-MR", 625), rep("E-FU", 1250), rep("D-FR", 1250))
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
param <- param.net(death.rate.gen = 1/(45*52), death.rate.aids = 1/40,
migration.rate = 1/30, birth.rate = 8*625/(45*52),
late.cutoff = 512, inf.prob = c(rep(1-(1-.0007*26)^3, 12), rep(1-(1-.0007)^3, 500), rep(1-(1-.0007*7)^3, 40)))
save.image(file='estimation_restricted.RData')
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
nsim <- 1
#control <- control.net(type="SI", nsims = 1, nsteps = 200, initialize.FUN = initialize.net.mig, deaths.FUN = deaths, births.FUN = births, edges_correct.FUN = edges_correct, infection.FUN = infection, migration.FUN = migration, get_prev.FUN = get_prev, verbose.FUN = verbose, depend = TRUE)
control <- control.net(type = "SI", nsims = 1, nsteps = 5000, initialize.FUN = initialize.net.mig,
deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
migrahetion.FUN = migration,
get_prev.FUN = get_prev_example_4,
verbose.FUN = verbose_example, depend = TRUE)
#status.vector <- rbinom(5000, 1, 0.1)
#status.vector <- replace(status.vector, status.vector == 0, "s")
#status.vector <- replace(status.vector, status.vector == 1, "i")
net <- network.initialize(50, directed=FALSE)
net
net%v%"att" <- c(rep("A", 25), rep("B", 25))
net%v%"att"
formation <- ~kstar(2, "att")
formation <- net~kstar(2, "att")
formation
?ergm
ergm(net ~ kstar(2, "att"))
ergm(net ~ kstar(2, "att"), target.stats=c(20, 20))
ergm(net ~ kstar(2, "att"), target.stats=c(20))
ergm(net ~ edges+kstar(2, "att"), target.stats=c(50, 20))
mod1 <- ergm(net ~ edges+kstar(2, "att"), target.stats=c(50, 20))
sim1 <- simulate(mod1)
sim1
mixingmatrix(sim1, "att")
sim2
sim2 <- simulate(mod1)
mixingmatrix(sim2, "att")
sim3 <- simulate(mod1)
mixingmatrix(sim3, "att")
mod1 <- ergm(net ~ edges+kstar(2, "att"), target.stats=c(50, 20))
q("no")
library(ergm)
net <- ergm(50, directed=FALSE
)
n0 <- initialize(50, directed=FALSE)
n0%v%"att" <- c(rep("A", 25), rep("B", 25))
n0 <- network.initialize(50, directed=FALSE)
n0%v%"att" <- c(rep("A", 25), rep("B", 25))
mod1 <- ergm(net ~ edges+kstar(2, "att"), target.stats=c(50, 0))
mod1 <- ergm(n0 ~ edges+kstar(2, "att"), target.stats=c(50, 0))
sim1 <- simulate(mod1); mixingmatrix(sim1, “att”)
sim1 <- simulate(mod1)
mixingmatrix(sim1, “att”)
mixingmatrix(sim1, "att")
sim2 <- simulate(mod1)
mixingmatrix(sim2, "att")
as.edgelist(sim1)
history()
library(statnet)
install.packages("statnet")
?simulate
?networkDynamic:simulate
?simulate.stergm
libtary(tergm)
library(tergm)
?simulate.stergm
simulate.stergm
?simulate.network
simulate.network
?stergm.getMCMCsample
stergm.getMCMCsample()
stergm.getMCMCsample
q("no")
source('~/Projects/circularmigrations/epimodel/network/Restricted_Model/estimation_restricted_khanna.R', echo=TRUE)
getwd()
setwd("epimodel/network/Restricted_Model/)
")"
""
setwd("epimodel/network/Restricted_Model/")
source('~/Projects/circularmigrations/epimodel/network/Restricted_Model/estimation_restricted_khanna.R', echo=TRUE)
rm(list=ls())
library(EpiModel)
source("migration_restricted.R")
source("infection_restricted.R")
source("births_restricted.R")
source("deaths_restricted.R")
source("get_prev_restricted_example_4.R")
source("initialize_net_restricted.R")
source("verbose_restricted_example.R")
load("estimation_restricted.RData")
nsim <- 1
#control <- control.net(type="SI", nsims = 1, nsteps = 200, initialize.FUN = initialize.net.mig, deaths.FUN = deaths, births.FUN = births, edges_correct.FUN = edges_correct, infection.FUN = infection, migration.FUN = migration, get_prev.FUN = get_prev, verbose.FUN = verbose, depend = TRUE)
control <- control.net(type = "SI", nsims = 1, nsteps = 5000, initialize.FUN = initialize.net.mig,
deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
migrahetion.FUN = migration,
get_prev.FUN = get_prev_example_4,
verbose.FUN = verbose_example, depend = TRUE)
init <- init.net(i.num = 0)
control <- control.net(type = "SI", nsims = 1, nsteps = 1, initialize.FUN = initialize.net.mig,
deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
migrahetion.FUN = migration,
get_prev.FUN = get_prev_example_4,
verbose.FUN = verbose_example, depend = TRUE)
#status.vector <- rbinom(5000, 1, 0.1)
#status.vector <- replace(status.vector, status.vector == 0, "s")
#status.vector <- replace(status.vector, status.vector == 1, "i")
init <- init.net(i.num = 0)
mod <- netsim(est_migrations, param, init, control)
control <- control.net(type = "SI", nsims = 1, nsteps = 3, initialize.FUN = initialize.net.mig,
deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
migrahetion.FUN = migration,
get_prev.FUN = get_prev_example_4,
verbose.FUN = verbose_example, depend = TRUE)
#status.vector <- rbinom(5000, 1, 0.1)
#status.vector <- replace(status.vector, status.vector == 0, "s")
#status.vector <- replace(status.vector, status.vector == 1, "i")
init <- init.net(i.num = 0)
mod <- netsim(est_migrations, param, init, control)
get_prev_example_4
source("migration_restricted.R")
source("infection_restricted.R")
source("births_restricted.R")
source("deaths_restricted.R")
source("get_prev_restricted_example_4.R")
source("initialize_net_restricted.R")
source("verbose_restricted_example.R")
get_prev_example_4
nsim <- 1
#control <- control.net(type="SI", nsims = 1, nsteps = 200, initialize.FUN = initialize.net.mig, deaths.FUN = deaths, births.FUN = births, edges_correct.FUN = edges_correct, infection.FUN = infection, migration.FUN = migration, get_prev.FUN = get_prev, verbose.FUN = verbose, depend = TRUE)
control <- control.net(type = "SI", nsims = 1, nsteps = 5000, initialize.FUN = initialize.net.mig,
deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
migration.FUN = migration,
get_prev.FUN = get_prev_example_4,
verbose.FUN = verbose_example, depend = TRUE)
control <- control.net(type = "SI", nsims = 1, nsteps = 3, initialize.FUN = initialize.net.mig,
deaths.FUN = deaths, births.FUN = births, infection.FUN = infection,
migration.FUN = migration,
get_prev.FUN = get_prev_example_4,
verbose.FUN = verbose_example, depend = TRUE)
init <- init.net(i.num = 0)
mod <- netsim(est_migrations, param, init, control)
ls()
mod
sim1
names(mod)
mod$network
net <- network.collapse(mod, at=1)
net <- network.collapse(mod$network, at=1)
mod$network
mod$network$sim1
nd <- mod$network$sim1
net <- network.collapse(nd, at=1)
net
list.vertex.attributes(net)
mixingmatrix(net, "loc")
825+880+778
825/2483
mixingmatrix(net, "sex")
mixingmatrix(net, "type")
list.vertex.attributes(net)
mixingmatrix(net, "migmalemix.rural")
mixingmatrix(net, "migmalemix.urban")
as.edgelist(net)
net.el <- as.edgelist(net)
head(net.el)
(net%v%"type")[600:651]
(net.el)[600:650,]
(net%v%"type")[768]
(net%v%"type")[3561]
(net%v%"type")[3521]
(net%v%"type")[4261]
history
history()
savehistory("~/Projects/circularmigrations/epimodel/network/Restricted_Model/testing-history.R")
