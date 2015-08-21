rm(list=ls())
set.seed(Sys.time())
library(ergm)

## 2 March 2012: Organized in folder for final runs

## 30 Jan 2012: distribution of infecteds has to be redefined
## because we late-stage has been reduced to 9 months

## 27 Jan 2012:  Add new ergm.term to have a 'restricted' partnership
## structure for migrant-men.

## 20 Jan 2012: Have 'time.infected' (time since infection)
## for infected people be uniformly
## distributed between 0 and life-span of infected individuals (600 weeks)

### 13 Jan 2012: Fix mean degree after discussion with Steve
### Still have 10 urban and 10 rural infected at outset

## 10 Jan 2012: Fix dyn interval
## after discussion with Steve

## Need to move away from
## bipartite graphs I think to
## do births and deaths correctly

### Set up population
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

## 15th element in mixing matrix not specified to avoind
## collinearity and other errors later

### Other Network Information
diss.network <- ~edges
gamma.network <- log(duration-1)
## constraints.network <- ~.
constraints.network <- ~.
constraints.degree <- ~bd(maxout=1) # not used
##dyninterval <- 1e4
dyninterval <- 1e6
maxit=50
burnin.sim.network=25e3
#burnin=25e3
burnin=25e4 # after discssion with steve

## Initiate Network
# san.obj <- n0 ~  edges+nodemix("type", base=c(15))+twopath
san.obj <- n0 ~  edges+nodemix("type", base=c(15))+
  kstar(2:9, "migmalemix.rural")+
  kstar(2:9, "migmalemix.urban")

meanstats  <- c(n.edges, MF.type.mix, 
                      c(#kstar.mmr.mix
                        rep(kstar.mmr.mix.2star,8)#,
                        #kstar.mmr.mix.2star
                        ),
                        c(#kstar.mmu.mix
                        rep(kstar.mmu.mix.2star, 8)
                          )
                       )


init.network <- san(san.obj, meanstats=meanstats,
                    interval=1000, # from ampham
                    MCMC.samplesize=1e4, # from ampham
                    burnin=burnin, # from ampham
                    maxit=maxit,
                    constraints=constraints.network)

## estimate  network
formula.network <- init.network ~ edges+nodemix("type", base=c(15))+
  kstar(2:9, "migmalemix.rural")+
  kstar(2:9, "migmalemix.urban")

migration.fit.network <- ergm(formula=formula.network,
                              control=control.ergm(dyninterval=
                                dyninterval), # 10 jan '12: from steve
                              meanstats=meanstats,
                              interval=1000, # from ampham
                              MCMC.samplesize=1e4, # from ampham
                              burnin=burnin, # from ampham
                              maxit=maxit, # from ampham
                              verbose=T) # from ampham


## simulate from cross-sectional network (to test if mean-stats make sense)
## sim.networks.take1<-simulate(migration.fit.network,
##                              nsim=100,statsonly=FALSE,sequential=TRUE)

## apply(sim.networks.take1$stats,2,mean)
## apply nicole's approximation
theta0 <- migration.fit.network$coef
theta0[1] <- theta0[1]-gamma.network
migration.fit.network$coef <- theta0
                              
## Make offset
offset.edges.network <- -log(network.size(init.network))
theta.network <- c((migration.fit.network$coef[1]-# error was right here!!!!
                    offset.edges.network),
                   migration.fit.network$coef[-1], 
                   offset.edges.network
                   )

names(theta.network) <- c(names(ergm(formula.network, MPLEonly=T)$coef),
                          "offset.edges"
                          )

formula.network.w.offset <- update.formula(formula.network,
                                           ~.+offset(edges)
                                           )

## check attributes in estimated network via
## net <- migration.fit.network$network
## get.vertex.attribut(net, "")

save.image(file="est_dur100_fixinterval_newterm.RData")
