rm(list=ls())
library(odesolve)
#library(simecol)


####################################################################
### Progression of Versions
####################################################################

## 9 Sep 2013:
## a. Change to 200 women initially infected to match network models
#  b. Change package to desolve
#  c. didn't actually change to "desolve" because of typo

## 21 Nov 2012: Fixed Bug in naming of new-infections compartment
## for NMR -- and no bug here in its computation

## 29 Oct 2012: This version works to calculate new incidence.
## Used Dobromir's suggestion of adding new compartments to track
## new infections. For prevalence information see "contactasact_recode.R"

## 28 Oct 2012: Add modification to compute incidence.

## 3 May 2012: After discussions with Dobromir, contact-rate for
## females split up and code cleaned up

## 24 April 2012: After conversation with Dobromir, changed the
## t_MM. and t_NM.  because both migrant and non-migrant men should be
## having 3 sex acts per week.

## 2 April 2012: possible error found and change made to
## tMM. tNM. tNF.???? -- no, was correct because average has to be 3
## per week.

## 30 March 2012: Add migration for late-stage individuals

### 6 March 2012: Change initial distribution of infectives to be in
### accordance with length of stage, instead of having all infectives
### be acute at the outset.

### 20 Jan 2012: Make new changes to make compatible with ODE model.
## 1. increase population sizes.
## 2. revise mean-contact rate.

### 16 Sep 2011: Version changed for manuscript
### This version  runs the deterministic model for
### contact defined as sexual act, as was in the MS thesis.

### 15 Sep 2011:
### Changes:
## 1. We add differential beta,
## for migrant and non-migrant men,
## as suggested by Steve
## 2. length of acute window changed to 12 weeks, according
## to Hollingsworth (2008): which says this stage lasts
## for 2.9 months.
## reflected in parameter gamma.
## added aids stage
## need better estimates for per coital act transmission
## in each of the three stages: acute (A), chronic (C), late-stage(L)

## ## 16 aug 2011:
## ## also test using the initial conditions in the thesis:
## ## only one rural female infected at the outset.
## ## with these alternate initial conditions, i.e.
## ## only one rural female infected at the outset,
## ## we see that all migration rates have about the same effect,
## ## but in the no migration case the disease equilibrates
## ## at a lower prevalence level.


# BIG DIFFERENCE HERE - ONE URBAN FEMALE AND ONE RURAL FEMALE INFECTED AT THE THE START OF THE SIMULATION. SO# DISEASE IS SYMMMETRICALLY INTRODUCED, AS OPPOSED TO MASTERS WHERE THE DISEASE WAS INTRODUCED ONLY IN ONE REGION. SO NOW EVEN IN THE CASE OF NO MIGRATION TOTAL PREVALENCE REACHED IS 100%. Every subpopulation now should have one infected individual.

####################################################################

####################################################################
### Function for modeling
####################################################################

mortality.take1.0b <- function(t, compartments, parms, filetosave){

### Urban State Variables

## Category MMU
s.mmu <- compartments[1]
a.mmu  <- compartments[2]
c.mmu <- compartments[3]
l.mmu <- compartments[4]
n.mmu <- s.mmu+a.mmu+c.mmu+l.mmu

##Category NFU
s.nfu <-compartments[5]
a.nfu <- compartments[6]
c.nfu <- compartments[7]
l.nfu <- compartments[8]
n.nfu <- s.nfu+a.nfu+c.nfu+l.nfu

##Category NMU
s.nmu <- compartments[9]
a.nmu <- compartments[10]
c.nmu <- compartments[11]
l.nmu <- compartments[12]
n.nmu <- s.nmu+a.nmu+c.nmu+l.nmu

### Rural State Variables

## Category MMR
s.mmr<-compartments[13]
a.mmr  <- compartments[14]
c.mmr <- compartments[15]
l.mmr <- compartments[16]
n.mmr <- s.mmr+a.mmr+c.mmr+l.mmr

## Category NFR 
s.nfr <-compartments[17]
a.nfr <- compartments[18]
c.nfr <- compartments[19]
l.nfr <- compartments[20]
n.nfr <- s.nfr+a.nfr+c.nfr+l.nfr

## Category NMR
s.nmr <- compartments[21]
a.nmr <- compartments[22]
c.nmr <- compartments[23]
l.nmr <- compartments[24]
n.nmr <- s.nmr+a.nmr+c.nmr+l.nmr

## to record new infections
new.inf.mmu <- compartments[25]
new.inf.nfu <- compartments[26]
new.inf.nmu <- compartments[27]
new.inf.mmr <- compartments[28]
new.inf.nfr <- compartments[29]
new.inf.nmr <- compartments[30]

with (as.list(parms),{

## contact-rate for females
## split by location and migration-status of men
t.nfu.m <- (t.mmu*n.mmu)/n.nfu
t.nfu.nm <- (t.nmu*n.nmu)/n.nfu
t.nfr.m <- (t.mmr*n.mmr)/n.nfr
t.nfr.nm <- (t.nmr*n.nmr)/n.nfr 
####################################################################
### Differential Equations
####################################################################

###For Urban Area
  ## migrant males urban (mmu)
  ds.mmu <- nu-(s.mmu*t.mmu*(a.nfu/n.nfu)*beta.a.m) -
    (s.mmu*t.mmu*(c.nfu/n.nfu)*beta.c.m) -
    (s.mmu*t.mmu*(l.nfu/n.nfu)*beta.l.m) - # for late-stage
      (delta*s.mmu) + (delta*s.mmr) - (mu*s.mmu) 

new.inf.mmu <- (s.mmu*t.mmu*(a.nfu/n.nfu)*beta.a.m) + 
    (s.mmu*t.mmu*(c.nfu/n.nfu)*beta.c.m) + # new infections in MMU
      (s.mmu*t.mmu*(l.nfu/n.nfu)*beta.l.m)  

  da.mmu <- (s.mmu*t.mmu*(a.nfu/n.nfu)*beta.a.m) +
    (s.mmu*t.mmu*(c.nfu/n.nfu)*beta.c.m) +
      (s.mmu*t.mmu*(l.nfu/n.nfu)*beta.l.m) - # late-stage 
      (delta*a.mmu) + (delta*a.mmr) - (gamma*a.mmu) - (mu*a.mmu)

  dc.mmu <- (gamma*a.mmu) - (delta*c.mmu) + (delta*c.mmr) - (mu*c.mmu) -
    (eta*c.mmu) # transition to late-stage
  dl.mmu <- (eta*c.mmu) - (mu*l.mmu) -
    (mu.d*l.mmu) - (delta*l.mmu) + (delta*l.mmr) # late-stag

  ## nfu : non-migrant female urban
  ds.nfu <- (2*nu) - (s.nfu*t.nfu.m*(a.mmu/n.mmu)*beta.a.m) -
    (s.nfu*t.nfu.m*(c.mmu/n.mmu)*beta.c.m) -
    (s.nfu*t.nfu.m*(l.mmu/n.mmu)*beta.l.m) - # late stage
      (s.nfu*t.nfu.nm*(a.nmu/n.nmu)*beta.a.n) -
        (s.nfu*t.nfu.nm*(c.nmu/n.nmu)*beta.c.n) -
        (s.nfu*t.nfu.nm*(l.nmu/n.nmu)*beta.l.n) - # late-stage
          (mu*s.nfu)  

new.inf.nfu <- (s.nfu*t.nfu.m/2*(a.mmu/n.mmu)*beta.a.m) +
  (s.nfu*t.nfu.m*(c.mmu/n.mmu)*beta.c.m) + # new infections in NFU
    (s.nfu*t.nfu.m*(l.mmu/n.mmu)*beta.l.m) + 
      (s.nfu*t.nfu.nm*(a.nmu/n.nmu)*beta.a.n) +
        (s.nfu*t.nfu.nm*(c.nmu/n.nmu)*beta.c.n) +
        (s.nfu*t.nfu.nm*(l.nmu/n.nmu)*beta.l.n) 

  da.nfu <- (s.nfu*t.nfu.m/2*(a.mmu/n.mmu)*beta.a.m) +
    (s.nfu*t.nfu.m*(c.mmu/n.mmu)*beta.c.m) +
    (s.nfu*t.nfu.m*(l.mmu/n.mmu)*beta.l.m) + # late-stage
      (s.nfu*t.nfu.nm*(a.nmu/n.nmu)*beta.a.n) +
        (s.nfu*t.nfu.nm*(c.nmu/n.nmu)*beta.c.n) +
        (s.nfu*t.nfu.nm*(l.nmu/n.nmu)*beta.l.n) - # late-stage
          (gamma*a.nfu) - (mu*a.nfu)

  dc.nfu <- (gamma*a.nfu) - (mu*c.nfu) -
    (eta*c.nfu) # transition to late-stage

  dl.nfu <- (eta*c.nfu)-(mu*l.nfu)-(mu.d*l.nfu) # late-stage compartment

  ### nmu: nonmigrant male urban
  ds.nmu <- nu - (s.nmu*t.nmu*(a.nfu/n.nfu)*beta.a.n) -
    (s.nmu*t.nmu*(c.nfu/n.nfu)*beta.c.n) -
      (s.nmu*t.nmu*(l.nfu/n.nfu)*beta.l.n) - #late-stage
        (mu*s.nmu)

new.inf.nmu <- (s.nmu*t.nmu*(a.nfu/n.nfu)*beta.a.n) + 
  (s.nmu*t.nmu*(c.nfu/n.nfu)*beta.c.n) + # new infections in NMU
    (s.nmu*t.nmu*(l.nfu/n.nfu)*beta.l.n) 

  da.nmu <- (s.nmu*t.nmu*(a.nfu/n.nfu)*beta.a.n) +
    (s.nmu*t.nmu*(c.nfu/n.nfu)*beta.c.n) +
      (s.nmu*t.nmu*(l.nfu/n.nfu)*beta.l.n)- # late-stage
        (gamma*a.nmu) - (mu*a.nmu)

  dc.nmu <- (gamma*a.nmu) - (mu*c.nmu) - (eta*c.nmu) # transition to late-stagexs

  dl.nmu <- (eta*c.nmu) - (mu*l.nmu) - (mu.d*l.nmu)  # late-stage compartment

#For Rural Area

# represents (mmr)
  ds.mmr <- nu-(s.mmr*t.mmr*(a.nfr/n.nfr)*beta.a.m) -
    (s.mmr*t.mmr*(c.nfr/n.nfr)*beta.c.m) -
    (s.mmr*t.mmr*(l.nfr/n.nfr)*beta.l.m) - # for late-stage
      (delta*s.mmr) + (delta*s.mmu) - (mu*s.mmr)

new.inf.mmr <- (s.mmr*t.mmr*(a.nfr/n.nfr)*beta.a.m) +
    (s.mmr*t.mmr*(c.nfr/n.nfr)*beta.c.m) +
    (s.mmr*t.mmr*(l.nfr/n.nfr)*beta.l.m) # new infections in MMR

  da.mmr <- (s.mmr*t.mmr*(a.nfr/n.nfr)*beta.a.m) +
    (s.mmr*t.mmr*(c.nfr/n.nfr)*beta.c.m) +
    (s.mmr*t.mmr*(l.nfr/n.nfr)*beta.l.m) - # late-stage 
      (delta*a.mmr) + (delta*a.mmu) - (gamma*a.mmr) - (mu*a.mmr)

  dc.mmr <- (gamma*a.mmr) - (delta*c.mmr) + (delta*c.mmu) - (mu*c.mmr)-
    (eta*c.mmr) # transition to late-stage
  
  dl.mmr <- (eta*c.mmr) - (mu*l.mmr) - (mu.d*l.mmr) - (delta*l.mmr) +
    (delta*l.mmu)
                                        # late-stage compartment,

  ## nfr : non-migrant female rural
  ds.nfr <- (2*nu) - (s.nfr*t.nfr.m*(a.mmr/n.mmr)*beta.a.m) -
    (s.nfr*t.nfr.m*(c.mmr/n.mmr)*beta.c.m) -
    (s.nfr*t.nfr.m*(l.mmr/n.nmr)*beta.l.m) - # late stage
      (s.nfr*t.nfr.nm*(a.nmr/n.nmr)*beta.a.n) -
        (s.nfr*t.nfr.nm*(c.nmr/n.nmr)*beta.c.n) -
        (s.nfr*t.nfr.nm*(l.nmr/n.nmr)*beta.l.n) - # late-stage
          (mu*s.nfr)

new.inf.nfr <- (s.nfr*t.nfr.m*(a.mmr/n.mmr)*beta.a.m) +
    (s.nfr*t.nfr.m*(c.mmr/n.mmr)*beta.c.m) +
    (s.nfr*t.nfr.m*(l.mmr/n.mmr)*beta.l.m) + # new infections in NFR
      (s.nfr*t.nfr.nm*(a.nmr/n.nmr)*beta.a.n) +
        (s.nfr*t.nfr.nm*(c.nmr/n.nmr)*beta.c.n) +
        (s.nfr*t.nfr.nm*(l.nmr/n.nmr)*beta.l.n)

  da.nfr <- (s.nfr*t.nfr.m*(a.mmr/n.mmr)*beta.a.m) +
    (s.nfr*t.nfr.m*(c.mmr/n.mmr)*beta.c.m) +
    (s.nfr*t.nfr.m*(l.mmr/n.mmr)*beta.l.m) + # late-stage
      (s.nfr*t.nfr.nm*(a.nmr/n.nmr)*beta.a.n) +
        (s.nfr*t.nfr.nm*(c.nmr/n.nmr)*beta.c.n) +
        (s.nfr*t.nfr.nm*(l.nmr/n.nmr)*beta.l.n) - # late-stage
          (gamma*a.nfr) - (mu*a.nfr)

  dc.nfr <- (gamma*a.nfr) - (mu*c.nfr) - (eta*c.nfr)
                                        # transition to late-stage

  dl.nfr <- (eta*c.nfr) - (mu*l.nfr)-(mu.d*l.nfr) # late-stage compartment
                                      # 30 March 2012: add migration
                                        # of late-stage people


    ## nmr: nonmigrant male rural
  ds.nmr <- nu - (s.nmr*t.nmr*(a.nfr/n.nfr)*beta.a.n) -
    (s.nmr*t.nmr*(c.nfr/n.nfr)*beta.c.n) -
      (s.nmr*t.nmr*(l.nfr/n.nfr)*beta.l.n) - #late-stage
        (mu*s.nmr)

new.inf.nmr <- (s.nmr*t.nmr*(a.nfr/n.nfr)*beta.a.n) +
  (s.nmr*t.nmr*(c.nfr/n.nfr)*beta.c.n) + # new infections in NMR
    (s.nmr*t.nmr*(l.nfr/n.nfr)*beta.l.n)

  da.nmr <- (s.nmr*t.nmr*(a.nfr/n.nfr)*beta.a.n) +
    (s.nmr*t.nmr*(c.nfr/n.nfr)*beta.c.n) +
      (s.nmr*t.nmr*(l.nfr/n.nfr)*beta.l.n)- # late-stage
        (gamma*a.nmr) - (mu*a.nmr)

  dc.nmr <- (gamma*a.nmr) - (mu*c.nmr) - (eta*c.nmr)
                                        # transition to late-stage

  dl.nmr <- (eta*c.nmr) - (mu*l.nmr) - (mu.d*l.nmr)
                                        # late-stage compartment

### Collect State Variables
list(c(ds.mmu,da.mmu,dc.mmu,dl.mmu,
       ds.nfu,da.nfu,dc.nfu,dl.nfu,
       ds.nmu, da.nmu, dc.nmu,dl.nmu,
	ds.mmr,da.mmr,dc.mmr,dl.mmr,
       ds.nfr,da.nfr,dc.nfr,dl.nfr,
       ds.nmr, da.nmr, dc.nmr, dl.nmr,
       new.inf.mmu, new.inf.nfu, new.inf.nmu,
       new.inf.mmr, new.inf.nfr, new.inf.nmr
       ))

}
)
}

####################################################################
### Time of Simulation
####################################################################
##times <- seq(0,10000,1)
times <- seq(0, 5e3, 1)
####################################################################

####################################################################
### Initial Population Sizes
####################################################################
n.s.mmu=625; n.a.mmu=0; n.c.mmu=0; n.l.mmu=0
# n.s.nfu=1200; n.a.nfu=50; n.c.nfu=0; n.l.nfu=0;
n.s.nfu=1150; n.a.nfu=(12/552*100); # nfu
n.c.nfu=(500/552*100); n.l.nfu=(40/552*100); # nfu
n.s.nmu=625; n.a.nmu=0; n.c.nmu=0; n.l.nmu=0;
n.s.mmr=625; n.a.mmr=0; n.c.mmr=0; n.l.mmr=0;
# n.s.nfr=1200; n.a.nfr=50; n.c.nfr=0; n.l.nfr=0;
n.s.nfr=1150; n.a.nfr=(12/552*100); # nfr
n.c.nfr=(500/552*100); n.l.nfr=(40/552*100); # nfr
n.s.nmr=625; n.a.nmr=0; n.c.nmr=0; n.l.nmr=0;
####################################################################

####################################################################
### Starting Values
####################################################################
xstart.manuscript <- c(s.mmu=n.s.mmu, a.mmu=n.a.mmu, c.mmu=n.c.mmu, l.mmu=n.l.mmu,
                       s.nfu=n.s.nfu,a.nfu=n.a.nfu,c.nfu=n.c.nfu,l.nfu=n.l.nfu,
                       s.nmu=n.s.nmu,a.nmu=n.a.nmu,c.nmu=n.c.nmu,l.nmu=n.l.nmu,
                       s.mmr=n.s.mmr,a.mmr=n.a.mmr,c.mmr=n.c.mmr,l.mmr=n.l.mmr,
                       s.nfr=n.s.nfr,a.nfr=n.a.nfr,c.nfr=n.c.nfr,l.nfr=n.l.nfr,
                       s.nmr=n.s.nmr,a.nmr=n.a.nmr,c.nmr=n.c.nmr, l.nmr=n.l.nmr,
                       new.inf.mmu=0, new.inf.nfu=0, new.inf.nmu=0,
                       new.inf.mmr=0, new.inf.nfr=0, new.inf.nmr=0)

####################################################################

####################################################################
### Biology
####################################################################

### Transmission Probabilites
chronic.prb <- 0.0007
acute.prb <- 0.0007*26
late.prb <- 0.0007*7

acute.dur <- 12
late.dur <- 40
chronic.dur <- 500

dur.100 <- 100
c <- 3 # number of sex-acts per week
n.100 <- 100 # duration of partnerships

## appropriately account for coital dilution -- very important
beta.a.m <- 1-(((1-acute.prb)^(acute.dur*c/2))*
               ((1-chronic.prb)^((n.100-acute.dur)*c/2))
               ) # migrant men, acute
beta.a.n <- 1-(((1-acute.prb)^(n.100*c))*
              ((1-chronic.prb)^((n.100-acute.dur)*c))
               )# non-migrant men, acute
beta.c.m <- 1-((1-chronic.prb)^(n.100*c/2)) # migrant men, chronic
beta.c.n <- 1-((1-chronic.prb)^(n.100*c)) # non-migrant men, chronic
beta.l.m <- 1-((1-late.prb)^(late.dur*c/2)) # migrant men, late 
beta.l.n <- 1-((1-late.prb)^(late.dur*c)) # non-migrant men, late

## beta.a.m <- 1-(1-acute.prb)^(dur.100*c)
## beta.a.n <- 1-(1-acute.prb)^(dur.100*c)
## beta.c.m <- 1-(1-chronic.prb)^(dur.100*c)
## beta.c.n <- 1-(1-chronic.prb)^(dur.100*c)
## beta.l.m <- 1-(1-late.prb)^(dur.100*c)
## beta.l.n <- 1-(1-late.prb)^(dur.100*c)

### length of stages
gamma <- 1/12 # acute stage 12 years
eta <- 1/500 # chronic stage 10 yrs
####################################################################

####################################################################
### Behavior: Mixing Information
####################################################################
mean.deg <- 1.2
mean.edges <- 5000*mean.deg/2 # should be divided by 2

n.100 <- 100 # duration of partnerships
## t.mmu <- (1/3*mean.edges)/(2500*n.100)
## t.mmr <- t.mmu
## t.nmu <- (1/6*mean.edges)/(2500*n.100)
## t.nmr <- t.nmu
t.mmu <- 0.016 # set up to be consistent with total mean degree of 1.2
t.mmr <- t.mmu # and 3 sex acts per week on average with active partners
t.nmu <- 0.0080 # which gives an average of 2.4 sex acts per week
t.nmr <- t.nmu # with all partners

### Derivation of t.mmu and t.mmr.
## There are a total of 3000 partnerships (ties) in the population.
## 2000 are accounted for by migrant-men, and 1000 are accounted for by non-migrant men. Since there are 1250 total non-migrant men in the population, the active mean-degree of non-migrant men is 1000/1250=0.8.

## On average, only half of the partnerships of migrant-men are active at
## any given time -- therefore the 1250 migrant-men in our population
## account for 1000 active partnerships. The active mean degree for migrant-men, therefore, is 1000/1250=0.8. Since the active mean-degrees for both
## migrant and non-migrant men is 0.8 -- the active mean degree
## for the population is 0.8. Thus there are 2000 active partnerships
## in the population on average.

## Considering an average partnership duration of 100, the partner-change
## rate, per unit-time, is 0.8/100=0.008.

###because 0.02*1.5*100=3 sex acts per week (in total) for migrant men
## and 0.01*3*100=3 sex acts per week for non-migrant men
####################################################################

####################################################################
### Vital Dynamics
####################################################################
mu <- 1/(45*52) # Modified to represent sexual lifespan
mu.d <- 1/40 # change duration of late stage to 40 weeks
nu <- 625/(45*52) # for new population size
####################################################################


####################################################################
### Unchanging Parameters
####################################################################
parms.exceptmigration <- c(t.mmu<-t.mmu,t.mmr<-t.mmr,
                           t.nmu<-t.nmu,t.nmr<-t.nmr,
                           beta.a.m<-beta.a.m,beta.c.m<-beta.c.m,
                           beta.l.m<-beta.l.m,
                           beta.a.n<-beta.a.n,beta.c.n<-beta.c.n,
                           beta.l.n<-beta.l.n,
                           gamma=gamma,
                           eta=eta,
                           mu=mu, mu.d=mu.d,nu=nu)
####################################################################


####################################################################
### Simulating Migration-Rates
####################################################################



## No Migrations 
### Migration Rate 0
delta <- 0
parms.mig0 <- c(parms.exceptmigration, delta)

mig0 <- as.data.frame(lsoda(xstart.manuscript,times,
                            mortality.take1.0b,
                            parms.mig0,
                            )) 


### Migration-Rate 10
delta <- 10
parms.mig10 <- c(parms.exceptmigration, delta)
mig10 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig10)
                           ) 


### Migration-Rate 7
delta <- 7
parms.mig7 <- c(parms.exceptmigration, delta)
mig7 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig7)
                           ) 


### Migration-Rate 4
delta <- 4
parms.mig4 <- c(parms.exceptmigration, delta)
mig4 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig4)
                           ) 

### Migration-Rate 3
delta <- 3
parms.mig3 <- c(parms.exceptmigration, delta)
mig3 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig3)
                           ) 


### Migration-Rate 1
delta <- 1
parms.mig1 <- c(parms.exceptmigration, delta)
mig1 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig1)
                           ) 


################################################    
### Experiments -- Non-Acute Phase migrations
################################################          

### Migration-Rate 20
delta <- 20
parms.mig20 <- c(parms.exceptmigration, delta)
mig20 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig20)
                           ) 

### Migration-Rate 30
delta <- 30
parms.mig30 <- c(parms.exceptmigration, delta)
mig30 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig30)
                           ) 

### Migration-Rate 40
delta <- 40
parms.mig40 <- c(parms.exceptmigration, delta)
mig40 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig40)
                           ) 
### Migration-Rate 60
delta <- 60
parms.mig60 <- c(parms.exceptmigration, delta)
mig60 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig60)
                           ) 

### Migration-Rate 80
delta <- 80
parms.mig80 <- c(parms.exceptmigration, delta)
mig80 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig80)
                           ) 

### Migration-Rate 100
delta <- 100
parms.mig100 <- c(parms.exceptmigration, delta)
mig100 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig100)
                           ) 


############################################################
### Incidence Summaries (over 5e3 time steps)
############################################################
sum(mig3[nrow(mig0),26:31]) # 3 step migration
sum(mig30[nrow(mig0),26:31]) # 30 step migration
############################################################
