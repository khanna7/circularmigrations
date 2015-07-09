## 9 Sep 2013:
## a. Change to 200 women initially infected to match network models
#  b. Change package to desolve
#  c. didn't actually change to "desolve" because of typo

## 23 Dec 2012: Revised to get prevalence graph with title

## 28 Oct 2012:Version used for thesis/defense/manuscript.

rm(list=ls())
 library(deSolve)
# library(desolve) #9 Sep 2013

mortality.take1.0b <- function(t, compartments, parms){

#######################################  
### Urban State Variables
#######################################    

# Category MMU
s.mmu <- compartments[1]
a.mmu  <- compartments[2]
c.mmu <- compartments[3]
l.mmu <- compartments[4]
n.mmu <- s.mmu+a.mmu+c.mmu+l.mmu

#Category NFU
s.nfu <-compartments[5]
a.nfu <- compartments[6]
c.nfu <- compartments[7]
l.nfu <- compartments[8]
n.nfu <- s.nfu+a.nfu+c.nfu+l.nfu

#Category NMU
s.nmu <- compartments[9]
a.nmu <- compartments[10]
c.nmu <- compartments[11]
l.nmu <- compartments[12]
n.nmu <- s.nmu+a.nmu+c.nmu+l.nmu

# Rural State Variables

# Category MMR
s.mmr<-compartments[13]
a.mmr  <- compartments[14]
c.mmr <- compartments[15]
l.mmr <- compartments[16]
n.mmr <- s.mmr+a.mmr+c.mmr+l.mmr

# Category NFR 
s.nfr <-compartments[17]
a.nfr <- compartments[18]
c.nfr <- compartments[19]
l.nfr <- compartments[20]
n.nfr <- s.nfr+a.nfr+c.nfr+l.nfr

# Category NMR
s.nmr <- compartments[21]
a.nmr <- compartments[22]
c.nmr <- compartments[23]
l.nmr <- compartments[24]
n.nmr <- s.nmr+a.nmr+c.nmr+l.nmr

with (as.list(parms),{

  ## contact-rate for females
  t.nfu.m <- ((t.mmu*n.mmu))/n.nfu
  t.nfu.nm <- (t.nmu*n.nmu)/n.nfu
  t.nfr.m <- (t.mmr*n.mmr)/n.nfr 
  t.nfr.nm <- (t.nmr*n.nmr)/n.nfr

###For Urban Area
  ## migrant males urban (mmu)
  ds.mmu <- nu-(s.mmu*t.mmu*(a.nfu/n.nfu)*beta.a.m) -
    (s.mmu*t.mmu*(c.nfu/n.nfu)*beta.c.m) -
    (s.mmu*t.mmu*(l.nfu/n.nfu)*beta.l.m) - # for late-stage
      (delta*s.mmu) + (delta*s.mmr) - (mu*s.mmu) 

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
    (s.nfr*t.nfr.m*(l.mmr/n.mmr)*beta.l.m) - # late stage
      (s.nfr*t.nfr.nm*(a.nmr/n.nmr)*beta.a.n) -
        (s.nfr*t.nfr.nm*(c.nmr/n.nmr)*beta.c.n) -
        (s.nfr*t.nfr.nm*(l.nmr/n.nmr)*beta.l.n) - # late-stage
          (mu*s.nfr)

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
       ds.nmr, da.nmr, dc.nmr, dl.nmr
))

}
)
}

#######################################        
### Time of Simulation
#######################################
times <- seq(0, 1e4, 1) 
#######################################

#######################################        
### Initial population sizes
#######################################  
n.s.mmu=625; n.a.mmu=0; n.c.mmu=0; n.l.mmu=0 # mmu
## n.s.nfu=1200; n.a.nfu=(12/552*50); # nfu
## n.c.nfu=(500/552*50); n.l.nfu=(40/552*50); #9Sep13 -- changed initial prev # nfu
n.s.nfu=1150; n.a.nfu=(12/552*100); # nfu
n.c.nfu=(500/552*100); n.l.nfu=(40/552*100); # nfu
n.s.nmu=625; n.a.nmu=0; n.c.nmu=0; n.l.nmu=0; # nfu
n.s.mmr=625; n.a.mmr=0; n.c.mmr=0; n.l.mmr=0; # nfu
## n.s.nfr=1200; n.a.nfr=(12/552*50); #nfr 
## n.c.nfr=(500/552*50); n.l.nfr=(40/552*50); #nfr
n.s.nfr=1150; n.a.nfr=(12/552*100); # nfr
n.c.nfr=(500/552*100); n.l.nfr=(40/552*100); # nfr
n.s.nmr=625; n.a.nmr=0; n.c.nmr=0; n.l.nmr=0; #nmr

#######################################        
### Provide Starting Values
#######################################  

xstart.manuscript <- c(s.mmu=n.s.mmu, a.mmu=n.a.mmu,c.mmu=n.c.mmu, l.mmu=n.l.mmu,
                       s.nfu=n.s.nfu,a.nfu=n.a.nfu,c.nfu=n.c.nfu,l.nfu=n.l.nfu,
                       s.nmu=n.s.nmu,a.nmu=n.a.nmu,c.nmu=n.c.nmu,l.nmu=n.l.nmu,
                       s.mmr=n.s.mmr,a.mmr=n.a.mmr,c.mmr=n.c.mmr,l.mmr=n.l.mmr,
                       s.nfr=n.s.nfr,a.nfr=n.a.nfr,c.nfr=n.c.nfr,l.nfr=n.l.nfr,
                       s.nmr=n.s.nmr,a.nmr=n.a.nmr,c.nmr=n.c.nmr, l.nmr=n.l.nmr)

#######################################        
### Behavioral Informaion
#######################################  

##n.100 <- 100 # Average partnership duratioan = 100 weeks
##c <- 3
##c <- 2

#######################################        
### Biological Informaion
#######################################  
### per act transmission probabiliites
chronic.prb <- 0.0007
acute.prb <- chronic.prb*26 # from hollingsworth
late.prb <- chronic.prb*7 # from hollingsworth

### duration of various stages    
acute.dur <- 12
chronic.dur <- 500
late.dur <- 40 # 6 March 2012: changed in accordance with Hollingsworth

### True transmission values for contact as partnership      

beta.a.m <- acute.prb
beta.a.n <- acute.prb
beta.c.m <- chronic.prb
beta.c.n <- chronic.prb
beta.l.m <- late.prb
beta.l.n <- late.prb


beta.a.m <- acute.prb
beta.a.n <- acute.prb
beta.c.m <- chronic.prb
beta.c.n <- chronic.prb
beta.l.m <- late.prb
beta.l.n <- late.prb

################################################    
### Mixing Information
################################################          

## t.mmu <- 3
## t.mmr <- 3
## t.nmu <- 3
## t.nmr <- 3

t.mmu <- 2.4 # has to be changed to 2.4 to be comparable with 
t.mmr <- 2.4 # contact as partnership model
t.nmu <- 2.4
t.nmr <- 2.4

## t.mmu<-2
## t.mmr<-2
## t.nmu<-2
## t.nmr<-2
################################################

################################################    
### Infection-Transition Probabilities &
### Vital Dynamics
################################################          
gamma <- 1/12
eta <- 1/500 # chronic stage fixed at 10 yrs

mu <- 1/(45*52) # Modified to represent sexual lifespan
mu.d <- 1/40 # change late-stage to 40 weeks
nu <- 625/(45*52) # birth-rate

################################################          

################################################    
### Parameters (except migration rate)
################################################          

parms.exceptmigration <- c(t.mmu<-t.mmu,t.mmr<-t.mmr,
                           t.nmu<-t.nmu,t.nmr<-t.nmr,
                           beta.a.m<-beta.a.m,beta.c.m<-beta.c.m,
                           beta.l.m<-beta.l.m,
                           beta.a.n<-beta.a.n,beta.c.n<-beta.c.n,
                           beta.l.n<-beta.l.n,
                           gamma=gamma, eta=eta,
                           mu=mu, mu.d=mu.d,nu=nu)

################################################    
### Experiments -- Acute Phase migrations
################################################          

### No Migrations
### Migration-Rate 0 # does not apply
delta <- 0
parms.mig0 <- c(parms.exceptmigration, delta)
mig0 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig0)
                           ) 

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

################################################
### Acute Phase -- Combine Data
################################################

## mig-step 10
susceptibles.10 <- mig10$s.mmu+mig10$s.nfu+mig10$s.nmu+
  mig10$s.mmr+mig10$s.nfr+mig10$s.nmr
acute.10 <- mig10$a.mmu+mig10$a.nfu+mig10$a.nmu+
  mig10$a.mmr+mig10$a.nfr+mig10$a.nmr
chronic.10 <- mig10$c.mmu+mig10$c.nfu+mig10$c.nmu+
  mig10$c.mmr+mig10$c.nfr+mig10$c.nmr
late.10 <- mig10$l.mmu+mig10$l.nfu+mig10$l.nmu+
  mig10$l.mmr+mig10$l.nfr+mig10$l.nmr
infected.10 <- acute.10+chronic.10+late.10

## mig-step 7
susceptibles.7 <- mig7$s.mmu+mig7$s.nfu+mig7$s.nmu+
  mig7$s.mmr+mig7$s.nfr+mig7$s.nmr
acute.7 <- mig7$a.mmu+mig7$a.nfu+mig7$a.nmu+
  mig7$a.mmr+mig7$a.nfr+mig7$a.nmr
chronic.7 <- mig7$c.mmu+mig7$c.nfu+mig7$c.nmu+
  mig7$c.mmr+mig7$c.nfr+mig7$c.nmr
late.7 <- mig7$l.mmu+mig7$l.nfu+mig7$l.nmu+
  mig7$l.mmr+mig7$l.nfr+mig7$l.nmr
infected.7 <- acute.7+chronic.7+late.7

## mig-step 4
susceptibles.4 <- mig4$s.mmu+mig4$s.nfu+mig4$s.nmu+
  mig4$s.mmr+mig4$s.nfr+mig4$s.nmr
acute.4 <- mig4$a.mmu+mig4$a.nfu+mig4$a.nmu+
  mig4$a.mmr+mig4$a.nfr+mig4$a.nmr
chronic.4 <- mig4$c.mmu+mig4$c.nfu+mig4$c.nmu+
  mig4$c.mmr+mig4$c.nfr+mig4$c.nmr
late.4 <- mig4$l.mmu+mig4$l.nfu+mig4$l.nmu+
  mig4$l.mmr+mig4$l.nfr+mig4$l.nmr

infected.4 <- acute.4+chronic.4+late.4

## mig-step 3
susceptibles.3 <- mig3$s.mmu+mig3$s.nfu+mig3$s.nmu+
  mig3$s.mmr+mig3$s.nfr+mig3$s.nmr
acute.3 <- mig3$a.mmu+mig3$a.nfu+mig3$a.nmu+
  mig3$a.mmr+mig3$a.nfr+mig3$a.nmr
chronic.3 <- mig3$c.mmu+mig3$c.nfu+mig3$c.nmu+
  mig3$c.mmr+mig3$c.nfr+mig3$c.nmr
late.3 <- mig3$l.mmu+mig3$l.nfu+mig3$l.nmu+
  mig3$l.mmr+mig3$l.nfr+mig3$l.nmr

infected.3 <- acute.3+chronic.3+late.3


## mig-step 1
susceptibles.1 <- mig1$s.mmu+mig1$s.nfu+mig1$s.nmu+
  mig1$s.mmr+mig1$s.nfr+mig1$s.nmr
acute.1 <- mig1$a.mmu+mig1$a.nfu+mig1$a.nmu+
  mig1$a.mmr+mig1$a.nfr+mig1$a.nmr
chronic.1 <- mig1$c.mmu+mig1$c.nfu+mig1$c.nmu+
  mig1$c.mmr+mig1$c.nfr+mig1$c.nmr
late.1 <- mig1$l.mmu+mig1$l.nfu+mig1$l.nmu+
  mig1$l.mmr+mig1$l.nfr+mig1$l.nmr

infected.1 <- acute.1+chronic.1+late.1

################################################
### Non-Acute Phase -- Combine Data
################################################

## mig-step 20
susceptibles.20 <- mig20$s.mmu+mig20$s.nfu+mig20$s.nmu+
  mig20$s.mmr+mig20$s.nfr+mig20$s.nmr
acute.20 <- mig20$a.mmu+mig20$a.nfu+mig20$a.nmu+
  mig20$a.mmr+mig20$a.nfr+mig20$a.nmr
chronic.20 <- mig20$c.mmu+mig20$c.nfu+mig20$c.nmu+
  mig20$c.mmr+mig20$c.nfr+mig20$c.nmr
late.20 <- mig20$l.mmu+mig20$l.nfu+mig20$l.nmu+
  mig20$l.mmr+mig20$l.nfr+mig20$l.nmr

infected.20 <- acute.20+chronic.20+late.20

## mig-step 30
susceptibles.30 <- mig30$s.mmu+mig30$s.nfu+mig30$s.nmu+
  mig30$s.mmr+mig30$s.nfr+mig30$s.nmr
acute.30 <- mig30$a.mmu+mig30$a.nfu+mig30$a.nmu+
  mig30$a.mmr+mig30$a.nfr+mig30$a.nmr
chronic.30 <- mig30$c.mmu+mig30$c.nfu+mig30$c.nmu+
  mig30$c.mmr+mig30$c.nfr+mig30$c.nmr
late.30 <- mig30$l.mmu+mig30$l.nfu+mig30$l.nmu+
  mig30$l.mmr+mig30$l.nfr+mig30$l.nmr

infected.30 <- acute.30+chronic.30+late.30

## mig-step 40

susceptibles.40 <- mig40$s.mmu+mig40$s.nfu+mig40$s.nmu+
  mig40$s.mmr+mig40$s.nfr+mig40$s.nmr
acute.40 <- mig40$a.mmu+mig40$a.nfu+mig40$a.nmu+
  mig40$a.mmr+mig40$a.nfr+mig40$a.nmr
chronic.40 <- mig40$c.mmu+mig40$c.nfu+mig40$c.nmu+
  mig40$c.mmr+mig40$c.nfr+mig40$c.nmr
late.40 <- mig40$l.mmu+mig40$l.nfu+mig40$l.nmu+
  mig40$l.mmr+mig40$l.nfr+mig40$l.nmr

infected.40 <- acute.40+chronic.40+late.40

## mig-step 60
susceptibles.60 <- mig60$s.mmu+mig60$s.nfu+mig60$s.nmu+
  mig60$s.mmr+mig60$s.nfr+mig60$s.nmr
acute.60 <- mig60$a.mmu+mig60$a.nfu+mig60$a.nmu+
  mig60$a.mmr+mig60$a.nfr+mig60$a.nmr
chronic.60 <- mig60$c.mmu+mig60$c.nfu+mig60$c.nmu+
  mig60$c.mmr+mig60$c.nfr+mig60$c.nmr
late.60 <- mig60$l.mmu+mig60$l.nfu+mig60$l.nmu+
  mig60$l.mmr+mig60$l.nfr+mig60$l.nmr

infected.60 <- acute.60+chronic.60+late.60

## mig-step 80
susceptibles.80 <- mig80$s.mmu+mig80$s.nfu+mig80$s.nmu+
  mig80$s.mmr+mig80$s.nfr+mig80$s.nmr
acute.80 <- mig80$a.mmu+mig80$a.nfu+mig80$a.nmu+
  mig80$a.mmr+mig80$a.nfr+mig80$a.nmr
chronic.80 <- mig80$c.mmu+mig80$c.nfu+mig80$c.nmu+
  mig80$c.mmr+mig80$c.nfr+mig80$c.nmr
late.80 <- mig80$l.mmu+mig80$l.nfu+mig80$l.nmu+
  mig80$l.mmr+mig80$l.nfr+mig80$l.nmr

infected.80 <- acute.80+chronic.80+late.80

## mig-step 100

susceptibles.100 <- mig100$s.mmu+mig100$s.nfu+mig100$s.nmu+
  mig100$s.mmr+mig100$s.nfr+mig100$s.nmr
acute.100 <- mig100$a.mmu+mig100$a.nfu+mig100$a.nmu+
  mig100$a.mmr+mig100$a.nfr+mig100$a.nmr
chronic.100 <- mig100$c.mmu+mig100$c.nfu+mig100$c.nmu+
  mig100$c.mmr+mig100$c.nfr+mig100$c.nmr
late.100 <- mig100$l.mmu+mig100$l.nfu+mig100$l.nmu+
  mig100$l.mmr+mig100$l.nfr+mig100$l.nmr

infected.100 <- acute.100+chronic.100+late.100

################################################
### No Migrations -- Combine data
################################################
susceptibles.0 <- mig0$s.mmu+mig0$s.nfu+mig0$s.nmu+
  mig0$s.mmr+mig0$s.nfr+mig0$s.nmr
acute.0 <- mig0$a.mmu+mig0$a.nfu+mig0$a.nmu+
  mig0$a.mmr+mig0$a.nfr+mig0$a.nmr
chronic.0 <- mig0$c.mmu+mig0$c.nfu+mig0$c.nmu+
  mig0$c.mmr+mig0$c.nfr+mig0$c.nmr
late.0 <- mig0$l.mmu+mig0$l.nfu+mig0$l.nmu+
  mig0$l.mmr+mig0$l.nfr+mig0$l.nmr

infected.0 <- acute.0+chronic.0+late.0

################################################
### plots
################################################

pdf("contactasact_prevalences.pdf")

plot(infected.1/(infected.1+susceptibles.1),
     type="l", ylim=c(0,1),
     xlab="Time (Weeks)",
     ylab="Prevalence")
lines(infected.3/(infected.3+susceptibles.3))
lines(infected.4/(infected.4+susceptibles.4))
lines(infected.7/(infected.7+susceptibles.7))
lines(infected.10/(infected.10+susceptibles.10))

lines(infected.20/(infected.20+susceptibles.20), lty=2)
lines(infected.30/(infected.30+susceptibles.30), lty=2)
lines(infected.40/(infected.40+susceptibles.40), lty=2)
lines(infected.60/(infected.60+susceptibles.60), lty=2)
lines(infected.80/(infected.80+susceptibles.80), lty=2)
lines(infected.100/(infected.100+susceptibles.100), lty=2)

legend("topleft", c("Frequent Migrations",
                    "Infrequent Migrations"),
                    lty=1:2, inset=0.05)
dev.off()

pdf("nomigrations_contactasact_prevalences.pdf")

plot(infected.0/(infected.0+susceptibles.0),
     type="l", ylim=c(0,1),
     xlab="Time (Weeks)",
     ylab="Prevalence")

dev.off()

## for defense
pdf("contactasact_prevalences_fordefense.pdf")

plot(infected.3/(infected.3+susceptibles.3),
     type="l", ylim=c(0, 0.5),
     xlab="Time (Weeks)",
     lwd=2,
     ylab="Prevalence",
     cex.lab=1.35,
     cex.axis=1.25)

lines(infected.30/(infected.20+susceptibles.30),
      lty=2, lwd=2)

legend("topright", c("3 weeks",
                    "30 weeks"),
       title="Migration Interval",
       lty=1:2, inset=0.05)

dev.off()

save.image(file="contactasact.RData")

### 23 Dec 2012: Revised to get prevalence graph with tile
load(file="contactasact.RData")

pdf("contactasact_prevalences_wtitle.pdf")

plot(infected.3/(infected.3+susceptibles.3),
     type="l", ylim=c(0, 0.5),
     xlab="Time (Weeks)",
     lwd=2,
     ylab="Prevalence",
     cex=1.35,
     main="Contact as Act")

lines(infected.30/(infected.30+susceptibles.30),#9Sep13: error found: said
                                        # infected.20 -- but doesn't change results
      lty=2, lwd=2)

legend("topright", c("3 weeks",
                    "30 weeks"),
       title="Migration Interval",
       lty=1:2, inset=0.05)

dev.off()
