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
n.s.nfu=1150; n.a.nfu=(12/552*100); # nfu
n.c.nfu=(500/552*100); n.l.nfu=(40/552*100); # nfu
n.s.nmu=625; n.a.nmu=0; n.c.nmu=0; n.l.nmu=0; # nfu
n.s.mmr=625; n.a.mmr=0; n.c.mmr=0; n.l.mmr=0; # nfu
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


t.mmu <- 2.4 # has to be changed to 2.4 to be comparable with 
t.mmr <- 2.4 # contact as partnership model
t.nmu <- 2.4
t.nmr <- 2.4

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


### Migration-Rate 3
delta <- 1/3
parms.mig3 <- c(parms.exceptmigration, delta)
mig3 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig3)
                           ) 

################################################    
### Experiments -- Non-Acute Phase migrations
################################################          

### Migration-Rate 30
delta <- 1/30
parms.mig30 <- c(parms.exceptmigration, delta)
mig30 <- as.data.frame(lsoda(xstart.manuscript,
                            times,
                            mortality.take1.0b,
                            parms.mig30)
                           ) 

################################################
### Acute Phase -- Combine Data
################################################


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



################################################
### Non-Acute Phase -- Combine Data
################################################

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



################################################
### plots
################################################

pdf("contactasact_prevalences.pdf")

plot(infected.3/(infected.3+susceptibles.3),
     type="l", ylim=c(0,1),
     xlab="Time (Weeks)",
     ylab="Prevalence")
lines(infected.3/(infected.3+susceptibles.3))

legend("topleft", c("Frequent Migrations",
                    "Infrequent Migrations"),
                    lty=1:2, inset=0.05)
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

lines(infected.30/(infected.30+susceptibles.30),
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
