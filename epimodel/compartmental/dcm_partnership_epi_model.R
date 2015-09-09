infCaP <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    #derived totals
    n.mmu <- s.mmu + a.mmu + c.mmu + l.mmu
    n.nmu <- s.nmu + a.nmu + c.nmu + l.nmu
    n.mmr <- s.mmr + a.mmr + c.mmr + l.mmr
    n.nmr <- s.nmr + a.nmr + c.nmr + l.nmr
    n.nfu <- s.nfu + a.nfu + c.nfu + l.nfu
    n.nfr <- s.nfr + a.nfr + c.nfr + l.nfr
    num <- n.mmu + n.mmr + n.nmu + n.nmr + n.nfu + n.nfr
    
    #varying parameters
    t.nfu.m <- ((t.mmu * n.mmu)/n.nfu)
    t.nfu.nm <- ((t.nmu * n.nmu)/n.nfu)
    t.nfr.m <- ((t.mmr * n.mmr)/n.nfr)
    t.nfr.nm <- ((t.nmr * n.nmr)/n.nfr)
    
    #main ODEs
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
    
    #output
    
    list(c(ds.mmu, da.mmu, dc.mmu, dl.mmu, ds.nfu, 
           da.nfu, dc.nfu, dl.nfu, ds.nmu, da.nmu, 
           dc.nmu, dl.nmu, ds.mmr, da.mmr, dc.mmr, 
           dl.mmr, ds.nfr, da.nfr, dc.nfr, dl.nfr, 
           ds.nmr, da.nmr, dc.nmr, dl.nmr), 
         num = num, prevalence=(num-s.mmu-s.nfu-s.nmu-s.nfr-s.mmr-s.nmr)/num)
    
  })
}

#per act transmission probabilities
pC <- 0.0007
pA <- 0.0007*26
pL <- 0.0007*7
n=3
d=100

param <- param.dcm(t.mmu = 2*2.4/(n*d), t.mmr = 2*2.4/(n*d), t.nmu = 2.4/(n*d), t.nmr = 2.4/(n*d),
                   beta.a.m = 1-(1-pA)^(n*d/2), beta.a.n = 1-(1-pA)^(n*d), beta.c.m = 1-(1-pC)^(n*d/2),
                   beta.c.n = 1-(1-pC)^(n*d), beta.l.m = 1-(1-pL)^(n*d/2), beta.l.n = 1-(1-pL)^(n*d),
                   gamma = 1/12, eta = 1/500, mu = 1/(45*52), mu.d = 1/40,
                   nu = 625/(45*52), delta = 1/3)

init <- init.dcm(s.mmu=625, a.mmu=0, c.mmu=0, l.mmu=0,
                 s.nfu=1000, a.nfu=(10/2)/500*200, c.nfu=(460/2)/500*200, l.nfu=((30)/2)/500*200,
                 s.nmu=625, a.nmu=0, c.nmu=0, l.nmu=0, 
                 s.mmr=625, a.mmr=0, c.mmr=0, l.mmr=0,
                 s.nfr=1000, a.nfr=(10/2)/500*200, c.nfr=(460/2)/500*200, l.nfr=(30/2)/500*200,
                 s.nmr=625, a.nmr=0, c.nmr=0, l.nmr=0)

control <- control.dcm(nsteps=5000, new.mod=infCaP)

mod2_3 <- dcm(param, init, control)
