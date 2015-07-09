## Committing to GitHub

infCaA <- function(t, t0, parms) {
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
    #urban
    ds.mmu <- v/8 - s.mmu*t.mmu*(a.nfu/n.nfu)*beta.a.m 
            - s.mmu*t.mmu*(c.nfu/n.nfu)*beta.c.m 
            - s.mmu*t.mmu*(l.nfu/n.nfu)*beta.l.m - delta*s.mmu
            + delta*s.mmr - mu*s.mmu
    
    da.mmu <- s.mmu*t.mmu*(a.nfu/n.nfu)*beta.a.m + s.mmu*t.mmu*(c.nfu/n.nfu)*beta.c.m
            + s.mmu*t.mmu*(l.nfu/n.nfu)*beta.l.m - delta*a.mmu + delta*a.mmr
            - gamma1*a.mmu - mu*a.mmu
    
    dc.mmu <- gamma1*a.mmu - gamma2*c.mmu - delta*c.mmu
            + delta*c.mmr - mu*c.mmu
    
    dl.mmu <- gamma2*c.mmu - delta*l.mmu + delta*l.mmr - (mu + mu.d)*l.mmu
    
    ds.nfu <- v/4 - s.nfu*t.nfu.m*(a.mmu/n.mmu)*beta.a.m - s.nfu*t.nfu.m*(c.mmu/n.mmu)*beta.c.m
            - s.nfu*t.nfu.m*(l.mmu/n.mmu)*beta.l.n - s.nfu*t.nfu.nm*(a.nmu/n.nmu)*beta.a.n
            - s.nfu*t.nfu.nm*(c.nmu/n.nmu)*beta.c.n - s.nfu*t.nfu.nm*(l.nmu/n.nmu)*beta.l.n
            - mu*s.nfu
    
    da.nfu <- s.nfu*t.nfu.m*(a.mmu/n.mmu)*beta.a.m + s.nfu*t.nfu.m*(c.mmu/n.mmu)*beta.c.m
            + s.nfu*t.nfu.m*(l.mmu/n.mmu)*beta.l.n + s.nfu*t.nfu.nm*(a.nmu/n.nmu)*beta.a.n
            + s.nfu*t.nfu.nm*(c.nmu/n.nmu)*beta.c.n + s.nfu*t.nfu.nm*(l.nmu/n.nmu)*beta.l.n
            - gamma1*a.nfu - mu*a.nfu
    
    dc.nfu <- gamma1*a.nfu - gamma2*c.nfu - mu*c.nfu
    
    dl.nfu <- gamma2*c.nfu - (mu + mu.d)*l.nfu
    
    ds.nmu <- v/8 - s.nmu*t.nmu*(a.nfu/n.nfu)*beta.a.n - s.nmu*t.nmu*(c.nfu/n.nfu)*beta.c.n
            - s.nmu*t.nmu*(l.nfu/n.nfu)*beta.l.n - mu*s.nmu
    
    da.nmu <- s.nmu*t.nmu*(a.nfu/n.nfu)*beta.a.n + s.nmu*t.nmu*(c.nfu/n.nfu)*beta.c.n
            + s.nmu*t.nmu*(l.nfu/n.nfu)*beta.l.n - mu*a.nmu - gamma1*a.nmu
    
    dc.nmu <- gamma1*a.nmu - gamma2*c.nmu - mu*c.nmu
    
    dl.nmu <- gamma2*c.nmu - (mu + mu.d)*l.nmu
    
    #rural
    ds.mmr <- v/8 - s.mmr*t.mmr*(a.nfr/n.nfr)*beta.a.m 
    - s.mmr*t.mmr*(c.nfr/n.nfr)*beta.c.m 
    - s.mmr*t.mmr*(l.nfr/n.nfr)*beta.l.m - delta*s.mmr
    + delta*s.mmu - mu*s.mmr
    
    da.mmr <- s.mmr*t.mmr*(a.nfr/n.nfr)*beta.a.m + s.mmr*t.mmr*(c.nfr/n.nfr)*beta.c.m
    + s.mmr*t.mmr*(l.nfr/n.nfr)*beta.l.m - delta*a.mmr + delta*a.mmu
    - gamma1*a.mmr - mu*a.mmr
    
    dc.mmr <- gamma1*a.mmr - gamma2*c.mmr - delta*c.mmr
    + delta*c.mmu - mu*c.mmr
    
    dl.mmr <- gamma2*c.mmr - delta*l.mmr + delta*l.mmu - (mu + mu.d)*l.mmr
    
    ds.nfr <- v/4 - s.nfr*t.nfr.m*(a.mmr/n.mmr)*beta.a.m - s.nfr*t.nfr.m*(c.mmr/n.mmr)*beta.c.m
    - s.nfr*t.nfr.m*(l.mmr/n.mmr)*beta.l.m - s.nfr*t.nfr.nm*(a.nmr/n.nmr)*beta.a.n
    - s.nfr*t.nfr.nm*(c.nmr/n.nmr)*beta.c.n - s.nfr*t.nfr.nm*(l.nmr/n.nmr)*beta.l.n
    - mu*s.nfr
    
    da.nfr <- s.nfr*t.nfr.m*(a.mmr/n.mmr)*beta.a.m + s.nfr*t.nfr.m*(c.mmr/n.mmr)*beta.c.m
    + s.nfr*t.nfr.m*(l.mmr/n.mmr)*beta.l.m + s.nfr*t.nfr.nm*(a.nmr/n.nmr)*beta.a.n
    + s.nfr*t.nfr.nm*(c.nmr/n.nmr)*beta.c.n + s.nfr*t.nfr.nm*(l.nmr/n.nmr)*beta.l.n
    - gamma1*a.nfr - mu*a.nfr
    
    dc.nfr <- gamma1*a.nfr - gamma2*c.nfr - mu*c.nfr
    
    dl.nfr <- gamma2*c.nfr - (mu + mu.d)*l.nfr
    
    ds.nmr <- v/8 - s.nmr*t.nmr*(a.nfr/n.nfr)*beta.a.n - s.nmr*t.nmr*(c.nfr/n.nfr)*beta.c.n
    - s.nmr*t.nmr*(l.nfr/n.nfr)*beta.l.n - mu*s.nmr
    
    da.nmr <- s.nmr*t.nmr*(a.nfr/n.nfr)*beta.a.n + s.nmr*t.nmr*(c.nfr/n.nfr)*beta.c.n
    + s.nmr*t.nmr*(l.nfr/n.nfr)*beta.l.n - mu*a.nmr - gamma1*a.nmr
    
    dc.nmr <- gamma1*a.nmr - gamma2*c.nmr - mu*c.nmr
    
    dl.nmr <- gamma2*c.nmr - (mu + mu.d)*l.nmr
    
    #output
    
    list(c(ds.mmu, da.mmu, dc.mmu, dl.mmu, ds.nfu, 
           da.nfu, dc.nfu, dl.nfu, ds.nmu, da.nmu, 
           dc.nmu, dl.nmu, ds.mmr, da.mmr, dc.mmr, 
           dl.mmr, ds.nfr, da.nfr, dc.nfr, dl.nfr, 
           ds.nmr, da.nmr, dc.nmr, dl.nmr), 
           num = num, prevalence=(num-s.mmu-s.nfu-s.nmu-s.nfr-s.mmr-s.nmr)/num)
    
  })
}

param <- param.dcm(t.mmu = 2.4, t.mmr = 2.4, t.nmu = 2.4, t.nmr = 2.4, 
                   beta.a.m = 0.0007*26, beta.a.n = 0.0007*26, beta.c.m = 0.0007,
                   beta.c.n = 0.0007, beta.l.m = 0.0007*7, beta.l.n = 0.0007*7,
                   gamma1 = 1/12, gamma2 = 1/500, mu = 1/(45*52), mu.d = 1/40,
                   v = 625*8/(45*52), delta = 1/30)

init <- init.dcm(s.mmu=625, a.mmu=0, c.mmu=0, l.mmu=0,
                 s.nfu=1000, a.nfu=(10/500*100)/2, c.nfu=(460/500*100)/2, l.nfu=(30/500*100)/2,
                 s.nmu=625, a.nmu=0, c.nmu=0, l.nmu=0, 
                 s.mmr=625, a.mmr=0, c.mmr=0, l.mmr=0,
                 s.nfr=1000, a.nfr=(10/500*100)/2, c.nfr=(460/500*100)/2, l.nfr=(30/500*100)/2,
                 s.nmr=625, a.nmr=0, c.nmr=0, l.nmr=0)

control <- control.dcm(nsteps=1500, new.mod=infCaA)

mod <- dcm(param, init, control)
