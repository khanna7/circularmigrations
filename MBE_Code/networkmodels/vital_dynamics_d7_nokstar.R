## 27 Jan 2012: Added new attributes to account for kstart terms
## in model

## 17 Jan 2012: Dying due to AIDS was not working because
## the 'time.infected' attribute was not getting updated in the
## proper place. This has now been fixed.

## 13 Jan 2012: After discussion with Steve,
## time.infected.fatal has been re-introduced.
## Total life-span after infection is 616 weeks.

## 3 Jan 2012: Following discussion with Steve, decided
## to reduce life-expectancy only dor late-stage people
## and remove the Wawer fatal death completely.
## See email to Steve on 3 Jan 2012 for details

## 2 Jan 2012: To model infection related reduced life expectancy
## Version 6 of this file: vital_dynamics_d6a.R not working

mig.vital.dynamics <- function(net,
                               curr.time,
                               pop.rate.of.birth,
                               birth.type.prob,
                               time.infected.fatal,
                               std.prob.death,
                               late.prob.death,
                               late.stage.begin,
                               filename=filename){

  time.since.infection <- get.vertex.attribute(net, "time.infected") # moved up top
  set.vertex.attribute(net, "time.infected", c(time.since.infection+1),
                       v=1:network.size(net)) # should happen in vital dynamics

  newpop <- net
  n.cum <- network.size(net)
  active.v <- (1:n.cum)[is.active(net,
                                  onset=prev.time,v=1:n.cum)]
  n.active <- length(active.v)

  ## update time since infection

  ## time.since.infection <- get.vertex.attribute(net, "time.infected")
  ## set.vertex.attribute(net, "time.infected", c(time.since.infection+1),
  ##                      v=1:network.size(net)) # should happen in vital dynamics
                                        # should happen up at the top of the code,
                                        # because 'net' is called 'newpop' henceforth 

  ## deaths
  #inf.ind <- which(!is.na(time.since.infection)) # there is a problem here
  inf.ind <- which(get.vertex.attribute(net,"inf.status")==1) 
  late.ind <- which(time.since.infection > late.stage.begin)
  cat("Number of total infected individuals is",
      length(inf.ind), ".") # for monitoring
  
  inf.alive.ind <- intersect(inf.ind, active.v)
  cat("Number of alive infected  individuals is",
      length(inf.alive.ind), ".\n") # for monitoring and calculating

  late.alive.ind <- intersect(late.ind, active.v)
  
  write.table(cbind(
                curr.time, length(active.v), length(inf.alive.ind)),
              file=filename,
              append=T,
              col.names=FALSE,
              row.names=FALSE
              )
  
  prob.of.death <- rep(NA, network.size(newpop))
  prob.of.death[-active.v] <- 1         # dead stay dead
  prob.of.death[active.v] <- std.prob.death
  #prob.of.death[inf.alive.ind] <- inf.prob.death
  #prob.of.death[late.alive.ind] <- late.prob.death+std.prob.death
      
  coin.flip.death <- runif(length(prob.of.death))

  ## dying.nonaids <- which(coin.flip.death <= prob.of.death
  ##                         # earlier had written this as prob
  ##                         # std.prob.death, which is a const, not
  ##                         # a vector
  ##                         )
  ##                                       # non-AIDS deaths
  dying.nonaids <- which(coin.flip.death <= prob.of.death)
  dying.aids <- which((newpop %v% "time.infected" > time.infected.fatal))
  ## browser()
  dying <- union(dying.nonaids, dying.aids)
  
                                        # 3 Jan 2012: update foll.
                                        # discussion with Steve
  
  new.dying <- intersect(dying, active.v)

  cat("Number of new deaths is ", length(new.dying), ".", 
      "Number of total deaths is ", length(dying), ".\n"
      )
   ## newpop <- set.vertex.attribute(newpop,"depart.time.bp", #(didn't do this)
        ## ##                        v=new.dying,0)
        ## ##                                set depart.time.bp to 0 for those newly dying

   survivors <- setdiff(active.v, dying)
   newpop <- activate.vertices(newpop,
                               onset=curr.time,
                               terminus=(curr.time+1),
                               v=survivors
                               )
                                        # extend surviving inds
    active.v <- survivors

    nedge.cum <- network.edgecount(newpop)
    if (length(dying)>0) {
      edges.with.dying.node <- sort(unique(unlist(sapply(dying, function(x)get.edgeIDs.active(newpop, prev.time, v=x)))))
    } else {
      edges.with.dying.node <- NULL
    } ## what is this?

    active.edgeIDs <- (1:nedge.cum)[is.active(newpop,prev.time,e=(1:nedge.cum))]
                                        # XXXXXXXXXXXX

     surviving.edges <- setdiff(active.edgeIDs,edges.with.dying.node)
     newpop <- activate.edges(newpop, onset=curr.time, terminus=(curr.time+1),
                              e=surviving.edges
  )
   ##                                      ## extend ties between two survivors

  ## add here, once the file is debugged from version d_3

  pop.mean.births <- pop.rate.of.birth # think about this
  ## if "rate" of births goes up, i.e. more number
  ## of births per unit time, then the mean number of
  ## births in a given time step should also go up
  ## because in this sense pop.rate.of.birth is mean number of 
  ## births per minute
  
  nbirths <- rpois(1, pop.mean.births)
  cat("Number of new births is ", nbirths, ".\n")
  
  if (nbirths>0){
    birthIDs <- (n.cum+1):(n.cum+nbirths)
    newpop <- add.vertices(newpop,nbirths)
    newpop <- set.vertex.attribute(newpop, "vertex.names", v=birthIDs, birthIDs)
    newpop <- activate.vertices(newpop,onset=curr.time,terminus=curr.time+1,v=birthIDs)
    newtypes <- rmult.sg(nbirths, birth.type.prob)
    newpop <- set.vertex.attribute(newpop, "type", v=birthIDs, newtypes)

    male.rural <- which(get.vertex.attribute(newpop, "type")==1)
    male.migrant <- which(get.vertex.attribute(newpop, "type")==2)
    male.urban <- which(get.vertex.attribute(newpop, "type")==3)
    female.rural <- which(get.vertex.attribute(newpop, "type")==4)
    female.urban <- which(get.vertex.attribute(newpop, "type")==5)
    male.rural.migrant <- as.integer(length(male.migrant)/2)

    ## set "type" attributes for new people
    set.vertex.attribute(newpop, "type", "A-MR", male.rural)
    set.vertex.attribute(newpop, "type", "B-MM", male.migrant)
    set.vertex.attribute(newpop, "type", "C-MU", male.urban)
    set.vertex.attribute(newpop, "type", "D-FR", female.rural)
    set.vertex.attribute(newpop, "type", "F-FU", female.urban)

    ## set "location" attributes for new people
    set.vertex.attribute(newpop, "location", 0, male.rural)
    set.vertex.attribute(newpop, "location", 1, male.urban)
    set.vertex.attribute(newpop, "location", 0, female.rural)
    set.vertex.attribute(newpop, "location", 1, female.urban)
    if (length(male.migrant) >= 1){
    set.vertex.attribute(newpop, "location", (rbinom(length(male.migrant),1,0.5)),
                                                    male.migrant[1:length(male.migrant)]
                         )


   #set.vertex.attribute(newpop, "location", 1, male.migrant[-c(1:male.rural.migrant)])
  }

    ## attributes to account for new k-star terms
    ## set.vertex.attribute(newpop, "migmalemix.rural",
    ##                      "B-MMR.MIX", male.migrant)
    ## set.vertex.attribute(newpop, "migmalemix.rural",
    ##                      "B-MMR.MIX", female.rural)

    ## set.vertex.attribute(newpop, "migmalemix.urban",
    ##                      "2MMU.MIX", male.migrant)
    ## ## set.vertex.attribute(newpop, "migmalemix.urban",
    ## ##                      "2MMU.MIX", female.rural)
    ## set.vertex.attribute(newpop, "migmalemix.urban",
    ##                      "2MMU.MIX", female.urban)

    ## set other attributes for new people
    set.vertex.attribute(newpop, "time.infected", NA, birthIDs)
    set.vertex.attribute(newpop, "inf.time", NA, birthIDs)
    set.vertex.attribute(newpop, "inf.status", 0, birthIDs) # they should be susceptible
                                        # not NA
    set.vertex.attribute(newpop, "infector.ID", NA, birthIDs)
    set.vertex.attribute(newpop, "infector.location", NA, birthIDs)
    
    active.v <- c(active.v, birthIDs) ## will only work after deaths are built in
  }

cat("Number alive is", length(active.v), ".\n")

return (newpop);
  
}


