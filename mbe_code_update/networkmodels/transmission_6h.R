## writing transmision function to
## start with and finally output full-network
## memoryless migrations

## 31 Jan 2012: Incidence data

## 27 Jan 2012: Rewrite migration bit to not be in C.

## 22 Jan 2012: introduce conditional
## for mig-step 0

## 13 jan 2012: try again
## fixed death-rate

## 9 Jan 2012: C-speed-up of migration-step is working
## now will vectorize disease transmission-step to speed it up
## seems to be working, much faster
## 6d is working but possible error in infection transmission
## check here

## transmission_d6a.R contains the 
## unvectorized transmission function that works

## 7 Jan 2012: To coordinate with C speed-up functions.


## 2 Jan 2012: This data is directly written to a csv file.
## Possibly I am using information on "all infected indviduals"
## instead of "alive infected individuals" in computing prevalence
## rates. need to correct that.

## 13 November 2011
## set up to write to a file

## 18 November 2011 record infector.type

transmission <- function(net,
                         mig_step,
                         acute.prob,
                         chronic.prob,
                         late.prob,
                         acute.chronic.transition,
                         chronic.late.transition,
                         curr.time,
                         filename=filename
                         ){
  ## pop <- net

  net <- net
  ##pop <-  network.collapse(net, at=curr.time)
  pop <- net
  
  migrant_males <- which(get.vertex.attribute(pop, "type")==
                        "B-MM")
  num_migrant_males <- length(migrant_males)
  current_location <- get.vertex.attribute(pop, "location")
  current_location <- current_location[migrant_males] # restrict to only migrant-males
  
  current_time  <- curr.time # for compatability with C

  #cat("Migration-Probability is", prob_migrating, "\n");
  #prob_migrating <- rep(prob_migrating, length(migrant_males))
  ##cat("Initial locations of migrants are", current_location,"\n")
  
  #if ((as.integer(curr.time/mig_step) %% 2) == 1){

  if(mig_step > 0){
  ## new_location <-
  ##     migration_algorithm(#coin_flips=coin_flips,
  ##                         #prob_migrating=prob_migrating,
  ##                         current_time=current_time,
  ##                         current_location=current_location,
  ##                         migrant_males=migrant_males,
  ##                         num_migrant_males=num_migrant_males,
  ##                         mig_step=mig_step
  ##                         )
    new_location <- current_location
    prob_mig <- 1/mig_step
    coin_flips <- runif(num_migrant_males, 0, 1)
    move <- which(coin_flips <= prob_mig)
    new_location[move] <- 1-new_location[move]

    ## cat("Moving Individuals are",
    ##     (network.vertex.names(pop))[migrant_males[move]], "\n")

    ## cat("And their types are",
    ##      (get.vertex.attribute(pop, "type"))[migrant_males[move]], "\n")

    ##browser()
    
    set.vertex.attribute(pop,"location", new_location, migrant_males)
    set.vertex.attribute(net,"location", new_location,
                         ((network.vertex.names(pop))[migrant_males])
                         )
    ## cat("and the locations of all migrant-males are",
    ##     get.vertex.attribute(net,"location")[((network.vertex.names(pop))[migrant_males])], "\n")
                                             
} # conditional for mig_step 0
  
  ## transmit disease
  c <- 3 # from ODE model
  acute.prob.week <- 1-(1-acute.prob)^c
  chronic.prob.week <- 1-(1-chronic.prob)^c
  late.prob.week <- 1-(1-late.prob)^c

  ## p1 infected, p2 uninfected
  pop <- net
  
  el <- as.edgelist(pop)
  p1.inf <- intersect(which((get.vertex.attribute(pop, "inf.status")
                             [el[,1]]==1)),
                      which((get.vertex.attribute(pop, "inf.status")
                             [el[,2]]==0))
                      )
                      
  ##browser()
  
  cat("Infected p1's with uninfected p2's are ",
      (network.vertex.names(pop))[el[p1.inf,1]], "\n")
  same.location.p1 <- which((pop%v%"location")[el[,1]]==
                            (pop%v%"location")[el[,2]]
                            )
  ##p1 infected, p2 unifected same location
  p1.inf.samel.p2 <- intersect(p1.inf, same.location.p1)
  ##browser()
  cat("Infected p1's in same location with uninfected p2's are ",
      (network.vertex.names(pop))
      [el[p1.inf.samel.p2,1]], "\n")
  cat("and those p2's are ",
      (network.vertex.names(pop))
      [el[p1.inf.samel.p2,2]], "\n")
  ##browser()


  total.new.infections <- 0  
  if(length(p1.inf.samel.p2)>0){
    ## define p1's in same locations with p2's
    
    time.since.infection <- get.vertex.attribute(pop, "time.infected")

    acute.p1 <- which(time.since.infection[el[p1.inf.samel.p2,1]] <=
                      acute.chronic.transition)

    late.p1 <- which(time.since.infection[el[p1.inf.samel.p2,1]] >=
                      chronic.late.transition)

    chronic.p1 <- intersect( 
                            which(time.since.infection[el[p1.inf.samel.p2,
                                                          1]] >
                                  acute.chronic.transition),
                            which(time.since.infection[el[p1.inf.samel.p2,
                                                          1]] <
                                  chronic.late.transition)
                            )

    newinf.acute.p2 <- which(rbinom(length(acute.p1),1,acute.prob.week)==1)
    newinf.late.p2 <- which(rbinom(length(late.p1), 1,late.prob.week)==1)
    newinf.chronic.p2 <- which(rbinom(length(chronic.p1), 1,chronic.prob.week)==1)

    ##browser()
    ## update inf.status and time.infected in pop
    set.vertex.attribute(pop, "inf.status",
                           1,
                           el[acute.p1[newinf.acute.p2],2]
                           )

    set.vertex.attribute(pop, "inf.status",
                           1,
                           el[late.p1[newinf.late.p2],2]
                           )

    set.vertex.attribute(pop, "inf.status",
                           1, #5 March 2012: found this error
                           el[chronic.p1[newinf.chronic.p2],2]
                           )

    set.vertex.attribute(pop, "time.infected",
                           0,
                           el[acute.p1[newinf.acute.p2],2]
                           )

    set.vertex.attribute(pop, "time.infected",
                           0,
                           el[late.p1[newinf.late.p2],2]
                           )

    set.vertex.attribute(pop, "time.infected",
                           0,
                           el[chronic.p1[newinf.chronic.p2],2]
                           )

    ## update in net -- acute
    set.vertex.attribute(net, "inf.status", 1,
                        ((network.vertex.names(pop) # check this -- should be pop
                          )[el[p1.inf.samel.p2[newinf.acute.p2],2]]
                        ))

    set.vertex.attribute(net, "time.infected", 0,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p1.inf.samel.p2[newinf.acute.p2],2]]
                        ))

    set.vertex.attribute(net, "inf.time", curr.time,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p1.inf.samel.p2[newinf.acute.p2],2]]
                        ))

    ## update in net -- late
    set.vertex.attribute(net, "inf.status", 1,
                        ((network.vertex.names(pop) # check this -- should be pop
                          )[el[p1.inf.samel.p2[newinf.late.p2],2]]
                        ))

    set.vertex.attribute(net, "time.infected", 0,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p1.inf.samel.p2[newinf.late.p2],2]]
                        ))

    set.vertex.attribute(net, "inf.time", curr.time,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p1.inf.samel.p2[newinf.late.p2],2]]
                        ))

    ## update in net -- chronic
    set.vertex.attribute(net, "inf.status", 1,
                        ((network.vertex.names(pop) # check this -- should be pop
                          )[el[p1.inf.samel.p2[newinf.chronic.p2],2]]
                        ))

    set.vertex.attribute(net, "time.infected", 0,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p1.inf.samel.p2[newinf.chronic.p2],2]]
                        ))

    set.vertex.attribute(net, "inf.time", curr.time,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p1.inf.samel.p2[newinf.chronic.p2],2]]
                        ))
    
    newinfp2 <- union(el[p1.inf.samel.p2[newinf.acute.p2],2],
                      union(el[p1.inf.samel.p2[newinf.chronic.p2],2],
                            el[p1.inf.samel.p2[newinf.late.p2],2])
                      ) 

    length.newinfp2 <- length(newinfp2)
    total.new.infections <- length.newinfp2

      newinfp2.ptns <- union(el[p1.inf.samel.p2[newinf.acute.p2],1],
                         union(el[p1.inf.samel.p2[newinf.chronic.p2],1],
                               el[p1.inf.samel.p2[newinf.late.p2],1])
                      )

  cat("Newly infected p2's are ",
      network.vertex.names(pop)[newinfp2],
      "\n",
      "partnered with ",
      network.vertex.names(pop)[newinfp2.ptns],
      ".\n",
      "NUMBER OF NEWLY INFECTED P2:",
      length(c(newinf.acute.p2, newinf.late.p2, newinf.chronic.p2)),
      "\n"
      )

  }

  
  ##assign("length.newinfp2", length(newinfp2), envir=.GlobalEnv)

  ##browser();

  ## may need to redefine 'pop' here since attributes in
  ## net have been updated
  ##pop <- network.collapse(net, at=curr.time)
  
  ## where p2 is the infected one
  p2.inf <- intersect(which((get.vertex.attribute(pop, "inf.status")
                             [el[,1]]==0)),
                      which((get.vertex.attribute(pop, "inf.status")
                             [el[,2]]==1))
                      )
                      

  
  cat("Infected p2's with uninfected p1's are ",
      (network.vertex.names(pop))[el[p2.inf,2]],
      "\n")

  same.location.p1 <- which((pop%v%"location")[el[,1]]==
                            (pop%v%"location")[el[,2]]
                            )

  p2.inf.samel.p1 <- intersect(p2.inf, same.location.p1)

  cat("Infected p2's in same location with uninfected p1's are ",
      (network.vertex.names(pop))[el[p2.inf.samel.p1,2]],
      "\n")
    
  cat("and those p1's are ",
      (network.vertex.names(pop))
      [el[p2.inf.samel.p1,1]], "\n"
      )

  
  ##p2 infected, p1 unifected, same location

  
  if(length(p2.inf.samel.p1)>0){
    ## define p1's in same locations with p2's
    
    time.since.infection <- get.vertex.attribute(pop, "time.infected")

    acute.p2 <- which(time.since.infection[el[p2.inf.samel.p1,2]] <=
                      acute.chronic.transition)

    late.p2 <- which(time.since.infection[el[p2.inf.samel.p1,2]] >=
                      chronic.late.transition)

    chronic.p2 <- intersect( 
                            which(time.since.infection[el[p2.inf.samel.p1,
                                                          2]] >
                                  acute.chronic.transition),
                            which(time.since.infection[el[p2.inf.samel.p1,
                                                          2]] <
                                  chronic.late.transition)
                            )

    newinf.acute.p1 <- which(rbinom(length(acute.p2),1,acute.prob.week)==1)
    newinf.late.p1 <- which(rbinom(length(late.p2), 1,late.prob.week)==1)
    newinf.chronic.p1 <- which(rbinom(length(chronic.p2), 1,chronic.prob.week)==1)

    ## update inf.status and time.infected in pop
    set.vertex.attribute(pop, "inf.status",
                           1,
                           el[acute.p2[newinf.acute.p1],1]
                           )

    set.vertex.attribute(pop, "inf.status",
                           1,
                           el[late.p2[newinf.late.p1],1]
                           )

    set.vertex.attribute(pop, "inf.status",
                           1, # 5 March 2012: found this error
                           el[chronic.p2[newinf.chronic.p1],1]
                           )

    ## update in net -- acute
    set.vertex.attribute(net, "inf.status", 1,
                        ((network.vertex.names(pop) # check this -- should be pop
                          )[el[p2.inf.samel.p1[newinf.acute.p1],1]]
                        ))

    set.vertex.attribute(net, "time.infected", 0,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p2.inf.samel.p1[newinf.acute.p1],1]]
                        ))

    set.vertex.attribute(net, "inf.time", curr.time,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p2.inf.samel.p1[newinf.acute.p1],1]]
                        ))

    ## update in net -- late
    set.vertex.attribute(net, "inf.status", 1,
                        ((network.vertex.names(pop) # check this -- should be pop
                          )[el[p2.inf.samel.p1[newinf.late.p1],1]]
                        ))

    set.vertex.attribute(net, "time.infected", 0,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p2.inf.samel.p1[newinf.late.p1],2]]
                        ))

    set.vertex.attribute(net, "inf.time", curr.time,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p2.inf.samel.p1[newinf.late.p1],1]]
                        ))

    ## update in net -- chronic
    set.vertex.attribute(net, "inf.status", 1,
                        ((network.vertex.names(pop) # check this -- should be pop
                          )[el[p2.inf.samel.p1[newinf.chronic.p1],1]]
                        ))

    set.vertex.attribute(net, "time.infected", 0,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p2.inf.samel.p1[newinf.chronic.p1],1]]
                        ))

    set.vertex.attribute(net, "inf.time", curr.time,
                       ((network.vertex.names(pop) # check this -- should be pop
                         )[el[p2.inf.samel.p1[newinf.chronic.p1],1]]
                        ))
  ##} # parallel to newinfp2

    newinfp1 <- union(el[p2.inf.samel.p1[newinf.acute.p1],1],
                    union(el[p2.inf.samel.p1[newinf.chronic.p1],1],
                          el[p2.inf.samel.p1[newinf.late.p1],1])
                    )

    length.newinfp1 <<- length(newinfp1)
    
    ##assign("length.newinfp1", length(newinfp1), envir=.GlobalEnv)
    
    newinfp1.ptns <- union(el[p2.inf.samel.p1[newinf.acute.p1],2],
                         union(el[p2.inf.samel.p1[newinf.chronic.p1],2],
                               el[p2.inf.samel.p1[newinf.late.p1],2])
                      )

    cat("Newly infected p1's are ",
        network.vertex.names(pop)[newinfp1],
        "\n",
        "partnered with ",
        network.vertex.names(pop)[newinfp1.ptns],
        ".\n",
        "NUMBER OF NEWLY INFECTED P1:",
        length(c(newinf.acute.p1, newinf.late.p1, newinf.chronic.p1)),
        "\n"
        )

      ##browser();



  total.new.infections <- total.new.infections + length(newinfp1)
     
  }

  incidence_data <- paste("incidence_data", filename, 
                          sep="")

  
  write.table(cbind(curr.time, total.new.infections),
              file=incidence_data,
              append=TRUE,
              col.names=FALSE,
              row.names=FALSE
              )

 
  return(net) #return full network
} #close function
