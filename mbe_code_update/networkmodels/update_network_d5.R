update.network <- function(net,
                           curr.time,
                           theta.network,
                           burnin.sim.network,
                           diss.network,
                           gamma.network,
                           formula.network.w.offset,
                           constraints.network,
                           activenet=activenet,
                           dyninterval){

#  curr.time <- timestep

  ##activepop <<- activepop <- network.copy(activenet) # update this
  activepop <<- activepop <- network.copy(activenet) # update this

  n.active <- network.size(activepop) # update this

  theta.network[length(theta.network)] <- -log(n.active)
  formula.network.with.activepop <- update.formula(
                                              formula.network.w.offset,activepop~.)		
  
  changes <- simulate.formula.stergm(formula.network.with.activepop,
                                   theta.form=theta.network,
                                        # should be theta.network or theta.st
                                   dissolution=diss.network,
                                   theta.diss=gamma.network,
                                   MH.burnin=dyninterval,
                                   constraints=constraints.network,
                                   nsim=1)


		

  active.changelist <- matrix(changes$changed[,2:3], ncol=2) # changelist

  new.pop <- network.copy(net) # cumulative pop
  
  edgevector <- 1:length(new.pop$mel)

  if(dim(active.changelist)[1] > 0){
    tiepresence <- sapply(1:dim(active.changelist)[1], function(x)
                        # decide which additions, deletions
                          activepop[active.changelist[x,1],active.changelist[x,2]]
                                        # whats this
                
                          )

    active.addlist <- active.changelist[tiepresence==0,]
                                        # Create edgelist of starting edges
    active.dellist <- active.changelist[tiepresence==1,]
                                        # Create edgelist of ending edges

    full.addlist <- matrix((activepop %v% "vertex.names")[active.addlist],ncol=2)
                                        #  Convert back to full network nodeIDs
    full.dellist <- matrix((activepop %v% "vertex.names")[active.dellist],ncol=2)
                                        #  Convert back to full network nodeIDs

    if (dim(full.dellist)[1]>0) {								# ENDING EDGE DE-ACTIVATION
      del.edgeIDs <- get.undir.dyad.edge.IDs.from.edgelist(
                                        # Get all edgeIDs for the
                                        #terminating dyads (NB: dyads 
                                                           new.pop,full.dellist)
                                        # can have more than one edgeID,
                                        # if they have had 
                                        # multiple relationships through time
      active.del.edgeIDs <- del.edgeIDs[is.active(# Determine which are active
                                                  new.pop,onset=curr.time,
                                                  e=del.edgeIDs)]					

      new.pop <- deactivate.edges(new.pop,onset=curr.time,terminus=curr.time+1,
                                  e=active.del.edgeIDs)							# De-active them
    }



    ##browser()
  if (dim(full.addlist)[1]>0) {
                                        # STARTING EDGE ACTIVATION
    nedge.cum <- network.edgecount(new.pop)
                                        # Figure out cum number of edges
                                        # to determine new edge IDs
    new.pop <- add.edges(new.pop,
                         tail=full.addlist[,1],head=full.addlist[,2])				# Add the new edges to the full pop		
    new.pop <- activate.edges(new.pop,onset=curr.time,terminus=curr.time+1,
                              e=(nedge.cum+1):(nedge.cum+dim(full.addlist)[1]))			# Activate them
  }
  }
  #cat("Number of active edges is ", nrow(as.edgelist(new.pop)), ".\n")
  return(new.pop)
}
