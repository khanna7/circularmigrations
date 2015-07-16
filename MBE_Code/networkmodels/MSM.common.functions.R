
#########################################

freshen <- function() {
  rm(list=ls())
  load("puma.model.100.RData")
}

#########################################
### logit
logit <- function(p)log(p/(1-p))


#########################################
### expit
expit <- function(p) exp(p)/(1+exp(p))


### dyn.diag = Dynamic diagnostics ######
# Basic checks for whether the network.series that comes out of a call to simulatedyn
# has the right average duration for relationships and has the right structure for the
# final network
#########################################

dyn.diag <- function(object) {

  dm <- duration.matrix(object)

  init.dur <- dm[dm$Start==0,]$duration
  mean.init.dur <- mean(init.dur)
  init.censored <- sum(dm[dm$Start==0,]$Noncensored==0)
  init.n.edges <- sum(dm$Start==0)
  
  final.edgelist <-  as.matrix(dm[dm$Noncensored==0,][,1:2])
  final.net <- network.initialize(n=network.size(object$network),directed=object$network$gal$directed)
  final.net <- add.edges(final.net,final.edgelist[,1],final.edgelist[,2])
  new.form <-  update.formula(object$formula,final.net~.)
  final.net.summary <- summary(new.form)


  cat("Duration of initial edges = ",mean.init.dur," with ",init.censored," out of ",init.n.edges," censored.\n")
  cat("Structure of final network:\n")
  print(final.net.summary)
}

##########################################

get.undir.dyad.edge.IDs.from.edgelist <- function(net,edgelist) {

  if (!is.network(net)) 		stop("Net argument must be a network.")
  if (!is.matrix(edgelist)) 	stop("Edgelist argument must be an edgelist.")
  if (dim(edgelist)[2]!=2) 	stop("Edgelist argument must be an edgelist.")

  n.edges <- dim(edgelist)[1]

  result <- sapply(1:n.edges, function(x) intersect(get.edgeIDs(net,edgelist[x,1]),get.edgeIDs(net,edgelist[x,2])))
  result <- unlist(result)
  return(result)
}


##########################################

update.formula.wo.simplifying <- function (old, new, ...) {
    tmp <- .Internal(update.formula(as.formula(old), as.formula(new)))
    out <- formula(terms.formula(tmp, simplify = FALSE))
    return(out)
}

##########################################

dead.nodes.with.ties <- function(network,timestep) {

  ncs <- network.crosssection(network,timestep)
  ncsv <- ncs %v% "vertex.names"
  a <- setdiff(1:network.size(network), ncsv)
  b <- ncsv[as.vector(as.edgelist(ncs))]
  return(intersect(a,b))
}


#############################################

# is.active(msm.pop,timestep,e=get.undir.dyad.edge.IDs.from.edgelist(msm.pop,as.edgelist(msm.pop.steady.coital)))

act.v <- function(network,timestep) {
	(1:network.size(network))[is.active(network,timestep,v=1:network.size(network))]
	}

act.e <- function(network,timestep) {
	(1:length(network$mel))[is.active(network,timestep,e=1:length(network$mel))]
	}

act.el <- function(network,timestep) {
	act <- act.e(network,timestep)
	as.edgelist(network)[act,]
	}

##############################################
# Subtract undirected networks

dir.net.diff <- function(e1,e2) {
    e1 <- as.sociomatrix(e1)
    e2 <- as.sociomatrix(e2)
    network(e1 - e2,directed=F)
}

##############################################
# Three-way set intersection

intersect3 <- function(x,y,z) {
	intersect(intersect(x,y),z)
}


###############################################
# Return a vector of numbers indicating bin based on probabilities
# What you might think rmultinom does, but it doesn't!

rmult.sg <- function(n,p) {
	aaa <- rmultinom(n,1,p)
	bbb <- matrix(1:length(p),length(p),n)
	colSums(aaa*bbb)
	}

###########

get.all.edgeIDs <- function(network) {
	result <- vector()
	for (v in 1:network.size(network)) {
	  result <- c(result,get.edgeIDs(network,v=v))
	}
	result <- sort(unique(result))
	return(result)
}

###############

simulate.formula.with.known.null.stats <- function (object, nsim = 1, seed = NULL, theta0, burnin = 1000, 
#simulate.formula.with.known.null.stats <- function (object, model, casual.model.change, nsim = 1, seed = NULL, theta0, burnin = 1000, 	# SPEED
    interval = 1000, basis = NULL, statsonly = FALSE, sequential = TRUE, 
    constraints = ~., control = control.simulate.formula(), verbose = FALSE, 
    ...) 
{
    if (!is.null(seed)) {
        set.seed(as.integer(seed))
    }
    if (is.null(nw <- basis)) {
        nw <- ergm.getnetwork(object)
    }
    if (class(nw) == "network.series") {
        nw <- nw$networks[[1]]
    }
    nw <- as.network(nw)
    if (!is.network(nw)) {
        stop("A network object on the LHS of the formula or via", 
            " the 'basis' argument must be given")
    }
    formula <- ergm.update.formula(object, nw ~ .)

#    if(casual.model.change) {						# SPEED
	    m <- ergm.getmodel(formula, nw, drop = FALSE)
#    } else {
#	    source("~/AmPham/code/ergm.update.casual.model.R")
#	    m <- new.ergmCasualUpdateModel(nw, model)
#    }

    Clist <- ergm.Cprepare(nw, m)

    MHproposal <- MHproposal(constraints, arguments = control$prop.args, 
        nw = nw, model = m, weights = control$prop.weights, class = "c")
    if (missing(theta0)) {
        theta0 <- rep(0, Clist$nstats)
        warning("No parameter values given, using Bernouli network\n\t")
    }
    if (any(is.infinite(theta0) | is.nan(theta0) | is.na(theta0))) 
        stop("Illegal value of theta0 passed to simulate.formula")
#   curstats <- summary(formula)
    curstats <- rep(0, Clist$nstats)

    names(curstats) <- m$coef.names
    MCMCparams <- list(samplesize = 1, maxedges = 1 + max(20000, 
        Clist$nedges), burnin = burnin, interval = interval, 
        parallel = control$parallel, packagenames = control$packagenames, 
        Clist.miss = ergm.design(nw, m, verbose = verbose))
    if (verbose) {
        cat(paste("Starting ", nsim, " MCMC iteration", ifelse(nsim > 
            1, "s", ""), " of ", burnin + interval * (MCMCparams$samplesize - 
            1), " steps", ifelse(nsim > 1, " each", ""), ".\n", 
            sep = ""))
    }
    if (sequential && statsonly) {
        MCMCparams$samplesize <- nsim
        MCMCparams$nmatrixentries <- nsim * length(curstats)
        z <- ergm.getMCMCsample(Clist, MHproposal, theta0, MCMCparams, 
            verbose = verbose)
        colnames(z$statsmatrix) <- m$coef.names
        return(sweep(z$statsmatrix, 2, curstats, "+"))
    }
    if (!statsonly) {
        nw.list <- list()
    }
    out.mat <- matrix(nrow = nsim, ncol = length(curstats), dimnames = list(NULL, 
        m$coef.names))
    MCMCparams$nmatrixentries <- length(curstats)
    for (i in 1:nsim) {
        MCMCparams$burnin <- ifelse(i == 1 || !sequential, burnin, 
            interval)
        z <- ergm.getMCMCsample(Clist, MHproposal, theta0, MCMCparams, 
            verbose)
        if (!statsonly) {
            nw.list[[i]] <- network.update(nw, z$newedgelist, 
                matrix.type = "edgelist")
        }
        out.mat[i, ] <- curstats + z$statsmatrix
        if (sequential) {
            nw <- nw.list[[i]]
            Clist <- ergm.Cprepare(nw, m)
            curstats <- curstats + z$statsmatrix
        }
    }
    if (statsonly) 
        return(out.mat[1:nsim, ])
    if (nsim == 1) {
        return(nw.list[[1]])
    }
    else {
        out.list <- list(formula = object, networks = nw.list, 
            stats = out.mat, coef = theta0)
        class(out.list) <- "network.series"
        return(out.list)
    }
}

#############################################
steves.debugger <-function() {}


