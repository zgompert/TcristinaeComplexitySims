## functions for simulating eco-evolutionary dynamics of color and pattern (stripe)
## in Timema cristinae
## various processes are included, which can be turned on or off



## main function
## ngen = number of generations post burnin
## nburn = number of burnin generations
## ndeme = number of demes
## dprob = probability of dispersing
## dsd = sd of normal kernel for dispersal
## acfunc = method for autocorrelation in host, 1 = random with equal freqs
## mnN = mean of Poisson carrying capacity per deme
eesim<-function(ngen=100,nburn=100,ndeme=30,dprob=0.01,dsd=0.1,
		acfunc=1,mnN=25){

	## initialize the network
	network<-init(ndeme,afunc,mnN)

	## loop
	Nt = ngen+nburn ## total number of iterations
	for(i in 1:Nt){
		## one generation of evolution
		network<-ongen(nw=network,ndeme=ndeme,dprob=dprob,dsd=dsd) 		
		if(i<=nburn){
			if((i %% 20) == 0){cat("Burnin gen ",i,"\n")}
		} else {
			if(((i-nburn) %% 20) == 0){cat("Post-burnin gen ",i-nburn,"\n")}
			## need to decide what to save or print
		}
	}
}

## initialization function
## network is always on a x-y grid with locations from a bivariate standard normal
## most initialization is "random", not tied to some level of adaptation
## ndeme = number of demes
## acfunc = method for autocorrelation in host, 1 = random with equal freqs
## mnN = mean of Poisson carrying capacity per deme
init<-function(ndeme=30,acfunc=1,mnN=25){
	## sample x and y coordinates independently
	xcoords<-rnorm(n=ndeme,mean=0,sd=1)
	ycoords<-rnorm(n=ndeme,mean=0,sd=1)

	## determine host for each deme
	hosts<-samhosts(x=xcoords,y=ycoords,ac=acfunc)

	## set max population size of each deme
	maxN<-rpois(n=ndeme,lambda=mnN)

	## set initial deme N
	N<-round(runif(ndeme,min=0,max=maxN),0)

	## set inital morph frequencies
	q<-rbeta(n=ndeme,shape1=15,shape2=85) ## melanic, expectation 15%
	## stripe, expectation 50% of non-melanic
	p<-(1-q) * rbeta(n=ndeme,shape1=20,shape2=20) 
	## green is thus 1 - p - q
	## we assume random mating, so genotypes do not need init here
	## they are created on the fly

	## initialize arthropod numbers
	arth<-rpois(ndeme,lambda=100) ## set arth number at 100 mean for now

	nw<-list(xcoords=xcoords,ycoords=ycoords,hosts=hosts,maxN=maxN,
		 N=N,p=p,q=q,arth=arth)
	return(nw)
}

## assign hosts
samhosts<-function(x=NA,y=NA,ac=1){
	hosts<-NA
	if(ac==1){ ## random and equal
		hosts<-sample(1:2,length(x),replace=TRUE)
	} else{
		warning("Undefined host autocorrelation function")
	}
	return(hosts)
}

## simulate one generation of eco-evo
## nw = network object, see init for components
## ndeme = number of demes
## dprob = probability of dispersing
## dsd = sd of normal kernel for dispersal
ongen<-function(nw=NA,ndeme=NA,dprob=NA,dsd=NA){
	
	## create genotypes, this sampling from allele frequencies
	## is where drift occurs
	G<-vector("list",ndeme)
	for(k in 1:ndeme){
		ppp<-c(nw$q[k],nw$p[k],(1-nw$p[k]-nw$q[k]))
		ppp[ppp<0]<-0 ## avoids flot imprecision
		## order is melanic, stripe, green
		G[[k]]<-rmultinom(n=nw$N[k],size=2,prob=ppp)
	}

	## dispersal
	for(k in 1:ndeme){ ## no double dispersal
		## sample dispersers, 1 = yes
		dd<-which(sample(0:1,nw$N[k],replace=TRUE,prob=c(1-dprob,dprob))==1)
		if(length(dd)>=1){ ## at least one disperser
			for(j in 1:length(dd)){ ## move the dispersers
				## need to precompute dispersal probs!!! ##
				dest<-sample(1:ndeme,1) ## equal dispersal, temp
				#dest<-sample(1:ndeme,1,prob=dkp[k,-k]) ## matrix of disperal probs
				G[[dest]]<-as.matrix(cbind(G[[dest]],G[[k]][,j]))
				nw$N[dest]<-nw$N[dest]+1
			}
			G[[k]]<-G[[k]][,-dd] ## remove from source
			nw$N[k]<-nw$N[[k]]-length(dd)
		}
	}

	## generate phenotypes assumes dominance: s > g > m
	P<-vector("list",ndeme)
	for(k in 1:ndeme){
		P[[k]]<-rep(NA,nw$N[k])
		## melanics, assign value 1
		P[[k]][G[[k]][1,]==2]<-1
		## stripes assign value 2
		P[[k]][G[[k]][2,]>=1]<-2
		## green assign value 3
		P[[k]][(G[[k]][3,]==2) | ((G[[k]][1,]+G[[k]][3,])==2)]<-3
	}


	## compute fitness
	wbar = 1 ## baseline mean fitness, to be modified
	w<-G ## relative fitness

	for(k in 1:ndeme){
		w[[k]]<-rep(wbar,nw$N[k]) ## start with 1
		## fluctuating selection
		if(sdmg > 0){ ## this means fluctuating
			mgk<-rnorm(1,mumg,sdmg) ## mgk is log( w(m/g)); 0 = w(m) = w(g)
		} else{ ## not fluctuating
			mgk<-mumg
		}
		
		## other multiplies
		## host and nfds
		## arthropod density
		## overdominance
	}
	
	## apply selection and update allele frequencies and N
	for(k in 1:ndeme){
		## take weighted average calculation
		asums<-G[[k]] %*% w[[k]]
		fr<-asums/sum(asums)
		nw$q[k]<-fr[1,1]
		nw$p[k]<-fr[2,1]
	}
	return(nw)
}
