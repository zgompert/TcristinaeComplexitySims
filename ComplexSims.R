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
## modHOST = boolean, include host and possible NFDS
## modOD = boolean, include overdominance
## alpha = base host effect on log scale, log(ws/wg) = alpha[host] + sfreq * beta[host]
## beta = NFDS host effect on log scale, log(ws/wg) = alpha[host] + sfreq * beta[host]
## mumg = mean log(wm/(wg | ws))
## sdmg = sd of log(wm/(wg | ws)) about mean, set to non-zero for fluctuating selection on color
## od = relative fitness (not log scale) of melanic hets
## note that all relative fitness values are multiplier
eesim<-function(ngen=100,nburn=100,ndeme=30,dprob=0.01,dsd=0.1,
		acfunc=1,mnN=25,modHOST=TRUE,modOD=FALSE,
		alpha=c(0,0),beta=c(0,0),mumg=0,sdmg=0,od=1){

	## collect model parameters
	smod<-list(modHOST=modHOST,modOD=modOD,alpha=alpha,beta=beta,
		   mumg=mumg,sdmg=sdmg,od=od)

	## initialize the network
	network<-init(ndeme,acfunc,mnN)

	## loop
	Nt = ngen+nburn ## total number of iterations
	
	## output
	## columns are gen, deme, N, q, p
	out<-matrix(NA,nrow=ngen*ndeme,ncol=5)

	for(i in 1:Nt){
		## one generation of evolution
		network<-ongen(nw=network,ndeme=ndeme,dprob=dprob,dsd=dsd,mm=smod) 		
		if(i<=nburn){
			if((i %% 20) == 0){cat("Burnin gen ",i,"\n")}
		} else {
			if(((i-nburn) %% 20) == 0){cat("Post-burnin gen ",i-nburn,"\n")}
			oneres<-as.matrix(cbind(rep(i,ndeme),1:ndeme,network$N,network$q,network$p))
			out[(((i-nburn-1) * ndeme) + 1):(((i-nburn-1) * ndeme) + ndeme),]<-oneres
		}
	}
	return(out)
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

## assign hosts, A = 1, C = 2
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
## mm = collection of model parameters
ongen<-function(nw=NA,ndeme=NA,dprob=NA,dsd=NA,mm=NA){
	
	## create genotypes, this sampling from allele frequencies
	## is where drift occurs
	G<-vector("list",ndeme)
	for(k in 1:ndeme){
		if(nw$N[k] >= 1){
			ppp<-c(nw$q[k],nw$p[k],(1-nw$p[k]-nw$q[k]))
			ppp[ppp<0]<-0 ## avoids flot imprecision
			ppp<-ppp/sum(ppp) 
			## order is melanic, stripe, green
			G[[k]]<-as.matrix(rmultinom(n=nw$N[k],size=2,prob=ppp))
		}
		else{ ## extripated
			G[[k]]<-NULL
		}
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

	## generate phenotypes assumes dominance: g > s > m, from Aaron's work
	P<-vector("list",ndeme)
	for(k in 1:ndeme){
		if(nw$N[k] >=1){
			G[[k]]<-as.matrix(G[[k]])
			P[[k]]<-rep(NA,nw$N[k])
			## melanics, assign value 1
			P[[k]][G[[k]][1,]==2]<-1
			## stripes assign value 2
			P[[k]][(G[[k]][2,]==2) | ((G[[k]][1,]+G[[k]][2,])==2)]<-2
			## green assign value 3
			P[[k]][G[[k]][3,]>=1]<-3
		}
	}


	## compute fitness
	wbar = 1 ## initial baseline mean fitness, to be modified
	w<-G ## relative fitness

	for(k in 1:ndeme){
		w[[k]]<-rep(wbar,nw$N[k]) ## start with 1
		## fluctuating selection
		if(mm$sdmg > 0){ ## this means fluctuating
			mgk<-rnorm(1,mm$mumg,mm$sdmg) ## mgk is log(wm/wg); 0 = wm = wg (wg = ws)
		} else{ ## not fluctuating
			mgk<-mm$mumg
		}
		## melanics
		w[[k]][P[[k]]==3]<- w[[k]][P[[k]]==3] * exp(mgk) ## multiply wbar * wm/wg, exp to reverse log
		## don't need to multiply green or stripe
		
		## host and nfds
		if(mm$modHOST==TRUE){
			sfreq<-sum(P[[k]]==2)/sum(P[[k]]==2 | P[[k]]==3) ## stripe freq relative to green + stripe	
			## index 1 = A, index 2 = C
			sgk<-mm$alpha[nw$host[k]] + sfreq * mm$beta[nw$host[k]] ## log(ws/wg) = alpha_host + sfreq * beta_host
			w[[k]][P[[k]]==2]<-w[[k]][P[[k]]==2]*log(sgk) ## multiply base ws * ws/wg
		}

		## overdominance
		if(mm$modOD==TRUE){ 
			## melanic hets
			w[[k]][G[[k]][1,]==1]<-w[[k]][G[[k]][1,]==1] * mm$od ## relative fitness of overdominant
		}

		## arthropod density
	}
	
	## apply selection and update allele frequencies and N
	for(k in 1:ndeme){
		## take weighted average calculation
		if(nw$N[k] >= 1){ ## population size is non-zero
			asums<-G[[k]] %*% w[[k]]
			fr<-asums/sum(asums)
			nw$q[k]<-fr[1,1]
			nw$p[k]<-fr[2,1]
		} else{
			nw$q[k]<-NA
			nw$p[k]<-NA
		}
	}
	return(nw)
}
