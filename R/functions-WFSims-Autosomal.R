###############
# DEPENDENCIES
###############


####################
# Functions for 
####################

#' Sample approximate stationary distribution for deleterious
#' mutations at n autosomal loci  
#'
#' @title Rejection Sampler for autosomal loci
#' @param n  Number of loci to sample.
#' @param Ne Effective population size.
#' @param u  Mutation rate (default value of u = 1e-6).
#' @param h  Dominance of deleterious mutations (default value of h = 0).
#' @param s  Selection against deleterious mutations.
#' @export
rejectionSampler  <-  function(n=100, Ne=100, u=1e-6, h=0, s=0.01) {
	
	# Empty vector for allele frequencies
	qi  <-  rep(0,times=n)
	
	# Rejection sampler loop
	for(i in 1:n) {
		accept  <-  FALSE
		while(accept == FALSE) {
			x  <-  rbeta(1, shape1=4*Ne*u, shape2=4*Ne*u)
			U  <-  runif(1)
			if(U < exp(2*Ne*(-2*(1-x)*x*h*s - (x^2)*s))) {
				qi[i]   <-  rbinom(1,Ne,x)/Ne
				accept  <-  TRUE
			}
		}
	}
	qi
}


#' Sample approximate stationary distribution for deleterious
#' mutations at n X-linked loci  
#'
#' @title Rejection Sampler for autosomal loci
#' @param n  Number of loci to sample.
#' @param Ne Effective population size.
#' @param u  Mutation rate (default value of u = 1e-6).
#' @param h  Dominance of deleterious mutations (default value of h = 0).
#' @param sf Selection against deleterious mutations in females.
#' @param sm Selection against deleterious mutations in males.
#' @export
rejectionSamplerX  <-  function(n=100, Ne=100, u=1e-6, h=0, sf=0.01, sm=0.01) {
	
	# Empty vector for allele frequencies
	qi  <-  rep(0,times=n)
	
	# Rejection sampler loop
	for(i in 1:n) {
		accept  <-  FALSE
		while(accept == FALSE) {
			x  <-  rbeta(1, shape1=4*Ne*u, shape2=4*Ne*u)
			U  <-  runif(1)
			if(U < exp((4/3)*Ne*x*(2*h*sf*(x-1)-x*sf-sm))) {
				qi[i]   <-  rbinom(1,Ne,x)/Ne
				accept  <-  TRUE
			}
		}
	}
	qi
}


#' Linkage Disequilibrium function for W-F recursions
#'
#' @title Linkage Disequilibrium (Dstar)
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @export
Dstar  <-  function(Fii=Fii, m=m, ...) {
	(((Fii[4] + Fii[16]) - (Fii[8] + Fii[12])) / 2)*(1 - m)
}

		# genotypes ordered:
		#	c(ABAB, ABAb, ABaB, ABab, ABba*,		c(x1x1, x1x2, x1x3, x1x4, x1x5, 
		#	  AbAB, AbAb, AbaB, Abab, Abba*,		  x2x1, x2x2, x2x3, x2x4, x2x5, 
		#	  aBAB, aBAb, aBaB, aBab, aBba*,		  x3x1, x3x2, x3x3, x3x4, x3x5, 
		#	  abAB, abAb, abaB, abab, abba*,		  x4x1, x4x2, x4x3, x4x4, x4x5, 
		#	  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)


#' Haplotype frequencies among gametes
#'
#' @title Haplotype frequencies among gametes
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @param r   Recombination rate
#' @export
x.1  <-  function(Fii=Fii, m=m, r=r) {
	((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
x.2  <-  function(Fii=Fii, m=m, r=r) {
	((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.3  <-  function(Fii=Fii, m=m, r=r) {
	((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.4  <-  function(Fii=Fii, m=m, r=r) {
	((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
x.5  <-  function(Fii=Fii, m=m, r=r) {
	((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}



#' Offspring frequencies after random mating
#'
#' @title Offspring frequencies after random mating
#' @param xi Vector of haplotype frequencies among gametes (of length = 5)
#' @export
offFreq  <-  function(xi) {
	O  <-  c(xi[1]^2, (xi[1]*xi[2]), (xi[1]*xi[3]), (xi[1]*xi[4]), (xi[1]*xi[5]),
			 (xi[2]*xi[1]), xi[2]^2, (xi[2]*xi[3]), (xi[2]*xi[4]), (xi[2]*xi[5]),
			 (xi[3]*xi[1]), (xi[3]*xi[2]), xi[3]^2, (xi[3]*xi[4]), (xi[3]*xi[5]),
			 (xi[4]*xi[1]), (xi[4]*xi[2]), (xi[4]*xi[3]), xi[4]^2, (xi[4]*xi[5]),
			 (xi[5]*xi[1]), (xi[5]*xi[2]), (xi[5]*xi[3]), (xi[5]*xi[4]), xi[5]^2
			 )
	O
}


#' Calculate deterministic equilibrium genotype frequencies 
#' in the absence of the inversion 
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1x1, x1x2, x1x3, x1x4, x1x5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
#'
#' @title Find the deterministic equilibeium genotype frequencies prior to introducing the inversion
#' @param W      Vector of fitness expressions for all 25 genotypes
#' @param m      Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param threshold The threshold change in genotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. 
#' @export
#' @seealso `offFreq`, `autoInvWrightFisherSim`
#' @author Colin Olito
findEqFreqs  <-  function(W, m, r, threshold = 1e-6) {
	
	Fii.init  <-  c(1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					   0,    0,    0,    0, 0)
	Fii    <-  Fii.init
	E.Fii  <-  Fii.init

	# Storage for gamete frequencies
	xi         <-  rep(0, times=5)
	names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')

	delta  <-  rep(1, times=16)
	while(any(delta > threshold)) {
		for (j in 1:length(xi)) {
			recFct  <-  get(names(xi)[j])
			xi[j]   <-  round(recFct(Fii = E.Fii, m = m, r = r), digits=3)
		}
	
	    # offspring genotype frequencies
	    O  <-  offFreq(xi)

		# mean fitness 
		Wbar      <-  sum(O*W)

	    # difference in expected frequency (has simulation reached equilibrium yet?)
		delta   <-  E.Fii[c(1:4,6:9,11:14,16:19)] - (O*W/Wbar)[c(1:4,6:9,11:14,16:19)]
		E.Fii   <-  O*W/Wbar
	}
	names(E.Fii)  <-  NULL
	E.Fii 
}



#' Run a single Wright-Fisher Forward simulation with introduced autosomal inversion
#' using multinomial sampling with linkage to deleterious mutations
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1x1, x1x2, x1x3, x1x4, x1x5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
#'
#' @title Find the deterministic equilibeium genotype frequencies prior to introducing the inversion
#' @param Fii.init  Vector of initial frequencies (deterministic eq. frequencies in absence of inversion)
#' @param N         Population size
#' @param W         Vector of fitness expressions for all 25 genotypes
#' @param m         Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r         Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? 
#' @export
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
autoInvFwdSim  <-  function(Fii.init = Fii.init, N = N, W = W, m = m, r = r, 
							saveTrajectories = FALSE, ...) {

	# Use deterministic eq. initial frequencies
	Fii  <-  Fii.init

	# Define threshold frequency for establishment of inversion
	pcrit  <-  2/(N*m)

	# Storage for gamete frequencies  
	xi         <-  rep(0, times=5)
	names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')


	if(saveTrajectories) {
		# Storage structures for individual simulation data
		InvFreq    <-  rep(0, times=(4*N+1))
		E.InvFreq  <-  rep(0, times=(4*N+1))
		W.mean     <-  rep(0, times=(4*N+1))
	
		# Initial inversion frequency 
		InvFreq[1]    <-  sum(Fii[c(5,10,15,20:24)]/2, Fii[25])
		E.InvFreq[1]  <-  InvFreq[1]
		
		## Start forward simulation with newly introduced inversion
		gen  <-  1
		while(gen < (4*N) & InvFreq[gen] != 0) {
		
			## Step through recursions:
			# 1) Calculate gamete frequencies
			for (j in 1:length(xi)) {
				recFct  <-  get(names(xi)[j])
				xi[j]   <-  round(recFct(Fii = Fii, m = m, r = r), digits=8)
			}
			# 2) Offspring genotype frequencies
			O      <-  offFreq(xi)
			# 3) Mean fitness 
			Wbar   <-  sum(O*W)
			# 4) Expected frequencies
			E.Fii  <-  O*W/Wbar
			# 5) Draw random frequencies in adults
			Fii    <-  as.vector(rmultinom(1, N/2, E.Fii)/(N/2))
	
			# Realized frequencies
			InvFreq[gen+1]    <-  sum(Fii[c(5,10,15,20:24)]/2, Fii[25])
			E.InvFreq[gen+1]  <-  sum(E.Fii[c(5,10,15,20:24)]/2, Fii[25])
			W.mean[gen+1]     <-  Wbar
	
			gen  <-  gen+1
		}
			
		# Has the inversion reached threshold frequency for establishment (pcrit)? 
		# When did it first reach pcrit?
		if(any(InvFreq >= pcrit)) {
			invEst      <-  1
			invEstTime  <-  gen[invFreq >= pcrit][1]
		}
	
		# Save  simulation data
		res  <-  list(
				  	"InvFreq"     =  InvFreq[1:gen-1],
				  	"E.InvFreq"   =  E.InvFreq[1:gen-1],
				  	"W.mean"      =  W.mean[1:gen-1],
					"nGen"        =  gen,
				  	"InvEst"      =  invEst,
				  	"InvEstTime"  =  invEstTime
 				 	)
	} 

	if(!saveTrajectories) {

		# Storage structures for individual simulation data
		InvFreq    <-  0
		E.InvFreq  <-  0
		W.mean     <-  0
	
		# Initial inversion frequency 
		InvFreq    <-  sum(Fii[c(5,10,15,20:24)]/2, Fii[25])
		E.InvFreq  <-  InvFreq[1]
		
		## Start forward simulation with newly introduced inversion
		gen  <-  1
#		while(InvFreq > 0 & InvFreq <= pcrit) {
		while(gen < (4*N) & InvFreq > 0 ) {
		
			## Step through recursions:
			# 1) Calculate gamete frequencies
			for (j in 1:length(xi)) {
				recFct  <-  get(names(xi)[j])
				xi[j]   <-  round(recFct(Fii = Fii, m = m, r = r), digits=8)
			}
			# 2) Offspring genotype frequencies
			O      <-  offFreq(xi)
			# 3) Mean fitness 
			Wbar   <-  sum(O*W)
			# 4) Expected frequencies
			E.Fii  <-  O*W/Wbar
			# 5) Draw random frequencies in adults
			Fii    <-  as.vector(rmultinom(1, N/2, E.Fii)/(N/2))
	
			# Realized frequencies
			InvFreq    <-  sum(Fii[c(5,10,15,20:24)]/2, Fii[25])
			E.InvFreq  <-  sum(E.Fii[c(5,10,15,20:24)]/2, Fii[25])
			W.mean     <-  Wbar

			gen  <-  gen+1
		}
			
		# Has the inversion reached threshold frequency for establishment (pcrit)? 
		# When did it first reach pcrit?
		if(InvFreq >= pcrit) {
			invEst      <-  1
			invEstTime  <-  gen
		} else {
			invEst      <-  0
			invEstTime  <-  NA
		}
	
		# Save  simulation data
		res  <-  list(
					  "InvFreq"     =  InvFreq,
					  "E.InvFreq"   =  E.InvFreq,
					  "W.mean"      =  sum(W.mean)/length(W.mean),
					  "nGen"        =  gen,
					  "invEst"      =  invEst,
					  "invEstTime"  =  invEstTime
 				 	  )
	}

return(res)

}


#' Wrapper function to run replicate forward simulations for invasion
#' of autosomal inversions in a Wright-Fisher population 
#'
#' @title Wright-Fisher forward simulation of genotypic frequencies (default parameter values in parentheses)
#' @param nReps  Numer of replicate simulations. With no deleterious mutations, and introducing 
#' 				 a single copy of the inversion, it takes 1,600,000 replicate simulations to 
#'				 get 10,000 where the inversion successfully establishes.
#' @param N      Effective population size
#' @param m      Migration rate for locally maladaptive alleles (m =  0.01)
#' @param s      Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02)
#' @param h      Dominance coefficient for locally adaptive alleles relative to migrant alleles (h = 0.5)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param n      Number of loci at which deleterious mutations may occur.
#' @param u      Mutation rate (default value of u = 1e-6).
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Defalut value of FALSE assumes that 
#' 				 selection against deleterious mutations is twice as strong as selection favouring
#' 				 the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				 there were no delterious mutations segregating in the population that are linked
#' 				 to the loci involved in local adaptation. 
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? Setting this 
#' 							to TRUE can become extremely memory intensive if you are running many
#' 							replicate simulations (see warning).
#' 							otherwise
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSim`
#' @export
#' @author Colin Olito.
runReplicateAutoInvSims  <-  function(nReps = 1000, N = 500, m = 0.01, s = 0.1, h = 1/2, r = 0.1, 
									  n = 100, u = 1e-5, h.del = 0, s.del = 1, noDel = FALSE,
									  saveTrajectories = FALSE) {

	##  Preemptive Warnings
	if(any(c(N,m,s,h,r,n,u,h.del) < 0) | any(c(m,s,h,r,u,h.del) > 1) | r > 0.5)
		stop('The chosen parameter values fall outside of reasonable parameter space')

	try({
		 if(m >= s )
			stop('Warning: migration is stronger than than selection, 
				  adaptive alleles will be swamped by maladaptive migrants')
	}, silent=FALSE)

	try({
		 if(nReps > 1000 & saveTrajectories)
			stop('Warning: You have chosen to save evolutionary trajectories 
				  for a large number of replicate simulations. Thiss will be
				  memory intensive. Consider setting saveTrajectories = FALSE')
	}, silent=FALSE)

	##  Define Fitness Expressions for determining eq. frequencies in absence of inversion
	W.init  <-  c(1,          (1 + h*s),         (1 + h*s),         (1 + h*s)^2,        0,
			     (1 + h*s),   (1 + s),           (1 + h*s)^2,       (1 + h*s)*(1 + s),  0,
			     (1 + h*s),   (1 + h*s)^2,       (1 + s),           (1 + s)*(1 + h*s),  0,
			     (1 + h*s)^2, (1 + h*s)*(1 + s), (1 + s)*(1 + h*s), (1 + s)^2,          0,
			     0,           0,                 0,                 0,                  0)
  		
 	## Find deterministic equilibrium frequencies in absence of inversion  
	Fii.init  <-  findEqFreqs(W=W.init, m=m, r=r, threshold=1e-7)

	# Use deterministic equilibrium frequencies of non-inversion genotypes
	# as initial conditions when introducing the inversion via a single
	# copy of the abba* genotype 
	Fii.init[19]  <-  Fii.init[19] - 1/N
	Fii.init[20]  <-  1/N
	
	# Storage structures for replicate simulation data
	finalInvFreq    <-  rep(0, times=nReps)
	finalE.InvFreq  <-  rep(0, times=nReps)
	finalW.mean     <-  rep(0, times=nReps)
	nGen            <-  rep(0, times=nReps)
	invEst          <-  rep(0, times=nReps)
	invEstTime      <-  rep(0, times=nReps)
	nDels           <-  rep(0, times=nReps)

	if(saveTrajectories) {
		replicateTraj  <-  c()
		InvFreqTraj    <-  c()
		E.InvFreqTraj  <-  c()
		W.meanTraj     <-  c()
	} 

	# Replicate simulation loop
#	print('Running Wright-Fisher Forward Simulations')
	pb   <-  txtProgressBar(min=0, max=nReps, style=3)
	setTxtProgressBar(pb, 0)
	for(i in 1:nReps) {

	## Sample stationary distribution of deleterious alleles
	delMutFreq  <-  rejectionSampler(n=n, Ne=N, u=u)
	n.del       <-  sum(delMutFreq > runif(n=n))

	# Define fitness expressions, including fitness effects of deleterious mutations
	if(noDel) {
		s.del  <-  0
	}
	W  <-  c(1,                                  (1 + h*s),                                 (1 + h*s),                                 (1 + h*s)^2,                       (1 + h*s)^2*(1 - h.del*s.del)^n.del,
			(1 + h*s),                           (1 + s),                                   (1 + h*s)^2,                               (1 + h*s)*(1 + s),                 (1 + h*s)*(1 + s)*(1 - h.del*s.del)^n.del,
			(1 + h*s),                           (1 + h*s)^2,                               (1 + s),                                   (1 + s)*(1 + h*s),                 (1 + s)*(1 + h*s)*(1 - h.del*s.del)^n.del,
			(1 + h*s)^2,                         (1 + h*s)*(1 + s),                         (1 + s)*(1 + h*s),                         (1 + s)^2,                         (1 + s)^2*(1 - h.del*s.del)^n.del,
			(1 + h*s)^2*(1 - h.del*s.del)^n.del, (1 + h*s)*(1 + s)*(1 - h.del*s.del)^n.del, (1 + s)*(1 + h*s)*(1 - h.del*s.del)^n.del, (1 + s)^2*(1 - h.del*s.del)^n.del, (1 + s)^2*(1 - s.del)^n.del)

		## RUN SIMULATION
		repRes  <-  autoInvFwdSim(Fii.init=Fii.init, N=N, W=W, m=m, r=r)

		# save results for each replicate
		finalInvFreq[i]    <-  repRes$InvFreq[length(repRes$InvFreq)]
		finalE.InvFreq[i]  <-  repRes$E.InvFreq[length(repRes$E.InvFreq)]
		finalW.mean[i]     <-  repRes$W.mean
		nGen[i]            <-  repRes$nGen
		invEst[i]          <-  repRes$invEst
		invEstTime[i]      <-  repRes$invEstTime
		nDels[i]           <-  n.del

		if(saveTrajectories) {
			replicateTraj  <-  c(replicateTraj, rep(i, times=length(repRes$InvFreq)))
			InvFreqTraj    <-  c(InvFreqTraj, repRes$InvFreq)
			E.InvFreqTraj  <-  c(E.InvFreqTraj, repRes$E.InvFreq)
			W.meanTraj     <-  c(W.meanTraj, repRes$W.mean)
			} 

	setTxtProgressBar(pb, i)
	}

	# Save results and return results as a list
	results.df  <-  data.frame(
							   "finalInvFreq"    =  finalInvFreq,
							   "finalE.InvFreq"  =  finalE.InvFreq,
							   "finalW.mean"     =  finalW.mean,
							   "nGen"            =  nGen,
							   "invEst"          =  invEst,
							   "invEstTime"      =  invEstTime,
							   "nDels"           =  nDels
							   )
	if(saveTrajectories) {
		traj.df  <-  data.frame(
								"replicateTraj"  =  replicateTraj,
								"InvFreqTraj"    =  InvFreqTraj,
								"E.InvFreqTraj"  =  E.InvFreqTraj,
								"W.meanTraj"     =  W.meanTraj
								)
	} else {
		traj.df  <-  NULL
	}
	res  <-  list(
				  "results.df"  =  results.df,
				  "traj.df"     =  traj.df
				  )
	return(res)
}


#' Wrapper function to run replicate forward simulations for invasion
#' of autosomal inversions in a Wright-Fisher population 
#' USING DESIGNATED PARAMETER VALUES
#'
#' @title Run replicate Wright-Fisher forward simulations for autosomal inversion under different parameter values 
#' @param N.vals Desired population sizes
#' @param m.vals desired migration rates for locally maladaptive alleles (m =  0.01)
#' @param s.del.vals  Desired selection coefficients for deleterious mutations (default value of s = 1).
#' @param s      Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02)
#' @param h      Dominance coefficient for locally adaptive alleles relative to migrant alleles (h = 0.5)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param n      Number of loci at which deleterious mutations may occur.
#' @param u      Mutation rate (default value of u = 1e-6).
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Defalut value of FALSE assumes that 
#' 				 selection against deleterious mutations is twice as strong as selection favouring
#' 				 the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				 there were no delterious mutations segregating in the population that are linked
#' 				 to the loci involved in local adaptation. 
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? Setting this 
#' 							to TRUE can become extremely memory intensive if you are running many
#' 							replicate simulations (see warning).
#' 							otherwise
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSim`
#' @export
#' @author Colin Olito.
makeReplicateAutoInvSimsData  <-  function(nReps = 1000, N.vals = c(500, 1000), m.vals = c(0.01, 0.05), 
										   s = 0.1, h = 1/2, r = 0.1, 
										   n = 100, u = 1e-5, h.del = 0) {

	# Simulate deleterious mutations that are either 
	# 1) recessive lethals OR
	# 2) recessive experiencing purifying selection
	#    that is twice as strong as the selective 
	#    advantage of the locally adaptive allels  
	s.del.vals = c(0, 1, 2*s)


	# create empty data frame with same structure as we are going to need
	data  <-  data.frame(matrix(ncol=10, nrow=0))

	# Convenience variables to monitor progress
	prog  <-  0
	tot   <-  length(N.vals)*length(m.vals)*length(s.del.vals)

	# Loop over parameter values we want to explore 
	for(j in 1:length(N.vals)) {
		for(k in 1:length(m.vals)) {
			for(l in 1:length(s.del.vals)) {

				# Display progress in terminal
				prog  <-  prog + 1
				cat("\n",paste('Running simulations for parameter set ', prog, "/", tot),"\n")

				# Run simulations  
				res  <-  runReplicateAutoInvSims(nReps = nReps, N = N.vals[j], m = m.vals[k], s = s, h = h, r = r, 
												 n = n, u = u, h.del = h.del, s.del = s.del.vals[l], 
												 noDel = FALSE, saveTrajectories = FALSE)

				# Save data 
				Ns      <-  rep(N.vals[j], times=nrow(res$results.df))
				ms      <-  rep(m.vals[k], times=nrow(res$results.df))
				s.dels  <-  rep(s.del.vals[l], times=nrow(res$results.df))

				# Append to data frame
				df      <-  cbind(res$results.df, Ns, ms, s.dels)
				data    <-  rbind(data, df)
				rm(df)

			}
		}
	}

	# Include constant variables in data frame
	ss      <-  rep(s, times=nrow(data))
	hs      <-  rep(h, times=nrow(data))
	rs      <-  rep(r, times=nrow(data))
	us      <-  rep(u, times=nrow(data))
	h.dels  <-  rep(h.del, times=nrow(data))
	data    <-  cbind(data, ss, hs, rs, us, h.dels)
	colnames(data)  <-  c("finalInvFreq","finalE.InvFreq","finalW.mean",
						  "nGen","invEst","invEstTime","nDels","N","m",
						  "s.dels","s","h","r","u","h.del")

	# create file name
	filename  <-  paste("./output/data/simResults/auto-InvSimsData", "_s", s, "_h", h, "_r", r, "_n", n, "_u", u, ".csv", sep="")

	# export data as .csv to ./output/data
	write.csv(data, file=filename, row.names = FALSE)

	#  Return results in case user wants it
	return(data)
	
}