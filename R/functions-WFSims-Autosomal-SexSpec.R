###############
# DEPENDENCIES
###############

######################
# Necessary functions  
######################

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


#' Linkage Disequilibrium function for W-F recursions
#'
#' @title Linkage Disequilibrium (Dstar)
#' @param Fii  Accepts a vector of adult genotypic frequencies in females (Fii.f; of length = 25)
#' 				OR a vector of adult genotypic frequencies in males (Fii.m; of length = 25)
#' @param m   Migration rate
#' @export

Dstar  <-  function(Fii, m, ...) {
	(((Fii[4] + Fii[16]) - (Fii[8] + Fii[12])) / 2)*(1 - m)
}

		# genotypes ordered:
		#Females
		#	c(ABAB, ABAb, ABaB, ABab, ABba*,		c(x1x1, x1x2, x1x3, x1x4, x1x5, 
		#	  AbAB, AbAb, AbaB, Abab, Abba*,		  x2x1, x2x2, x2x3, x2x4, x2x5, 
		#	  aBAB, aBAb, aBaB, aBab, aBba*,		  x3x1, x3x2, x3x3, x3x4, x3x5, 
		#	  abAB, abAb, abaB, abab, abba*,		  x4x1, x4x2, x4x3, x4x4, x4x5, 
		#	  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
		#Males		
		#	c(ABAB, ABAb, ABaB, ABab, ABba*,		c(y1y1, y1y2, y1y3, y1y4, y1y5, 
		#	  AbAB, AbAb, AbaB, Abab, Abba*,		  y2y1, y2y2, y2y3, y2y4, y2y5, 
		#	  aBAB, aBAb, aBaB, aBab, aBba*,		  y3y1, y3y2, y3y3, y3y4, y3y5, 
		#	  abAB, abAb, abaB, abab, abba*,		  y4y1, y4y2, y4y3, y4y4, y4y5, 
		#	  baAB*, baAb*, baaB*, baab*, baba*)	  y5y1, y5y2, y5y3, y5y4, y5y5)

#' Haplotype frequencies among gametes
#'
#' @title Haplotype frequencies among gametes
#' @param Fii  Accepts a vector of adult genotypic frequencies in females (Fii.f; of length = 25)
#' 				OR a vector of adult genotypic frequencies in males (Fii.m; of length = 25)
#' @param m   Migration rate
#' @param r   Recombination rate
#' @export
# Haplotype frequency equations for ovules
x.1  <-  function(Fii=Fii.f, m=mf, r=rf) { 
	((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
x.2  <-  function(Fii=Fii.f, m=mf, r=rf) {
	((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.3  <-  function(Fii=Fii.f, m=mf, r=rf) {
	((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.4  <-  function(Fii=Fii.f, m=mf, r=rf) {
	((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
x.5  <-  function(Fii=Fii.f, m=mf, r=rf) {
	((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}
# Haplotype frequency equations for sperm/pollen
y.1  <-  function(Fii=Fii.m, m=mm, r=rm) {
	((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
y.2  <-  function(Fii=Fii.m, m=mm, r=rm) {
	((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.3  <-  function(Fii=Fii.m, m=mm, r=rm) {
	((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.4  <-  function(Fii=Fii.m, m=mm, r=rm) {
	((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
y.5  <-  function(Fii=Fii.m, m=mm, r=rm) {
	((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}


#' Offspring frequencies after random mating
#'
#' @title Offspring frequencies after random mating, ( males and females have the same allele frequencies)
#' @param xi Vector of haplotype frequencies among female gametes (of length = 5)
#' @param yi Vector of haplotype frequencies among male gametes (of length = 5)
#' @export
offFreq  <-  function(xi,yi) {
	O  <-  c((xi[1]*yi[1]), (xi[1]*yi[2]), (xi[1]*yi[3]), (xi[1]*yi[4]), (xi[1]*yi[5]),
			 (xi[2]*yi[1]), (xi[2]*yi[2]), (xi[2]*yi[3]), (xi[2]*yi[4]), (xi[2]*yi[5]),
			 (xi[3]*yi[1]), (xi[3]*yi[2]), (xi[3]*yi[3]), (xi[3]*yi[4]), (xi[3]*yi[5]),
			 (xi[4]*yi[1]), (xi[4]*yi[2]), (xi[4]*yi[3]), (xi[4]*yi[4]), (xi[4]*yi[5]),
			 (xi[5]*yi[1]), (xi[5]*yi[2]), (xi[5]*yi[3]), (xi[5]*yi[4]), (xi[5]*yi[5])
			 )
	O
}


#' Calculate deterministic equilibrium genotype frequencies 
#' in the absence of the inversion 
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#' Females
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1x1, x1x2, x1x3, x1x4, x1x5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
#' Males
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(y1y1, y1y2, y1y3, y1y4, y1y5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  y2y1, y2y2, y2y3, y2y4, y2y5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  y3y1, y3y2, y3y3, y3y4, y3y5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  y4y1, y4y2, y4y3, y4y4, y4y5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  y5y1, y5y2, y5y3, y5y4, y5y5)
#' @title Find the deterministic equilibeium genotype frequencies prior to introducing the inversion
#' @param Wf     Vector of fitness expressions for all 25 genotypes in females
#' @param Wm     Vector of fitness expressions for all 25 genotypes in males
#' @param mf     Female migration rate for locally maladaptive alleles
#' @param mm     Male migration rate for locally maladaptive alleles
#' @param rf     Female recombination rate among the two loci involved in local adaptation
#' @param rm     Male recombination rate among the two loci involved in local adaptation
#' @param threshold The threshold change in genotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. 
#' @export
#' @seealso `offFreq`, `autoInvWrightFisherSim`
#' @author Ludovic Dutoit based on Colin Olito
findEqFreqs  <-  function(Wf,Wm, mf, mm, rf, rm, threshold = 1e-6, ...) {
	Fii.init  <-  c(1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					   0,    0,    0,    0, 0)
	Fii.f    <-  Fii.init
	E.Fii.f  <-  Fii.init

	Fii.m    <-  Fii.init
	E.Fii.m  <-  Fii.init

	# Storage for female gamete frequencies
	xi         <-  rep(0, times=5)
	names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')

	# Storage for female gamete frequencies
	yi         <-  rep(0, times=5)
	names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')
	
	# Storage for differences in frequencies across 1 generation
	deltaF  <-  rep(1, times=16)
	deltaM  <-  rep(1, times=16)
	delta   <-  append(deltaF,deltaM)
	
	# Initiate generation counter 
	gen  <-  0 
	
	# Simulation loop
	while(any(delta > threshold)) {
		gen  <-  gen+1
		for (j in 1:length(xi)) {
			recFct  <-  get(names(xi)[j])
			xi[j]   <-  round(recFct(Fii = E.Fii.f, m = mf, r = rf), digits=6)
			recFct  <-  get(names(yi)[j])
			yi[j]   <-  round(recFct(Fii = E.Fii.m, m = mm, r = rm), digits=6)
		}

	    # offspring genotype frequencies
	    O  <-  offFreq(xi,yi)

		# mean fitness 
		Wfbar      <-  sum(O*Wf)
		Wmbar      <-  sum(O*Wm)

	    # difference in expected frequency (has simulation reached equilibrium yet?)
		deltaF   <-  E.Fii.f[c(1:4,6:9,11:14,16:19)] - (O*Wf/Wfbar)[c(1:4,6:9,11:14,16:19)]
		deltaM   <-  E.Fii.m[c(1:4,6:9,11:14,16:19)] - (O*Wm/Wmbar)[c(1:4,6:9,11:14,16:19)]
		delta    <-  append(deltaF,deltaM)

		# Expected adult genotypic frequencies in next generation
		E.Fii.f   <-  O*Wf/Wfbar
		E.Fii.m   <-  O*Wm/Wmbar
	}

	names(E.Fii.f)  <-  NULL
	names(E.Fii.m)  <-  NULL
	return (rbind(E.Fii.f,E.Fii.m) )
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
#' @param Fii.f.init  Vector of initial frequencies (deterministic eq. frequencies in absence of inversion) in females
#' @param Fii.m.init  Vector of initial frequencies (deterministic eq. frequencies in absence of inversion) in males
#' @param N           Population size
#' @param Wf          Vector of fitness expressions for all 25 genotypes in females
#' @param Wm          Vector of fitness expressions for all 25 genotypes in males
#' @param mm          Migration rate for locally maladaptive alleles (m =  0.01) in males
#' @param mf          Migration rate for locally maladaptive alleles (m =  0.01) in females
#' @param rf          Female recombination rate among the two loci involved in local adaptation
#' @param rm          Male recombination rate among the two loci involved in local adaptation
#' @param fastSim     Logical. Use threshold frequency for establishment of inversion? 
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? 
#' @export
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Ludovic Dutoit, based on Colin Olito
autoInvFwdSimSexSpec  <-  function(Fii.f.init = Fii.f.init, Fii.m.init = Fii.m.init, 
									N = N, Wf = Wf, Wm = Wm, mm = mm, mf = mf, sm = sm, sf = sf, rf = rf, rm = rm, 
									fastSim = TRUE, saveTrajectories = FALSE, ...) {

	# Use deterministic eq. initial frequencies
	Fii.f  <-  Fii.f.init 
	Fii.m  <-  Fii.m.init

	# Define threshold frequency for establishment of inversion
	pcrit  <-  2/(N*((mm+mf)/2))

	# Storage for gamete frequencies  
	xi         <-  rep(0, times=5)
	yi         <-  rep(0, times=5)
	names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')
	names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')

	# use threshold frequency for establishment of inversion?
	if(fastSim) {

		# Storage structures for individual simulation data
		InvFreq    <-  0
		E.InvFreq  <-  0
		W.mean     <-  0
		# Initial inversion frequency 
		InvFreq    <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) +  sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))/2
		E.InvFreq  <-  InvFreq[1]
		
		## Start forward simulation with newly introduced inversion
		gen  <-  1

		# use threshold frequency for establishment of inversion?
		while(InvFreq > 0 & InvFreq < pcrit) {		
		
			## Step through recursions:
			# 1) Calculate gamete frequencies
			for (j in 1:length(xi)) {
				recFctx  <-  get(names(xi)[j])
				recFcty  <-  get(names(yi)[j])
				xi[j]    <-  round(recFctx(Fii = Fii.f, m = mf, r = rf), digits=8)
				yi[j]    <-  round(recFcty(Fii = Fii.m, m = mm, r = rm), digits=8)
			}
			# 2) Offspring genotype frequencies
			O          <-  offFreq(xi,yi) 
			# 3) Mean fitness 
			Wbar.f     <-  sum(O*Wf)
			Wbar.m     <-  sum(O*Wm)
			Wbar       <-  (Wbar.f + Wbar.m)/2
			# 4) Expected frequencies
			E.Fii.f    <-  (O*Wf)/Wbar.f # sex specific expected freq
			E.Fii.m    <-  (O*Wm)/Wbar.m # sex specific expected freq
			# 5) Draw random frequencies in adults
			Fii.f      <-  as.vector(rmultinom(1, N/2, E.Fii.f)/(N/2))
			Fii.m      <-  as.vector(rmultinom(1, N/2, E.Fii.m)/(N/2))
			# Realized frequencies
			InvFreq    <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))/2
			E.InvFreq  <-  (sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25]) + sum(E.Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))/2
			W.mean     <-  Wbar
			# next gen
			gen        <-  gen + 1		
		}

		# Store relevant data 
		Wf.mean      <-  Wbar.f
		Wm.mean      <-  Wbar.m
		InvFreq_f    <-  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25])
		InvFreq_m    <-  sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25])
		E.InvFreq_f  <-  sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25])
		E.InvFreq_m  <-  sum(E.Fii.m[c(5,10,15,20:24)]/2, Fii.m[25])

		# Has the inversion reached threshold frequency for establishment (pcrit)? 
		# When did it first reach pcrit?
		if(InvFreq >= pcrit) {
			invEst      <-  1
			invEstTime  <-  gen
		} else {	
			invEst      <-  0
			invEstTime  <-  NA
		}
		res  <-  list(
				  	"InvFreq"      =  InvFreq,
				  	"E.InvFreq"    =  E.InvFreq,
				  	"Wbar"         =  sum(W.mean)/length(W.mean),
					"Wbar_f"       =  sum(Wf.mean)/length(Wf.mean), 
					"Wbar_m"       =  sum(Wm.mean)/length(Wm.mean),
					"InvFreq_f"    =  InvFreq_f,
					"InvFreq_m"    =  InvFreq_m, 
					"E.InvFreq_f"  =  E.InvFreq_f,
					"E.InvFreq_m"  =  E.InvFreq_m,
					"nGen"         =  gen,
				  	"InvEst"       =  invEst,
				  	"InvEstTime"   =  invEstTime
			 )
	return(res)
	}
	else {	# Run simulations until 4N generations

		if(saveTrajectories) {
			# Storage structures for individual simulation data
			InvFreq      <-  rep(0, times=(4*N+1))
			E.InvFreq    <-  rep(0, times=(4*N+1))
			W.mean       <-  rep(0, times=(4*N+1))
			Wf.mean      <-  rep(0, times=(4*N+1)) 
			Wm.mean      <-  rep(0, times=(4*N+1))
			InvFreq_f    <-  rep(0, times=(4*N+1))
			InvFreq_m    <-  rep(0, times=(4*N+1)) 
			E.InvFreq_f  <-  rep(0, times=(4*N+1))
			E.InvFreq_m  <-  rep(0, times=(4*N+1))	
			
			# Initial inversion frequency 
			InvFreq[1]    <-  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25])
			E.InvFreq[1]  <-  InvFreq[1]
			
			## Start forward simulation with newly introduced inversion
			gen  <-  1
			while(gen < (4*N) & InvFreq[gen] > 0 ) {
	
				## Step through recursions:
				# 1) Calculate gamete frequencies
				for (j in 1:length(xi)) {
					recFctx  <-  get(names(xi)[j])
					recFcty  <-  get(names(yi)[j])
					xi[j]    <-  round(recFctx(Fii = Fii.f, m = mf, r = rf), digits=6)
					yi[j]    <-  round(recFcty(Fii = Fii.m, m = mm, r = rm), digits=6)
				}
				# 2) Offspring genotype frequencies
				O        <-  offFreq(xi,yi)
				# 3) Mean fitness 
				Wbar.f   <-  sum(O*Wf)
				Wbar.m   <-  sum(O*Wm)
				Wbar     <-  (sum(Wbar.f)+sum(Wbar.m))/2
				# 4) Expected frequencies
				E.Fii.f  <-  O*Wf/Wbar.f 
				E.Fii.m  <-  O*Wm/Wbar.m 
				# 5) Draw random frequencies in adults
				Fii.f    <-  as.vector(rmultinom(1, (N/2), E.Fii.f)/(N/2)) 
				Fii.m    <-  as.vector(rmultinom(1, (N/2), E.Fii.m)/(N/2))
				
	
				# Realized frequencies
				InvFreq[gen+1]     <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) +  sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))/2
				E.InvFreq[gen+1]   <-  (sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25]) +  sum(E.Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))/2
				W.mean[gen+1]      <-  Wbar
				
				#The variable below are stored specifically for the SexSpecific simulation 
				Wf.mean[gen+1]      <-  sum(O*Wf)
				Wm.mean[gen+1]      <-  sum(O*Wm)
				InvFreq_f[gen+1]    <-  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25])
				InvFreq_m[gen+1]    <-  sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25])
				E.InvFreq_f[gen+1]  <-  sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25])
				E.InvFreq_m[gen+1]  <-  sum(E.Fii.m[c(5,10,15,20:24)]/2, Fii.m[25])
				# next gen
				gen  <-  gen+1
			}
				
			# Has the inversion reached threshold frequency for establishment (pcrit)? 
			# When did it first reach pcrit?
			if(any(InvFreq >= pcrit)) {
				invEst      <-  1
				invEstTime  <-  gen[InvFreq >= pcrit][1]
			} else {	
				invEst      <-  0
				invEstTime  <-  NA
			}
		
			# Save  simulation data
			res  <-  list(
					  	"InvFreq"      =  InvFreq[1:gen-1],
					  	"E.InvFreq"    =  E.InvFreq[1:gen-1],
					  	"Wbar"         =  W.mean[1:gen-1],
						"Wbar_f" 	   =  Wf.mean[1:gen-1], 
						"Wbar_m" 	   =  Wm.mean[1:gen-1],
						"InvFreq_f"    =  InvFreq_f[1:gen-1],
						"InvFreq_m"    =  InvFreq_m[1:gen-1], 
						"E.InvFreq_f"  =  E.InvFreq_f[1:gen-1],
						"E.InvFreq_m"  =  E.InvFreq_m[1:gen-1],
						"nGen"         =  gen,
					  	"InvEst"       =  invEst,
					  	"InvEstTime"   =  invEstTime
	 				 	)
		} 
	
		if(!saveTrajectories) {
			# Storage structures for individual simulation data
			InvFreq    <-  0
			E.InvFreq  <-  0
			W.mean     <-  0
			# Initial inversion frequency 
			InvFreq    <-  0.5 * (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) +  sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))
			E.InvFreq  <-  InvFreq[1]
			
			## Start forward simulation with newly introduced inversion
			gen  <-  1
	
			while(gen < (4*N) & InvFreq[gen] > 0 ) {
		
				## Step through recursions:
				# 1) Calculate gamete frequencies
				for (j in 1:length(xi)) {
					recFctx  <-  get(names(xi)[j])
					recFcty  <-  get(names(yi)[j])
					xi[j]    <-  round(recFctx(Fii = Fii.f, m = mf, r = rf), digits=8)
					yi[j]    <-  round(recFcty(Fii = Fii.m, m = mm, r = rm), digits=8)
				}
				# 2) Offspring genotype frequencies
				O          <-  offFreq(xi,yi) # no need of two vectors, we just apply twice selection 
				# 3) Mean fitness 
				Wbar.f     <-  sum(O*Wf)
				Wbar.m     <-  sum(O*Wm)	
				# 4) Expected frequencies
				E.Fii.f    <-  O*Wf/Wbar.f # sex specific expected freq
				E.Fii.m    <-  O*Wm/Wbar.m # sex specific expected freq
				# 5) Draw random frequencies in adults
				Fii.f      <-  as.vector(rmultinom(1, N/2, E.Fii.f)/(N/2))
				Fii.m      <-  as.vector(rmultinom(1, N/2, E.Fii.m)/(N/2))
				# Realized frequencies
				InvFreq    <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) +  sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))/2
				E.InvFreq  <-  (sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25]) +  sum(E.Fii.m[c(5,10,15,20:24)]/2, Fii.m[25]))/2
				# next gen
				gen        <-  gen + 1		
				#The variable below are stored specifically for the SexSpecific simulation 
				
			}
			
			# Store relevant data 
			Wf.mean      <-  sum(O*Wf)
			Wm.mean      <-  sum(O*Wm)
			InvFreq_f    <-  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25])
			InvFreq_m    <-  sum(Fii.m[c(5,10,15,20:24)]/2, Fii.m[25])
			E.InvFreq_f  <-  sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25])
			E.InvFreq_m  <-  sum(E.Fii.m[c(5,10,15,20:24)]/2, Fii.m[25])
	
			# Has the inversion reached threshold frequency for establishment (pcrit)? 
			# When did it first reach pcrit?
			if(InvFreq >= pcrit) {
				invEst      <-  1
				invEstTime  <-  gen
			} else {	
				invEst      <-  0
				invEstTime  <-  NA
			}
			res  <-  list(
						"InvFreq"      =  InvFreq,
						"E.InvFreq"    =  E.InvFreq,
						"Wbar_f"       =  sum(Wf.mean)/length(Wf.mean), 
						"Wbar_m"       =  sum(Wm.mean)/length(Wm.mean),
						"InvFreq_f"    =  InvFreq_f,
						"InvFreq_m"    =  InvFreq_m, 
						"E.InvFreq_f"  =  E.InvFreq_f,
						"E.InvFreq_m"  =  E.InvFreq_m,
						"nGen"         =  gen,
						"InvEst"       =  invEst,
						"InvEstTime"   =  invEstTime
				 )
		}
	return(res)
	}
}



#' Introduce new mutant inverion genotype
#'
#' @title Introduce new mutant inverion genotype
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#' 					 chosen randomly, given initial genotypic frequencies (Fii.f.init and Fii.m.init). 
#'					 A 2 positions vector ("m"/"f"/"random", "numeric"/"random"). the first position specify wether the inversions come in males, 
#'					 in females, or randomnly. The second select either a numerical position in the haplotype vector for the inversion genotype to be created
#'					 or wether it is random.
#' 					 See params for runReplicateAutoInvSims().
#' @param Fii.f.init   Initial genotypic frequencies, from which to calculate probability of new 
#' 					 mutant inversion occurring
#' @param Fii.m.init   Initial genotypic frequencies, from which to calculate probability of new 
#' 					 mutant inversion occurring
#' @param N			 Population size
introduceInversion  <-  function(newMutant, Fii.f.init, Fii.m.init, N, ...) {
	
	#extract newMutant Parameters
	sex        <-  newMutant[1]
	specified  <-  is.numeric(newMutant[2])
	mutant     <-  newMutant[2]

	#Preemptive warnings
	if(sex!="m" &  sex!="f" & sex!="random") 
			stop('Warning: newMutant is of wrong format')
	if(specified & all(mutant != c(4,9,14,16:19)))
		stop('If specifying the genotype of new inversion mutants, newMutant must take 
			  one of the following values: 4,9,14,16:19')
	if(!specified & newMutant[2] != 'random')
		stop('If the genotype of new inversion mutants is being chosen randomly, 
			  the parameter newMutant must equal random')

	# First determine which sex the inversion will occur if it is random
	if (sex == "random"){
		sex  <-  sample(c("m","f"), 1)
	}

	# Choose mutant genotype randomly
	if(!specified) {
		if (sex == "f"){
			probNewMutant       <-  c(Fii.f.init[c(4,9,14,16:18)], Fii.f.init[19]*2)/sum(c(Fii.f.init[c(4,9,14,16:18)], Fii.f.init[19]*2))
			mutant              <-  c(4,9,14,16:19)[as.vector(rmultinom(1,1,probNewMutant)) == 1]

			# Subtract new mutant individual from frequency of old genotype
			Fii.f.init[mutant]  <-  Fii.f.init[mutant] - 1/(N/2)
		}
		if (sex == "m"){
			probNewMutant       <-  Fii.m.init[c(4,9,14,16:19)]/sum(Fii.m.init[c(4,9,14,16:19)])
			mutant              <-  c(4,9,14,16:19)[as.vector(rmultinom(1,1,probNewMutant)) == 1]

			# Subtract new mutant individual from frequency of old genotype
			Fii.m.init[mutant]  <-  Fii.m.init[mutant] - 1/(N/2)
		}
	}

	# Specify mutant genotype
	if(specified) {
		if (sex=="f"){
			# Subtract new mutant individual from frequency of old genotype
			Fii.f.init[mutant]  <-  Fii.f.init[mutant] - 1/(N/2)
	    }
		if (sex=="m"){
			# Subtract new mutant individual from frequency of old genotype
			Fii.m.init[mutant]  <-  Fii.m.init[mutant] - 1/(N/2)
	    }	    
	}

	# Add mutant individual to frequency of new inversion genotype
	if (sex=="f"){
		if(mutant == 4 | mutant == 9 | mutant == 14)
			Fii.f.init[mutant + 1]  <-  1/(N/2)
		if(mutant == 16 | mutant == 17 | mutant == 18)
			Fii.f.init[mutant + 5]  <-  1/(N/2)

		# if inversion occurs on abab genotype, choose randomly whether it occurs on
		# the maternally or paternally inherited chromosome 
		if(mutant == 19) {
			if(runif(1) >= 1/2) {
				Fii.f.init[mutant + 1]   <-  1/(N/2)
			}
			else Fii.f.init[mutant + 5]  <-  1/(N/2)
		}
	}
	if (sex=="m"){
		if(mutant == 4 | mutant == 9 | mutant == 14)
			Fii.m.init[mutant + 1]  <-  1/(N/2)
		if(mutant == 16 | mutant == 17 | mutant == 18)
			Fii.m.init[mutant + 5]  <-  1/(N/2)

		# if inversion occurs on abab genotype, choose randomly whether it occurs on
		# the maternally or paternally inherited chromosome 
		if(mutant == 19) {
			if(runif(1) >= 1/2) {
				Fii.m.init[mutant + 1]   <-  1/(N/2)
			}
			else Fii.m.init[mutant + 5]  <-  1/(N/2)
		}
	}
	return (rbind(Fii.f.init, Fii.m.init))
}


#' Wrapper function to run replicate forward simulations for invasion
#' of autosomal inversions in a Wright-Fisher population with sex specific selection 
#'
#' @title Wright-Fisher forward simulation of genotypic frequencies (default parameter values in parentheses)
#' @param nReps  Number of replicate simulations. With no deleterious mutations, and introducing 
#' 				 a single copy of the inversion, it takes 1,600,000 replicate simulations to 
#'				 get 10,000 where the inversion successfully establishes.
#' @param N      Effective population size
#' @param mf     Migration rate for locally maladaptive alleles (m =  0.01) in females
#' @param mm     Migration rate for locally maladaptive alleles (m =  0.01) in males
#' @param sf     Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02) in females
#' @param sm     Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02) in males
#' @param h      Dominance coefficient for locally adaptive alleles relative to migrant alleles (h = 0.5)
#' @param rf     Female recombination rate among the two loci involved in local adaptation
#' @param rm     Male recombination rate among the two loci involved in local adaptation
#' @param n      Number of loci at which deleterious mutations may occur.
#' @param u      Mutation rate (default value of u = 1e-6).
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Defalut value of FALSE assumes that 
#' 				 selection against deleterious mutations is twice as strong as selection favouring
#' 				 the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				 there were no delterious mutations segregating in the population that are linked
#' 		 		 to the loci involved in local adaptation. 
#' @param fastSim     Logical. Use threshold frequency for establishment of inversion? 
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? Setting this 
#' 							to TRUE can become extremely memory intensive if you are running many
#' 							replicate simulations (see warning).
#' 							otherwise
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#' 					 chosen randomly, given initial genotypic frequencies (Fii.f.init and Fii.m.init). 
#'					 A 2 positions vector ("m"/"f"/"random", "numeric"/"random"). the first position specify wether the inversions come in males, 
#'					 in females, or randomnly. The second select either a numerical position in the haplotype vector for the inversion genotype to be created
#'					 or wether it is random.
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSimSexSpec`
#' @export
#' @author Ludovic Dutoit based on Colin Olito.
runReplicateAutoInvSimsSexSpec  <-  function(nReps = 1000, N = 5000, mm = 0.01, mf = 0.01, sf = 0.1, sm = 0.1, h = 1/2, rf = 0.5, rm = 0.5, 
									  		 n = 100, u = 1e-5, h.del = 0, s.del = 1, noDel = FALSE,
									  		 fastSim = TRUE, saveTrajectories = FALSE, newMutant = c("random","random")) {

	##  Preemptive Warnings
	if(any(c(N,mm,mf,sf,sm,h,rf,rm,n,u,h.del) < 0) | any(c(mm,mf,sf,sm,h,rf,rm,u,h.del) > 1) | any(c(rf,rm) > 0.5))
		stop('The chosen parameter values fall outside of reasonable parameter space')

	try({
		 if(nReps > 1000 & saveTrajectories)
			stop('Warning: You have chosen to save evolutionary trajectories 
				  for a large number of replicate simulations. Thiss will be
				  memory intensive. Consider setting saveTrajectories = FALSE')
	}, silent=FALSE)

	##  Define Fitness Expressions for determining eq. frequencies in absence of inversion
	Wf.init  <-  c(1,           (1 + h*sf),          (1 + h*sf),          (1 + h*sf)^2,        0,
			      (1 + h*sf),   (1 + sf),            (1 + h*sf)^2,        (1 + h*sf)*(1 + sf), 0,
			      (1 + h*sf),   (1 + h*sf)^2,        (1 + sf),            (1 + sf)*(1 + h*sf), 0,
			      (1 + h*sf)^2, (1 + h*sf)*(1 + sf), (1 + sf)*(1 + h*sf), (1 + sf)^2,          0,
			       0,            0,                   0,                   0,                  0)
  	Wm.init  <-  c(1,           (1 + h*sm),          (1 + h*sm),          (1 + h*sm)^2,        0,
			      (1 + h*sm),   (1 + sm),            (1 + h*sm)^2,        (1 + h*sm)*(1 + sm), 0,
			      (1 + h*sm),   (1 + h*sm)^2,        (1 + sm),            (1 + sm)*(1 + h*sm), 0,
			      (1 + h*sm)^2, (1 + h*sm)*(1 + sm), (1 + sm)*(1 + h*sm), (1 + sm)^2,          0,
			       0,            0,                   0,                   0,                  0)

  	# Calculate initial frequencies
	Fii         <-  findEqFreqs(Wf.init, Wm.init, mm = mm, mf = mf, rf = rf, rm = rm, threshold = 1e-6)
	Fii.f.init  <- Fii[1,]
	Fii.m.init  <- Fii[2,]
	
	#Introduce inversion
	Fii         <-  introduceInversion(newMutant = newMutant, Fii.f.init = Fii.f.init, Fii.m.init = Fii.m.init, N = N)
	Fii.f.init  <-  Fii[1,]
	Fii.m.init  <-  Fii[2,]

	# Storage structures for replicate simulation data
	finalInvFreq      <-  rep(0, times=nReps)
	finalE.InvFreq    <-  rep(0, times=nReps)
	finalWbar_f       <-  rep(0, times=nReps)
	finalWbar_m       <-  rep(0, times=nReps)
	finalInvFreq_f    <-  rep(0, times=nReps)
	finalInvFreq_m    <-  rep(0, times=nReps)
	finalE.InvFreq_f  <-  rep(0, times=nReps)
	finalE.InvFreq_m  <-  rep(0, times=nReps)
	nGen              <-  rep(0, times=nReps)
	InvEst            <-  rep(0, times=nReps)
	InvEstTime        <-  rep(0, times=nReps)
	nDels             <-  rep(0, times=nReps)
	
	if(saveTrajectories) {
		replicateTraj    <-  c()
		InvFreqTraj      <-  c()
		E.InvFreqTraj    <-  c()
		Wbar_fTraj       <-  c() 
		Wbar_mTraj       <-  c()
		InvFreq_fTraj    <-  c()
		InvFreq_mTraj    <-  c() 
		E.InvFreq_fTraj  <-  c()
		E.InvFreq_mTraj  <-  c()
	} 

	# Replicate simulation loop
	pb   <-  txtProgressBar(min=0, max=nReps, style=3)
	setTxtProgressBar(pb, 0)
	for(i in 1:nReps) {

	## Sample stationary distribution of deleterious alleles
	delMutFreq  <-  rejectionSampler(n = n, Ne = N, u = u)
	n.del       <-  sum(delMutFreq > runif(n=n))

	# Define fitness expressions, including fitness effects of deleterious mutations
	if(noDel) {
		s.del  <-  0
	}
	Wf  <-  c(1,                                   (1 + h*sf),                                  (1 + h*sf),                                  (1 + h*sf)^2,                       (1 + h*sf)^2*(1 - h.del*s.del)^n.del,
			 (1 + h*sf),                           (1 + sf),                                    (1 + h*sf)^2,                                (1 + h*sf)*(1 + sf),                (1 + h*sf)*(1 + sf)*(1 - h.del*s.del)^n.del,
			 (1 + h*sf),                           (1 + h*sf)^2,                                (1 + sf),                                    (1 + sf)*(1 + h*sf),                (1 + sf)*(1 + h*sf)*(1 - h.del*s.del)^n.del,
			 (1 + h*sf)^2,                         (1 + h*sf)*(1 + sf),                         (1 + sf)*(1 + h*sf),                         (1 + sf)^2,                         (1 + sf)^2*(1 - h.del*s.del)^n.del,
			 (1 + h*sf)^2*(1 - h.del*s.del)^n.del, (1 + h*sf)*(1 + sf)*(1 - h.del*s.del)^n.del, (1 + sf)*(1 + h*sf)*(1 - h.del*s.del)^n.del, (1 + sf)^2*(1 - h.del*s.del)^n.del, (1 + sf)^2*(1 - s.del)^n.del)
	Wm  <-  c(1,                                   (1 + h*sm),                                  (1 + h*sm),                                  (1 + h*sm)^2,                       (1 + h*sm)^2*(1 - h.del*s.del)^n.del,
			 (1 + h*sm),                           (1 + sm),                                    (1 + h*sm)^2,                                (1 + h*sm)*(1 + sm),                (1 + h*sm)*(1 + sm)*(1 - h.del*s.del)^n.del,
			 (1 + h*sm),                           (1 + h*sm)^2,                                (1 + sm),                                    (1 + sm)*(1 + h*sm),                (1 + sm)*(1 + h*sm)*(1 - h.del*s.del)^n.del,
			 (1 + h*sm)^2,                         (1 + h*sm)*(1 + sm),                         (1 + sm)*(1 + h*sm),                         (1 + sm)^2,                         (1 + sm)^2*(1 - h.del*s.del)^n.del,
			 (1 + h*sm)^2*(1 - h.del*s.del)^n.del, (1 + h*sm)*(1 + sm)*(1 - h.del*s.del)^n.del, (1 + sm)*(1 + h*sm)*(1 - h.del*s.del)^n.del, (1 + sm)^2*(1 - h.del*s.del)^n.del, (1 + sm)^2*(1 - s.del)^n.del)
	
	## RUN SIMULATION
	repRes  <-  autoInvFwdSimSexSpec(Fii.f.init = Fii.f.init, Fii.m.init = Fii.m.init, N = N, Wf = Wf, Wm = Wm, mm = mm, mf = mf, rf = rf, rm = rm, saveTrajectories = saveTrajectories)

	# save results for each replicate
	finalInvFreq[i]      <-  repRes$InvFreq[length(repRes$InvFreq)]
  	finalE.InvFreq[i]    <-  repRes$E.InvFreq[length(repRes$E.InvFreq)]
	finalWbar_f[i]  	 <-  repRes$Wbar_f[length(repRes$Wbar_f)]
	finalWbar_m[i]  	 <-  repRes$Wbar_m[length(repRes$Wbar_m)]
	finalInvFreq_f[i]  	 <-  repRes$InvFreq_f[length(repRes$InvFreq_f)]
	finalInvFreq_m[i]  	 <-  repRes$InvFreq_m[length(repRes$InvFreq_m)]
	finalE.InvFreq_f[i]  <-  repRes$E.InvFreq_f[length(repRes$E.InvFreq_f)]
	finalE.InvFreq_m[i]  <-  repRes$E.InvFreq_m[length(repRes$E.InvFreq_m)]
	nGen[i]              <-  repRes$nGen
  	InvEst[i]            <-  repRes$InvEst
  	InvEstTime[i]        <-  repRes$InvEstTime
	nDels[i]             <-  n.del

		if(saveTrajectories) {
				
			replicateTraj    <-  c(replicateTraj,rep(i, times=length(repRes$InvFreq)))
			InvFreqTraj      <-  c(InvFreqTraj,repRes$InvFreq)
			E.InvFreqTraj    <-  c(E.InvFreqTraj,repRes$E.InvFreq)
			Wbar_fTraj       <-  c(Wbar_fTraj,repRes$Wbar_f)
			Wbar_mTraj       <-  c(Wbar_mTraj,repRes$Wbar_m)
			InvFreq_fTraj    <-  c(InvFreq_fTraj,repRes$InvFreq_f)
			InvFreq_mTraj    <-  c(InvFreq_mTraj,repRes$InvFreq_m) 
			E.InvFreq_fTraj  <-  c(E.InvFreq_fTraj,repRes$E.InvFreq_f)
			E.InvFreq_mTraj  <-  c(E.InvFreq_mTraj,repRes$E.InvFreq_m)

			} 	
	setTxtProgressBar(pb, i)
	}

	# Save results and return results as a list
	results.df  <-  data.frame(
							  	"finalInvFreq"     =  finalInvFreq,
								"finalE.InvFreq"   =  finalE.InvFreq,
								"finalWbar_f"      =  finalWbar_f,
								"finalWbar_m"      =  finalWbar_m,
								"finalInvFreq_f"   =  finalInvFreq_f,
								"finalInvFreq_m"   =  finalInvFreq_m,
								"finalE.InvFreq_f" =  finalE.InvFreq_f,
								"finalE.InvFreq_m" =  finalE.InvFreq_m,
								"nGen"             =  nGen,
								"InvEst"           =  InvEst,
								"InvEstTime"       =  InvEstTime,
								"nDels"            =  nDels
							   )
	if(saveTrajectories) {
		traj.df  <-  data.frame(
								"replicateTraj"   =  replicateTraj,
								"InvFreqTraj"     =  InvFreqTraj,
								"E.InvFreqTraj"   =  E.InvFreqTraj,
								"Wbar_fTraj"      =  Wbar_fTraj,
								"Wbar_mTraj"      =  Wbar_mTraj,
								"InvFreq_fTraj"   =  InvFreq_fTraj,
								"InvFreq_mTraj"   =  InvFreq_mTraj,
								"E.InvFreq_fTraj" =  E.InvFreq_fTraj, 
								"E.InvFreq_mTraj" =  E.InvFreq_mTraj
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



#' Generate corner case values to explore equal/female-limited/male-limited cases of a parameter
#' 
#' @title Generate corner case values to explore equal/female-limited/male-limited cases of a parameter
#' @param mu Vector of mean values about which corner cases are to be calucalted
#' @seealso `makeReplicateAutoInvSimsDataSexSpec`
#' @export
#' @author Colin Olito
makeCornerCaseVals  <-  function(mu = c(0.001, 0.002), delta = c(0.001, 0.002)) {
	fVals  <-  c()
	mVals  <-  c()

	for(i in 1:length(mu)) {
		tmp1  <-  c(mu[i], (mu[i] + delta[i]), (mu[i] - delta[i]))
		tmp2  <-  c(mu[i], (mu[i] - delta[i]), (mu[i] + delta[i]))
		fVals  <-  c(fVals,tmp1)
		mVals  <-  c(mVals,tmp2)
	}
	return(rbind(fVals,mVals))
}


#' Wrapper function to run replicate forward simulations for invasion
#' of autosomal inversions in a Wright-Fisher population 
#' 
#' @title Run replicate Wright-Fisher forward simulations for autosomal inversion under different parameter values 
#' @param nReps  Number of replicate simulations to run for each parameter set
#' @param N      Total population size
#' @param h      Dominance coefficient for locally adaptive alleles relative to migrant alleles (h = 0.5)
#' @param m.vals Vector of mean migration rates for locally maladaptive alleles 
#' @param m.deltas Vector of differences between mean and sex-biased migration rates (if m.deltas = NULL, default
#'                 behaviour is to explore equal, female-limited, and male-limited migration)
#' @param s.vals Vector of mean selection coefficients favouring locally adaptive alleles 
#' @param s.deltas Vector of differences between mean and sex-biased selection coefficients (if s.deltas = NULL, default
#'                 behaviour is to explore equal, female-limited, and male-limited expression of locally adaptive alleles)
#' @param r.vals Vector of desired recombination rates among the two loci involved in local adaptation (r = 0.1).
#'                (NOTE: see comments in function for instructions to explore sex-limited recombination).
#' @param n      Number of loci at which deleterious mutations may occur.
#' @param u      Mutation rate (default value of u = 1e-6).
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Defalut value of FALSE assumes that 
#' 				 selection against deleterious mutations is twice as strong as selection favouring
#' 				 the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				 there were no delterious mutations segregating in the population that are linked
#' 				 to the loci involved in local adaptation. 
#' @param fastSim     Logical. Use threshold frequency for establishment of inversion? 
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? Setting this 
#' 							to TRUE can become extremely memory intensive if you are running many
#' 							replicate simulations (see warning).
#' 							otherwise
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#' 					 chosen randomly, given initial genotypic frequencies (Fii.f.init and Fii.m.init). 
#'					 A 2 positions vector ("m"/"f"/"random", "numeric"/"random"). the first position specify wether the inversions come in males, 
#'					 in females, or randomnly. The second select either a numerical position in the haplotype vector for the inversion genotype to be created
#'					 or wether it is random.
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSimSexSpec` runReplicateAutoInvSimsSexSpec
#' @export
#' @author Ludovic Dutoit based on Colin Olito.
makeFastReplicateAutoSexSpecInvSimsData  <-  function(nReps = 1000, N = 20000, h = 1/2, 
													  m.vals = c(0.0005, 0.001), m.deltas = NULL,
													  s.vals = c(0.001, 0.05), s.deltas = NULL, 
													  r.vals = seq(from = 0, to = 0.5, by = 0.025),
													  n = 100, u = 1e-5, h.del = 0, noDel = FALSE, 
													  fastSim = TRUE, newMutant=c("random","random"), 
													  saveTrajectories = FALSE) {

	# create empty data frame with same structure as we are going to need
	data  <-  data.frame(matrix(ncol=9, nrow=0))

	# make parameter values for equal/female-limited/male-limited cases
	ms   <-  makeCornerCaseVals(mu = m.vals, delta = m.vals) # use delta = m.deltas for alternative sex-biased parameterizations
	ss   <-  makeCornerCaseVals(mu = s.vals, delta = s.vals) # use delta = s.deltas for alternative sex-biased parameterizations
	sMu  <-  colSums(ss)
# uncomment to explore sex-limited recombination rates
#	rfs  <-  c(r.vals, r.vals, rep(0, length=length(r.vals)))
#	rms  <-  c(r.vals, rep(0, length=length(r.vals)), r.vals)
	rfs  <-  r.vals
	rms  <-  r.vals
	rs   <-  rbind(rfs,rms)

	# Convenience variables to monitor progress
	prog  <-  0
	tot   <-  ncol(ms)*ncol(ss)*2*length(rfs)
	cat("\n",paste('Running simulations for parameter set ', 1, "/", tot),"\n")

		for (i in 1:ncol(ms)) {
			for (j in 1:ncol(ss)) {
				# Simulate deleterious mutations that are either 
				# 1) neutral
				# 2) lethals 
				# 3) strongly deleterious (twice the selective advantage of locally adaptive alleles)
# uncomment to explore effects of deleterious mutations
#				s.del.vals  <-  c(0, 1, 2*sMu[j]) 
				s.del.vals  <-  0
				
				for(k in 1:length(s.del.vals)) {
					for(l in 1:ncol(rs)) {

						# Display progress in terminal
						prog  <-  prog + 1
						if(prog %% 100 == 0) {
							cat("\n",paste('Running simulations for parameter set ', prog, "/", tot),"\n")
						}

						# Run simulations
						res  <-  runReplicateAutoInvSimsSexSpec(nReps = nReps, N = N, h = h,
																mf = ms[1,i], mm = ms[2,i], 
																sf = ss[1,j], sm = ss[2,j],
																n = n, u = u, h.del = h.del, s.del = s.del.vals[k], 
																rf = rfs[l], rm = rms[l],
																newMutant = newMutant)
						# Append data 
						dat   <-  c(ms[1,i], ms[2,i], ss[1,j], ss[2,j], s.del.vals[k], rfs[l], rms[l], mean(res$results$nDels), (sum(res$results.df$InvEst)/length(res$results.df$InvEst)))
						data  <-  rbind(data, dat)
						rm(dat,res)
					}
				}
			}
		}

	# add variable names to data frame
	colnames(data)  <-  c("mf",
						  "mm",
						  "sf",
						  "sm",
						  "s.dels",
						  "rf",
						  "rm",
						  "av.nDels",
						  "PropEst"
						  )
	
	# create file name
	filename  <-  paste("./output/data/simResults/testFast", "_N", N, "_h", h, "_n", n, "_u", u, "_nReps", nReps, ".csv", sep="")

	# export data as .csv to ./output/data
	write.csv(data, file=filename, row.names = FALSE)

	#  Return results in case user wants it
	return(data)
	
}


########################################################################
## To do:  write slow version of makeReplicateAutoSexSpecInvSimsData()
##         that saves data from all replicate simulations, and uses 
##         fastSim = FALSE. This will enable us to look at equilibrium 
##         frequencies for the inversion, as well as prob of establishment. 
##         Loops will need to be structured to explore more targeted pieces of 
##         parameter space. 
########################################################################
#makeReplicateAutoSexSpecInvSimsData  <-  function(nReps = 1000, N = 20000, h = 1/2, 
#													  m.vals = c(0.0005, 0.001), m.deltas = NULL,
#													  s.vals = c(0.001, 0.05), s.deltas = NULL, 
#													  r.vals = seq(from=0, to=0.5, by=0.025),
#										   			  n = 100, u = 1e-5, h.del = 0, noDel = FALSE, 
#										   			  fastSim = TRUE, newMutant=c("random","random"), saveTrajectories = FALSE) {

	# create empty data frame with same structure as we are going to need
#	data  <-  data.frame(matrix(ncol=13, nrow=0))

	# make parameter values for equal/female-limited/male-limited cases
#	ms          <-  makeCornerCaseVals(mu = m.vals, delta = m.vals) # use delta = m.deltas for alternative sex-biased parameterizations
#	ss          <-  makeCornerCaseVals(mu = s.vals, delta = s.vals) # use delta = s.deltas for alternative sex-biased parameterizations
#	rfs         <-  c(r.vals, r.vals, rep(0, length=length(r.vals)))
#	rms         <-  c(r.vals, rep(0, length=length(r.vals)), r.vals)
#	rs          <-  rbind(rfs,rms)


	# Convenience variables to monitor progress
#	prog  <-  0
#	tot   <-  length(ms)*length(ss)*3*length(rs)

#			for (i in length(ms)) {
#				for (j in length(ss)) {
					
					# Simulate deleterious mutations that are either 
					# 1) neutral
					# 2) lethals 
					# 3) strongly deleterious (twice the selective advantage of locally adaptive alleles)
#					s.del.vals  <-  c(0, 1, 2*ss[j])
				
#					for(k in length(s.del.vals)) {
#						for(l in length(rs)) {

							# Display progress in terminal
#							prog  <-  prog + 1
#							cat("\n",paste('Running simulations for parameter set ', prog, "/", tot),"\n")

							# Run simulations
#							res  <-  runReplicateAutoInvSimsSexSpec(nReps = nReps, N = N, h = h,
#																	mf = ms[1,i], mm = ms[2,i], 
#																	sf = ss[1,j], sm = ss[2,j],
#																	n = n, u = u, h.del = h.del, s.del = s.del.vals[k], noDel = noDel, 
#																	rf = rfs[l], rm = rms[l],
#																	newMutant = newMutant)
#							sum(res$results.df$InvEst)
							# Save data 
#							mms      <-  rep(mm, times=nrow(res$results.df))
#							mfs      <-  rep(mf, times=nrow(res$results.df))
#							s.delss  <-  rep(s.del.vals[k], times=nrow(res$results.df))
#							sfs		 <-  rep(sm, times=nrow(res$results.df))
#							sms  	 <-  rep(sf, times=nrow(res$results.df))
#							rfss     <-  rep(rfs[l], times=nrow(res$results.df))
#							rmss     <-  rep(rfs[l], times=nrow(res$results.df))

							# Append to data frame						
#							df      <-  cbind(res$results.df, N, mms, mfs, s.delss, sfs, sms, rfs, rms)
#							data    <-  rbind(data, df)
#							rm(df)
#						}
#					}
#				}
#			}
#		}
#	}

	# Include constant variables in data frame
#	hs      <-  rep(h, times=nrow(data))
#	us      <-  rep(u, times=nrow(data))
#	h.dels  <-  rep(h.del, times=nrow(data))
#	data    <-  cbind(data,  hs, rs, us, h.dels)
#	colnames(data)  <-  c("finalInvFreq",
#						  "finalE.InvFreq",
#						  "finalWf.mean",
#						  "finalWm.mean",
#						  "finalInvFreq_f",
#						  "finalInvFreq_m",
#						  "finalE.InvFreq_f",
#						  "finalE.InvFreq_m",
#						  "nGen",
#						  "InvEst",
#						  "InvEstTime",
#						  "nDels",
#						  "N",
#						  "mm",
#						  "mf",
#						  "s.dels",
#						  "sf",
#						  "sm",
#						  "rf",
#						  "rm",
#						  "h",
#						  "u",
#						  "h.del"
#						  )
	

	# create file name
#	filename  <-  paste("test",  "_h", h, "_n", n, "_u", u, ".csv", sep="")

	# export data as .csv to ./output/data
#	write.csv(data, file=filename, row.names = FALSE)

	#  Return results in case user wants it
#	return(data)
	
#}
