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
#rejectionSamplerX  <-  function(n=100, Ne=100, u=1e-6, h=0, sf=0.01, sm=0.01) {
#	
#	# Empty vector for allele frequencies
#	qi  <-  rep(0,times=n)
#	
#	# Rejection sampler loop
#	for(i in 1:n) {
#		accept  <-  FALSE
#		while(accept == FALSE) {
#			x  <-  rbeta(1, shape1=4*Ne*u, shape2=4*Ne*u)
#			U  <-  runif(1)
#			if(U < exp((4/3)*Ne*x*(2*h*sf*(x-1)-x*sf-sm))) {
#				qi[i]   <-  rbinom(1,Ne,x)/Ne
#				accept  <-  TRUE
#			}
#		}
#	}
#	qi
#}


# Whoeveer is going to do the Autosomal model with sex-specific selection
# should start writing their functions here (using the code provided in 
# './R/functions-WFSims-Autosomal.R' as a template)...
#
# just clone into the github repository, then run 
#


#' Linkage Disequilibrium function for W-F recursions
#'
#' @title Linkage Disequilibrium (Dstar)
#' @param Fiix Vector of adult genotypic frequencies in females (of length = 25)
#' @param Fiiy Vector of adult genotypic frequencies in males (of length = 25)
#' @param m   Migration rate
#' @export

Dstar  <-  function(Fii=Fii, m=m, ...) {
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
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @param r   Recombination rate
#' @export

#It takes Fiix or Fiiy as an argument

x.1  <-  function(Fii=Fiix, m=m, r=r) { 
	((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
x.2  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.3  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.4  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
x.5  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}


y.1  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
y.2  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.3  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.4  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
y.5  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}


#' Offspring frequencies after random mating
#'
#' @title Offspring frequencies after random mating, ( males and females have the same allele frequencies)
#' @param xi Vector of haplotype frequencies among female gametes (of length = 5)
#' @param yi Vector of haplotype frequencies among male gametes (of length = 5)

#' @export
offFreq  <-  function(xi,yi) {
	O  <-  c(xi[1] *yi[1], (xi[1]*yi[2]), (xi[1]*yi[3]), (xi[1]*yi[4]), (xi[1]*yi[5]),
			 (xi[2]*yi[1]), xi[2]*yi[2], (xi[2]*yi[3]), (xi[2]*yi[4]), (xi[2]*yi[5]),
			 (xi[3]*yi[1]), (xi[3]*yi[2]), xi[3]*yi[3], (xi[3]*yi[4]), (xi[3]*yi[5]),
			 (xi[4]*yi[1]), (xi[4]*yi[2]), (xi[4]*yi[3]), xi[4]*yi[4], (xi[4]*yi[5]),
			 (xi[5]*yi[1]), (xi[5]*yi[2]), (xi[5]*yi[3]), (xi[5]*yi[4]), xi[5]*yi[4]
			 )
	O
}


#' Calculate deterministic equilibrium genotype frequencies 
#' in the absence of the inversion 
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#.Females
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1x1, x1x2, x1x3, x1x4, x1x5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
#.<ales
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(y1y1, y1y2, y1y3, y1y4, y1y5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  y2y1, y2y2, y2y3, y2y4, y2y5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  y3y1, y3y2, y3y3, y3y4, y3y5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  y4y1, y4y2, y4y3, y4y4, y4y5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  y5y1, y5y2, y5y3, y5y4, y5y5)
#' @title Find the deterministic equilibeium genotype frequencies prior to introducing the inversion
#' @param Wf     Vector of fitness expressions for all 25 genotypes in females
#' @param Wm     Vector of fitness expressions for all 25 genotypes in males
#' @param m      Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param threshold The threshold change in genotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. 
#' @export
#' @seealso `offFreq`, `autoInvWrightFisherSim`
#' @author Ludovic Dutoit based on Colin Olito




findEqFreqs  <-  function(Wf,Wm, mm, mf r, threshold = 1e-6) {
	
	Fii.init  <-  c(1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					1/16, 1/16, 1/16, 1/16, 0, 
					   0,    0,    0,    0, 0)
	Fiix    <-  Fii.init
	E.Fiix  <-  Fii.init

	Fiiy    <-  Fii.init
	E.Fiiy  <-  Fii.init

	# Storage for female gamete frequencies
	xi         <-  rep(0, times=5)
	names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')

	# Storage for female gamete frequencies
	yi         <-  rep(0, times=5)
	names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')
	
	xdelta  <-  rep(1, times=16) # change of allele frequencies in 1gen in males
	ydelta  <-  rep(1, times=16)#  change of allele frequencies in 1gen in females
	delta<-append(xdelta,ydelta) # just to be able to check if ANY allele is not stable
	
	gen=0 # counter of generations
	while(any(delta > threshold)) {
		gen=gen+1
		for (j in 1:length(xi)) {
			recFct  <-  get(names(xi)[j])
			xi[j]   <-  round(recFct(Fii = E.Fiix, m = mf, r = r), digits=3)
			recFct  <-  get(names(yi)[j])
			yi[j]   <-  round(recFct(Fii = E.Fiiy, m = mm, r = r), digits=3)
		}
	    # offspring genotype frequencies
	    O  <-  offFreq(xi,yi)

		# mean fitness 
		Wfbar      <-  sum(O*Wf)
		Wmbar      <-  sum(O*Wm)

	    # difference in expected frequency (has simulation reached equilibrium yet?)
		xdelta   <-  E.Fiix[c(1:4,6:9,11:14,16:19)] - (O*Wf/Wfbar)[c(1:4,6:9,11:14,16:19)] # not the inversion freq which will not move
		ydelta   <-  E.Fiiy[c(1:4,6:9,11:14,16:19)] - (O*Wm/Wmbar)[c(1:4,6:9,11:14,16:19)]
		delta<-append(xdelta,ydelta)

		E.Fiix   <-  O*Wf/Wfbar
		E.Fiiy   <-  O*Wm/Wmbar

	}
	print (c("reach equilibrium of initial frequencies at gen",gen))

	names(E.Fiix)  <-  NULL
	names(E.Fiiy)  <-  NULL
	return (rbind(E.Fiix,E.Fiiy) )
}


##Test findEqFreqs 
#Wf<-rep(1,25)
#Wm<-rep(1,25)
#m=0
#r=1
#findEqFreqs(Wf,Wm, m, r, threshold = 1e-6)#

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
#' @param Fiix.init  Vector of initial frequencies (deterministic eq. frequencies in absence of inversion) in females
#' @param Fiiy.init  Vector of initial frequencies (deterministic eq. frequencies in absence of inversion) in males
#' @param N         Population size
#' @param Wf        Vector of fitness expressions for all 25 genotypes in females
#' @param Wm        Vector of fitness expressions for all 25 genotypes in males
#' @param mm         Migration rate for locally maladaptive alleles (m =  0.01) in males
#' @param mf         Migration rate for locally maladaptive alleles (m =  0.01) in females
#' @param r         Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? 
#' @export
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Ludovic Dutoit, based on Colin Olito

autoInvFwdSimSexSpec  <-  function(Fiix.init = Fiix.init  , Fiiy.init = Fiiy.init,   N = N, Wf = Wf, Wm = Wm, mm = sm,, mf=sf, r = r, 
							saveTrajectories = FALSE, ...) {

	# Use deterministic eq. initial frequencies
	Fiix  <-  Fiix.init 
	Fiiy  <-  Fiiy.init


	# Define threshold frequency for establishment of inversion
	pcrit  <-  2/(N*m)

	# Storage for gamete frequencies  
	xi         <-  rep(0, times=5)
	yi         <-  rep(0, times=5)
	names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')
	names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')


	if(saveTrajectories) {
		# Storage structures for individual simulation data
		InvFreq    <-  rep(0, times=(4*N+1))
		E.InvFreq  <-  rep(0, times=(4*N+1))
		W.mean     <-  rep(0, times=(4*N+1))
		Wf.mean <-  rep(0, times=(4*N+1)) 
		Wm.mean <-  rep(0, times=(4*N+1))
		InvFreq_f <-  rep(0, times=(4*N+1))
		InvFreq_m <-  rep(0, times=(4*N+1)) 
		E.InvFreq_f <-  rep(0, times=(4*N+1))
		E.InvFreq_m <-  rep(0, times=(4*N+1))	
				# Initial inversion frequency 
		InvFreq[1]    <-  sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25]) +  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])
		E.InvFreq[1]  <-  InvFreq[1]
		
		## Start forward simulation with newly introduced inversion
		gen  <-  1
#		while(InvFreq > 0 & InvFreq <= pcrit) {
		while(gen < (4*N) & InvFreq[gen] > 0 ) {
		
			## Step through recursions:
			# 1) Calculate gamete frequencies
			for (j in 1:length(xi)) {
				recFctx  <-  get(names(xi)[j])
				recFcty  <-  get(names(yi)[j])
				xi[j]   <-  round(recFctx(Fii = Fiix, m = mf, r = r), digits=8)
				yi[j]   <-  round(recFcty(Fii = Fiiy, m = mm, r = r), digits=8)
			}
			# 2) Offspring genotype frequencies
			O      <-  offFreq(xi,yi) # no need of two vectors, we just apply twice selection 
			# 3) Mean fitness 
			Wbar.f  <-  sum(O*Wf)
			Wbar.m  <-  sum(O*Wm)
			Wbar   <-  0.5 * (sum(Wbar.f)+sum(Wbar.m)) #Wbar is the average fitness across sexes, CHECK
			# 4) Expected frequencies
			E.Fiix  <-  O*Wf/Wbar.f # sex specific expected freq
			E.Fiiy  <-  O*Wm/Wbar.m # sex specific expected freq
			# 5) Draw random frequencies in adults
			Fiix    <-  as.vector(rmultinom(1, N/2, E.Fiix)/(N/2)) #IS THAT RIGHT?? or does the draw change
			Fiiy    <-  as.vector(rmultinom(1, N/2, E.Fiiy)/(N/2)) #IS THAT RIGHT??
			

			# Realized frequencies
			InvFreq[gen+1]     <-  0.5 * (sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25]) +  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			E.InvFreq[gen+1]   <-  0.5 * (sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25]) +  sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			W.mean[gen+1]     <-  Wbar
			
			#The variable below are stored specifically for the SexSpecific simulation 
			Wf.mean[gen+1]   =  sum(O*Wf)
			Wm.mean[gen+1]   =  sum(O*Wm)
			InvFreq_f[gen+1]  =  sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25])
			InvFreq_m[gen+1]  =  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])
			E.InvFreq_f[gen+1]=  sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25])
			E.InvFreq_m[gen+1]=  sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])

			gen  <-  gen+1
			print (gen)
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
				  	"InvFreq"     =  InvFreq[1:gen-1],
				  	"E.InvFreq"   =  E.InvFreq[1:gen-1],
				  	"W.mean"      =  W.mean[1:gen-1],
					"Wf.mean" 	  =  Wf.mean[1:gen-1], 
					"Wm.mean" 	  =	 Wm.mean[1:gen-1],
					"InvFreq_f" 	  =  InvFreq_f[1:gen-1],
					"InvFreq_m" 	  =  InvFreq_m[1:gen-1], 
					"E.InvFreq_f"  =  E.InvFreq_f[1:gen-1],
					"E.InvFreq_m"  =  E.InvFreq_m[1:gen-1],
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
		InvFreq    <-  0.5 * (sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25]) +  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
		E.InvFreq  <-  InvFreq[1]
		
		## Start forward simulation with newly introduced inversion
		gen  <-  1
#		while(InvFreq > 0 & InvFreq <= pcrit) {
		while(gen < (4*N) & InvFreq > 0 ) {
		
			## Step through recursions:
			# 1) Calculate gamete frequencies
			for (j in 1:length(xi)) {
				recFctx  <-  get(names(xi)[j])
				recFcty  <-  get(names(yi)[j])
				xi[j]   <-  round(recFctx(Fii = Fiix, m = m, r = r), digits=8)
				yi[j]   <-  round(recFcty(Fii = Fiiy, m = m, r = r), digits=8)
			}
			# 2) Offspring genotype frequencies
			O      <-  offFreq(xi,yi) # no need of two vectors, we just apply twice selection 
			# 3) Mean fitness 
			Wbar.f  <-  sum(O*Wf)
			Wbar.m  <-  sum(O*Wm)	
			Wbar   <-  0.5 * (sum(O*Wbar.f)+sum(O*Wbar.m))		
			# 4) Expected frequencies
			E.Fiix  <-  O*Wf/Wbar.f # sex specific expected freq
			E.Fiiy  <-  O*Wm/Wbar.m # sex specific expected freq

			# 5) Draw random frequencies in adults
			Fiix    <-  as.vector(rmultinom(1, N/2, E.Fiix)/(N/2)) #IS THAT RIGHT?? or does the draw change
			Fiiy    <-  as.vector(rmultinom(1, N/2, E.Fiiy)/(N/2)) #IS THAT RIGHT??
			
			# Realized frequencies
			InvFreq    <-  0.5 * (sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25]) +  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			E.InvFreq  <-  0.5 * (sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25]) +  sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			W.mean     <-  Wbar
			
			#The variable below are stored specifically for the SexSpecific simulation 
			Wf.mean =  sum(O*Wf)
			Wm.mean =  sum(O*Wm)
			InvFreq_f = sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25])
			InvFreq_m  = sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])
			E.InvFreq_f = sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25])
			E.InvFreq_m = sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])
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
		res  <-  list(
				  	"InvFreq"     =  InvFreq,
				  	"E.InvFreq"   =  E.InvFreq,
				  	"W.mean"      =  sum(W.mean)/length(W.mean),
					"Wf.mean" 	  =  sum(Wf.mean)/length(Wf.mean), 
					"Wm.mean" 	  =	 sum(Wm.mean)/length(Wm.mean),
					"InvFreq_f" 	  =  InvFreq_f,
					"InvFreq_m" 	  =  InvFreq_m, 
					"E.InvFreq_f"  =  E.InvFreq_f,
					"E.InvFreq_m"  =  E.InvFreq_m,
					"nGen"        =  gen,
				  	"InvEst"      =  invEst,
				  	"InvEstTime"  =  invEstTime
			 )
		}

		 	 
return(res)

}

#' Wrapper function to run replicate forward simulations for invasion
#' of autosomal inversions in a Wright-Fisher population with sex specific selection 
#'
#' @title Wright-Fisher forward simulation of genotypic frequencies (default parameter values in parentheses)
#' @param nReps  Numer of replicate simulations. With no deleterious mutations, and introducing 
#' 				 a single copy of the inversion, it takes 1,600,000 replicate simulations to 
#'				 get 10,000 where the inversion successfully establishes.
#' @param N      Effective population size
#' @param mf      Migration rate for locally maladaptive alleles (m =  0.01) in females
#' @param mm     Migration rate for locally maladaptive alleles (m =  0.01) in males
#' @param sf     Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02) in females
#' @param sm     Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02) in males
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
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSimSexSpec`
#' @export
#' @author Ludovic Dutoit based on Colin Olito.
runReplicateAutoInvSimsSexSpec  <-  function(nReps = 1000, N = 500, mm = 0.01, mf=0.01, sf = 0.1, sm=0.1, h = 1/2, r = 0.1, 
									  n = 100, u = 1e-5, h.del = 0, sf.del = 1, sm.del=1, noDel = FALSE,
									  saveTrajectories = FALSE) {

	##  Preemptive Warnings
	if(any(c(N,m,sf,sm,h,r,n,u,h.del) < 0) | any(c(m,sf,sm,h,r,u,h.del) > 1) | r > 0.5)
		stop('The chosen parameter values fall outside of reasonable parameter space')

	try({ 
		 if(m >= sf )
			stop('Warning: migration is stronger than than selection, 
				  adaptive alleles will be swamped by maladaptive migrants')
	}, silent=FALSE)

	try({ 
		 if(m >= sm )
			stop('Warning: migration is stronger than than selection, 
				  adaptive alleles will be swamped by maladaptive migrants')
	}, silent=FALSE)

	try({
		 if(nReps > 1000 & saveTrajectories)
			stop('Warning: You have chosen to save evolutionary trajectories 
				  for a large number of replicate simulations. Thiss will be
				  memory intensive. Consider setting saveTrajectories = FALSE')
	}, silent=FALSE)

	##  Define Fitness Expressions for males and females
	Wf.init  <-  c(1,          (1 + h*sf),         (1 + h*sf),         (1 + h*sf)^2,        0,
			     (1 + h*sf),   (1 + sf),           (1 + h*sf)^2,       (1 + h*sf)*(1 + sf),  0,
			     (1 + h*sf),   (1 + h*sf)^2,       (1 + sf),           (1 + sf)*(1 + h*sf),  0,
			     (1 + h*sf)^2, (1 + h*sf)*(1 + sf), (1 + sf)*(1 + h*sf), (1 + sf)^2,          0,
			     0,           0,                 0,                 0,                  0)
  	Wm.init  <-  c(1,          (1 + h*sm),         (1 + h*sm),         (1 + h*sm)^2,        0,
			     (1 + h*sm),   (1 + sm),           (1 + h*sm)^2,       (1 + h*sm)*(1 + sm),  0,
			     (1 + h*sm),   (1 + h*sm)^2,       (1 + sm),           (1 + sm)*(1 + h*sm),  0,
			     (1 + h*sm)^2, (1 + h*sm)*(1 + sm), (1 + sm)*(1 + h*sm), (1 + sm)^2,          0,
			     0,           0,                 0,                 0,                  0)
	##  Define Fitness Expressions for determining eq. frequencies in absence of inversion
	Fii<-findEqFreqs(Wf.init ,Wm.init, m, r, threshold = 1e-6)
	Fiix.init <- Fii[1,]
	Fiiy.init <- Fii[2,]
	# Use deterministic equilibrium frequencies of non-inversion genotypes
	# as initial conditions when introducing the inversion via a single
	# copy of the abba* genotype , just introduce in females?
	Fiix.init[19]  <-  Fiix.init[19] - 1/N #(IS IT HERE 1/(N/2))
	Fiix.init[20]  <-  1/N
	# Storage structures for replicate simulation data
	finalInvFreq     <-  rep(0, times=nReps)
  	finalE.InvFreq   <-  rep(0, times=nReps)
  	finalW.mean      <-  rep(0, times=nReps)
	finalWf.mean 	  <-  rep(0, times=nReps)
	finalWm.mean 	  <-  rep(0, times=nReps)
	finalInvFreq_f 	  <-  rep(0, times=nReps)
	finalInvFreq_m 	  <-  rep(0, times=nReps)
	finalE.InvFreq_f  <-  rep(0, times=nReps)
	finalE.InvFreq_m  <-  rep(0, times=nReps)
	nGen        <-  rep(0, times=nReps)
  	InvEst      <-  rep(0, times=nReps)
  	InvEstTime  <-  rep(0, times=nReps)
	nDels <-  rep(0, times=nReps)
	
	if(saveTrajectories) {
		replicateTraj  <-  c()
		InvFreqTraj    <-  c()
		E.InvFreqTraj  <-  c()
		W.meanTraj     <-  c()
		Wf.meanTraj   <-  c() 
		Wm.meanTraj   <-  c()
		InvFreq_fTraj   <-  c()
		InvFreq_mTraj   <-  c() 
		E.InvFreq_fTraj   <-  c()
		E.InvFreq_mTraj   <-  c()
	} 

	# Replicate simulation loop
	#print('Running Wright-Fisher Forward Simulations')
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
	Wf  <-  c(1,                                  (1 + h*sf),                                 (1 + h*sf),                                 (1 + h*sf)^2,                       (1 + h*sf)^2*(1 - h.del*sf.del)^n.del,
			(1 + h*sf),                           (1 + sf),                                   (1 + h*sf)^2,                               (1 + h*sf)*(1 + sf),                 (1 + h*sf)*(1 + sf)*(1 - h.del*sf.del)^n.del,
			(1 + h*sf),                           (1 + h*sf)^2,                               (1 + sf),                                   (1 + sf)*(1 + h*sf),                 (1 + sf)*(1 + h*sf)*(1 - h.del*sf.del)^n.del,
			(1 + h*sf)^2,                         (1 + h*sf)*(1 + sf),                         (1 + sf)*(1 + h*sf),                         (1 + sf)^2,                         (1 + sf)^2*(1 - h.del*sf.del)^n.del,
			(1 + h*sf)^2*(1 - h.del*sf.del)^n.del, (1 + h*sf)*(1 + sf)*(1 - h.del*sf.del)^n.del, (1 + sf)*(1 + h*sf)*(1 - h.del*sf.del)^n.del, (1 + sf)^2*(1 - h.del*sf.del)^n.del, (1 + sf)^2*(1 - sf.del)^n.del)
	Wm  <-  c(1,                                  (1 + h*sm),                                 (1 + h*sm),                                 (1 + h*sm)^2,                       (1 + h*sm)^2*(1 - h.del*sm.del)^n.del,
				(1 + h*sm),                           (1 + sm),                                   (1 + h*sm)^2,                               (1 + h*sm)*(1 + sm),                 (1 + h*sm)*(1 + sm)*(1 - h.del*sm.del)^n.del,
				(1 + h*sm),                           (1 + h*sm)^2,                               (1 + sm),                                   (1 + sm)*(1 + h*sm),                 (1 + sm)*(1 + h*sm)*(1 - h.del*sm.del)^n.del,
				(1 + h*sm)^2,                         (1 + h*sm)*(1 + sm),                         (1 + sm)*(1 + h*sm),                         (1 + sm)^2,                         (1 + sm)^2*(1 - h.del*sm.del)^n.del,
				(1 + h*sm)^2*(1 - h.del*sm.del)^n.del, (1 + h*sm)*(1 + sm)*(1 - h.del*sm.del)^n.del, (1 + sm)*(1 + h*sm)*(1 - h.del*sm.del)^n.del, (1 + sm)^2*(1 - h.del*sm.del)^n.del, (1 + sm)^2*(1 - sm.del)^n.del)
	## RUN SIMULATION
	repRes<-autoInvFwdSimSexSpec(Fiix.init =Fiix.init  , Fiiy.init = Fiiy.init,   N = N, Wf = Wf, Wm = Wm, mm = mm, mf =mf, r = r,saveTrajectories=saveTrajectories)
	# save results for each replicate
	finalInvFreq[i]      <-  repRes$InvFreq[length(repRes$InvFreq)]
  	finalE.InvFreq[i]    <-  repRes$E.InvFreq[length(repRes$E.InvFreq)]
  	finalW.mean[i]       <-  repRes$W.mean[length(repRes$W.mean)]
	finalWf.mean[i]  	  <-  repRes$Wf.mean[length(repRes$Wf.mean)]
	finalWm.mean[i]  	  <-  repRes$Wm.mean[length(repRes$Wm.mean)]
	finalInvFreq_f[i]  	  <-  repRes$InvFreq_f[length(repRes$InvFreq_f)]
	finalInvFreq_m[i]  	  <-  repRes$InvFreq_m[length(repRes$InvFreq_m)]
	finalE.InvFreq_f[i]   <-  repRes$E.InvFreq_f[length(repRes$E.InvFreq_f)]
	finalE.InvFreq_m[i]   <-  repRes$E.InvFreq_m[length(repRes$E.InvFreq_m)]
	nGen[i]         <-  repRes$nGen
  	InvEst[i]       <-  repRes$InvEst
  	InvEstTime[i]   <-  repRes$InvEstTime
	nDels[i]           <-  n.del

		if(saveTrajectories) {
				
			replicateTraj  <-  c(replicateTraj,rep(i, times=length(repRes$InvFreq)))
			InvFreqTraj    <-  c(InvFreqTraj,repRes$InvFreq)
			E.InvFreqTraj  <-  c(E.InvFreqTraj,repRes$E.InvFreq)
			W.meanTraj     <-  c(W.meanTraj,repRes$W.mean)
			Wf.meanTraj   <-  c(Wf.meanTraj,repRes$Wf.mean)
			Wm.meanTraj   <-  c(Wm.meanTraj,repRes$Wm.mean)
			InvFreq_fTraj   <-  c(InvFreq_fTraj,repRes$InvFreq_f)
			InvFreq_mTraj   <-  c(InvFreq_mTraj,repRes$InvFreq_m) 
			E.InvFreq_fTraj   <-  c(E.InvFreq_fTraj,repRes$E.InvFreq_f)
			E.InvFreq_mTraj   <-  c(E.InvFreq_mTraj,repRes$E.InvFreq_m)

			} 	
	setTxtProgressBar(pb, i)
	}

	# Save results and return results as a list
	results.df  <-  data.frame(
							  	"finalInvFreq" = finalInvFreq,
								"finalE.InvFreq" = 	finalE.InvFreq,
								"finalW.mean" = finalW.mean,
								"finalWf.mean" = finalWf.mean,
								"finalWm.mean" = finalWm.mean,
								"finalInvFreq_f" = finalInvFreq_f,
								"finalInvFreq_m" = finalInvFreq_m,
								"finalE.InvFreq_f" = finalE.InvFreq_f,
								"finalE.InvFreq_m" = finalE.InvFreq_m,
								"nGen" =  nGen,
								"InvEst" = InvEst,
								"InvEstTime" = InvEstTime,
								"nDels" = nDels
							   )
	if(saveTrajectories) {
		traj.df  <-  data.frame(
								"replicateTraj" = replicateTraj,
								"InvFreqTraj" = InvFreqTraj,
								"E.InvFreqTraj" = E.InvFreqTraj,
								"W.meanTraj" = W.meanTraj,
								"Wf.meanTraj" = Wf.meanTraj,
								"Wm.meanTraj" = Wm.meanTraj,
								"InvFreq_fTraj" = InvFreq_fTraj,
								"InvFreq_mTraj" = InvFreq_mTraj,
								"E.InvFreq_fTraj" = E.InvFreq_fTraj, 
								"E.InvFreq_mTraj" = E.InvFreq_mTraj
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
#' @param N.vals Vector of desired population sizes
#' @param m.vals Vector ofdesired migration rates for locally maladaptive alleles (m =  0.01)
#' @param sf.vals      Desired selectives advantage of locally adaptive alleles over migrant alleles in females (s = 0.02)
#' @param sm.vals      Desired selectives advantage of locally adaptive alleles over migrant alleles in males (s = 0.02)
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
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSimSexSpec` runReplicateAutoInvSimsSexSpec
#' @export
#' @author Ludovic Dutoit based on Colin Olito.
makeReplicateAutoInvSimsDataSexSpec  <-  function(nReps = 1000, N.vals = c(500, 1000), mf.vals = c(0.01, 0.05), mm.vals,
										   sf.vals = 0.1, sm.vals=0.1,h = 1/2, r = 0.1, 
										   n = 100, u = 1e-5, h.del = 0, noDel = FALSE,saveTrajectories = FALSE) {



	# create empty data frame with same structure as we are going to need
	data  <-  data.frame(matrix(ncol=14, nrow=0))

	# Convenience variables to monitor progress
	prog  <-  0
	s.del.vals <- rep(NA,3) # now defined for every combination than sf and sm
	tot   <-  length(N.vals)*length(m.vals)*length(s.del.vals)*length(sf.vals)*length(sm.vals)

	# Loop over parameter values we want to explore 
	CHANGE THAT FOR R AND FOR mm mf
	(Think of having in pair, just here)
	for(j in 1:length(N.vals)) {
		for(k in 1:length(m.vals)) {
				for (m in 1:length(sf.vals))   { 
					for (n in 1:length(sm.vals)) {

						# Simulate deleterious mutations that are either 
						# 1) recessive lethals OR
						# 2) recessive experiencing purifying selection
						#    that is twice as strong as the selective 
						#    advantage of the locally adaptive allels  


						s.del.vals = c(0, 1, 2*(0.5*(sf.vals[m]+sm.vals[n])))#Not sure Colin 
						for(l in 1:length(s.del.vals)) {

							# Display progress in terminal
							prog  <-  prog + 1
							cat("\n",paste('Running simulations for parameter set ', prog, "/", tot),"\n")

							# Run simulations  
							res  <-  runReplicateAutoInvSimsSexSpec(nReps = nReps, N = N.vals[j], mm = m.vals[k],mf =  m.vals[k] , sf = sf.vals[m], sm=sm.vals[n],  h = h, r = r, 
															 n = n, u = u, h.del = h.del, sf.del = s.del.vals[l], sm.del=s.del.vals[l], 
															 noDel = noDel, saveTrajectories = FALSE)

							# Save data 
							Ns      <-  rep(N.vals[j], times=nrow(res$results.df))
							ms      <-  rep(m.vals[k], times=nrow(res$results.df))
							s.dels  <-  rep(s.del.vals[l], times=nrow(res$results.df))
							sf	<-      rep(sf.vals[m], times=nrow(res$results.df))
							sm 		<-  rep(sf.vals[n], times=nrow(res$results.df))
							# Append to data frame
							df      <-  cbind(res$results.df, Ns, ms, s.dels,sf,sm)
							data    <-  rbind(data, df)
							rm(df)
					}
				}
			}
		}
	}

	# Include constant variables in data frame
	hs      <-  rep(h, times=nrow(data))
	rs      <-  rep(r, times=nrow(data))
	us      <-  rep(u, times=nrow(data))
	h.dels  <-  rep(h.del, times=nrow(data))
	data    <-  cbind(data,  hs, rs, us, h.dels)
	colnames(data)  <-  c("finalInvFreq" ,"finalE.InvFreq" ,"finalW.mean" ,"finalWf.mean" ,"finalWm.mean" ,"finalInvFreq_f" ,"finalInvFreq_m" ,"finalE.InvFreq_f" ,"finalE.InvFreq_m" ,"nGen" ,"InvEst" ,"InvEstTime" ,"nDels","N","m", "s.dels","sf","sm","h","r","u","h.del")
	

	# create file name
	filename  <-  paste("test",  "_h", h, "_r", r, "_n", n, "_u", u, ".csv", sep="")

	# export data as .csv to ./output/data
	write.csv(data, file=filename, row.names = FALSE)

	#  Return results in case user wants it
	return(data)
	
}
 
###TEST the whole thing
#source("R/functions-WFSims-Autosomal-SexSpec.R")# from repos home folder
#a<-makeReplicateAutoInvSimsDataSexSpec(nReps = 2, N.vals = c(500, 1000), m.vals = c(0.01, 0.05), sf.vals = c(0.1), sm.vals=c(0.1),h = 1/2, r = 0.1,  n = 100, u = 1e-5, h.del = 0)

