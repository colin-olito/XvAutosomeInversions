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




#' Forward simulation of genotypic frequencies in a Wright-Fisher
#' population using multinomial sampling with linkage to deleterious
#' mutations
#'
#' @title Wright-Fisher forward simulation of genotypic frequencies (default parameter values in parentheses)
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
#' @export
autoInvWrightFisherSim  <-  function(N = 1000, m = 0.01, s = 0.1, h = 1/2, r = 0.1, 
									 n = 100, u = 1e-6, h.del = 0, noDel = FALSE) {

	##  Preemptive Warnings
	if(any(c(N,m,s,h,r,n,u,h.del) < 0) | any(c(m,s,h,r,u,h.del) > 1) | r > 0.5)
		stop('The chosen parameter values fall outside of reasonable parameter space')

	## Sample stationary distribution of deleterious alleles
	delMutFreq  <-  rejectionSampler(n=n, Ne=N, u=u)
	n.del       <-  sum(delMutFreq > runif(n=n))

	##  Define Fitness Expressions
	if(noDel) {
		s.del  <-  0
	} else {
		s.del  <-  2*s
	}
	W   <-  c(1,                                 (1 + h*s),                               (1 + h*s),                               (1 + h*s)^2,                     (1 + h*s)^2 - (n.del*h.del*s.del),
			  (1 + h*s),                         (1 + s),                                 (1 + h*s)^2,                             (1 + h*s)*(1 + s),               (1 + h*s)*(1 + s) - (n.del*h.del*s.del),
			  (1 + h*s),                         (1 + h*s)^2,                             (1 + s),                                 (1 + s)*(1 + h*s),               (1 + s)*(1 + h*s) - (n.del*h.del*s.del),
			  (1 + h*s)^2,                       (1 + h*s)*(1 + s),                       (1 + s)*(1 + h*s),                       (1 + s)^2,                       (1 + s)^2 - (n.del*h.del*s.del),
			  (1 + h*s)^2 - (n.del*h.del*s.del), (1 + h*s)*(1 + s) - (n.del*h.del*s.del), (1 + s)*(1 + h*s) - (n.del*h.del*s.del), (1 + s)^2 - (n.del*h.del*s.del), (1 + s)^2 - (n.del*s.del))
  	
	## Calculate equilibrium frequencies in absence of inversion  
	# genotypes ordered:
	#	c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1x1, x1x2, x1x3, x1x4, x1x5, 
	#	  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
	#	  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
	#	  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
	#	  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
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

	
	# Run simulation to get equilibrium genotypic fequencies
	delta  <-  rep(1, times=16)
	while(any(delta > 1e-6)) {
		for (j in 1:length(xi)) {
			recFct  <-  get(names(xi)[j])
			xi[j]   <-  round(recFct(Fii = Fii, m = m, r = r), digits=3)
		}
	
	    # offspring genotype frequencies
	    O  <-  offFreq(xi)

		# mean fitness 
		Wbar      <-  sum(O*W)

	    # difference in expected frequency (has simulation reached equilibrium yet?)
		delta   <-  E.Fii[c(1:4,6:9,11:14,16:19)] - (O*W/Wbar)[c(1:4,6:9,11:14,16:19)]
		E.Fii   <-  O*W/Wbar
		Fii     <-  E.Fii
	} 

	# Use deterministic equilibrium frequencies of non-inversion genotypes
	# as initial conditions when introducing the inversion via genotype abba* 
	# at low frequency 
	Fii[19]  <-  Fii[19] - 0.1*N/N
	Fii[20]  <-  0.1*N/N
	
#Fii2  <-  Fii
Fii  <-  Fii2
	# Initiate storage structures
	InvFreq    <-  rep(0, times=(4*N+1))
	E.InvFreq  <-  rep(0, times=(4*N+1))
	W.mean     <-  rep(0, times=(4*N+1))
	rm(xi)
	xi         <-  rep(0, times=5)
	names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')

	# Initial frequencies 
	InvFreq[1]    <-  sum(Fii[c(5,10,15,20:24)]/2, Fii[25])
	E.InvFreq[1]  <-  sum(E.Fii[c(5,10,15,20:24)]/2,Fii[25])

	## Start forward simulation with newly introduced inversion
	gen  <-  1
	while(gen < (4*N) & InvFreq[gen] != 0){

		# Calculate gamete frequencies
		for (j in 1:length(xi)) {
			recFct  <-  get(names(xi)[j])
			xi[j]   <-  round(recFct(Fii = Fii, m = m, r = r), digits=8)
		}
	
	    # Offspring genotype frequencies
	    O  <-  offFreq(xi)

		# Mean fitness 
		Wbar  <-  sum(O*W)

	    # Expected frequencies
	    E.Fii  <-  O*W/Wbar
		
		#Frequencies in adults
		Fii  <-  as.vector(rmultinom(1, N/2, E.Fii)/(N/2))
    
		# Realized frequencies
		InvFreq[gen+1]    <-  sum(Fii[c(5,10,15,20:24)]/2, Fii[25])
		E.InvFreq[gen+1]  <-  sum(E.Fii[c(5,10,15,20:24)]/2, Fii[25])
		W.mean[gen+1]     <-  Wbar
print(gen)
		gen  <-  gen+1
	}
	
	# clear unused storage
	InvFreq    <-  InvFreq[1:gen]
	E.InvFreq  <-  E.InvFreq[1:gen]
	W.mean     <-  W.mean[1:gen+1]
print(c(gen, InvFreq[gen], sum(W.mean[-1])/length(W.mean[-1])))
plot(InvFreq, ylim=c(0,1),type='l', lwd=2)
lines(E.InvFreq, lwd=2, col=2)
}
  	
