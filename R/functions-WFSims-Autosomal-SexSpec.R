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
#' @param Wx     Vector of fitness expressions for all 25 genotypes in females
#' @param Wy     Vector of fitness expressions for all 25 genotypes in males
#' @param m      Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param threshold The threshold change in genotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. 
#' @export
#' @seealso `offFreq`, `autoInvWrightFisherSim`
#' @author Colin Olito




findEqFreqs  <-  function(Wx,Wy, m, r, threshold = 1e-6) {
	
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
			xi[j]   <-  round(recFct(Fii = E.Fiix, m = m, r = r), digits=3)
			recFct  <-  get(names(yi)[j])
			yi[j]   <-  round(recFct(Fii = E.Fiiy, m = m, r = r), digits=3)
		}
	    # offspring genotype frequencies
	    O  <-  offFreq(xi,yi)

		# mean fitness 
		Wxbar      <-  sum(O*Wx)
		Wybar      <-  sum(O*Wy)

	    # difference in expected frequency (has simulation reached equilibrium yet?)
		xdelta   <-  E.Fiix[c(1:4,6:9,11:14,16:19)] - (O*Wx/Wxbar)[c(1:4,6:9,11:14,16:19)] # not the inversion freq which will not move
		ydelta   <-  E.Fiiy[c(1:4,6:9,11:14,16:19)] - (O*Wy/Wybar)[c(1:4,6:9,11:14,16:19)]
		delta<-append(xdelta,ydelta)

		E.Fiix   <-  O*Wx/Wxbar
		E.Fiiy   <-  O*Wy/Wybar

	}
	print (c("reach equilibrium of initial frequencies at gen",gen))

	names(E.Fiix)  <-  NULL
	names(E.Fiiy)  <-  NULL
	return (rbind(E.Fiix,E.Fiiy) )
}


##Test findEqFreqs 
#Wx<-rep(1,25)
#Wy<-rep(1,25)
#m=0
#r=1
#findEqFreqs(Wx,Wy, m, r, threshold = 1e-6)#

#ADD HEADER HERE
autoInvFwdSimSexSpec  <-  function(Fiix.init = Fiix.init  , Fiiy.init = Fiiy.init,   N = N, Wx = Wx, Wy = Wy, m = m, r = r, 
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
		Wx.mean <-  rep(0, times=(4*N+1)) 
		Wy.mean <-  rep(0, times=(4*N+1))
		InvFreqx <-  rep(0, times=(4*N+1))
		InvFreqy <-  rep(0, times=(4*N+1)) 
		E.InvFreqx <-  rep(0, times=(4*N+1))
		E.InvFreqy <-  rep(0, times=(4*N+1))	
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
				recFctx  <-  get(names(xi)[j])
				recFcty  <-  get(names(yi)[j])
				xi[j]   <-  round(recFctx(Fii = Fiix, m = m, r = r), digits=8)
				yi[j]   <-  round(recFcty(Fii = Fiiy, m = m, r = r), digits=8)
			}
			# 2) Offspring genotype frequencies
			O      <-  offFreq(xi,yi) # no need of two vectors, we just apply twice selection 
			# 3) Mean fitness 
			Wbar   <-  0.5 * (sum(O*Wx)+sum(0*Wy)) #Wbar is the average fitness across sexes, CHECK
			# 4) Expected frequencies
			E.Fiix  <-  O*Wx/Wbar # sex specific expected freq
			E.Fiiy  <-  O*Wy/Wbar # sex specific expected freq
			# 5) Draw random frequencies in adults
			Fiix    <-  as.vector(rmultinom(1, N/2, E.Fiix)/(N/2)) #IS THAT RIGHT?? or does the draw change
			Fiiy    <-  as.vector(rmultinom(1, N/2, E.Fiiy)/(N/2)) #IS THAT RIGHT??
			

			# Realized frequencies
			InvFreq    <-  0.5 * (sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25]) +  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			E.InvFreq  <-  0.5 * (sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25]) +  sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			W.mean     <-  Wbar
			
			#The variable below are stored specifically for the SexSpecific simulation 
			Wx.mean[gen+1]   =  sum(O*Wx)
			Wy.mean[gen+1]   =  sum(O*Wy)
			InvFreqx[gen+1]  =  sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25])
			InvFreqy[gen+1]  =  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])
			E.InvFreqx[gen+1]=  sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25])
			E.InvFreqy[gen+1]=  sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])

			gen  <-  gen+1
		}
			
		# Has the inversion reached threshold frequency for establishment (pcrit)? 
		# When did it first reach pcrit?
		if(any(InvFreq >= pcrit)) {
			invEst      <-  1
			invEstTime  <-  gen[invFreq >= pcrit][1]
		} else {	
			invEst      <-  0
			invEstTime  <-  NA
		}
	
		# Save  simulation data
		res  <-  list(
				  	"InvFreq"     =  InvFreq[1:gen-1],
				  	"E.InvFreq"   =  E.InvFreq[1:gen-1],
				  	"W.mean"      =  W.mean[1:gen-1],
					"Wx.mean" 	  =  Wx.mean[1:gen-1], 
					"Wy.mean" 	  =	 Wy.mean[1:gen-1],
					"InvFreqx" 	  =  InvFreqx[1:gen-1],
					"InvFreqy" 	  =  InvFreqy[1:gen-1], 
					"E.InvFreqx"  =  E.InvFreqx[1:gen-1],
					"E.InvFreqy"  =  E.InvFreqy[1:gen-1],
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
				recFctx  <-  get(names(xi)[j])
				recFcty  <-  get(names(yi)[j])
				xi[j]   <-  round(recFctx(Fii = Fiix, m = m, r = r), digits=8)
				yi[j]   <-  round(recFcty(Fii = Fiiy, m = m, r = r), digits=8)
			}
			# 2) Offspring genotype frequencies
			O      <-  offFreq(xi,yi) # no need of two vectors, we just apply twice selection 
			# 3) Mean fitness 
			Wbar   <-  0.5 * (sum(O*Wx)+sum(0*Wy)) #Wbar is the average fitness across sexes, CHECK
			# 4) Expected frequencies
			E.Fiix  <-  O*Wx/Wbar # sex specific expected freq
			E.Fiiy  <-  O*Wy/Wbar # sex specific expected freq
			# 5) Draw random frequencies in adults
			Fiix    <-  as.vector(rmultinom(1, N/2, E.Fiix)/(N/2)) #IS THAT RIGHT?? or does the draw change
			Fiiy    <-  as.vector(rmultinom(1, N/2, E.Fiiy)/(N/2)) #IS THAT RIGHT??
			
			# Realized frequencies
			InvFreq    <-  0.5 * (sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25]) +  sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			E.InvFreq  <-  0.5 * (sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25]) +  sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25]))
			W.mean     <-  Wbar
			
			#The variable below are stored specifically for the SexSpecific simulation 
			Wx.mean =  sum(O*Wx)
			Wy.mean =  sum(O*Wy)
			InvFreqx = sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25])
			InvFreqy  = sum(Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])
			E.InvFreqx = sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25])
			E.InvFreqy = sum(E.Fiiy[c(5,10,15,20:24)]/2, Fiiy[25])
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
					"Wx.mean" 	  =  sum(Wx.mean)/length(Wx.mean), 
					"Wy.mean" 	  =	 sum(Wy.mean)/length(Wy.mean),
					"InvFreqx" 	  =  InvFreqx,
					"InvFreqy" 	  =  InvFreqy, 
					"E.InvFreqx"  =  E.InvFreqx,
					"E.InvFreqy"  =  E.InvFreqy,
					"nGen"        =  gen,
				  	"InvEst"      =  invEst,
				  	"InvEstTime"  =  invEstTime
			 )
		}

		 	 
return(res)

}


#source("FILE LOCATION")#

##ridiculously fit mutation without sex difference
#Wx<-c(rep(1,25))
#Wy<-c(rep(1,25))
#Wx[c(5,10,15,20:24)]<-1.25
#Wy[c(5,10,15,20:24)]<-1.25
#Wx[25]<-1.5
#Wy[25]<-1.5#

##other parameters
#N = 1000
#m=0
#r=1
#Fii<-findEqFreqs(Wx,Wy, m, r, threshold = 1e-6)
#Fiix.init <- Fii[1,]
#Fiiy.init <- Fii[2,]##

##introduce a clumsy mutation
#Fiix.init[19]  <-  Fiix.init[19] - 1/N
#Fiix.init[20]  <-  1/N#

#test<- autoInvFwdSimSexSpec(Fiix.init = Fiix.init  , Fiiy.init = Fiiy.init  , N = N, Wx = Wx, Wy = Wy, m = m, r = r, 
#							saveTrajectories = T)
#test



