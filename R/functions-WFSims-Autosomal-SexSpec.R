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

!!!TOMODIFY Dstar  <-  function(Fii=Fii, m=m, ...) {
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


x.1  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fiix[1] + (Fiix[2] + Fiix[6]) + (Fiix[3] + Fiix[11]) + (Fiix[4] + Fiix[16]) + (Fiix[5] + Fiix[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
x.2  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fiix[7] + (Fiix[2] + Fiix[6]) + (Fiix[8] + Fiix[12]) + (Fiix[9] + Fiix[17]) + (Fiix[10] + Fiix[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.3  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fiix[13] + (Fiix[3] + Fiix[11]) + (Fiix[8] + Fiix[12]) + (Fiix[14] + Fiix[18]) + (Fiix[15] + Fiix[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.4  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fiix[19] + (Fiix[4] + Fiix[16]) + (Fiix[9] + Fiix[17]) + (Fiix[14] + Fiix[18]) + (Fiix[20] + Fiix[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
x.5  <-  function(Fii=Fiix, m=m, r=r) {
	((2*Fiix[25] + (Fiix[5] + Fiix[21]) + (Fiix[10] + Fiix[22]) + (Fiix[15] + Fiix[23]) + (Fiix[20] + Fiix[24])) / 2)*(1 - m)
}


y.1  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fiiy[1] + (Fiiy[2] + Fiiy[6]) + (Fiiy[3] + Fiiy[11]) + (Fiiy[4] + Fiiy[16]) + (Fiiy[5] + Fiiy[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
y.2  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fiiy[7] + (Fiiy[2] + Fiiy[6]) + (Fiiy[8] + Fiiy[12]) + (Fiiy[9] + Fiiy[17]) + (Fiiy[10] + Fiiy[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.3  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fiiy[13] + (Fiiy[3] + Fiiy[11]) + (Fiiy[8] + Fiiy[12]) + (Fiiy[14] + Fiiy[18]) + (Fiiy[15] + Fiiy[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.4  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fiiy[19] + (Fiiy[4] + Fiiy[16]) + (Fiiy[9] + Fiiy[17]) + (Fiiy[14] + Fiiy[18]) + (Fiiy[20] + Fiiy[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
y.5  <-  function(Fii=Fiiy, m=m, r=r) {
	((2*Fiiy[25] + (Fiiy[5] + Fiiy[21]) + (Fiiy[10] + Fiiy[22]) + (Fiiy[15] + Fiiy[23]) + (Fiiy[20] + Fiiy[24])) / 2)*(1 - m)
}

