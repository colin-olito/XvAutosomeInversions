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


# Whoeveer is going to do the X-linked model
# should start writing their functions here (using the code provided in 
# './R/functions-WFSims-Autosomal.R' as a template)...
#
#' Linkage Disequilibrium function for W-F recursions
#'
#' @title Linkage Disequilibrium (Dstar)
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @export
#' Change Dstar for X chromosome:
Dstar  <-  function(Fii=Fii, m=m, ...) {
  (((Fii[4] + Fii[16])/2 - (Fii[8] + Fii[12])/2) / 2)*(1 - m)
}

# genotypes ordered:
#	c(ABAB, ABAb, ABaB, ABab, ABba*,		c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#	  AbAB, AbAb, AbaB, Abab, Abba*,		  x2y1, x2y2, x2y3, x2y4, x2y5, 
#	  aBAB, aBAb, aBaB, aBab, aBba*,		  x3y1, x3y2, x3y3, x3y4, x3y5, 
#	  abAB, abAb, abaB, abab, abba*,		  x4y1, x4y2, x4y3, x4y4, x4y5, 
#	  baAB*, baAb*, baaB*, baab*, baba*)	  x5y1, x5y2, x5y3, x5y4, x5y5)


#' Haplotype frequencies among gametes
#'
#' @title Haplotype frequencies among gametes
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @param r   Recombination rate
#' @export

#' Haplotype frequencies for X in females (XX) is a function of their frequencies in females and males in the previous generation
x.1  <-  function(Fii=Fii, m=m, r=r) {
  (((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21]))+(Fii[1] + Fii[6] + Fii[11] + Fii[16] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
x.2  <-  function(Fii=Fii, m=m, r=r) {
  (((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22]))+(Fii[7] + Fii[2] + Fii[12] + Fii[17] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.3  <-  function(Fii=Fii, m=m, r=r) {
  (((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23]))+(Fii[13] + Fii[3] + Fii[8] + Fii[18] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
x.4  <-  function(Fii=Fii, m=m, r=r) {
  (((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24]))+(Fii[19] + Fii[4] + Fii[9] + Fii[14] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
x.5  <-  function(Fii=Fii, m=m, r=r) {
  (((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24]))+(Fii[25] + Fii[5] + Fii[10] + Fii[15] + Fii[20])) / 2)*(1 - m)
}

# Haplotype frequencies for X in males (XY) is a function of only female frequencies in the previous generation
y.1 <- function(Fii=Fii, m=m, r=r) {
  ((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
y.2 <- function(Fii=Fii, m=m, r=r) {
  ((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.3 <- function(Fii=Fii, m=m, r=r) {
  ((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=m)
}
y.4 <- function(Fii=Fii, m=m, r=r) {
  ((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m)
}
y.5 <- function(Fii=Fii, m=m, r=r) {
  ((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}

#' Offspring frequencies after random mating for X linked haplotypes
#' @title Offspring frequencies after random mating
#' @param xi Vector of haplotype frequencies among gametes (of length = 5)
#' @export
#' Females
offFreq_females  <-  function(xi, yi) {
  O  <-  c((xi[1]*yi[1]), (xi[1]*yi[2]), (xi[1]*yi[3]), (xi[1]*yi[4]), (xi[1]*yi[5]),
           (xi[2]*yi[1]), (xi[2]*yi[2]), (xi[2]*yi[3]), (xi[2]*yi[4]), (xi[2]*yi[5]),
           (xi[3]*yi[1]), (xi[3]*yi[2]), (xi[3]*yi[3]), (xi[3]*yi[4]), (xi[3]*yi[5]),
           (xi[4]*yi[1]), (xi[4]*yi[2]), (xi[4]*yi[3]), (xi[4]*yi[4]), (xi[4]*yi[5]),
           (xi[5]*yi[1]), (xi[5]*yi[2]), (xi[5]*yi[3]), (xi[5]*yi[4]), (xi[5]*yi[5])
  )
  O
}

# Males
offFreq_males  <-  function(xi) {
  O  <-  c(xi[1], xi[2], xi[3], xi[4], xi[5]
  )
  O
}

#***** Change genotype frequency functions Fii to correspond to males and females and fitness functions

#' Calculate deterministic equilibrium genotype frequencies 
#' in the absence of the inversion 
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2y1, x2y2, x2y3, x2y4, x2y5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3y1, x3y2, x3y3, x3y4, x3y5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4y1, x4y2, x4y3, x4y4, x4y5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5y1, x5y2, x5y3, x5y4, x5y5)
#'
#' @title Find the deterministic equilibrium genotype frequencies prior to introducing the inversion
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
  xi        <-  rep(0, times=5)
  yi        <-  rep(0, times=5)
  names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')
  names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')
  
  # *************From here
  delta  <-  rep(1, times=16)
  while(any(delta > threshold)) {
    for (j in 1:length(xi)) {
      recFct_x  <-  get(names(xi)[j])
      recFct_y <- get(names(yi)[j])
      xi[j]   <-  round(recFct_x(Fii = E.Fii, m = m, r = r), digits=3)
      yi[j]   <-  round(recFct_y(Fii = E.Fii, m = m, r = r), digits=3)
    }
    
    # offspring genotype frequencies
    O_females  <-  offFreq_females(xi, yi)
    O_males <- offFreq_males(xi)
    
    # mean fitness 
    Wbar      <-  sum(O*W)
    
    # difference in expected frequency (has simulation reached equilibrium yet?)
    delta   <-  E.Fii[c(1:4,6:9,11:14,16:19)] - (O*W/Wbar)[c(1:4,6:9,11:14,16:19)]
    E.Fii   <-  O*W/Wbar
  }
  names(E.Fii)  <-  NULL
  E.Fii 
}


# Make some changes here...###HOMA
# just clone into the github repository, then run 
#
# git checkout -b 'homa' OR git checkout -b 'ludo'
#
# to create your branch. then maybe make a couple trivial edits and push 
# the changes to get you branch started.
