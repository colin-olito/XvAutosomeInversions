###############
# DEPENDENCIES
###############

######################
# Necessary functions  
######################

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
#********************************************************************************************
# Run rejectionSamplerX
# rejectionSamplerX(100,100,1e-6,0,0.01,0.01)
# It returns a vector of length 100 of 0 and 1.
#********************************************************************************************

#********************************************************************************************
# For X-linked loci:
# Female genotypes ordered:
# c(ABAB, ABAb, ABaB, ABab, ABba*,    c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#   AbAB, AbAb, AbaB, Abab, Abba*,      x2y1, x2y2, x2y3, x2y4, x2y5, 
#   aBAB, aBAb, aBaB, aBab, aBba*,      x3y1, x3y2, x3y3, x3y4, x3y5, 
#   abAB, abAb, abaB, abab, abba*,      x4y1, x4y2, x4y3, x4y4, x4y5, 
#   baAB*, baAb*, baaB*, baab*, baba*)    x5y1, x5y2, x5y3, x5y4, x5y5)
#
# Male genotypes ordered:
# c(AB, Ab, aB, ab, ba*)              c(y1, y2, y3, y4, y5)
#********************************************************************************************


#' Linkage Disequilibrium function for W-F recursions
#'
#' @title Linkage Disequilibrium (Dstar)
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @export
#' Change Dstar for X chromosome
Dstar  <-  function(Fii=Fii, m=m, ...) {
  (((Fii[4] + Fii[16]) - (Fii[8] + Fii[12])) / 2)*(1 - m)
}
#********************************************************************************************
# Run Dstar
# Returns 0 with the Fii.f vector of x1x4=1/16 and x2x3=1/16
# Change Fii.f as follows to run D.
# Fii.f[4] = 0.0625*2
# Fii.f[8] = 0
# Dstar(Fii.f, m=0) returns 0.0625
# Dstar(Fii.f, m=0.9) return 0.00625
# This means that migration reduces linkage disequilibrium. 
#********************************************************************************************


#' Haplotype frequencies among gametes
#' 
#' @title Haplotype frequencies among gametes
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @param r   Recombination rate
#' @export
#' @author Homa Papoli based on Colin Olito.

#' Haplotype frequencies for X in females (XX)
x.1  <-  function(Fii=Fii, m=mf, r=rf) {
  ((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=m) + m
} 
x.2  <-  function(Fii=Fii, m=mf, r=rf) {
  ((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=mf)
}
x.3  <-  function(Fii=Fii, m=mf, r=rf) {
  ((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fii, m=mf)
}
x.4  <-  function(Fii=Fii, m=mf, r=rf) {
  ((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fii, m=mf)
}
x.5  <-  function(Fii=Fii, m=mf, r=rf) {
  ((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}
# Haplotype frequencies for X in males (XY)
y.1 <- function(Fii=Fii, m=mm) {
  Fii[1]*(1 - m) + m
} 
y.2 <- function(Fii=Fii, m=mm) {
  Fii[2]*(1 - m)
}
y.3 <- function(Fii=Fii, m=mm) {
  Fii[3]*(1 - m)
}
y.4 <- function(Fii=Fii, m=mm) {
  Fii[4]*(1 - m)
}
y.5 <- function(Fii=Fii, m=mm) {
  Fii[5]*(1 - m)
}


#' Offspring frequencies after random mating for X linked haplotypes
#' @title Offspring frequencies after random mating
#' @param xi Vector of haplotype frequencies among gametes (of length = 5)
#' @export
#' Females
offFreq_f  <-  function(xi, yi) {
  O  <-  c((xi[1]*yi[1]), (xi[1]*yi[2]), (xi[1]*yi[3]), (xi[1]*yi[4]), (xi[1]*yi[5]),
           (xi[2]*yi[1]), (xi[2]*yi[2]), (xi[2]*yi[3]), (xi[2]*yi[4]), (xi[2]*yi[5]),
           (xi[3]*yi[1]), (xi[3]*yi[2]), (xi[3]*yi[3]), (xi[3]*yi[4]), (xi[3]*yi[5]),
           (xi[4]*yi[1]), (xi[4]*yi[2]), (xi[4]*yi[3]), (xi[4]*yi[4]), (xi[4]*yi[5]),
           (xi[5]*yi[1]), (xi[5]*yi[2]), (xi[5]*yi[3]), (xi[5]*yi[4]), (xi[5]*yi[5])
  )
  O
}

# Males
offFreq_m  <-  function(xi) {
  O  <-  c(xi[1], xi[2], xi[3], xi[4], xi[5]
  )
  O
}

#' Calculate deterministic equilibrium genotype frequencies 
#' in the absence of the inversion 
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#' Females:
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2y1, x2y2, x2y3, x2y4, x2y5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3y1, x3y2, x3y3, x3y4, x3y5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4y1, x4y2, x4y3, x4y4, x4y5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5y1, x5y2, x5y3, x5y4, x5y5)
#' Males:
#'    c(AB, Ab, aB, ab, ba*)                c(y1, y2, y3, y4, y5)
#' @title Find the deterministic equilibrium genotype frequencies prior to introducing the inversion
#' @param Wf     Vector of fitness expressions for all 25 female genotypes
#' @param Wm     Vector of fitness expressions for all 5 male genotypes 
#' @param mf     Female migration rate for locally maladaptive alleles (mm = mf = 0.01)
#' @param mm     Male migration rate for locally maladaptive alleles (mm = mf = 0.01)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param threshold The threshold change in genotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. 
#' @export
#' @seealso `offFreq`, `autoInvWrightFisherSim`
#' @author Colin Olito - Modified for X-linked by Homa Papoli
findEqFreqsX  <-  function(Wf, Wm, mm, mf, r, threshold = 1e-6) {
  
  # Arbitrarily initiate with equal frequences 
  # for all non-inversion genotypes
  Fii.f.init  <-  c(1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  0,    0,    0,    0, 0)

  Fii.m.init <- c(1/4, 1/4, 1/4, 1/4, 0)
  
  Fii.f    <-  Fii.f.init
  E.Fii.f  <-  Fii.f.init 
  
  Fii.m    <-  Fii.m.init
  E.Fii.m  <-  Fii.m.init 
  
  # Storage for female gamete frequencies
  xi         <-  rep(0, times=5)
  names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')

  # Storage for male gamete frequencies
  yi         <-  rep(0, times=5)
  names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')

  # Storage for differences in frequencies across 1 generation
  deltaF  <-  rep(1, times=16) # change of haplotype frequencies in 1 generation in females
  deltaM  <-  rep(1, times=4) # change of haplotye frequencies in 1 generation in males
  delta   <- append(deltaF, deltaM) # putting deltaF and deltaM in one vector
  
  # Initiate generation counter
  gen  <-  0

  # Simulation loop
  while(any(delta > threshold)) { # It stops when all of them has reached the equilibrium.
    gen  <-  gen+1
    for (j in 1:length(xi)) {
      recFct_x  <-  get(names(xi)[j])
      recFct_y  <- get(names(yi)[j])
      xi[j]     <-  round(recFct_x(Fii = E.Fii.f, m = mf, r = r), digits=3)
      yi[j]     <-  round(recFct_y(Fii = E.Fii.m, m = mm), digits=3)
    }
    
    # offspring genotype frequencies
    O_f  <-  offFreq_f(xi, yi)
    O_m  <- offFreq_m(xi)
    
    # mean fitness 
    Wfbar  <-  sum(O_f*Wf) 
    Wmbar  <-  sum(O_m*Wm)
    
    # difference in expected frequency (has simulation reached equilibrium yet?)
    deltaF   <-  E.Fii.f[c(1:4,6:9,11:14,16:19)] - (O_f*Wf/Wfbar)[c(1:4,6:9,11:14,16:19)]
    deltaM   <-  E.Fii.m[c(1:4)] - (O_m*Wm/Wmbar)[c(1:4)]
    delta    <- append(deltaF, deltaM)
    E.Fii.f  <-  O_f*Wf/Wfbar
    E.Fii.m  <-  O_m*Wm/Wmbar
  }
  print (c("reach equilibrium of initial frequencies at gen", gen))
  
  names(E.Fii.f)  <-  NULL
  names(E.Fii.m)  <-  NULL
  return (list(E.Fii.f, E.Fii.m)) # It returns two vectors, one eq. vector for females and one eq. vector for males. I can't use rbind() because the two vectors have different lengths.
}
#*********************************************************************************
# Run findEqFreqsX
# mm = mf = 0
# r = 0.5
# sf = sm = 0.01
# findEqFreqsX(Wf=Wf.init, Wm=Wm.init, mf=mf, mm=mm, r, threshold = 1e-6)
# It is deterministic so everytime it must give the same result. Everytime, with the above setting, the equilibrium is reached at generation 439. 
# If mm=mf=0.8, equilibrium is reached at generation 5.
# sum(findEqFreqsX(Wf, Wm, 0, 0, 0.5)[[1]]) is genotype frequencies in females and is equal to 1.
# sum(findEqFreqsX(Wf, Wm, 0, 0, 0.5)[[2]]) is genotype frequencies in males and is equal to 1. 
#*********************************************************************************


#' Run a single Wright-Fisher Forward simulation with introduced inversion on the X chromosome
#' using multinomial sampling with linkage to deleterious mutations
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#' Females:
#'  	c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2y1, x2y2, x2y3, x2y4, x2y5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3y1, x3y2, x3y3, x3y4, x3y5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4y1, x4y2, x4y3, x4y4, x4y5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5y1, x5y2, x5y3, x5y4, x5y5)
#'
#'Males:
#'    c(AB, Ab, aB, ab, ba*)                c(y1, y2, y3, y4, y5)
#'
#' @title Find the deterministic equilibrium genotype frequencies prior to introducing the inversion
#' @param Fii.init  Vector of initial frequencies (deterministic eq. frequencies in absence of inversion)
#' @param N         Population size
#' @param Wf        Vector of fitness expressions for all 25 female genotypes
#' @param Wm        Vector of fitness expressions for all 5 male genotypes
#' @param mm, mf         Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r         Recombination rate among the two loci involved in local adaptation in females (r = 0.1).
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? 
#' @export
#' @seealso `offFreq`, `findEqFreqsX`, `x.1`, ...
#' @author Colin Olito - Modified for X-linked by Homa Papoli

InvFwdSimXlinked  <-  function(Fii.f.init = Fii.f.init, Fii.m.init = Fii.m.init, N = N, Wf = Wf, Wm = Wm, mm = mm, mf = mf, sm = sm, sf = sf, r = r, 
                            saveTrajectories = FALSE, ...) {
  
  # Use deterministic eq. initial frequencies
  Fii.f  <-  Fii.f.init
  Fii.m  <-  Fii.m.init
  
  # Define threshold frequency for establishment of inversion
  mbar = 1/3*(2*mf + mm) #!!! Is this correct?
  pcrit  <-  2/(N*mbar) #!!! Is this correct?
  
  # Storage for female gamete frequencies  
  xi         <-  rep(0, times=5)
  names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')
  # Storage for male gamete frequencies  
  yi         <-  rep(0, times=5)
  names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')
    
  if(saveTrajectories) {
    # Storage structures for individual simulation data
    InvFreq    <-  rep(0, times=(4*N+1))
    E.InvFreq  <-  rep(0, times=(4*N+1))
    InvFreq_f    <-  rep(0, times=(4*N+1))
    InvFreq_m    <-  rep(0, times=(4*N+1))
    E.InvFreq_f    <-  rep(0, times=(4*N+1))
    E.InvFreq_m    <-  rep(0, times=(4*N+1))
    W.mean     <-  rep(0, times=(4*N+1)) # Overall fitness of the population regardless of sex
    Wf.mean    <-  rep(0, times=(4*N+1))
    Wm.mean    <-  rep(0, times=(4*N+1))
    
    # Initial inversion frequency
    InvFreq[1]    <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + Fii.m[5])/2 #!!! Sum of all inverted genotypes in females and males.
    E.InvFreq[1]  <-  InvFreq[1]
    InvFreq_f[1]  <-  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25])
    E.InvFreq_f[1] <- InvFreq_f[1]
    InvFreq_m[1]  <-  Fii.m[5]
    E.InvFreq_m[1] <- InvFreq_m[1]
    
    ## Start forward simulation with newly introduced inversion
    gen  <-  1
    while(gen < (4*N) & InvFreq[gen] != 0) { # There is no mutation rate for inversions so loop stops when inversion freq. goes to zero.
      
      ## Step through recursions:
      # 1) Calculate gamete frequencies
      for (j in 1:length(xi)) {
        recFctx  <-  get(names(xi)[j])
        recFcty  <-  get(names(yi)[j])
        xi[j]   <-  round(recFctx(Fii = Fii.f, m = mf, r = r), digits=8)
        yi[j]   <-  round(recFcty(Fii = Fii.m, m = mm), digits=8)
      }
      # 2) Offspring genotype frequencies
      O_f     <-  offFreq_f(xi, yi)
      O_m       <-  offFreq_m(xi)
      # 3) Mean fitness 
      Wfbar   <-  sum(O_f*Wf)
      Wmbar   <-  sum(O_m*Wm)
      Wbar    <-  (Wfbar+Wmbar)/2 #!!! Mean fitness of the population. If Nm=Nf, then it seems division by 2 would be right.
      # 4) Expected frequencies
      E.Fii.f  <-  O_f*Wf/Wfbar
      E.Fii.m   <-  O_m*Wm/Wmbar
      # In deterministic case, E.Fii.f, that is the expected frequency would become the Fii.f for the next round. Here, however, we draw from the multinomial distribution
      # using the expected frequency to get the Fii.f for the new round. In this way, we generate stochasticity in the sampling of genotypes. 
      # 5) Draw random frequencies in adults
      Fii.f    <-  as.vector(rmultinom(1, N/2, E.Fii.f)/(N/2)) 
      Fii.m    <-  as.vector(rmultinom(1, N/2, E.Fii.m)/(N/2))
      # Explanation for rmultinom: 
      # Considering equal sex ratio = Nm = Nf
      # Females: We have 25 genotypes, each has an expected frequency stored in E.Fii.f. rmultinom draws 1 value for genotype ABAB out of N/2 (this means 1 to 50 if N is 100)
      # with a probability which is the frequency of ABAB. This number can be 10. We then divide 10 by 50 which means the new frequency of genotype ABAB will be 0.2. If a genotype
      # has a higher frequency, it is more likely that a larger value will be chosen from 1:50. For example, if genotype ABAB has a frequency of 0.9 and genotype abab has a frequency of
      # 0.01, it is more likely that the number from multinomial draw for ABAB will be around 50. 
      # If Nm = Nf = 20, if 10 males are AB, then frequency of AB males will be 0.5. Dividing by N/2 is correct as we are dealing with genotypes not with haplotypes here. 
      
      # Realized frequencies
      InvFreq[gen+1]    <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + Fii.m[5])/2 #!!!
      E.InvFreq[gen+1]  <-  (sum(E.Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + Fii.m[5])/2 #!!!
      W.mean[gen+1]     <-  Wbar
      
      #The variable below are stored specifically for the X linked simulation
      Wf.mean[gen+1]   =  Wfbar
      Wm.mean[gen+1]   =  Wmbar
      InvFreq_f[gen+1]  =  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25])
      InvFreq_m[gen+1]  =  Fii.m[5]
      E.InvFreq_f[gen+1]=  sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25])
      E.InvFreq_m[gen+1]= Fii.m[5]
      
      gen  <-  gen+1
      print (gen)
    }
   
    # Has the inversion reached threshold frequency for establishment (pcrit)? 
    # When did it first reach pcrit?
    if(any(InvFreq >= pcrit)) { # inversion got established
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
      "InvFreq_f"   =  InvFreq_f[1:gen-1],
      "E.InvFreq_f" =  E.InvFreq_f[1:gen-1],
      "InvFreq_m"   =  InvFreq_m[1:gen-1],
      "E.InvFreq_m" =  E.InvFreq_m[1:gen-1],
      "InvEst"      =  invEst,
      "InvEstTime"  =  invEstTime,
      "W.mean"      =  W.mean[1:gen-1],
      "Wf.mean"     =  Wf.mean[1:gen-1],
      "Wm.mean"     =  Wm.mean[1:gen-1],
      "nGen"        =  gen
    )
  } 
  
  if(!saveTrajectories) { # change as above just not save the trajectories.
    
    # Storage structures for individual simulation data
    InvFreq    <-  0
    E.InvFreq  <-  0
    InvFreq_f    <-  0
    InvFreq_m    <-  0
    E.InvFreq_f    <-  0
    E.InvFreq_m    <-  0
    W.mean     <-  0
    Wf.mean    <-  0
    Wm.mean    <-  0
        
    # Initial inversion frequency 
    InvFreq    <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + Fii.m[5])/2
    E.InvFreq  <-  InvFreq[1]
    InvFreq_f  <-  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25])
    E.InvFreq_f <- InvFreq_f[1]
    InvFreq_m  <-  Fii.m[5]
    E.InvFreq_m <- InvFreq_m[1]
    
    ## Start forward simulation with newly introduced inversion
    gen  <-  1
    #		while(InvFreq > 0 & InvFreq <= pcrit) {
    while(gen < (4*N) & InvFreq > 0 ) {
      
      ## Step through recursions:
      # 1) Calculate gamete frequencies
      for (j in 1:length(xi)) {
        recFctx  <-  get(names(xi)[j])
        recFcty  <-  get(names(yi)[j])
        xi[j]   <-  round(recFctx(Fii = Fii.f, m = mf, r = r), digits=8)
        yi[j]   <-  round(recFcty(Fii = Fii.m, m = mm), digits=8)
      }
      # 2) Offspring genotype frequencies
      O_f     <-  offFreq_f(xi, yi)
      O_m       <-  offFreq_m(xi)
      # 3) Mean fitness 
      Wfbar   <-  sum(O_f*Wf)
      Wmbar   <-  sum(O_m*Wm)
      Wbar    <-  (Wfbar+Wmbar)/2 #!!! Mean fitness of the population. If Nm=Nf, then it seems division by 2 would be right.
      # 4) Expected frequencies
      E.Fii.f  <-  O_f*Wf/Wfbar
      E.Fii.m   <-  O_m*Wm/Wmbar
      # 5) Draw random frequencies in adults
      Fii.f    <-  as.vector(rmultinom(1, N/2, E.Fii.f)/(N/2))
      Fii.m    <-  as.vector(rmultinom(1, N/2, E.Fii.m)/(N/2))
      
      # Realized frequencies
      InvFreq    <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + Fii.m[5])/2
      E.InvFreq  <-  (sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25]) + Fii.m[5])/2
      W.mean     <-  Wbar
      
      #The variable below are stored specifically for the X linked simulation
      Wf.mean   =  sum(O_f*Wf)
      Wm.mean   =  sum(O_m*Wm)
      InvFreq_f  =  sum(Fii.f[c(5,10,15,20:24)]/2, Fii.f[25])
      InvFreq_m  =  Fii.m[5]
      E.InvFreq_f =  sum(E.Fii.f[c(5,10,15,20:24)]/2, E.Fii.f[25])
      E.InvFreq_m = Fii.m[5]
      
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
      "InvFreq_f"   =  InvFreq_f,
      "E.InvFreq_f" =  E.InvFreq_f,
      "InvFreq_m"   =  InvFreq_m,
      "E.InvFreq_m" =  E.InvFreq_m,
      "InvEst"      =  invEst,
      "InvEstTime"  =  invEstTime,
      "Wf.mean"     =  Wf.mean, #!!! Here the code differs from your autosomal version. There you have sum(W.mean)/length(W.mean), but here, as I understand, we're not storing it in a vector so I couldn't understand how sum() & length() would work.
      "Wm.mean"     =  Wm.mean,
      "W.mean"      =  W.mean,
      "nGen"        =  gen
    )
  }
  
  return(res)
}
#*********************************************************************************
# Run InvFwdSimXlinked
# InvFwdSimXlinked(Fii.f, Fii.m, 100, Wf, Wm, 0, 0, 0.01, 0.01, 0.5, saveTrajectories = TRUE) : returns a list of 12 elements
# InvFwdSimXlinked(Fii.f, Fii.m, 100, Wf, Wm, 0, 0, 0.01, 0.01, 0.5, saveTrajectories = FALSE) : returns a list of 12 elements
#*********************************************************************************

#' Introduce new mutant inversion genotype
#'
#' @title Introduce new mutant inversion genotype
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#'   				 chosen randomly, given initial genotypic frequencies (Fii.init). 
#' 					 See params for runReplicateInvSimsXlinked().
#' @param Fii.init   Initial genotypic frequencies, from which to calculate probability of new 
#' 					 mutant inversion occurring
#' @param N			 Population size

# added m, which is true if the inversion is introduced into the male genotype.
# if inversion needs to be introduced in males, set m = TRUE.
introduceInversion  <-  function(newMutant, male = FALSE, Fii.f.init, Fii.m.init, N) {
#' Females:
#'    c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2y1, x2y2, x2y3, x2y4, x2y5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3y1, x3y2, x3y3, x3y4, x3y5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4y1, x4y2, x4y3, x4y4, x4y5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5y1, x5y2, x5y3, x5y4, x5y5)
#'
#'Males:
#'    c(AB, Ab, aB, ab, ba*)                c(y1, y2, y3, y4, y5)
  
  # Toggle
  specifyNewMutant  <-  is.numeric(newMutant) # it returns TRUE if newMutant is the number of a specific genotype given. 
                                              # In this case, sex and genotype are pre-defined and are not random.
  
  # Choose mutant genotype randomly
  if(!specifyNewMutant) {
    
      # Probability of new mutant occuring on X in a genotype containing ab in females and males
      # 7 genotypes in females contain ab haplotype, only one contains 2 ab, in this case, the frequency of abab is multiplied by 2. For the other genotypes, it's only 1 ab within each. 
      probNewMutantX     <-  c(Fii.f.init[c(4,9,14,16:18)], Fii.f.init[19]*2, Fii.m.init[4])/sum(Fii.f.init[c(4,9,14,16:18)], Fii.f.init[19]*2 , sum(Fii.m.init[4]))
      newMutX            <-  c(4,9,14,16:19,100)[as.vector(rmultinom(1,1,probNewMutantX)) == 1] # In the vector c(4,9,14,16:19,100), 100 is a number that indicates the male genotype. If this is chosen by rmultinom, 
      if(newMutX == 100) {                                                                      # then the inversion is introduced in male which is indicated by the condition if(newMutX == 100)
        # Subtract new mutant individual from frequency of old genotype
        Fii.m.init[4]  <-  Fii.m.init[4] - 1/N
        Fii.m.init[5] <- 1/N
      }
      else {
        Fii.f.init[newMutX]  <-  Fii.f.init[newMutX] - 1/N
      }
    }
  
  # Specify mutant genotype
  if(specifyNewMutant) {
    # Subtract new mutant individual from frequency of old genotype
    if(male) {
      newMutY  <-  newMutant
      Fii.m.init[newMutY]  <-  Fii.m.init[newMutY] - 1/N
        Fii.m.init[5] <- 1/N
    }
    else {
      newMutX  <-  newMutant
      Fii.f.init[newMutX]  <-  Fii.f.init[newMutX] - 1/N
    }
  }

  # Add mutant individual to frequency of new inversion genotype
  if(newMutX == 4 | newMutX == 9 | newMutX == 14)
    Fii.f.init[newMutX + 1]  <-  1/N
  if(newMutX == 16 | newMutX == 17 | newMutX == 18)
    Fii.f.init[newMutX + 5]  <-  1/N
  
  
  # if inversion occurs on abab genotype, choose randomly whether it occurs on
  # the maternally or paternally inherited chromosome 
  if(newMutX == 19) {
    if(runif(1) >= 1/2) {
      Fii.f.init[newMutX + 1]  <-  1/N
    }
    else Fii.f.init[newMutX + 5]  <-  1/N
  } 
  return (list(Fii.f.init, Fii.m.init))
}

#*********************************************************************************
# Run introduceInversion
# Returns a list. First element is Fii.f.init and second element is Fii.m.init after
# introducing the inversion.
# Introduce inversion on the 4th genotype of males
# introduceInversion(4, male = TRUE, Fii.f.init, Fii.m.init, 100)
# Introduce inversion on the 4th genotype of females
# introduceInversion(4, male = FALSE, Fii.f.init, Fii.m.init, 100)
# Introduce inversion in a random genotype and sex
# introduceInversion("random", male = FALSE, Fii.f.init, Fii.m.init, N)
#*********************************************************************************

#' Wrapper function to run replicate forward simulations for invasion
#' of X-linked inversions in a Wright-Fisher population 
#'
#' @title Wright-Fisher forward simulation of genotypic frequencies (default parameter values in parentheses)
#' @param nReps  Number of replicate simulations. With no deleterious mutations, and introducing 
#' 				 a single copy of the inversion, it takes 1,600,000 replicate simulations to 
#'				 get 10,000 where the inversion successfully establishes.
#' @param N      Effective population size
#' @param mm, mf Male and female migration rate for locally maladaptive alleles (mm = mf =  0.01)
#' @param sm, sf  Male and female selective advantage of locally adaptive alleles over migrant alleles (sm = sf = 0.02)
#' @param h      Dominance coefficient for locally adaptive alleles relative to migrant alleles (h = 0.5)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1)
#' @param n      Number of loci at which deleterious mutations may occur
#' @param u      Mutation rate (default value of u = 1e-6)
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Default value of FALSE assumes that 
#' 				 selection against deleterious mutations is twice as strong as selection favouring
#' 				 the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				 there were no delterious mutations segregating in the population that are linked
#' 				 to the loci involved in local adaptation. 
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? Setting this 
#' 							to TRUE can become extremely memory intensive if you are running many
#' 							replicate simulations (see warning).
#' 							otherwise
#' @seealso `offFreq`, `findEqFreqsX`, `InvFwdSimXlinked`
#' @export
#' @author Colin Olito - Modified for X-linked by Homa Papoli

runReplicateInvSimsXlinked  <-  function(nReps = 1000, N = 500, mm = 0.01, mf = 0.01, sf = 0.02, sm = 0.02,  h = 1/2, r = 0.1, 
                                      n = 100, u = 1e-5, h.del = 0, sf.del = 1, sm.del = 1, noDel = FALSE,
                                      saveTrajectories = FALSE, newMutant = "random", male = FALSE) { # male = FALSE is for the introduceInversion function
                                                                                                      # if new inversion should be introduced in male, set male = TRUE, otherwise, it will be introduced in female.
                                                                                                      # if newMutant="random", set male = FALSE
  ##  Preemptive Warnings
  if(any(c(N,mm,mf,sf,sm,h,r,n,u,h.del) < 0) | any(c(mm,mf,sf,sm,h,r,u,h.del) > 1) | r > 0.5)
    stop('The chosen parameter values fall outside of reasonable parameter space')
  
  specifyNewMutant <- is.numeric(newMutant)
  if(specifyNewMutant & all(newMutant != c(4,9,14,16:19)))
    stop('If specifying the genotype of new inversion mutants, newMutant must take one of the following values: 4,9,14,16:19')
  
  if(!specifyNewMutant & newMutant != 'random')
    stop('If the genotype of new inversion mutants is being chosen randomly, the parameter newMutant must equal random')
  
  #!!! I set saveTrajectories = TRUE but I couldn't get it issue warning. It also returned error for the first try about migration.
#   try({
#     if(m >= s )
#       stop('Warning: migration is stronger than selection, 
#            adaptive alleles will be swamped by maladaptive migrants')
#   }#, silent=FALSE)
#   
#   try({
#     if(nReps > 1000 & saveTrajectories)
#       stop('Warning: You have chosen to save evolutionary trajectories 
#            for a large number of replicate simulations. This will be
#            memory intensive. Consider setting saveTrajectories = FALSE')
#   }#, silent=FALSE)
#   

  ##  Define Fitness Expressions for females for determining eq. frequencies in absence of inversion
  Wf.init  <-  c(1,          (1 + h*sf),         (1 + h*sf),         (1 + h*sf)^2,        0,
                (1 + h*sf),   (1 + sf),           (1 + h*sf)^2,       (1 + h*sf)*(1 + sf),  0,
                (1 + h*sf),   (1 + h*sf)^2,       (1 + sf),           (1 + sf)*(1 + h*sf),  0,
                (1 + h*sf)^2, (1 + h*sf)*(1 + sf), (1 + sf)*(1 + h*sf), (1 + sf)^2,          0,
                0,           0,                 0,                 0,                  0)
  ## Define Fitness Expressions for males for determining eq. frequencies in absence of inversion
  Wm.init  <-  c(1, (1 + sm), (1 + sm), (1 + sm)^2, 0)
  
  ## Find deterministic equilibrium frequencies in absence of inversion in females and males
  eq_freq <- findEqFreqsX(Wf=Wf.init, Wm=Wm.init, mm=mm, mf=mf, r=r, threshold = 1e-7)
  Fii.f.init <- eq_freq[[1]] # First element of the list for females  
  Fii.m.init <- eq_freq[[2]] # Second element of the list for for males
  
  ## Introduce rare inversion mutations
  new_inversion <- introduceInversion(newMutant=newMutant, male = male, Fii.f.init=Fii.f.init, Fii.m.init=Fii.m.init, N = N)
  Fii.f.init <- new_inversion[[1]] # First element of the list for females 
  Fii.m.init <- new_inversion[[2]] # Second element of the list for for males
  
  # Storage structures for replicate simulation data
  finalInvFreq    <-  rep(0, times=nReps)
  finalE.InvFreq  <-  rep(0, times=nReps)
  finalInvFreq_f   <-  rep(0, times=nReps)
  finalE.InvFreq_f   <-  rep(0, times=nReps)
  finalInvFreq_m   <-  rep(0, times=nReps)
  finalE.InvFreq_m   <-  rep(0, times=nReps)
  invEst          <-  rep(0, times=nReps)
  invEstTime      <-  rep(0, times=nReps)
  finalW.mean     <-  rep(0, times=nReps)
  finalWf.mean     <-  rep(0, times=nReps)
  finalWm.mean     <-  rep(0, times=nReps)
  nGen            <-  rep(0, times=nReps)
  nDels           <-  rep(0, times=nReps)
  
  if(saveTrajectories) {
    replicateTraj  <-  c()
    InvFreqTraj    <-  c()
    E.InvFreqTraj  <-  c()
    InvFreq_fTraj   <-  c()
    E.InvFreq_fTraj  <- c()
    InvFreq_mTraj   <-  c()
    E.InvFreq_mTraj  <- c()
    W.meanTraj     <-  c()
    Wf.meanTraj     <-  c()
    Wm.meanTraj     <-  c()
  } 
  


  # Replicate simulation loop
# 	print('Running Wright-Fisher Forward Simulations')
  pb   <-  txtProgressBar(min=0, max=nReps, style=3)
  setTxtProgressBar(pb, 0)
  for(i in 1:nReps) {
    
    ## Sample stationary distribution of deleterious alleles
    delMutFreq  <-  rejectionSamplerX(n=n, Ne=N, u=u)
    n.del       <-  sum(delMutFreq > runif(n=n))
    
    # Define fitness expressions, including fitness effects of deleterious mutations for females and males
    if(noDel) {
      sf.del  <-  0
      sm.del  <-  0
    }
    Wf  <-  c(1,                                  (1 + h*sf),                                 (1 + h*sf),                                 (1 + h*sf)^2,                       (1 + h*sf)^2*(1 - h.del*sf.del)^n.del,
             (1 + h*sf),                           (1 + sf),                                   (1 + h*sf)^2,                               (1 + h*sf)*(1 + sf),                 (1 + h*sf)*(1 + sf)*(1 - h.del*sf.del)^n.del,
             (1 + h*sf),                           (1 + h*sf)^2,                               (1 + sf),                                   (1 + sf)*(1 + h*sf),                 (1 + sf)*(1 + h*sf)*(1 - h.del*sf.del)^n.del,
             (1 + h*sf)^2,                         (1 + h*sf)*(1 + sf),                         (1 + sf)*(1 + h*sf),                         (1 + sf)^2,                         (1 + sf)^2*(1 - h.del*sf.del)^n.del,
             (1 + h*sf)^2*(1 - h.del*sf.del)^n.del, (1 + h*sf)*(1 + sf)*(1 - h.del*sf.del)^n.del, (1 + sf)*(1 + h*sf)*(1 - h.del*sf.del)^n.del, (1 + sf)^2*(1 - h.del*sf.del)^n.del, (1 + sf)^2*(1 - sf.del)^n.del)
    
    Wm  <-  c(1, (1 + sm), (1 + sm), (1 + sm)^2, (1 + sm)^2*(1 - sm.del)^n.del)

    ## RUN SIMULATION
    repRes  <-  InvFwdSimXlinked(Fii.f.init=Fii.f.init, Fii.m.init=Fii.m.init, N=N, Wf=Wf, Wm=Wm, mm=mm, mf=mf, r=r)
    
    # save results for each replicate
    finalInvFreq[i]    <-  repRes$InvFreq[length(repRes$InvFreq)]
    finalE.InvFreq[i]  <-  repRes$E.InvFreq[length(repRes$E.InvFreq)]
    finalInvFreq_f[i]  <-  repRes$InvFreq_f[length(repRes$InvFreq_f)]
    finalE.InvFreq_f[i]  <- repRes$E.InvFreq_f[length(repRes$E.InvFreq_f)]
    finalInvFreq_m[i]  <-  repRes$InvFreq_m[length(repRes$InvFreq_m)]
    finalE.InvFreq_m[i]  <- repRes$E.InvFreq_m[length(repRes$E.InvFreq_m)]
    finalW.mean[i]     <-  repRes$W.mean
    finalWf.mean[i]     <-  repRes$Wf.mean
    finalWm.mean[i]     <-  repRes$Wm.mean
    nGen[i]            <-  repRes$nGen
    invEst[i]          <-  repRes$InvEst
    invEstTime[i]      <-  repRes$InvEstTime
    nDels[i]           <-  n.del
    
    
    
    if(saveTrajectories) {
      replicateTraj  <-  c(replicateTraj, rep(i, times=length(repRes$InvFreq)))
      InvFreqTraj    <-  c(InvFreqTraj, repRes$InvFreq)
      E.InvFreqTraj  <-  c(E.InvFreqTraj, repRes$E.InvFreq)
      InvFreq_fTraj  <-  c(InvFreq_fTraj, repRes$InvFreq_f)
      InvFreq_mTraj  <-  c(InvFreq_mTraj, repRes$InvFreq_m)
      E.InvFreq_fTraj <- c(E.InvFreq_fTraj, repRes$E.InvFreq_f)
      E.InvFreq_mTraj <- c(E.InvFreq_mTraj, repRes$E.InvFreq_m)
      W.meanTraj     <-  c(W.meanTraj, repRes$W.mean)
      Wf.meanTraj     <-  c(Wf.meanTraj, repRes$Wf.mean)
      Wm.meanTraj     <-  c(Wm.meanTraj, repRes$Wm.mean)
    } 
    
    setTxtProgressBar(pb, i)
  }
  
  # Save results and return results as a list
  results.df  <-  data.frame(
    "finalInvFreq"    =  finalInvFreq,
    "finalE.InvFreq"  =  finalE.InvFreq,
    "finalInvFreq_f"  =  finalInvFreq_f,
    "finalE.InvFreq_f" = finalE.InvFreq_f,
    "finalInvFreq_m"  =  finalInvFreq_m,
    "finalE.InvFreq_m" = finalE.InvFreq_m,
    "finalW.mean"      = finalW.mean,  
    "finalWf.mean"     =  finalWf.mean,
    "finalWm.mean"     =  finalWm.mean,
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
      "InvFreq_fTraj"  = InvFreq_fTraj,
      "E.InvFreq_fTraj" = E.InvFreq_fTraj,
      "InvFreq_mTraj"  = InvFreq_mTraj,
      "E.InvFreq_mTraj" = E.InvFreq_mTraj,
      "W.meanTraj"      = W.meanTraj,
      "Wf.meanTraj"     =  Wf.meanTraj,
      "Wm.meanTraj"     =  Wm.meanTraj
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
#' of X-linked inversions in a Wright-Fisher population 
#' USING DESIGNATED PARAMETER VALUES
#'
#' @title Run replicate Wright-Fisher forward simulations for autosomal inversion under different parameter values 
#' @param N.vals Desired population sizes
#' @param mm.vals desired male migration rates for locally maladaptive alleles (mm =  0.01)
#' @param mf.vals desired female migration rates for locally maladaptive alleles (mf =  0.01)
#' @param s.del.vals  Desired selection coefficients for deleterious mutations (default value of s = 1).
#' @param sf      Female selective advantage of locally adaptive alleles over migrant alleles (sf = 0.02)
#' @param sm      Male selective advantage of locally adaptive alleles over migrant alleles (sm = 0.02)
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
#' @seealso `offFreq`, `findEqFreqsX`, `InvFwdSimXlinked`
#' @export
#' @author Colin Olito - Modified for X-linked by Homa Papoli
makeReplicateInvSimsDataXlinked  <-  function(nReps = 1000, N = 5000, mm.vals = c(0.01, 0.05), mf.vals = c(0.01, 0.05), 
                                           sf = 0.1, sm = 0.1, h = 1/2, r.vals = c(0.5, 0.1), n = 100, u = 1e-5, h.del = 0, newMutant="random", male = FALSE, saveTrajectories = FALSE) {
  
  # Simulate deleterious mutations that are either 
  # 1) recessive lethals OR
  # 2) recessive experiencing purifying selection
  #    that is twice as strong as the selective 
  #    advantage of the locally adaptive allels  
  sf.del.vals = c(0, 1, 2*sf)
  sm.del.vals = c(0, 1, 2*sm)
  
  # create empty data frame with same structure as we are going to need
  data  <-  data.frame(matrix(ncol=13, nrow=0))
  
  # Convenience variables to monitor progress
  prog  <-  0
  tot   <-  length(r.vals)*length(mm.vals)*length(sf.del.vals)
  # For single population size. Loop of population size will change with recombination rate. 
  # Loop over parameter values we want to explore 
  for(j in 1:length(r.vals)) {
    for(k in 1:length(mm.vals)) {
      for(l in 1:length(sf.del.vals)) {
        
        # Display progress in terminal
        prog  <-  prog + 1
        cat("\n",paste('Running simulations for parameter set ', prog, "/", tot),"\n")
        
        # Run simulations  
        res  <-  runReplicateInvSimsXlinked(nReps = nReps, N = N, mm = mm.vals[k], mf = mf.vals[k], sf = sf, sm = sm, h = h, r = r.vals[j], 
                                         n = n, u = u, h.del = h.del, sf.del = sf.del.vals[l], sm.del = sm.del.vals[l],
                                         noDel = FALSE, saveTrajectories = FALSE)
        
        # Save data 
        rs      <-  rep(r.vals[j], times=nrow(res$results.df))
        mms      <-  rep(mm.vals[k], times=nrow(res$results.df))
        mfs      <-  rep(mf.vals[k], times=nrow(res$results.df))
        sf.dels  <-  rep(sf.del.vals[l], times=nrow(res$results.df))
        sm.dels  <-  rep(sm.del.vals[l], times=nrow(res$results.df))
        
        # Append to data frame
        df      <-  cbind(res$results.df, rs, mms, mfs, sf.dels, sm.dels)
        data    <-  rbind(data, df)
        rm(df)
        
      }
    }
  }
  
  # Include constant variables in data frame
  sfs      <-  rep(sf, times=nrow(data))
  sms      <-  rep(sm, times=nrow(data))
  hs      <-  rep(h, times=nrow(data))
  Ns      <-  rep(N, times=nrow(data))
  us      <-  rep(u, times=nrow(data))
  h.dels  <-  rep(h.del, times=nrow(data))
  data    <-  cbind(data, sfs, sms, hs, Ns, us, h.dels)
  colnames(data)  <-  c("finalInvFreq", "finalE.InvFreq", "finalInvFreq_f", "finalE.InvFreq_f", "finalInvFreq_m", "finalE.InvFreq_m", 
                        "finalW.mean", "finalWf.mean", "finalWm.mean", "nGen", "invEst", "invEstTime", "nDels", "rs", "mms", "mfs",
                        "sf.dels", "sm.dels", "sfs", "sms", "hs", "Ns", "us", "h.dels")
  
  # create file name
  filename  <-  paste("./output/data/simResults/Xlinked-InvSimsData", "_sf", sf, "_sm", sm, "_h", h, "_N", N, "_n", n, "_u", u, "_hdel", h.del, ".csv", sep="")

  # export data as .csv to ./output/data
  write.csv(data, file=filename, row.names = FALSE)
  
  #  Return results in case user wants it
  return(data)
  
}    
