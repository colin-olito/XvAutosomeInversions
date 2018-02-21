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

# female genotypes ordered:
#	c(ABAB, ABAb, ABaB, ABab, ABba*,		c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#	  AbAB, AbAb, AbaB, Abab, Abba*,		  x2y1, x2y2, x2y3, x2y4, x2y5, 
#	  aBAB, aBAb, aBaB, aBab, aBba*,		  x3y1, x3y2, x3y3, x3y4, x3y5, 
#	  abAB, abAb, abaB, abab, abba*,		  x4y1, x4y2, x4y3, x4y4, x4y5, 
#	  baAB*, baAb*, baaB*, baab*, baba*)	  x5y1, x5y2, x5y3, x5y4, x5y5)

# male genotypes ordered:
# c(AB, Ab, aB, ab, ba*)              c(y1, y2, y3, y4, y5)

#' Haplotype frequencies among gametes
#'
#' @title Haplotype frequencies among gametes
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @param r   Recombination rate
#' @export

#' Haplotype frequencies for X in females (XX) is a function of their frequencies in females and males in the previous generation
x.1  <-  function(Fiix=Fiix, m=m, r=r) {
  (((2*Fiix[1] + (Fiix[2] + Fiix[6]) + (Fiix[3] + Fiix[11]) + (Fiix[4] + Fiix[16]) + (Fiix[5] + Fiix[21])) / 2)*(1 - m) - r*Dstar(Fiix=Fiix, m=m) + m
} 
x.2  <-  function(Fiix=Fiix, m=m, r=r) {
  (((2*Fiix[7] + (Fiix[2] + Fiix[6]) + (Fiix[8] + Fiix[12]) + (Fiix[9] + Fiix[17]) + (Fiix[10] + Fiix[22])) / 2)*(1 - m) + r*Dstar(Fiix=Fiix, m=m)
}
x.3  <-  function(Fiix=Fiix, m=m, r=r) {
  (((2*Fiix[13] + (Fiix[3] + Fiix[11]) + (Fiix[8] + Fiix[12]) + (Fiix[14] + Fiix[18]) + (Fiix[15] + Fiix[23])) / 2)*(1 - m) + r*Dstar(Fiix=Fiix, m=m)
}
x.4  <-  function(Fiix=Fiix, m=m, r=r) {
  (((2*Fiix[19] + (Fiix[4] + Fiix[16]) + (Fiix[9] + Fiix[17]) + (Fiix[14] + Fiix[18]) + (Fiix[20] + Fiix[24])) / 2)*(1 - m) - r*Dstar(Fiix=Fiix, m=m)
}
x.5  <-  function(Fiix=Fiix, m=m, r=r) {
  (((2*Fiix[25] + (Fiix[5] + Fiix[21]) + (Fiix[10] + Fiix[22]) + (Fiix[15] + Fiix[23]) + (Fiix[20] + Fiix[24])) / 2)*(1 - m)
}

# Haplotype frequencies for X in males (XY) is a function of only female frequencies in the previous generation
y.1 <- function(Fiiy=Fiiy, m=m) {
  Fiiy[1]*(1-m)+m
} 
y.2 <- function(Fiiy=Fiiy, m=m) {
  Fiiy[2]*(1-m)
}
y.3 <- function(Fiiy=Fiiy, m=m) {
  Fiiy[3]*(1-m)
}
y.4 <- function(Fiiy=Fiiy, m=m) {
  Fiiy[4]*(1-m)
}
y.5 <- function(Fiiy=Fiiy, m=m) {
  Fiiy[5]*(1-m)
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
#' Females:
#'		c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1y1, x1y2, x1y3, x1y4, x1y5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2y1, x2y2, x2y3, x2y4, x2y5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3y1, x3y2, x3y3, x3y4, x3y5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4y1, x4y2, x4y3, x4y4, x4y5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5y1, x5y2, x5y3, x5y4, x5y5)
#'Males:
#'    c(AB, Ab, aB, ab, ba*)                c(y1, y2, y3, y4, y5)
#' @title Find the deterministic equilibrium genotype frequencies prior to introducing the inversion
#' @param Wf     Vector of fitness expressions for all 25 female genotypes
#' @param Wm     Vector of fitness expressions for all 5 male genotypes 
#' @param m      Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param threshold The threshold change in genotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. 
#' @export
#' @seealso `offFreq`, `autoInvWrightFisherSim`
#' @author Colin Olito

findEqFreqs  <-  function(W, m, r, threshold = 1e-6) {
  
  Fiix.init  <-  c(1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  0,    0,    0,    0, 0)
  Fiiy.init <- c(1/4, 1/4, 1/4, 1/4, 0)
  
  Fiix    <-  Fiix.init
  E.Fiix  <-  Fiix.init
  
  Fiiy    <-  Fiiy.init
  E.Fiiy  <-  Fiiy.init
  
  # Storage for female gamete frequencies
  xi        <-  rep(0, times=5)
  names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')
  # Storage for male gamete frequencies
  yi        <-  rep(0, times=5)
  names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')
  # *************From here
  xdelta  <-  rep(1, times=16) # change of haplotype frequencies in 1 generation in females
  ydelta  <-  rep(1, times=4) # change of haplotye frequencies in 1 generation in males
  delta <- append(xdelta, ydelta)
  
  gen = 0 # counter of generations
  while(any(delta > threshold)) {
    gen=gen+1
    for (j in 1:length(xi)) {
      recFct_x  <-  get(names(xi)[j])
      recFct_y <- get(names(yi)[j])
      xi[j]   <-  round(recFct_x(Fii = E.Fii, m = m, r = r), digits=3)
      yi[j]   <-  round(recFct_y(Fii = E.Fii, m = m), digits=3)
    }
    
    # offspring genotype frequencies
    O_females  <-  offFreq_females(xi, yi)
    O_males <- offFreq_males(xi)
    
    # mean fitness 
    Wfbar      <-  sum(O*Wf)
    Wmbar      <-  sum(O*Wm)
    
    # difference in expected frequency (has simulation reached equilibrium yet?)
    xdelta  <-  E.Fiix[c(1:4,6:9,11:14,16:19)] - (O*Wf/Wfbar)[c(1:4,6:9,11:14,16:19)]
    ydelta   <-  E.Fiiy[c(1:4,6:9,11:14,16:19)] - (O*Wm/Wmbar)[c(1:4,6:9,11:14,16:19)]
    delta <- append(xdelta, ydelta)
    E.Fiix   <-  O*Wf/Wfbar
    E.Fiiy   <-  O*Wm/Wmbar
  }
  print (c("reacg equilibrium of initial frequencies at gen", gen))
  
  names(E.Fiix)  <-  NULL
  names(E.Fiiy)  <-  NULL
  return (rbind(E.Fiix, E.Fiiy)) 
}

#' Run a single Wright-Fisher Forward simulation with introduced autosomal inversion
#' using multinomial sampling with linkage to deleterious mutations
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#'  	c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1x1, x1x2, x1x3, x1x4, x1x5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
#'
#' @title Find the deterministic equilibeium genotype frequencies prior to introducing the inversion
#' @param Fii.init  Vector of initial frequencies (deterministic eq. frequencies in absence of inversion)
#' @param N         Population size
#' @param Wf        Vector of fitness expressions for all 25 female genotypes
#' @param Wm        Vector of fitness expressions for all 5 male genotypes
#' @param m         Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r         Recombination rate among the two loci involved in local adaptation in females (r = 0.1).
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? 
#' @export
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
autoInvFwdSim  <-  function(Fiix.init = Fiix.init, Fiiy.init = Fiiy.init, N = N, Wf = Wf, Wm = Wm, m = m, r = r, 
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
  
############# WORKED UP UNTIL HERE 21 FEBRUARY  
  
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
  
  ##  Define Fitness Expressions for females for determining eq. frequencies in absence of inversion
  Wf.init  <-  c(1,          (1 + h*sf),         (1 + h*sf),         (1 + h*sf)^2,        0,
                (1 + h*sf),   (1 + sf),           (1 + h*sf)^2,       (1 + h*sf)*(1 + sf),  0,
                (1 + h*sf),   (1 + h*sf)^2,       (1 + sf),           (1 + sf)*(1 + h*sf),  0,
                (1 + h*sf)^2, (1 + h*sf)*(1 + sf), (1 + sf)*(1 + h*s), (1 + sf)^2,          0,
                0,           0,                 0,                 0,                  0)
  ## Define Fitness Expressions for males for determining eq. frequencies in absence of inversion
  Wm.init  <-  c(1,(1 + sm),(1 + sm),(1 + sm)^2,(1 + sm)^2)
  
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
    