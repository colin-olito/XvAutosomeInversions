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

#' Linkage Disequilibrium function for W-F recursions
#'
#' @title Linkage Disequilibrium (Dstar)
#' @param Fii Vector of adult genotypic frequencies (of length = 25)
#' @param m   Migration rate
#' @export
#' Change Dstar for X chromosome:
Dstar  <-  function(Fii=Fiix, m=m, ...) {
  (((Fii[4] + Fii[16]) - (Fii[8] + Fii[12])) / 2)*(1 - m)
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

#' Haplotype frequencies for X in females (XX)
x.1  <-  function(Fii=Fiix, m=m, r=r) {
  ((2*Fii[1] + (Fii[2] + Fii[6]) + (Fii[3] + Fii[11]) + (Fii[4] + Fii[16]) + (Fii[5] + Fii[21])) / 2)*(1 - m) - r*Dstar(Fii=Fiix, m=m) + m
} 
x.2  <-  function(Fii=Fiix, m=m, r=r) {
  ((2*Fii[7] + (Fii[2] + Fii[6]) + (Fii[8] + Fii[12]) + (Fii[9] + Fii[17]) + (Fii[10] + Fii[22])) / 2)*(1 - m) + r*Dstar(Fii=Fiix, m=m)
}
x.3  <-  function(Fii=Fiix, m=m, r=r) {
  ((2*Fii[13] + (Fii[3] + Fii[11]) + (Fii[8] + Fii[12]) + (Fii[14] + Fii[18]) + (Fii[15] + Fii[23])) / 2)*(1 - m) + r*Dstar(Fii=Fiix, m=m)
}
x.4  <-  function(Fii=Fiix, m=m, r=r) {
  ((2*Fii[19] + (Fii[4] + Fii[16]) + (Fii[9] + Fii[17]) + (Fii[14] + Fii[18]) + (Fii[20] + Fii[24])) / 2)*(1 - m) - r*Dstar(Fii=Fiix, m=m)
}
x.5  <-  function(Fii=Fiix, m=m, r=r) {
  ((2*Fii[25] + (Fii[5] + Fii[21]) + (Fii[10] + Fii[22]) + (Fii[15] + Fii[23]) + (Fii[20] + Fii[24])) / 2)*(1 - m)
}

# Haplotype frequencies for X in males (XY)
y.1 <- function(Fii=Fiiy, m=m) {
  Fii[1]*(1-m)+m
} 
y.2 <- function(Fii=Fiiy, m=m) {
  Fii[2]*(1-m)
}
y.3 <- function(Fii=Fiiy, m=m) {
  Fii[3]*(1-m)
}
y.4 <- function(Fii=Fiiy, m=m) {
  Fii[4]*(1-m)
}
y.5 <- function(Fii=Fiiy, m=m) {
  Fii[5]*(1-m)
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
#' @param mm, mf      Migration rate for locally maladaptive alleles (m =  0.01)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param threshold The threshold change in genotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. 
#' @export
#' @seealso `offFreq`, `autoInvWrightFisherSim`
#' @author Colin Olito

findEqFreqs  <-  function(Wf, Wm, mm, mf, r, threshold = 1e-6) {
  
  Fiix.init  <-  c(1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  1/16, 1/16, 1/16, 1/16, 0, 
                  0,    0,    0,    0, 0)
  
  Fiix    <-  Fiix.init
  E.Fiix  <-  Fiix.init 
  
  Fiiy.init <- c(1/4, 1/4, 1/4, 1/4, 0)
  
  Fiiy    <-  Fiiy.init
  E.Fiiy  <-  Fiiy.init 
  
  # Storage for female gamete frequencies
  xi        <-  rep(0, times=5)
  names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')
  # Storage for male gamete frequencies
  yi        <-  rep(0, times=5)
  names(yi)  <-  c('y.1', 'y.2', 'y.3', 'y.4', 'y.5')

  xdelta  <-  rep(1, times=16) # change of haplotype frequencies in 1 generation in females
  ydelta  <-  rep(1, times=4) # change of haplotye frequencies in 1 generation in males
  delta <- append(xdelta, ydelta) # putting xdelta and ydelta in one vector
  
  gen = 0 # counter of generations
  while(any(delta > threshold)) { # It stops when all of them has reached the equilibrium.
    gen=gen+1
    for (j in 1:length(xi)) {
      recFct_x  <-  get(names(xi)[j])
      xi[j]   <-  round(recFct_x(Fii = E.Fiix, m = mf, r = r), digits=3)
      recFct_y <- get(names(yi)[j])
      yi[j]   <-  round(recFct_y(Fii = E.Fiiy, m = mm), digits=3)
    }
    
    # offspring genotype frequencies
    O_females  <-  offFreq_females(xi, yi)
    O_males <- offFreq_males(xi)
    
    # mean fitness 
    Wfbar      <-  sum(O_females*Wf)
    Wmbar      <-  sum(O_males*Wm)
    
    # difference in expected frequency (has simulation reached equilibrium yet?)
    xdelta  <-  E.Fiix[c(1:4,6:9,11:14,16:19)] - (O_females*Wf/Wfbar)[c(1:4,6:9,11:14,16:19)]
    ydelta   <-  E.Fiiy[c(1:4)] - (O_males*Wm/Wmbar)[c(1:4)]
    delta <- append(xdelta, ydelta)
    E.Fiix   <-  O_females*Wf/Wfbar
    E.Fiiy   <-  O_males*Wm/Wmbar
  }
  print (c("reach equilibrium of initial frequencies at gen", gen))
  
  names(E.Fiix)  <-  NULL
  names(E.Fiiy)  <-  NULL
  return (list(E.Fiix, E.Fiiy)) # It returns two vectors, one eq. vector for females and one eq. vector for males.
}



#' Run a single Wright-Fisher Forward simulation with introduced inversion on the X chromosome
#' using multinomial sampling with linkage to deleterious mutations
#'
#' Genotypes are organized as numeric vectors of length = 25. For consistency
#' they are always ordered as follows:
#' Females:
#'  	c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,	c(x1x1, x1x2, x1x3, x1x4, x1x5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
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
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
autoInvFwdSim  <-  function(Fiix.init = Fiix.init, Fiiy.init = Fiiy.init, N = N, Wf = Wf, Wm = Wm, mm = mm, mf = mf, sm = sm, sf = sf, r = r, 
                            saveTrajectories = FALSE, ...) {
  
  # Use deterministic eq. initial frequencies
  Fiix  <-  Fiix.init
  Fiiy  <-  Fiiy.init
  
  # Define threshold frequency for establishment of inversion
  pcrit  <-  2/(N*(mm+mf)/2)
  
  # Storage for gamete frequencies  
  xi         <-  rep(0, times=5)
  names(xi)  <-  c('x.1', 'x.2', 'x.3', 'x.4', 'x.5')
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
    
    # Initial inversion frequency #!InvFreq for females and males separate or as is it?!
    InvFreq[1]    <-  sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25]) + Fiiy[5] # Sum of all inverted genotypes in females and males.
    E.InvFreq[1]  <-  InvFreq[1]
    
    ## Start forward simulation with newly introduced inversion
    gen  <-  1
    while(gen < (4*N) & InvFreq[gen] != 0) { # no mutation rate for inversions so loop stops when inversion freq. goes to zero.
      
      ## Step through recursions:
      # 1) Calculate gamete frequencies
      for (j in 1:length(xi)) {
        recFctx  <-  get(names(xi)[j])
        recFcty  <-  get(names(yi)[j])
        xi[j]   <-  round(recFctx(Fii = Fiix, m = mf, r = r), digits=8)
        yi[j]   <-  round(recFcty(Fii = Fiiy, m = mm, r = r), digits=8)
      }
      # 2) Offspring genotype frequencies
      O_females     <-  offFreq_females(xi, yi)
      O_males       <-  offFreq_males(xi)
      # 3) Mean fitness 
      Wfbar   <-  sum(O_females*Wf)
      Wmbar   <-  sum(O_males*Wm)
      # 4) Expected frequencies
      E.Fiix  <-  O*Wf/Wfbar
      E.Fiiy   <-  O*Wm/Wmbar
      # 5) Draw random frequencies in adults
      Fiix    <-  as.vector(rmultinom(1, N/2, E.Fiix)/(N/2)) # N/2 also applies for the X chromosome?
      Fiiy    <-  as.vector(rmultinom(1, N/2, E.Fiiy)/(N/2))
      
      # Realized frequencies
      InvFreq[gen+1]    <-  (sum(Fiix[c(5,10,15,20:24)]/2, Fii[25]) + Fiiy[5])
      E.InvFreq[gen+1]  <-  (sum(E.Fiix[c(5,10,15,20:24)]/2, Fii[25]) + Fiiy[5])
      W.mean[gen+1]     <-  Wbar
      
      #The variable below are stored specifically for the X linked simulation
      Wf.mean[gen+1]   =  sum(O_females*Wf)
      Wm.mean[gen+1]   =  sum(O_males*Wm)
      InvFreq_f[gen+1]  =  sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25])
      InvFreq_m[gen+1]  =  Fiiy[5]
      E.InvFreq_f[gen+1]=  sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25])
      E.InvFreq_m[gen+1]= Fiiy[5]
      
      gen  <-  gen+1
      print (gen)
    }
   
    # Has the inversion reached threshold frequency for establishment (pcrit)? 
    # When did it first reach pcrit?
    if(any(InvFreq >= pcrit)) { # inversion got established
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
      "InvFreq_f"   =  InvFreq_f(1:gen-1),
      "InvFreq_m"   =  InvFreq_m(1:gen-1),
      "E.InvFreq_f" =  E.InvFreq_f(1:gen-1),
      "E.InvFreq_m" =  E.InvFreq_m(1:gen-1),
      "InvEst"      =  invEst,
      "InvEstTime"  =  invEstTime,
      "W.mean"      =  W.mean[1:gen-1],
      "Wf.mean"     =  Wf.mean(1:gen-1),
      "Wm.mean"     =  Wm.mean(1:gen-1),
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
    Wf.mean    <-  rep(0, times=(4*N+1))
    Wm.mean    <-  rep(0, times=(4*N+1))
        
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
        xi[j]   <-  round(recFctx(Fii = Fiix, m = mf, r = r), digits=8)
        recFcty  <-  get(names(yi)[j])
        yi[j]   <-  round(recFcty(Fii = Fiiy, m = mm, r = r), digits=8)
      }
      # 2) Offspring genotype frequencies
      O_females     <-  offFreq_females(xi, yi)
      O_males       <-  offFreq_males(xi)
      # 3) Mean fitness 
      Wfbar   <-  sum(O_females*Wf)
      Wmbar   <-  sum(O_males*Wm)
      # 4) Expected frequencies
      E.Fiix  <-  O*Wf/Wfbar
      E.Fiiy   <-  O*Wm/Wmbar
      # 5) Draw random frequencies in adults
      Fiix    <-  as.vector(rmultinom(1, N/2, E.Fiix)/(N/2)) # N/2 also applies for the X chromosome?
      Fiiy    <-  as.vector(rmultinom(1, N/2, E.Fiiy)/(N/2))
      
      # Realized frequencies
      InvFreq    <-  sum(Fii[c(5,10,15,20:24)]/2, Fii[25])
      E.InvFreq  <-  sum(E.Fii[c(5,10,15,20:24)]/2, Fii[25])
      W.mean     <-  Wbar
      
      #The variable below are stored specifically for the X linked simulation
      Wf.mean[gen+1]   =  sum(O_females*Wf)
      Wm.mean[gen+1]   =  sum(O_males*Wm)
      InvFreq_f[gen+1]  =  sum(Fiix[c(5,10,15,20:24)]/2, Fiix[25])
      InvFreq_m[gen+1]  =  Fiiy[5]
      E.InvFreq_f[gen+1]=  sum(E.Fiix[c(5,10,15,20:24)]/2, E.Fiix[25])
      E.InvFreq_m[gen+1]= Fiiy[5]
      
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
      "InvFreq"     =  InvFreq[1:gen-1],
      "E.InvFreq"   =  E.InvFreq[1:gen-1],
      "InvFreq_f"   =  InvFreq_f(1:gen-1),
      "InvFreq_m"   =  InvFreq_m(1:gen-1),
      "E.InvFreq_f" =  E.InvFreq_f(1:gen-1),
      "E.InvFreq_m" =  E.InvFreq_m(1:gen-1),
      "InvEst"      =  invEst,
      "InvEstTime"  =  invEstTime,
      "Wf.mean"     =  Wf.mean(1:gen-1),
      "Wm.mean"     =  Wm.mean(1:gen-1),
      "W.mean"      =  W.mean[1:gen-1],
      "nGen"        =  gen
    )
  }
  
  return(res)
}

#' Introduce new mutant inversion genotype
#'
#' @title Introduce new mutant inversion genotype
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#'   				 chosen randomly, given initial genotypic frequencies (Fii.init). 
#' 					 See params for runReplicateAutoInvSims().
#' @param Fii.init   Initial genotypic frequencies, from which to calculate probability of new 
#' 					 mutant inversion occurring
#' @param N			 Population size

# add m, which is true if the inversion is introduced into the male genotype.
introduceInversion  <-  function(newMutant, m = FALSE, Fiix.init, Fiiy.init, N) {
#' Females:
#'    c(ABAB,  ABAb,  ABaB,  ABab,  ABba*,  c(x1x1, x1x2, x1x3, x1x4, x1x5, 
#'		  AbAB,  AbAb,  AbaB,  Abab,  Abba*,	  x2x1, x2x2, x2x3, x2x4, x2x5, 
#'		  aBAB,  aBAb,  aBaB,  aBab,  aBba*,	  x3x1, x3x2, x3x3, x3x4, x3x5, 
#'		  abAB,  abAb,  abaB,  abab,  abba*,	  x4x1, x4x2, x4x3, x4x4, x4x5, 
#'		  baAB*, baAb*, baaB*, baab*, baba*)	  x5x1, x5x2, x5x3, x5x4, x5x5)
#'
#'Males:
#'    c(AB, Ab, aB, ab, ba*)                c(y1, y2, y3, y4, y5)
  
  # Toggle
  specifyNewMutant  <-  is.numeric(newMutant)
  
  # Choose mutant genotype randomly
  if(!specifyNewMutant) {
    
      # Probability of new mutant occuring on X in females
      probNewMutantX     <-  c(Fiix.init[c(4,9,14,16:18)], Fiix.init[19]*2, Fiiy.init[4])/sum(Fiix.init[c(4,9,14,16:18)], Fiix.init[19]*2 , sum(Fiiy.init[4]))
      newMutX            <-  c(4,9,14,16:19,100)[as.vector(rmultinom(1,1,probNewMutantX)) == 1]
      if(newMutX == 100) {
        # Subtract new mutant individual from frequency of old genotype
        Fiiy.init[4]  <-  Fiiy.init[4] - 1/N
      }
      else {
        Fiix.init[newMutX]  <-  Fiix.init[newMutX] - 1/N
      }
    }
  
  # Specify mutant genotype
  if(specifyNewMutant) {
    # Subtract new mutant individual from frequency of old genotype
    if(m) {
      newMutY  <-  newMutant
      Fiiy.init[newMutY]  <-  Fiiy.init[newMutY] - 1/N
    }
    else {
      newMutX  <-  newMutant
      Fiix.init[newMutX]  <-  Fiix.init[newMutX] - 1/N
    }
  }

  # Add mutant individual to frequency of new inversion genotype
  if(newMutX == 4 | newMutX == 9 | newMutX == 14)
    Fiix.init[newMutX + 1]  <-  1/N # add up mutations?
  if(newMutX == 16 | newMutX == 17 | newMutX == 18)
    Fiix.init[newMutX + 5]  <-  1/N
  
  # For males
  Fiiy.init[5] <- 1/N # add up mutations?
  
  # if inversion occurs on abab genotype, choose randomly whether it occurs on
  # the maternally or paternally inherited chromosome 
  # for an X linked inversion, it is more probable to have occurred on X in females
  # than X in males. 
  if(newMutX == 19) {
    if(runif(1) >= 1/2) {
      Fiix.init[newMutX + 1]  <-  1/N
    }
    else Fiix.init[newMutX + 5]  <-  1/N
  } 
  return (list(Fiix.init, Fiiy.init))
}


#' Wrapper function to run replicate forward simulations for invasion
#' of X-linked inversions in a Wright-Fisher population 
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
runReplicateAutoInvSims  <-  function(nReps = 1000, N = 500, mm = 0.01, mf = 0.01, sf = 0.1, sm = 0.1,  h = 1/2, r = 0.1, 
                                      n = 100, u = 1e-5, h.del = 0, s.del = 1, noDel = FALSE,
                                      saveTrajectories = FALSE, newMutant=c("random", "random")) {
  
  ##  Preemptive Warnings
  if(any(c(N,mm,mf,sf,sm,h,r,n,u,h.del) < 0) | any(c(mm,mf,sf,sm,h,r,u,h.del) > 1) | r > 0.5)
    stop('The chosen parameter values fall outside of reasonable parameter space')
  
  specifyNewMutant <- is.numeric(newMutant)
  if(specifyNewMutant & all(newMutant != c(4,9,14,16:19)))
    stop('If specifying the genotype of new inversion mutants, newMutant must take one of the following values: 4,9,14,16:19')
  
  if(!specifyNewMutant & newMutant != 'random')
    stop('If the genotype of new inversion mutants is being chosen randomly, the parameter newMutant must equal random')
  
  try({
    if(m >= s )
      stop('Warning: migration is stronger than than selection, 
           adaptive alleles will be swamped by maladaptive migrants')
  }, silent=FALSE)
  
  try({
    if(nReps > 1000 & saveTrajectories)
      stop('Warning: You have chosen to save evolutionary trajectories 
           for a large number of replicate simulations. This will be
           memory intensive. Consider setting saveTrajectories = FALSE')
  }, silent=FALSE)
  
  ##  Define Fitness Expressions for females for determining eq. frequencies in absence of inversion
  Wf.init  <-  c(1,          (1 + h*sf),         (1 + h*sf),         (1 + h*sf)^2,        0,
                (1 + h*sf),   (1 + sf),           (1 + h*sf)^2,       (1 + h*sf)*(1 + sf),  0,
                (1 + h*sf),   (1 + h*sf)^2,       (1 + sf),           (1 + sf)*(1 + h*sf),  0,
                (1 + h*sf)^2, (1 + h*sf)*(1 + sf), (1 + sf)*(1 + h*s), (1 + sf)^2,          0,
                0,           0,                 0,                 0,                  0)
  ## Define Fitness Expressions for males for determining eq. frequencies in absence of inversion
  Wm.init  <-  c(1, (1 + sm), (1 + sm), (1 + sm)^2, 0)
  
  ## Find deterministic equilibrium frequencies in absence of inversion in females and males
  Fiix.init <- findEqFreqs(Wf=Wf.init, Wm=Wm.init, mm=mm, mf=mf, r=r, threshold = 1e-7)[1] # First row for females  
  Fiiy.init <- findEqFreqs(Wf=Wf.init, Wm=Wm.init, mm=mm, mf=mf, r=r, threshold = 1e-7)[2] # Second row for males
  
  # Introduce rare inversion mutations
  Fiix.init <- introduceInversion(newMutant=newMutant, Fiix.init=Fiix.init, N = N)[1]
  Fiiy.init <- introduceInversion(newMutant=newMutant, Fiiy.init=Fiiy.init, N = N)[2]
    
  # Use deterministic equilibrium frequencies of non-inversion genotypes
  # as initial conditions when introducing the inversion via a single
  # copy of the abba* genotype
  # Females
  # Fiix.init[19]  <-  Fiix.init[19] - 1/N
  # Fiix.init[20]  <-  1/N
  
  # Males 
  # Fiiy.init[4] <- Fiiy.init[4] - 1/N
  # Fiiy.init[5] <- 1/N
  
  # Storage structures for replicate simulation data
  finalInvFreq    <-  rep(0, times=nReps)
  finalE.InvFreq  <-  rep(0, times=nReps)
  finalW.mean     <-  rep(0, times=nReps)
  finalWf.mean     <-  rep(0, times=nReps)
  finalWm.mean     <-  rep(0, times=nReps)
  finalInvFreq_f   <-  rep(0, times=nReps)
  finalInvFreq_m   <-  rep(0, times=nReps)
  finalE.InvFreq_f   <-  rep(0, times=nReps)
  finalE.InvFreq_m   <-  rep(0, times=nReps)
  nGen            <-  rep(0, times=nReps)
  invEst          <-  rep(0, times=nReps)
  invEstTime      <-  rep(0, times=nReps)
  nDels           <-  rep(0, times=nReps)
  
  if(saveTrajectories) {
    replicateTraj  <-  c()
    InvFreqTraj    <-  c()
    E.InvFreqTraj  <-  c()
    W.meanTraj     <-  c()
    Wf.meanTraj     <-  c()
    Wm.meanTraj     <-  c()
    InvFreq_fTraj   <-  c()
    InvFreq_mTraj   <-  c()
    E.InvFreq_fTraj  <- c()
    E.InvFreq_mTraj  <- c()
  } 
  
  # Replicate simulation loop
  #	print('Running Wright-Fisher Forward Simulations')
  pb   <-  txtProgressBar(min=0, max=nReps, style=3)
  setTxtProgressBar(pb, 0)
  for(i in 1:nReps) {
    
    ## Sample stationary distribution of deleterious alleles
    delMutFreq  <-  rejectionSamplerX(n=n, Ne=N, u=u)
    n.del       <-  sum(delMutFreq > runif(n=n))
    
    # Define fitness expressions, including fitness effects of deleterious mutations for females and males
    if(noDel) {
      s.del  <-  0
    }
    Wf  <-  c(1,                                  (1 + h*sf),                                 (1 + h*sf),                                 (1 + h*sf)^2,                       (1 + h*sf)^2*(1 - h.del*sf.del)^n.del,
             (1 + h*sf),                           (1 + sf),                                   (1 + h*sf)^2,                               (1 + h*sf)*(1 + sf),                 (1 + h*sf)*(1 + sf)*(1 - h.del*sf.del)^n.del,
             (1 + h*sf),                           (1 + h*sf)^2,                               (1 + sf),                                   (1 + sf)*(1 + h*sf),                 (1 + sf)*(1 + h*sf)*(1 - h.del*sf.del)^n.del,
             (1 + h*sf)^2,                         (1 + h*sf)*(1 + sf),                         (1 + sf)*(1 + h*sf),                         (1 + sf)^2,                         (1 + sf)^2*(1 - h.del*sf.del)^n.del,
             (1 + h*sf)^2*(1 - h.del*sf.del)^n.del, (1 + h*sf)*(1 + sf)*(1 - h.del*sf.del)^n.del, (1 + sf)*(1 + h*sf)*(1 - h.del*sf.del)^n.del, (1 + sf)^2*(1 - h.del*sf.del)^n.del, (1 + sf)^2*(1 - sf.del)^n.del)
    
    Wm  <-  c(1, (1 + sm), (1 + sm), (1 + sm)^2, (1 + sm)^2*(1 - sm.del)^n.del)

    ## RUN SIMULATION
    repRes  <-  autoInvFwdSim(Fiix.init=Fiix.init, Fiiy.init=Fiiy.init, N=N, Wf=Wf, Wm=Wm, mm=mm, mf=mf, r=r)
    
    # save results for each replicate
    finalInvFreq[i]    <-  repRes$InvFreq[length(repRes$InvFreq)]
    finalE.InvFreq[i]  <-  repRes$E.InvFreq[length(repRes$E.InvFreq)]
    finalInvFreq_f[i]  <-  repRes$InvFreq_f[length(repRes$InvFreq_f)]
    finalInvFreq_m[i]  <-  repRes$InvFreq_m[length(repRes$InvFreq_m)]
    finalE.InvFreq_f[i]  <- repRes$E.InvFreq_f[length(repRes$E.InvFreq_f)]
    finalE.InvFreq_m[i]  <- repRes$E.InvFreq_m[length(repRes$E.InvFreq_m)]
    finalW.mean[i]     <-  repRes$W.mean
    finalWf.mean[i]     <-  repRes$Wf.mean
    finalWm.mean[i]     <-  repRes$Wm.mean
    nGen[i]            <-  repRes$nGen
    invEst[i]          <-  repRes$invEst
    invEstTime[i]      <-  repRes$invEstTime
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
    "finalInvFreq_m"  =  finalInvFreq_m,
    "finalE.InvFreq_f" = finalE.InvFreq_f,
    "finalE.InvFreq_m" = finalE.InvFreq_m,
    "finalW.mean"      = finalW.mean  
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
      "InvFreq_mTraj"  = InvFreq_mTraj,
      "E.InvFreq_fTraj" = E.InvFreq_fTraj,
      "E.InvFreq_mTraj" = E.InvFreq_mTraj,
      "W.meanTraj"      = W.meanTraj,
      "Wf.meanTraj"     =  Wf.meanTraj
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
#****************************************************************************************************************    
# 27 Feb - Check the functions above and modify the code below    

#' Wrapper function to run replicate forward simulations for invasion
#' of X-linked inversions in a Wright-Fisher population 
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
                                           sf = 0.1, sm = 0.1, h = 1/2, r = 0.1, 
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
        res  <-  runReplicateAutoInvSims(nReps = nReps, N = N.vals[j], m = m.vals[k], sf = sf, sm = sm, h = h, r = r, 
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
  colnames(data)  <-  c("finalInvFreq","finalE.InvFreq","finalWf.mean", "finalWm.mean",
                        "nGen","invEst","invEstTime","nDels","N","m",
                        "s.dels","sf", "m", "h","r","u","h.del")
  
  # create file name
  filename  <-  paste("./output/data/simResults/auto-InvSimsData", "_s", s, "_h", h, "_r", r, "_n", n, "_u", u, ".csv", sep="")
  
  # export data as .csv to ./output/data
  write.csv(data, file=filename, row.names = FALSE)
  
  #  Return results in case user wants it
  return(data)
  
}    
