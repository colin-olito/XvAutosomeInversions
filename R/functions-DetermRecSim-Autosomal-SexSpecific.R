#####################################################
#  Functions for deterministic simulation of haplotype
#  frequency recursions: 
#
#  Autosomal Inversions w/ sex-specific selection
#
#
#  Author: Colin Olito
#
#  NOTES:  
#          - Two Loci *A* and *B*, with alleles A,a and B,b.
#          - Alleles a,b locally adaptive
#          - Continent --> Island migration scheme at rate m
#          - Continent fixed for locally maladaptive alleles (A,B)
#          - Fitness parameterized as:
#              wf_AA = w_BB = 1
#              wf_Aa = w_Bb = 1 + hf sf
#              wf_aa = w_bb = 1 + sf
#
#              wm_AA = w_BB = 1
#              wm_Aa = w_Bb = 1 + hm sm
#              wm_aa = w_bb = 1 + sm
#          

###################################
##  Key to xi haplotype subscripts
#  Female gametes
#  x1  =  xi[1]  =  AABB
#  x2  =  xi[2]  =  AABb
#  x3  =  xi[3]  =  AaBB
#  x4  =  xi[4]  =  AaBb
#  x5  =  xi[5]  =  AbBa *Inv
#
#  Male Gametes
#  y1  =  yi[1]  =  AABB
#  y2  =  yi[2]  =  AABb
#  y3  =  yi[3]  =  AaBB
#  y4  =  yi[4]  =  AaBb
#  y5  =  yi[5]  =  AbBa *Inv

#######################################################
## Necessary Functions
#######################################################

# Average Fitness in Females
wfbar  <-  function(xi, yi, Wf.mat, ...) { 
	xi[1]*yi[1]*Wf.mat[1,1] + (xi[1]*yi[2] + xi[2]*yi[1])*(Wf.mat[1,2]) + (xi[1]*yi[3] + xi[3]*yi[1])*(Wf.mat[1,3]) + (xi[1]*yi[4] + xi[4]*yi[1])*(Wf.mat[1,4]) + (xi[1]*yi[5] + xi[5]*yi[1])*(Wf.mat[1,5]) + 
	xi[2]*yi[2]*Wf.mat[2,2] + (xi[2]*yi[3] + xi[3]*yi[2])*(Wf.mat[2,3]) + (xi[2]*yi[4] + xi[4]*yi[2])*(Wf.mat[2,4]) + (xi[2]*yi[5] + xi[5]*yi[2])*(Wf.mat[2,5]) + 
	xi[3]*yi[3]*Wf.mat[3,3] + (xi[3]*yi[4] + xi[4]*yi[3])*(Wf.mat[3,4]) + (xi[3]*yi[5] + xi[5]*yi[3])*(Wf.mat[3,5]) +  
	xi[4]*yi[4]*Wf.mat[4,4] + (xi[4]*yi[5] + xi[5]*yi[4])*(Wf.mat[4,5]) + 
	xi[5]*yi[5]*Wf.mat[5,5]
}

# Average Fitness in Males
wmbar  <-  function(xi, yi, Wm.mat, ...) { 
	yi[1]*yi[1]*Wm.mat[1,1] + (yi[1]*yi[2] + yi[2]*yi[1])*(Wm.mat[1,2]) + (yi[1]*yi[3] + yi[3]*yi[1])*(Wm.mat[1,3]) + (yi[1]*yi[4] + yi[4]*yi[1])*(Wm.mat[1,4]) + (yi[1]*yi[5] + yi[5]*yi[1])*(Wm.mat[1,5]) + 
	yi[2]*yi[2]*Wm.mat[2,2] + (yi[2]*yi[3] + yi[3]*yi[2])*(Wm.mat[2,3]) + (yi[2]*yi[4] + yi[4]*yi[2])*(Wm.mat[2,4]) + (yi[2]*yi[5] + yi[5]*yi[2])*(Wm.mat[2,5]) + 
	yi[3]*yi[3]*Wm.mat[3,3] + (yi[3]*yi[4] + yi[4]*yi[3])*(Wm.mat[3,4]) + (yi[3]*yi[5] + yi[5]*yi[3])*(Wm.mat[3,5]) +  
	yi[4]*yi[4]*Wm.mat[4,4] + (yi[4]*yi[5] + yi[5]*yi[4])*(Wm.mat[4,5]) + 
	yi[5]*yi[5]*Wm.mat[5,5]
}

# Linkage Disequilibrium in females
Dfstar  <-  function(xi, yi, Wf.mat, par.list, ...) {
	(((xi[1]*yi[4] + xi[4]*yi[1])*(Wf.mat[1,4]) - (xi[2]*yi[3] + xi[3]*yi[2])*(Wf.mat[2,3])) / (2*wfbar(xi=xi, yi=yi, Wf.mat=Wf.mat)))*(1 - par.list$mf)
}

# Linkage Disequilibrium in males
Dmstar  <-  function(xi, yi, Wm.mat, par.list, ...) {
	(((yi[1]*yi[4] + yi[4]*yi[1])*(Wm.mat[1,4]) - (yi[2]*yi[3] + yi[3]*yi[2])*(Wm.mat[2,3])) / (2*wmbar(xi=xi, yi=yi, Wm.mat=Wm.mat)))*(1 - par.list$mm)
}


# Haplotype Recursions
xPr.1  <-  function(xi = xi, yi = yi, Wf.mat = Wf.mat, par.list = par.list, ...) {
	((2*(xi[1]*yi[1]*(Wf.mat[1,1])) + (xi[1]*yi[2] + xi[2]*yi[1])*(Wf.mat[1,2]) + (xi[1]*yi[3] + xi[3]*yi[1])*(Wf.mat[1,3]) + (xi[1]*yi[4] + xi[4]*yi[1])*(Wf.mat[1,4]) + (xi[1]*yi[5] + xi[5]*yi[1])*(Wf.mat[1,5])) / (2*wfbar(xi=xi, yi=yi, Wf.mat=Wf.mat)))*(1 - par.list$mf) + par.list$mf - par.list$r*Dfstar(xi=xi, yi=yi, Wf.mat=Wf.mat, par.list=par.list)
} 
xPr.2  <-  function(xi = xi, yi = yi, Wf.mat = Wf.mat, par.list = par.list, ...) {
	((2*(xi[2]*yi[2]*(Wf.mat[2,2])) + (xi[1]*yi[2] + xi[2]*yi[1])*(Wf.mat[1,2]) + (xi[2]*yi[3] + xi[3]*yi[2])*(Wf.mat[2,3]) + (xi[2]*yi[4] + xi[4]*yi[2])*(Wf.mat[2,4]) + (xi[2]*yi[5] + xi[5]*yi[2])*(Wf.mat[2,5])) / (2*wfbar(xi=xi, yi=yi, Wf.mat=Wf.mat)))*(1 - par.list$mf) + par.list$r*Dfstar(xi=xi, yi=yi, Wf.mat=Wf.mat, par.list=par.list)
}
xPr.3  <-  function(xi = xi, yi = yi, Wf.mat = Wf.mat, par.list = par.list, ...) {
	((2*(xi[3]*yi[3]*(Wf.mat[3,3])) + (xi[1]*yi[3] + xi[3]*yi[1])*(Wf.mat[1,3]) + (xi[2]*yi[3] + xi[3]*yi[2])*(Wf.mat[2,3]) + (xi[3]*yi[4] + xi[4]*yi[3])*(Wf.mat[3,4]) + (xi[3]*yi[5] + xi[5]*yi[3])*(Wf.mat[3,5])) / (2*wfbar(xi=xi, yi=yi, Wf.mat=Wf.mat)))*(1 - par.list$mf) + par.list$r*Dfstar(xi=xi, yi=yi, Wf.mat=Wf.mat, par.list=par.list)
}
xPr.4  <-  function(xi = xi, yi = yi, Wf.mat = Wf.mat, par.list = par.list, ...) {
	((2*(xi[4]*yi[4]*(Wf.mat[4,4])) + (xi[1]*yi[4] + xi[4]*yi[1])*(Wf.mat[1,4]) + (xi[2]*yi[4] + xi[4]*yi[2])*(Wf.mat[2,4]) + (xi[3]*yi[4] + xi[4]*yi[3])*(Wf.mat[3,4]) + (xi[4]*yi[5] + xi[5]*yi[4])*(Wf.mat[4,5])) / (2*wfbar(xi=xi, yi=yi, Wf.mat=Wf.mat)))*(1 - par.list$mf) - par.list$r*Dfstar(xi=xi, yi=yi, Wf.mat=Wf.mat, par.list=par.list)
}
xPr.5  <-  function(xi = xi, yi = yi, Wf.mat = Wf.mat, par.list = par.list, ...) {
	((2*(xi[5]*yi[5]*(Wf.mat[5,5])) + (xi[1]*yi[5] + xi[5]*yi[1])*(Wf.mat[1,5]) + (xi[2]*yi[5] + xi[5]*yi[2])*(Wf.mat[2,5]) + (xi[3]*yi[5] + xi[5]*yi[3])*(Wf.mat[3,5]) + (xi[4]*yi[5] + xi[5]*yi[4])*(Wf.mat[4,5])) / (2*wfbar(xi=xi, yi=yi, Wf.mat=Wf.mat)))*(1 - par.list$mf)
}



yPr.1  <-  function(xi = xi, yi = yi, Wm.mat = Wm.mat, par.list = par.list, ...) {
	((2*(xi[1]*yi[1]*(Wm.mat[1,1])) + (xi[1]*yi[2] + xi[2]*yi[1])*(Wm.mat[1,2]) + (xi[1]*yi[3] + xi[3]*yi[1])*(Wm.mat[1,3]) + (xi[1]*yi[4] + xi[4]*yi[1])*(Wm.mat[1,4]) + (xi[1]*yi[5] + xi[5]*yi[1])*(Wm.mat[1,5])) / (2*wmbar(xi=xi, yi=yi, Wm.mat=Wm.mat)))*(1 - par.list$mm) + par.list$mm - par.list$r*Dmstar(xi=xi, yi=yi, Wm.mat=Wm.mat, par.list=par.list)
} 
yPr.2  <-  function(xi = xi, yi = yi, Wm.mat = Wm.mat, par.list = par.list, ...) {
	((2*(xi[2]*yi[2]*(Wm.mat[2,2])) + (xi[1]*yi[2] + xi[2]*yi[1])*(Wm.mat[1,2]) + (xi[2]*yi[3] + xi[3]*yi[2])*(Wm.mat[2,3]) + (xi[2]*yi[4] + xi[4]*yi[2])*(Wm.mat[2,4]) + (xi[2]*yi[5] + xi[5]*yi[2])*(Wm.mat[2,5])) / (2*wmbar(xi=xi, yi=yi, Wm.mat=Wm.mat)))*(1 - par.list$mm) + par.list$r*Dmstar(xi=xi, yi=yi, Wm.mat=Wm.mat, par.list=par.list)
}
yPr.3  <-  function(xi = xi, yi = yi, Wm.mat = Wm.mat, par.list = par.list, ...) {
	((2*(xi[3]*yi[3]*(Wm.mat[3,3])) + (xi[1]*yi[3] + xi[3]*yi[1])*(Wm.mat[1,3]) + (xi[2]*yi[3] + xi[3]*yi[2])*(Wm.mat[2,3]) + (xi[3]*yi[4] + xi[4]*yi[3])*(Wm.mat[3,4]) + (xi[3]*yi[5] + xi[5]*yi[3])*(Wm.mat[3,5])) / (2*wmbar(xi=xi, yi=yi, Wm.mat=Wm.mat)))*(1 - par.list$mm) + par.list$r*Dmstar(xi=xi, yi=yi, Wm.mat=Wm.mat, par.list=par.list)
}
yPr.4  <-  function(xi = xi, yi = yi, Wm.mat = Wm.mat, par.list = par.list, ...) {
	((2*(xi[4]*yi[4]*(Wm.mat[4,4])) + (xi[1]*yi[4] + xi[4]*yi[1])*(Wm.mat[1,4]) + (xi[2]*yi[4] + xi[4]*yi[2])*(Wm.mat[2,4]) + (xi[3]*yi[4] + xi[4]*yi[3])*(Wm.mat[3,4]) + (xi[4]*yi[5] + xi[5]*yi[4])*(Wm.mat[4,5])) / (2*wmbar(xi=xi, yi=yi, Wm.mat=Wm.mat)))*(1 - par.list$mm) - par.list$r*Dmstar(xi=xi, yi=yi, Wm.mat=Wm.mat, par.list=par.list)
}
yPr.5  <-  function(xi = xi, yi = yi, Wm.mat = Wm.mat, par.list = par.list, ...) {
	((2*(xi[5]*yi[5]*(Wm.mat[5,5])) + (xi[1]*yi[5] + xi[5]*yi[1])*(Wm.mat[1,5]) + (xi[2]*yi[5] + xi[5]*yi[2])*(Wm.mat[2,5]) + (xi[3]*yi[5] + xi[5]*yi[3])*(Wm.mat[3,5]) + (xi[4]*yi[5] + xi[5]*yi[4])*(Wm.mat[4,5])) / (2*wmbar(xi=xi, yi=yi, Wm.mat=Wm.mat)))*(1 - par.list$mm)
}



#################################
##  Simulation Wrapper Functions
#################################

#' Forward deterministic simulation of haplotype recursions for
#' 2-locus local adaptation w/ autosomal inversion
#'
#' @title Forward deterministic simulation of haplotype recursions
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen  =  25000,
#'				   mf    =  0.01,
#'				   mm    =  0.01,
#'				   sf    =  0.02,
#'				   sm    =  0.02,
#'				   hf    =  0.5,
#'				   hm    =  0.5,
#'				   r    =  0.5
#'				   )
#' @param xi.init A vector of initial haplotype frequencies among female gametes (must have length = 5).
#'                For example, c(0,m/s,m/s,1-(2*(m/s)),0) to initiate simulation at approximate 
#'                single-locus equilibrium frequencies of [A],[a],[B],[b].
#' @param yi.init A vector of initial haplotype frequencies among male gametes (must have length = 5).
#'                For example, c(0,m/s,m/s,1-(2*(m/s)),0) to initiate simulation at approximate 
#'                single-locus equilibrium frequencies of [A],[a],[B],[b].
#' @param threshold The threshold change in haplotype frequency at which the simulation
#'                  will accept the current state as having reached equilibrium. Values
#'                  of 1e-6 or 1e-7 have worked pretty well in the past. Be changing this
#'                  parameter, as it will influence both how accurate the results are, 
#'                  and how long the simulations take to run.
#' @param silent Should the simulation loop print a warning if m > s? 
#' @return Returns a list with timeseries for each haplotype, equilibrium frequencies, and a numeric (0,1) for whether the 
#' equilibrium was polymorphic (with tolerance 1E-6).
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSim(par.list, xi.init, threshold = 1e-6, verbose=FALSE) 
recursionFwdSim  <-  function(par.list, xi.init, yi.init, threshold = 1e-6, silent = FALSE, ...) {

	##  Preemptive Warnings
	if(any(par.list[1:8] < 0) | any(par.list[2:7] > 1) | par.list$r > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	try({
		 if(par.list$mf  >=  par.list$sf | par.list$mm  >=  par.list$sm | par.list$mf  >=  par.list$sm | par.list$mm  >=  par.list$sf)
			stop('Warning: migration is stronger than than selection, 
				  adaptive alleles will be swamped by maladaptive migrants')
	}, silent=silent)

	if(round(sum(xi.init), digits=4) != 1)
		stop('Initial frequencies among female gametes must sum to 1')
	if(round(sum(yi.init), digits=2) != 1)
		stop('Initial frequencies among male gametes must sum to 1')

	##  Fitness Expression Matrices
	Wf.mat   <-  matrix(
					   c( 1,                              (1 + par.list$hf*par.list$sf),                   (1 + par.list$hf*par.list$sf),                   (1 + par.list$hf*par.list$sf)^2,                 (1 + par.list$hf*par.list$sf)^2,
                         (1 + par.list$hf*par.list$sf),   (1 + par.list$sf),                               (1 + par.list$hf*par.list$sf)^2,                 (1 + par.list$hf*par.list$sf)*(1 + par.list$sf), (1 + par.list$hf*par.list$sf)*(1 + par.list$sf),
                         (1 + par.list$hf*par.list$sf),   (1 + par.list$hf*par.list$sf)^2,                 (1 + par.list$sf),                               (1 + par.list$sf)*(1 + par.list$hf*par.list$sf), (1 + par.list$sf)*(1 + par.list$hf*par.list$sf),
                         (1 + par.list$hf*par.list$sf)^2, (1 + par.list$hf*par.list$sf)*(1 + par.list$sf), (1 + par.list$sf)*(1 + par.list$hf*par.list$sf), (1 + par.list$sf)^2,                             (1 + par.list$sf)^2,
                         (1 + par.list$hf*par.list$sf)^2, (1 + par.list$hf*par.list$sf)*(1 + par.list$sf), (1 + par.list$sf)*(1 + par.list$hf*par.list$sf), (1 + par.list$sf)^2,                             (1 + par.list$sf)^2),
						 nrow=5, byrow=TRUE
					    )
	Wm.mat   <-  matrix(
					   c( 1,                              (1 + par.list$hm*par.list$sm),                   (1 + par.list$hm*par.list$sm),                   (1 + par.list$hm*par.list$sm)^2,                 (1 + par.list$hm*par.list$sm)^2,
                         (1 + par.list$hm*par.list$sm),   (1 + par.list$sm),                               (1 + par.list$hm*par.list$sm)^2,                 (1 + par.list$hm*par.list$sm)*(1 + par.list$sm), (1 + par.list$hm*par.list$sm)*(1 + par.list$sm),
                         (1 + par.list$hm*par.list$sm),   (1 + par.list$hm*par.list$sm)^2,                 (1 + par.list$sm),                               (1 + par.list$sm)*(1 + par.list$hm*par.list$sm), (1 + par.list$sm)*(1 + par.list$hm*par.list$sm),
                         (1 + par.list$hm*par.list$sm)^2, (1 + par.list$hm*par.list$sm)*(1 + par.list$sm), (1 + par.list$sm)*(1 + par.list$hm*par.list$sm), (1 + par.list$sm)^2,                             (1 + par.list$sm)^2,
                         (1 + par.list$hm*par.list$sm)^2, (1 + par.list$hm*par.list$sm)*(1 + par.list$sm), (1 + par.list$sm)*(1 + par.list$hm*par.list$sm), (1 + par.list$sm)^2,                             (1 + par.list$sm)^2),
						 nrow=5, byrow=TRUE
					    )

	##  Initilize data storage structures
	xi.gen  <-  matrix(0, ncol=5, nrow=par.list$gen)
	yi.gen  <-  matrix(0, ncol=5, nrow=par.list$gen)
	colnames(xi.gen)  <-  c('xPr.1', 'xPr.2', 'xPr.3', 'xPr.4', 'xPr.5')
	colnames(yi.gen)  <-  c('yPr.1', 'yPr.2', 'yPr.3', 'yPr.4', 'yPr.5')
	
	##  Generation Loop
		# initialize
	for (j in 1:ncol(xi.gen)) {
		recFctX          <-  get(colnames(xi.gen)[j])
		recFctY          <-  get(colnames(yi.gen)[j])
		xi.gen[1, j]     <-  round(recFctX(xi = xi.init, yi = yi.init, par.list = par.list, Wf.mat = Wf.mat), digits=8)
		yi.gen[1, j]     <-  round(recFctY(xi = xi.init, yi = yi.init, par.list = par.list, Wm.mat = Wm.mat), digits=8)
	}

	# Start simulation
	i      <-  2
	diffs  <-  rep(1,ncol(xi.gen))
	while (i < par.list$gen & any(diffs[diffs != 0] > threshold)) {
		for (j in 1:ncol(xi.gen)) {
			recFctX          <-  get(colnames(xi.gen)[j])
			recFctY          <-  get(colnames(yi.gen)[j])
			xi.gen[i, j]   <-  round(recFctX(xi = xi.gen[i-1, ], yi = yi.gen[i-1, ], par.list = par.list, Wf.mat = Wf.mat), digits=8)
			yi.gen[i, j]   <-  round(recFctY(xi = xi.gen[i-1, ], yi = yi.gen[i-1, ], par.list = par.list, Wm.mat = Wm.mat), digits=8)
		}

		diffs  <-  c(abs(xi.gen[i,] - xi.gen[i-1,]), abs(yi.gen[i,] - yi.gen[i-1,]))
		i      <-  i+1
	}
	if(i == par.list$gen) {
		print('Warning: maximum runtime reached. Results may not represent equilibrium frequencies')
	}

	##  Output
	LDx   <-  xi.gen[i-1,][1]*xi.gen[i-1,][4] - xi.gen[i-1,][2]*xi.gen[i-1,][3]
	LDy   <-  yi.gen[i-1,][1]*yi.gen[i-1,][4] - yi.gen[i-1,][2]*yi.gen[i-1,][3]
	res  <-  list(
				  "par.list"  =  par.list,
				  "xi.gen"    =  xi.gen[1:i-1,],
				  "yi.gen"    =  yi.gen[1:i-1,],
				  "Fi.gen"    =  (xi.gen[1:i-1,] + yi.gen[1:i-1,])/2,
				  "EQxy.freq" =  c(xi.gen[i-1,], yi.gen[i-1,]),
				  "EQ.freq"   =  (xi.gen[i-1,] + yi.gen[i-1,])/2,
				  "LDx"       =  LDx,
				  "LDy"       =  LDy,
				  "LD"        =  (LDx + LDy)/2
 				 )
	colnames(res$Fi.gen)  <-  c("F.1","F.2","F.3","F.4","F.5")
	names(res$EQ.freq)    <-  c("F.1","F.2","F.3","F.4","F.5")
	return(res)
}





#' Simulation loop to repeat forward deterministic simulations 
#' of haplotype recursions for different parameter values 
#' 
#' @title Forward deterministic simulation of genotypic recursions.
#' @param gen        Maximum number of generations for each simulation (as in par.list).
#' @param sf.vals    Vector of selection coefficients for females (for use in par.list).
#' @param sm.vals    Vector of selection coefficients for males (for use in par.list).
#' @param h          Dominance coefficient (for use in par.list). Character  
#' @param resolution 'by' arg for sequence of selection coefficient values. 
#'                    Determines length of innermost loop of simulation. 
#'                    Recommend 0.01 for exploratory analyses, 0.005 for plotting.
#' @param m.vals Values of migration rate to explore. Determines length of second simulation loop.
#' @param r.vals Values of recombination rate to explore(as in par.list). Determines length of outermost simulation loop.
#' @param threshold Threshold difference between genotypic frequencies before simulation cuts off.
#' @return Returns a data frame with parameter values, final frequencies of the inversion,
#'         a variable describing whether the final state of the simulation was polymorphic, 
#'         whether evaluating the eigenvalues predicts polymorphism, and whether these two 
#'         methods agree with one another.
#' @seealso `recursionFwdSim`
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSimLoop(n = 5000, gen = 5000, sRange = c(0,1), C = 0, delta = 0, 
#'                     hf = 0.5, hm = 0.5, r.vals = c(0.0, 0.01, 0.02, 0.1, 0.2, 0.5), 
#'                     seed = 3497016, threshold = 1e-7)
recursionFwdSimLoop  <-  function(gen = 25000, hf = 0.5, hm = 0.5, threshold = 1e-7,
								  sf.vals = c(0.01, 0.05), sm.vals = c(0.01, 0.05),
								  mf.vals = c(0.01, 0.05), mm.vals = c(0.01, 0.05),
								  r.vals = c(0.0, 0.01, 0.1)) {

	## Warnings
	if(length(sf.vals)  !=  length(sm.vals))
		stop('Vectors of male and female selection coefficients must be of equal length')

	if(length(mf.vals)  !=  length(mm.vals))
		stop('Vectors of male and female migration rates must be of equal length')

	if(any(c(gen,hf,hm,sf.vals,sm.vals,mf.vals,mf.vals,r.vals) < 0) | any(c(hf,hm) > 1) | any(r.vals > 0.5) | any(c(sf.vals, sm.vals,mf.vals,mm.vals) > 0.2))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(threshold > 1e-6)
		stop('Carefully consider whether you want to change this threshold, 
			  as it will affect how many generations are required to reach
			  convergence, and thus how long the simulations take')

	# initialize storage structures
	eqxyFreqs  <-  matrix(0, nrow=length(r.vals)*length(mf.vals)*length(sf.vals), ncol=10)
	LDx        <-  rep(0, length(r.vals)*length(mf.vals)*length(sf.vals))
	LDy        <-  rep(0, length(r.vals)*length(mf.vals)*length(sf.vals))
	
	##  Simulation Loop over values of r, m, s 
	print('Running Deterministic Recursion Simulations')
	pb   <-  txtProgressBar(min=0, max=length(r.vals)*length(mf.vals), style=3)
	setTxtProgressBar(pb, 0)

	for (i in 1:length(r.vals)) {
		for (j in 1:length(mf.vals)) {
			for (k in 1:length(sf.vals)) {
				
				par.list  <-  list(
				   				gen  =  gen,
 				   				mf   =  mf.vals[j],
 				   				mm   =  mm.vals[j],
 				   				sf   =  sf.vals[k],
 				   				sm   =  sm.vals[k],
 				   				hf   =  hf,
 				   				hm   =  hm,
 				   				r    =  r.vals[i]
 				   				)

				##  Calculate equilibrium frequencies prior to invasion of the inversion
				inits       <-  recursionFwdSim(par.list=par.list, xi.init = c(0.25,0.25,0.25,0.25,0), yi.init = c(0.25,0.25,0.25,0.25,0), threshold = threshold, silent=TRUE)
				x.inits     <-  inits[[5]][1:5]
				x.inits[x.inits == max(x.inits)]  <-  x.inits[x.inits == max(x.inits)] - 0.01
				x.inits[5]  <-  0.01
				x.inits
				y.inits     <-  inits[[5]][6:10]
				y.inits[y.inits == max(y.inits)]  <-  y.inits[y.inits == max(y.inits)] - 0.01
				y.inits[5]  <-  0.01
				y.inits

		 		# Run simulation for given parameter values
				res  <-  recursionFwdSim(par.list = par.list, xi.init = x.inits, yi.init = y.inits, threshold = threshold, silent=TRUE)

				# Store equilibrium frequencies
				eqxyFreqs[((i-1)*length(mf.vals)*length(sf.vals)) + ((j-1)*length(sf.vals)) + k,]  <-  res$EQxy.freq
				LDx[((i-1)*length(mf.vals)*length(sf.vals))    + ((j-1)*length(sf.vals)) + k]      <-  res$LDx
				LDy[((i-1)*length(mf.vals)*length(sf.vals))    + ((j-1)*length(sf.vals)) + k]      <-  res$LDy
			}
		setTxtProgressBar(pb, ((i-1)*length(mf.vals) + j))
		}
	}
	setTxtProgressBar(pb, length(r.vals)*length(mf.vals))

	#  Compile results as data.frame
	rs   <-  c()
	mfs  <-  c()
	mms  <-  c()
	for(i in 1:length(r.vals)) {
		rs  <-  c(rs, rep(r.vals[i], length(mf.vals)*length(sf.vals)))
	}
	for(i in 1:length(r.vals)) {
		for(j in 1:length(mf.vals)) {
			mfs  <-  c(mfs, rep(mf.vals[j], length(sf.vals)))
			mms  <-  c(mms, rep(mm.vals[j], length(sf.vals)))
		}
	}
	results.df  <-  data.frame("hf"  =  rep(hf, length(r.vals)*length(mf.vals)*length(sf.vals)),
							   "hm"  =  rep(hm, length(r.vals)*length(mf.vals)*length(sf.vals)),
							   "sf"  =  rep(sf.vals, length(r.vals)*length(mf.vals)),
							   "sm"  =  rep(sm.vals, length(r.vals)*length(mf.vals)),
							   "r"   =  rs,
							   "mf"  =  mfs,
							   "mm"  =  mms
							   )
	results.df  <-  cbind(results.df, eqxyFreqs, LDx, LDy)
	colnames(results.df)  <-  c('hf','hm','sf','sm','r','mf','mm', 
								'x1', 'x2', 'x3', 'x4', 'x5', 'y1', 'y2', 'y3', 'y4', 'y5', 'LDx','LDy')
	
	#  Write results.df to .txt file
	filename  <-  paste("./output/data/simResults/determRecSim-SexSpecific", "_hf", hf, "_hm", hm, ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}


