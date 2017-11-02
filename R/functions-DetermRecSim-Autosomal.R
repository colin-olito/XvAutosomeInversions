#####################################################
#  Functions for deterministic simulation of haplotype
#  frequency recursions: Autosomal Inversions 
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
#              w_AA = w_BB = 1
#              w_Aa = w_Bb = 1 + h s
#              w_aa = w_bb = 1 + s
#          

###################################
##  Key to xi haplotype subscripts
#  x1  =  xi[1]  =  AABB
#  x2  =  xi[2]  =  AABb
#  x3  =  xi[3]  =  AaBB
#  x4  =  xi[4]  =  AaBb
#  x5  =  xi[5]  =  AbBa *Inv

#######################################################
## Necessary Functions
#######################################################

# Average Fitness
wbar  <-  function(xi, W.mat, ...) { 
	xi[1]*xi[1]*W.mat[1,1] + (xi[1]*xi[2] + xi[2]*xi[1])*(W.mat[1,2]) + (xi[1]*xi[3] + xi[3]*xi[1])*(W.mat[1,3]) + (xi[1]*xi[4] + xi[4]*xi[1])*(W.mat[1,4]) + (xi[1]*xi[5] + xi[5]*xi[1])*(W.mat[1,5]) + 
	xi[2]*xi[2]*W.mat[2,2] + (xi[2]*xi[3] + xi[3]*xi[2])*(W.mat[2,3]) + (xi[2]*xi[4] + xi[4]*xi[2])*(W.mat[2,4]) + (xi[2]*xi[5] + xi[5]*xi[2])*(W.mat[2,5]) + 
	xi[3]*xi[3]*W.mat[3,3] + (xi[3]*xi[4] + xi[4]*xi[3])*(W.mat[3,4]) + (xi[3]*xi[5] + xi[5]*xi[3])*(W.mat[3,5]) +  
	xi[4]*xi[4]*W.mat[4,4] + (xi[4]*xi[5] + xi[5]*xi[4])*(W.mat[4,5]) + 
	xi[5]*xi[5]*W.mat[5,5]
}


# Linkage Disequilibrium
Dstar  <-  function(xi, W.mat, par.list, ...) {
	(((xi[1]*xi[4] + xi[4]*xi[1])*(W.mat[1,4]) - (xi[2]*xi[3] + xi[3]*xi[2])*(W.mat[2,3])) / (2*wbar(xi=xi, W.mat=W.mat)))*(1 - par.list$m)
}


# Haplotype Recursions
xPr.1  <-  function(xi = xi, W.mat = W.mat, par.list = par.list, ...) {
	((2*(xi[1]*xi[1]*(W.mat[1,1])) + (xi[1]*xi[2] + xi[2]*xi[1])*(W.mat[1,2]) + (xi[1]*xi[3] + xi[3]*xi[1])*(W.mat[1,3]) + (xi[1]*xi[4] + xi[4]*xi[1])*(W.mat[1,4]) + (xi[1]*xi[5] + xi[5]*xi[1])*(W.mat[1,5])) / (2*wbar(xi=xi,W.mat=W.mat)))*(1 - par.list$m) + par.list$m - par.list$r*Dstar(xi=xi, W.mat=W.mat, par.list=par.list)
} 
xPr.2  <-  function(xi = xi, W.mat = W.mat, par.list = par.list, ...) {
	((2*(xi[2]*xi[2]*(W.mat[2,2])) + (xi[1]*xi[2] + xi[2]*xi[1])*(W.mat[1,2]) + (xi[2]*xi[3] + xi[3]*xi[2])*(W.mat[2,3]) + (xi[2]*xi[4] + xi[4]*xi[2])*(W.mat[2,4]) + (xi[2]*xi[5] + xi[5]*xi[2])*(W.mat[2,5])) / (2*wbar(xi=xi,W.mat=W.mat)))*(1 - par.list$m) + par.list$r*Dstar(xi=xi, W.mat=W.mat, par.list=par.list)
}
xPr.3  <-  function(xi = xi, W.mat = W.mat, par.list = par.list, ...) {
	((2*(xi[3]*xi[3]*(W.mat[3,3])) + (xi[1]*xi[3] + xi[3]*xi[1])*(W.mat[1,3]) + (xi[2]*xi[3] + xi[3]*xi[2])*(W.mat[2,3]) + (xi[3]*xi[4] + xi[4]*xi[3])*(W.mat[3,4]) + (xi[3]*xi[5] + xi[5]*xi[3])*(W.mat[3,5])) / (2*wbar(xi=xi,W.mat=W.mat)))*(1 - par.list$m) + par.list$r*Dstar(xi=xi, W.mat=W.mat, par.list=par.list)
}
xPr.4  <-  function(xi = xi, W.mat = W.mat, par.list = par.list, ...) {
	((2*(xi[4]*xi[4]*(W.mat[4,4])) + (xi[1]*xi[4] + xi[4]*xi[1])*(W.mat[1,4]) + (xi[2]*xi[4] + xi[4]*xi[2])*(W.mat[2,4]) + (xi[3]*xi[4] + xi[4]*xi[3])*(W.mat[3,4]) + (xi[4]*xi[5] + xi[5]*xi[4])*(W.mat[4,5])) / (2*wbar(xi=xi,W.mat=W.mat)))*(1 - par.list$m) - par.list$r*Dstar(xi=xi, W.mat=W.mat, par.list=par.list)
}
xPr.5  <-  function(xi = xi, W.mat = W.mat, par.list = par.list, ...) {
	((2*(xi[5]*xi[5]*(W.mat[5,5])) + (xi[1]*xi[5] + xi[5]*xi[1])*(W.mat[1,5]) + (xi[2]*xi[5] + xi[5]*xi[2])*(W.mat[2,5]) + (xi[3]*xi[5] + xi[5]*xi[3])*(W.mat[3,5]) + (xi[4]*xi[5] + xi[5]*xi[4])*(W.mat[4,5])) / (2*wbar(xi=xi,W.mat=W.mat)))*(1 - par.list$m)
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
#'				   gen  =  5000,
#'				   m    =  0.01,
#'				   s    =  0.02,
#'				   h    =  0.5,
#'				   r    =  0.5
#'				   )
#' @param xi.init A vector of initial haplotype frequencies (must have length = 5).
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
recursionFwdSim  <-  function(par.list, xi.init, threshold = 1e-6, silent = FALSE, ...) {

	##  Preemptive Warnings
	if(any(par.list[1:5] < 0) | any(par.list[2:5] > 1) | par.list$r > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	try({
		 if(par.list$m  >=  par.list$s)
			stop('Warning: migration rate is larger than the selection coefficient,
				  adaptive alleles will be swamped by maladaptive migrants')
	}, silent=silent)

	if(round(sum(xi.init), digits=4) != 1)
		stop('Initial frequencies must sum to 1')

	##  Fitness Expression Matrix
	W.mat   <-  matrix(
					   c( 1,                            (1 + par.list$h*par.list$s),                  (1 + par.list$h*par.list$s),                  (1 + par.list$h*par.list$s)^2,                (1 + par.list$h*par.list$s)^2,
                         (1 + par.list$h*par.list$s),   (1 + par.list$s),                             (1 + par.list$h*par.list$s)^2,                (1 + par.list$h*par.list$s)*(1 + par.list$s), (1 + par.list$h*par.list$s)*(1 + par.list$s),
                         (1 + par.list$h*par.list$s),   (1 + par.list$h*par.list$s)^2,                (1 + par.list$s),                             (1 + par.list$s)*(1 + par.list$h*par.list$s), (1 + par.list$s)*(1 + par.list$h*par.list$s),
                         (1 + par.list$h*par.list$s)^2, (1 + par.list$h*par.list$s)*(1 + par.list$s), (1 + par.list$s)*(1 + par.list$h*par.list$s), (1 + par.list$s)^2,                           (1 + par.list$s)^2,
                         (1 + par.list$h*par.list$s)^2, (1 + par.list$h*par.list$s)*(1 + par.list$s), (1 + par.list$s)*(1 + par.list$h*par.list$s), (1 + par.list$s)^2,                           (1 + par.list$s)^2),
						nrow=5, byrow=TRUE
					   )

	##  Initilize data storage structures
	xi.gen  <-  matrix(0, ncol=5, nrow=par.list$gen)
	colnames(xi.gen)  <-  c('xPr.1', 'xPr.2', 'xPr.3', 'xPr.4', 'xPr.5')
	
	##  Generation Loop
		# initialize
	for (j in 1:ncol(xi.gen)) {
		recFct          <-  get(colnames(xi.gen)[j])
		xi.gen[1, j]   <-  round(recFct(xi = xi.init, par.list = par.list, W.mat = W.mat), digits=8)
	}

	# Start simulation
	i      <-  2
	diffs  <-  rep(1,ncol(xi.gen))
	while (i < par.list$gen & any(diffs[diffs != 0] > threshold)) {
		for (j in 1:ncol(xi.gen)) {
			recFct          <-  get(colnames(xi.gen)[j])
			xi.gen[i, j]   <-  round(recFct(xi = xi.gen[i-1, ], par.list = par.list, W.mat = W.mat), digits=8)
		}

		diffs  <-  abs(xi.gen[i,] - xi.gen[i-1,])
		i      <-  i+1
	}
	if(i == par.list$gen) {
		print('Warning: maximum runtime reached. Results may not represent equilibrium frequencies')
	}

	##  Output
	LD   <-  xi.gen[i-1,][1]*xi.gen[i-1,][4] - xi.gen[i-1,][2]*xi.gen[i-1,][3]
	res  <-  list(
				  "par.list"  =  par.list,
				  "xi.gen"    =  xi.gen[1:i-1,],
				  "EQ.freq"   =  xi.gen[i-1,],
				  "LD"        =  LD
 				 )
	return(res)
}





#' Simulation loop to repeate forward deterministic simulations 
#' of haplotype recursions for different parameter values 
#' 
#' @title Forward deterministic simulation of genotypic recursions.
#' @param gen        Maximum number of generations for each simulation (as in par.list).
#' @param s.vals     Selection coefficient (for use in par.list).
#' @param h          Dominance coefficient (for use in par.list). Character  
#' @param m.vals     Values of migration rate to explore. Determines length of second simulation loop.
#' @param r.vals     Values of recombination rate to explore(as in par.list). Determines length of outermost simulation loop.
#' @param threshold  Threshold difference between genotypic frequencies before simulation cuts off.
#' @return Returns a data frame with parameter values, final frequencies of the inversion,
#'         a variable describing whether the final state of the simulation was polymorphic, 
#'         whether evaluating the eigenvalues predicts polymorphism, and whether these two 
#'         methods agree with one another.
#' @seealso `recursionFwdSim`
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSimLoop  <-  function(gen = 25000, h = 0.5, resolution = 0.005,
#'	                              m.vals = c(0.01, 0.05), r.vals = c(0.0, 0.01, 0.1), 
#'	                              threshold = 1e-7) {
recursionFwdSimLoop  <-  function(gen = 25000, h = 0.5, s.vals = c(0.01, 0.05),
	                              m.vals = c(0.01, 0.05), r.vals = c(0.0, 0.01, 0.1), 
	                              threshold = 1e-7) {

	## Warnings
	if(any(c(gen,h,m.vals,r.vals) < 0) | any(c(h) > 1) | any(r.vals > 0.5) | any(m.vals >= 0.2))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(threshold > 1e-6)
		stop('Carefully consider whether you want to change this threshold, 
			  as it will affect how many generations are required to reach
			  convergence, and thus how long the simulations take')

	# initialize storage structures
	eqFreqs  <-  matrix(0, nrow=length(r.vals)*length(m.vals)*length(s.vals), ncol=5)
	LD       <-  rep(0, length(r.vals)*length(m.vals)*length(s.vals))
	
	##  Simulation Loop over values of r, m, s 
	print('Running Deterministic Recursion Simulations')
	pb   <-  txtProgressBar(min=0, max=length(r.vals)*length(m.vals), style=3)
	setTxtProgressBar(pb, 0)

	for (i in 1:length(r.vals)) {
		for (j in 1:length(m.vals)) {
			for (k in 1:length(s.vals)) {
				
                par.list  <-  list(
                				   gen  =  gen,
                				   m    =  m.vals[j],
                				   s    =  s.vals[k],
                				   h    =  h,
                				   r    =  r.vals[i]
                				   )

				##  Calculate equilibrium frequencies prior to invasion of the inversion
				xi.init     <-  recursionFwdSim(par.list, xi.init = c(0.25,0.25,0.25,0.25,0), threshold = threshold, silent=TRUE)
				xi.init     <-  xi.init[[3]]
				xi.init[xi.init == max(xi.init)]  <-  xi.init[xi.init == max(xi.init)] - 0.01
				xi.init[5]  <-  0.01

		 		# Run simulation for given parameter values
				res  <-  recursionFwdSim(par.list = par.list, xi.init = xi.init, threshold = threshold, silent=TRUE)
#if(m==3) browser()

				# Store equilibrium frequencies
				eqFreqs[((i-1)*length(m.vals)*length(s.vals)) + ((j-1)*length(s.vals)) + k,]  <-  res$EQ.freq
				LD[((i-1)*length(m.vals)*length(s.vals))    + ((j-1)*length(s.vals)) + k]     <-  res$LD
			}
		setTxtProgressBar(pb, ((i-1)*length(m.vals) + j))
		}
	}
	setTxtProgressBar(pb, length(r.vals)*length(m.vals))

	#  Compile results as data.frame
	rs  <-  c()
	ms  <-  c()
	for(i in 1:length(r.vals)) {
		rs  <-  c(rs, rep(r.vals[i], length(m.vals)*length(s.vals)))
	}
	for(i in 1:length(r.vals)) {
		for(j in 1:length(m.vals)) {
			ms  <-  c(ms, rep(m.vals[j], length(s.vals)))
		}
	}
	results.df  <-  data.frame("h"  =  rep(h, length(r.vals)*length(m.vals)*length(s.vals)),
							   "s"  =  rep(s.vals, length(r.vals)*length(m.vals)),
							   "r"  =  rs,
							   "m"  =  ms
							   )
	results.df  <-  cbind(results.df, eqFreqs, LD)
	colnames(results.df)  <-  c('h','s','r','m', 'x1', 'x2', 'x3', 'x4', 'x5', 'LD')
	
	#  Write results.df to .txt file
	filename  <-  paste("./output/data/simResults/determRecSim-Autosomal", "_h", h, ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}


