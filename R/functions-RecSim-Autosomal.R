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
#' For example, c(0,m/s,m/s,1-(2*(m/s)),0) to initiate simulation at approximate 
#' single-locus equilibrium frequencies of [A],[a],[B],[b].
#' @return Returns a list with timeseries for each haplotype, equilibrium frequencies, and a numeric (0,1) for whether the 
#' equilibrium was polymorphic (with tolerance 1E-6).
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSim(par.list, xi.init, threshold = 1e-6) 
recursionFwdSim  <-  function(par.list, xi.init, threshold = 1e-6, ...) {

	##  Preemptive Warnings
	if(any(par.list[1:5] < 0) | any(par.list[2:5] > 1) | par.list$r > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	try({
		 if(par.list$m  >=  par.list$s)
			stop('Warning: migration rate is larger than the seleciton coefficient, \n
				  adaptive alleles will be swamped by maladaptive migrants')
	})

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
	res  <-  list(
				  "par.list" =  par.list,
				  "xi.gen"  =  xi.gen[1:i-1,],
				  "EQ.freq"  =  xi.gen[i-1,]
 				 )
	return(res)
}







