##############################################################
#  Local adaptation and the evolution of autosomal vs.
#  X-linked inversions
#
#  R code to run deterministic simulations of the
#  haplotype frequency recursions displayed in Figure X
#
#
#  Author: Colin Olito
#
#  NOTES:  
#          


###################################################
##  Autosomal 
#     -- no sex-specific selection
#     -- no sex-specific migration
#     -- Additive fitness effects (h = 1/2)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-DetermRecSim-Autosomal.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.01)
gen        <-  30000
h          <-  0.5
m.vals     <-  c(0.01, 0.05)
r.vals     <-  c(0, 0.1)
threshold  <-  1e-7

# Run simulation loop
test  <-  recursionFwdSimLoop(gen = gen, h = h, s.vals = s.vals,
                    		  m.vals = m.vals, r.vals = r.vals, threshold = threshold)


###################################################
##  Autosomal Sex-Specific
##    -- sf = sm = s (no sex-specific selection)
##    -- mf = mm = m (no sex-specific migration)
##    -- hf = hm = h = 1/2 (additive fitness effects)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-Autosomal-SexSpecific.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.01)
gen        <-  30000
h          <-  0.5
m.vals     <-  c(0.01, 0.05)
r.vals     <-  c(0, 0.1)
threshold  <-  1e-7

# Run simulation loop
test  <-  recursionFwdSimLoop(gen = gen, hf = h, hm = h, threshold = threshold,
							  sf.vals = s.vals, sm.vals = s.vals,
							  mf.vals = m.vals, mm.vals = m.vals,
							  r.vals  = r.vals)


###################################################
##  X-linked 
##    -- sf = sm = s (no sex-specific selection)
##    -- mf = mm = m (no sex-specific migration)
##    -- Additive fitness effects (h = 1/2)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-DetermRecSim-X-linked.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.01)
gen        <-  30000
h          <-  0.5
m.vals     <-  c(0.01, 0.05)
r.vals     <-  c(0, 0.1)
threshold  <-  1e-7

# Run simulation loop
test  <-  recursionFwdSimLoop(gen = gen, h = h, threshold = threshold,
							  sf.vals = s.vals, sm.vals = s.vals,
							  mf.vals = m.vals, mm.vals = m.vals,
							  r.vals  = r.vals)
