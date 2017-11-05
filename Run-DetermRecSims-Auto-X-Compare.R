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


###############################
## Low recombination (r = 0.01)
###############################

###################################################
##  Autosomal 
#     -- no sex-specific selection
#     -- no sex-specific migration
#     -- low recombination (r = 0.01)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-DetermRecSim-Autosomal.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.005)
gen        <-  35000
h.vals     <-  c(0, 0.5, 1)
m.vals     <-  c(0.01, 0.05)
r          <-  0.01
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, r = r, s.vals = s.vals,
					m.vals = m.vals, h.vals = h.vals, threshold = threshold)


###################################################
##  X-linked 
##    -- sf = sm = s (no sex-specific selection)
##    -- mf = mm = m (no sex-specific migration)
#     -- low recombination (r = 0.01)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-DetermRecSim-X-linked.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.005)
gen        <-  35000
h.vals     <-  c(0, 0.5, 1)
m.vals     <-  c(0.01, 0.05)
r          <-  0.01
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, r = r, threshold = threshold,
					sf.vals = s.vals, sm.vals = s.vals,
					mf.vals = m.vals, mm.vals = m.vals,
					h.vals  = h.vals)



###################################################
##  Autosomal Sex-Specific
##    -- sf = sm = s (no sex-specific selection)
##    -- mf = mm = m (no sex-specific migration)
#     -- low recombination (r = 0.01)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-Autosomal-SexSpecific.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.01)
gen        <-  35000
h.vals     <-  c(0, 0.5, 1)
m.vals     <-  c(0.01, 0.05)
r          <-  0.01
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, hf = h, hm = h, threshold = threshold,
				    sf.vals = s.vals, sm.vals = s.vals,
					mf.vals = m.vals, mm.vals = m.vals,
					h.vals  = h.vals)



################################
## High recombination (r = 0.1)
################################

###################################################
##  Autosomal 
#     -- no sex-specific selection
#     -- no sex-specific migration
#     -- ## High recombination (r = 0.1)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-DetermRecSim-Autosomal.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.005, 0.2, by=0.005)
gen        <-  50000
h.vals     <-  c(0, 0.5, 1)
m.vals     <-  c(0.01, 0.05)
r          <-  0.1
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, r = r, s.vals = s.vals,
					m.vals = m.vals, h.vals = h.vals, threshold = threshold)


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
s.vals     <-  seq(0.005, 0.2, by=0.005)
gen        <-  50000
h.vals     <-  c(0, 0.5, 1)
m.vals     <-  c(0.01, 0.05)
r          <-  0.1
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, r = r, threshold = threshold,
					sf.vals = s.vals, sm.vals = s.vals,
					mf.vals = m.vals, mm.vals = m.vals,
					h.vals  = h.vals)



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
gen        <-  50000
h.vals     <-  c(0, 0.5, 1)
m.vals     <-  c(0.01, 0.05)
r          <-  0.1
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, hf = h, hm = h, threshold = threshold,
				    sf.vals = s.vals, sm.vals = s.vals,
					mf.vals = m.vals, mm.vals = m.vals,
					h.vals  = h.vals)



########################################################
## Locally adaptive alleles completely dominant (h = 1)
########################################################

###################################################
##  Autosomal 
#     -- no sex-specific selection
#     -- no sex-specific migration
#     -- Locally adaptive alleles recessive (h = 1)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-DetermRecSim-Autosomal.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.005)
gen        <-  35000
h          <-  1
m.vals     <-  c(0.01, 0.05)
r          <-  
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, r = r, s.vals = s.vals,
					m.vals = m.vals, r.vals = r.vals, threshold = threshold)


###################################################
##  X-linked 
##    -- sf = sm = s (no sex-specific selection)
##    -- mf = mm = m (no sex-specific migration)
#     -- Locally adaptive alleles recessive (h = 1)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-DetermRecSim-X-linked.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.005)
gen        <-  35000
h          <-  1
m.vals     <-  c(0.01, 0.05)
r          <-  
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, r = r, threshold = threshold,
					sf.vals = s.vals, sm.vals = s.vals,
					mf.vals = m.vals, mm.vals = m.vals,
					h.vals  = h.vals)



###################################################
##  Autosomal Sex-Specific
##    -- sf = sm = s (no sex-specific selection)
##    -- mf = mm = m (no sex-specific migration)
##    -- Locally adaptive alleles recessive (hf = hm = h = 1)

# Clear existing objects/functions
rm(list=ls())

# Dependencies
source('R/functions-figures.R')
source('R/functions-Autosomal-SexSpecific.R')

# Set s.vals, other parameters for consistent plotting 
# Use consistent s.vals for pretty plots
s.vals     <-  seq(0.01, 0.2, by=0.01)
gen        <-  35000
h          <-  1
m.vals     <-  c(0.01, 0.05)
r          <-  
threshold  <-  1e-7

# Run simulation loop
recursionFwdSimLoop(gen = gen, hf = h, hm = h, threshold = threshold,
				    sf.vals = s.vals, sm.vals = s.vals,
					mf.vals = m.vals, mm.vals = m.vals,
					h.vals  = h.vals)

