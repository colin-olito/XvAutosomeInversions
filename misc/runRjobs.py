#!/usr/bin/env python2


# This file launch individual rscript called job1.R,job2.R ...jobN.R to the slurm cluster
# 1 job per node, adapted for Homa Papoli and Ludo Dutoit launch on Milou. See the end of this file for an example of Rscript 


###usage 
runRjobs.py numberofjobfiles 

###example usage for 24 jobs
unRjobs.py 24

import os

numberofjobs  = sys.argv[2]

for i in range(1,numberofjobs+1):
	output = open("run"+str(i)+".sh","w")
	output.write('#!/bin/sh\nR CMD BATCH job'+str(i)+'.R')
	output.close()
	os.system("mkdir -p output/data/simResults/") # just to make sure it exsits
	os.system("sbatch -A b2010010 -t 10-00:00:00 -p node -J run%i run%i.sh"  %(i,i) )



	###Example job1.R#

#rm(list=ls())
#source('R/functions-WFSims-Autosomal-SexSpec.R')#

#nReps  <-  1000000
#N      <-  30000
#m      <-  0.002
#s      <-  0.05
## Locally adaptive alleles recessive
#        makeFigReplicateAutoSexSpecInvSimsData(nReps = nReps, N = N, h = 0,
# 			m.vals = m, m.deltas = m,
#			s.vals = s, s.deltas = s,
#			r.vals = seq(from = 0, to = 0.5, by = 0.05), r.deltas = TRUE,
#			s.del.opt = "none", n = 100, u = 1e-5, h.del = 0,
#			newMutant=c("random","random"))
