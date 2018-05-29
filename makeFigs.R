#  Functions to generate figures for: 
#    
#  Title: 	Local adaptation and the evolution of  
#         	X-linked vs. autosomal inversions 
#
#			A Contributed paper for Phil. Trans.
#			Roy. Soc. Theme Issue put together by 
#		   	the ESEB Special Topics Network: linking 
#			local adaptation with the evolution of 
#			sex-differences.
#
#
#  Author: Colin Olito
#
#
#  NOTES: Run this file, either from terminal using Rscript,
#		  or interactively in R. This should create all the 
#		  figures needed to correctly compile the mansucript
#		  LaTeX file.  
#          

rm(list=ls())
###############
# Dependencies
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)

source('R/functions-figures.R')


########################
# Figures for the paper
########################

# toPdf(Fig.1(), figPath(name='Fig1.pdf'), width=5, height=7.75)
# embed_fonts(figPath(name='Fig1.pdf'))

# toPdf(Fig.2(), figPath(name='Fig2.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig2.pdf'))

# toPdf(Fig3Alt(), 
#             figPath(name='Fig3Alt.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig3Alt.pdf'))

toPdf(Fig4(), 
            figPath(name='Fig4.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig4.pdf'))

toPdf(propEstSuppFigs(), 
            figPath(name='ProportionEstablishSuppFig.pdf'), width=18, height=10)
embed_fonts(figPath(name='ProportionEstablishSuppFig.pdf'))

toPdf(finalFreqSuppFig(), 
            figPath(name='FinalInvFreqSuppFig.pdf'), width=7, height=15)
embed_fonts(figPath(name='FinalInvFreqSuppFig.pdf'))

########################
# Supplementary Figures
########################






######################
# Exploratory Figures
######################

source('R/functions-figures.R')
toPdf(autoVsXEqFreqs(), figPath(name='XvAuto-EqFreqs.pdf'), width=7, height=7)
embed_fonts(figPath(name='XvAuto-EqFreqs.pdf'))

