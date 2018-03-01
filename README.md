# The evolution of X-linked vs. autosomal inversions during local adaptation

## Overview

This is a GitHub repository for the development of a collaborative research project begun at the 2017 ESEB Special Topics Network: *linking local adaptation with the evolution of sex-differences*. A major goal of the workshop was to put together a proposal for a Theme Issue of Phil Trans. Roy. Soc. B on sex-specific local adaptation. This project is intended to result in a contributed paper to this project. Here you can find all of the necessary code to reproduce the analyses and the pre-publication manuscript. The *Mathematica* `.nb` files where the models are developed are not tracked by git, but can be downloaded as an online appendix from the publisher website above (or contact me and I will be happy to provide if you don't have access to the appendices). Aside from the `.nb` files, all necessary code for creating the figures can be found in the `./R/functions-*.R` files. 


## Abstract for the proposed Phil. Trans. Theme Issue

Local selection can favour the evolution of tightly linked clusters of locally adapted alleles – so-called genomic islands of differentiation. While there has been much theoretical emphasis linking the presence of inversions with the evolution of clusters of locally adapted alleles, there has been little attention devoted to the potential effects of sex-linkage on the accumulation of locally co-adapted gene combinations. In principle, it is possible that locally adapted gene clusters may accumulate at different rates on autosomes and sex chromosomes as a result of their unique patterns of recombination, differential susceptibility to invasion of inversion polymorphisms, and distinct responses to natural selection. Here, we model the evolution of inversions that capture locally adapted gene combinations and consider the opportunity for inversion establishment and polymorphism on the X chromosome relative to the autosomes. We consider both the initial invasion of inversions, as well as their fixation or persistence as inversion polymorphisms. Our results show that inversions should differentially accumulate and segregate on sex chromosomes and autosomes as a consequence of selection for local adaptation. We review the empirical data on fixed and polymorphic inversions in animal genomes with an emphasis on autosome/X differences. We discuss the new theoretical results in light of this data.


## Citing information

Citing information for the resulting paper will be provided when it is made [available through the publisher](http://XXXXX). You can also contact me directly if you would like a reprint. 


## Reproducing the manuscript

The easiest way to reproduce the manuscript is to simply clone the repo, run `createFigs.R`, and then compile the manuscript file `doc/XvAuto.tex` using whatever default LaTeX editor/engine you have. 


## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the XvAutosomeInversions github [issues page](https://github.com/colin-olito/XvAutosomeInversions/issues). If you would like to report a bug/issue, and do not have a github account (and don't want to get one), please send a brief email detailing the problem you encountered to colin.olito at monash dot edu.


## Ludo branch notes

I added a file that track the progress and the current issues of understanding