# ihMTRAGE-optimize

This repository contains Matlab code for simulation of the inhomogeneous Magnetization Transfer (ihMT) experiment with the ability to include an MPRAGE like acquisition. The ihMT experiment is based on two Magnetization Transer (MT) experiments and the ihMT signal is calculated by data relating to MT at a single offset frequency from which, data relating to MT dual offset frequency, is subtracted. The Matlab code makes use of numerical integration of differential equations from the two pool or binary spin bath model for MT with and without a dipolar reservoir for the ihMT signal (see [*Varma et al.*](https://doi.org/10.1016/j.jmr.2015.08.024)).

The Matlab code in **ihMTRAGEoptimize.m** was developed to produce simulations related to optimization for ihMT within an MPRAGE acquisition, i.e. ihMTRAGE, and to produce Figure 2 in our paper describing 3D ihMTRAGE MRI ((https://doi.org/10.1002/mrm.28324)). Type **ihMTRAGEoptimize** in Matlab (with all the relevant files available) to create plots of the results from the simulation. Thus far, the code has been tested in Matlab R2017b and R2018a.

Please contact Gopal Varma ([@gvarma617](https://twitter.com/gvarma617)) to report any bugs.
