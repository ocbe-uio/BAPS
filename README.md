# BAPS: Bayesian Analysis of Population Structure

BAPS is a MATLAB package for Bayesian inference of the genetic structure in a population. BAPS treats both the allele frequencies of the molecular markers (or nucleotide frequencies for DNA sequence data) and the number of genetically diverged groups in population as random variables. However, analyses and model comparisons can also be performed using a fixed number of genetically diverged groups or pre-specified population structures.

For installation and usage instructions, see the following sections and [the PDF manual](https://github.com/ocbe-uio/BAPS/blob/develop/BAPS6manual.pdf). For methodological information, see reference below.

Corander, J. and Marttinen, P. (2006), Bayesian identification of admixture events using multilocus molecular markers. Molecular Ecology, 15: 2833-2843. [doi:10.1111/j.1365-294X.2006.02994.x](https://doi.org/10.1111/j.1365-294X.2006.02994.x)

The BAPS version available in this repository is BAPS 6.0, authored by Jukka Corander, Pekka Marttinen, Jukka Siren, Jing Tang and Lu Cheng. Copyright 2005-2012, all rights reserved.

# Installation

## Pre-requisites

To install and run BAPS, you will probably need:

1. A local installation of MATLAB (i.e., not the web version). You can download the latest version [here](https://se.mathworks.com/downloads/) (assuming you have a license)
2. MATLAB Compiler and/or MATLAB Compiler SDK (may be an available option as you install MATLAB)
3. A Bash shell

## Shortcut with Make

If you have [Make](https://www.gnu.org/software/make/) installed, you can simply type `make` from the root folder of this repository to install BAPS. Otherwise, please follow the instructions below.

## Running MATLAB from the command line

Open a terminal. From the root folder of this repository, open MATLAB. For example, if MATLAB is installed on `/usr/local/MATLAB/R2023b`, this can be done by typing:

```bash
/usr/local/MATLAB/R2023b/bin/matlab -nodesktop -nosplash
```

If you don't know where MATLAB is installed, you can find out by opening the MATLAB GUI from your OS's application's menu and typing `matlabroot` in the MATLAB command window.

## Compiling BAPS

Once MATLAB is running, compile BAPS by typing:

```MATLAB
run add_BAPS_to_path.m
run compileBaps6.m
```

This should create a `BAPS_package` folder in the root directory of this repository.

# Running BAPS

After installation, run BAPS by typing the following from a terminal, replacing `/usr/local/MATLAB/R2023b` with the path to MATLAB on your system:

```bash
bash BAPS_package/run_baps6.sh /usr/local/MATLAB/R2023b/
```

If you have [Make](https://www.gnu.org/software/make/) installed, you can simply type `make run` from the root folder of this repository to install BAPS.

Eventually, a graphical interface like the one below should appear:

![baps home screen](/aux/home_screen.png)

For more information about usage, please read [the BAPS6 manual](https://github.com/ocbe-uio/BAPS/blob/develop/BAPS6manual.pdf).

# Additional (legacy) documentation

For preservation, a [Wiki page](https://github.com/ocbe-uio/BAPS/wiki) has been setup containing a slightly-modified version of the [original webpage for BAPS](http://www.helsinki.fi/bsg/software/BAPS/), which is no longer available. Please note that the content available in those pages may not be applicable to the current version of the software available on this repository, but might still provide some useful information to its users.
