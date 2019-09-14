# IRGA

Repository with the code used for the manuscript in preparation titled "Approximating high-dimensional posteriors with nuisance parameters via integrated rotated Gaussian approximation".


## Description of files

* [`BVS_IRGA.R`](BVS_IRGA.R) contains functions for Bayesian variable selection (BVS) using integrated rotated Gaussian approximation (IRGA). It loads [`VAMP.R`](VAMP.R) which should thus be present in the working directory. It also uses parallel computing which I have only tested on MacOS and Linux and might therefore not work on Windows operating systems.

* [`BVS_other.R`](BVS_other.R) contains functions for Bayesian variable selection using variational Bayes, expectation propagation, and Gibbs sampling. It loads [`epBVS.R`](epBVS.R) which should thus be present in the working directory.

* [`epBVS.R`](epBVS.R) contains code for Bayesian variable selection using expectation propagation by José Miguel Hernández-Lobato downloaded from [Bitbucket](https://bitbucket.org/jmh233/ep-lrmssp/src/e85e3170757ec1a4c6032c0aeca22b0ee57a6c67/methods/ep/epBVS.R?at=default&fileviewer=file-view-default).

* [`setup.R`](setup.R) sets up the R environment for the [Applications](Applications) and [Simulations](Simulations). It loads [`BVS_IRGA.R`](BVS_IRGA.R) and [`BVS_other.R`](BVS_other.R) which should thus be present in the working directory.

* [`VAMP.R`](VAMP.R) implements vector approximate message passing (VAMP).


## Folders

[Applications](Applications) and [Simulations](Simulations) contain the scripts that produce the applications and simulations of the manuscript.
