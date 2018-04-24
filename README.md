CoMM
===
Collaborative mixed model to dissecting genetic contributions to complex traits by leveraging regulatory information.

Installation 
===========

To install the development version of CoMM, it's easiest to use the 'devtools' package. Note that REMI depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("gordonliu810822/CoMM")
```

Usage
===========
[The 'CoMM' vignette](https://github.com/gordonliu810822/CoMM/blob/master/vignettes/CoMM.pdf) will provide a good start point for the genetic analysis using CoMM package. The following help page will also provide quick references for CoMM package and the example command lines:

```
library(CoMM)
package?CoMM
```

Reproducing Results of Yang et al. (2018)
===========
All the simulation results can be reproduced by using the code at [simulation](https://github.com/gordonliu810822/CoMM/tree/master/simulation). Before running simulation to reproduce the results, please familarize yourself with CoMM using ['CoMM' vignette](https://github.com/gordonliu810822/CoMM/blob/master/vignettes/CoMM.pdf). Simulation results can be reproduced using [simulation.R](https://github.com/gordonliu810822/CoMM/blob/master/simulation/simulation.R) with a batch script [nscc_sim.txt](https://github.com/gordonliu810822/CoMM/blob/master/simulation/nscc_sim.txt). 


References
==========
Can Yang$$^*$$, Xiang Wan$$^*$$, Xinyi Lin, Mengjie Chen, Xiang Zhou, Jin Liu$$^+$$. (2018) CoMM: a collaborative mixed model to dissecting genetic contributions to complex traits by leveraging regulatory information.


Development
===========

This package is developed and maintained by Jin Liu (jin.liu@duke-nus.edu.sg).
