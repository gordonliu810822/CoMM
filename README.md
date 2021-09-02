CoMM
===
Collaborative mixed model to dissecting genetic contributions to complex traits by leveraging regulatory information.

Installation 
===========

To install the development version of CoMM, it's easiest to use the 'devtools' package. Note that REMI depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("Shufeyangyi2015310117/CoMM")
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

Reproducing Results of Yang et al. (2019)
===========
The simulation results of CoMM-S<sup>2</sup> can be reproduced by the code of file in the [simulation](https://github.com/gordonliu810822/CoMM/tree/master/simulation). The introduction of usage of CoMM-S<sup>2</sup> is in the ['CoMM' vignette](https://github.com/gordonliu810822/CoMM/blob/master/vignettes/CoMM.pdf), please refer to the code file of CoMM-S<sup>2</sup> [CoMM_ S2_ simulation_power.R](https://github.com/gordonliu810822/CoMM/blob/master/simulation/CoMM_S2_simulation_power.R) and [CoMM_ S2_ simulation_t1e.R](https://github.com/gordonliu810822/CoMM/blob/master/simulation/CoMM_S2_simulation_t1e.R) in the vignette to reproduce the simulation results. The excel tables for the results of analysis of 14 GWAS summary statistics can be found at  [gene_significance4_centering_excel.zip](https://github.com/gordonliu810822/CoMM/blob/master/Results/gene_significance4_centering_excel.zip)

Reproducing Results of Yang et al. (2021)
===========
The simulation results of CoMM-S<sup>4</sup> can be reproduced by the code of file in the [simulation](https://github.com/Shufeyangyi2015310117/CoMM/tree/master/simulation). The introduction of usage of CoMM-S<sup>4</sup> is in the ['CoMM' vignette](https://github.com/Shufeyangyi2015310117/CoMM/blob/master/vignettes/CoMM.pdf), please refer to the code file of CoMM-S<sup>4</sup> [CoMM_ S4_ simulation_power.R](https://github.com/Shufeyangyi2015310117/CoMM/blob/master/simulation/CoMM_S4_simulation_power.R) and [CoMM_ S4_ simulation_t1e.R](https://github.com/Shufeyangyi2015310117/CoMM/blob/master/simulation/CoMM_S4_simulation_t1e.R) in the vignette to reproduce the simulation results. 


References
==========
1. Can Yang<sup>c</sup>,
Xiang Wan<sup>c</sup>, 
Xinyi Lin, Mengjie Chen, Xiang Zhou, 
Jin Liu<sup>+</sup>. (2018) [CoMM: a collaborative mixed model to dissecting genetic contributions to complex traits by leveraging regulatory information](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty865/5123355), Bioinformatics, 35(10), 1644–1652.
2. Yi Yang, Xingjie Shi, Yuling Jiao, Jian Huang, Min Chen, Xiang Zhou, Lei Sun, Xinyi Lin, Can Yang, Jin Liu<sup>+</sup>. (2019) [CoMM-S<sup>2</sup>: a collaborative mixed model using summary statistics in transcriptome-wide association studies](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz880/5637753), Bioinformatics, btz880.
3. Yi Yang †, Kar-Fu Yeung † and Jin Liu<sup>*</sup>. (2021) CoMM-S<sup>4</sup>: a collaborative mixed model using summary-level eQTL and GWAS datasets in transcriptome-wide association studies.

Development
===========

This package is developed and maintained by Yi Yang (gmsyany@nus.edu.sg) and Jin Liu (jin.liu@duke-nus.edu.sg).

# Demonstration

For an example of typical CoMM usage, please see our [Package vignette](https://shufeyangyi2015310117.github.io/CoMM/) for a demonstration and overview of the functions included in CoMM.


