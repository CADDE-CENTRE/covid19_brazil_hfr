# covid19_brazil_hfr
code to quantify age specific HFR of SARS-CoV-2 P1.

## Installation on HPC

The code is structured as an R package. A ```yml``` file is provided and can be used to build a conda virtual environment containing most dependencies. Create the environment using:
```bash
$ cd covid19_brazil_hfr
$ conda env create -f covid19_brazil_hfr.yml
$ source activate covid19_brazil_hfr
$ export TBB_CXX_TYPE=gcc
$ export CXXFLAGS+=-fPIE
```
or alternatively
```bash
$ cd covid19_brazil_hfr
$ conda create --name covid19_brazil_hfr
$ conda install -c conda-forge r-base r-rcpp r-ggplot2 r-data.table r-viridis r-gtools r-bayesplot r-mglm r-fitdistrplus r-actuar r-abind r-knitr r-rmarkdown r-yaml r-stringi r-codetools r-ggsci  r-bh r-matrix r-inline r-gridextra r-rcppparallel r-loo r-pkgbuild r-withr r-v8
$ source activate covid19_brazil_hfr
$ export TBB_CXX_TYPE=gcc
$ export CXXFLAGS+=-fPIE
```

Then fire up R and install cmdstanr [as described here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html#installing-cmdstan-1) 
```R
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
check_cmdstan_toolchain()
install_cmdstan(cores = 4)
```

Finally on the command line in the parent directory of ```covid19_brazil_hfr```
```bash
$ R CMD build covid19_brazil_hfr
$ R CMD INSTALL covid19_brazil_hfr_1.0-0.tar.gz
```
## Usage on HPC

To submit jobs to Imperial's HPC, edit the script ```covid19_brazil_hfr/inst/scripts/startme.HPC.R``` as needed, then do
```bash
$ cd covid19_brazil_hfr/inst/scripts
$ Rscript startme.HPC.R
```

Just a few things to note:

- the jobs are generally run using 4 processors, and current runtimes are about 1 hour
- to write to the  default ```BrazilP1``` folder, you need access to the RDS project allocation ```ratmann_covid19```. Email olli with your Imperial user ID, or alternatively change ```args$out_dir``` to write to your individual allocation
- you can use the private queue ```pqcovid19c``` if you have obtained access to it, or otherwise set ```hpc.q=NA```

