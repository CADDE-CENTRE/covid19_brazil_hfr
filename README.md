# covid19_brazil_hfr 

## Overview

This repository contains the code and the data necessary to run the model and associated Bayesian analysis via [Stan](https://mc-stan.org/users/documentation/)  in:

-  A Brizzi, C Whittaker,  LMS Servo et al. "Factors driving extensive spatial and temporal fluctuations in COVID-19 fatality rates in Brazilian hospitals". Imperial College London (06-10-2021), doi  https://doi.org/10.25561/91875.

## Data

The directory  ```inst/data``` contains:

- The hospital admission data in ``` SIVEP_hospital_20-09-2021-all.RDS ``` obtained from the [SIVEP-Gripe Platform](https://opendatasus.saude.gov.br/dataset/bd-srag-2021) is processed and redistributed under [CC BY-4.0](https://creativecommons.org/licenses/by/4.0/).
- Brazilian registry death data in  ```registry_covid_detailed_02-08-2021.csv```   available at https://transparencia.registrocivil.org.br/  were downloaded from [Marcelo Oliveira's github](https://github.com/capyvara/brazil-civil-registry-data). These are preprocessed and redistributed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). processed deaths data from the [Brazilian Civil registry](https://transparencia.registrocivil.org.br/).
- Hospital resource data in ```IPEA_ICUbeds_physicians_210928.csv```  obtained from the [National Register of Health Facilities](cnes.datasus.gov.br) is processed and redistributed under [CC BY-4.0](https://creativecommons.org/licenses/by/4.0/). 
- Vaccination data in ```aggregated_vaccinations_210805.rds ``` obtained from the [Brazilian Ministry of Health](https://opendatasus.saude.gov.br/dataset/covid-19-vacinacao/resource/ef3bd0b8-b605-474b-9ae5-c97390c197a8) is processed and redistributed under [CC BY-4.0](https://creativecommons.org/licenses/by/4.0/).
- Genomic data in``` genomic_data_210702.csv```   from [GISAID](https://www.gisaid.org/) and [Fiocruz](http://www.genomahcov.ﬁocruz.br). Gisaid data were obtained on June 28, 2021 (See ``` acknowledgments_GISAID_Tables```). Fiocruz data were manually scraped from the website on June 15, 2021.
- ```PNADc_populationpyramids_210617.csv ```: Population projections from the [Continuous National
  Household Sample Survey](https://www.ibge.gov.br/estatisticas/sociais/populacao/9171-pesquisa-nacional-por-amostra-de-domicilios-continua-mensal.html?=&t=o-que-e)
- `private_public_burials_manaus_210706.csv` : data on daily burials from [Manaus Mayor's office](https://www.manaus.am.gov.br/noticia/) obtained from [here](https://drive.google.com/drive/folders/1DYmuzzOJwHrB3LtXMrmguc1mIauJH4V_).
- Other helper datasets.

## License

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png



The code in this repository is licensed under [CC BY-4.0](https://creativecommons.org/licenses/by/4.0/).

[![CC BY 4.0][cc-by-image]][cc-by]

## Warranty

Imperial makes no representation or warranty about the accuracy or  completeness of the data nor that the results will not constitute in  infringement of third-party rights. Imperial accepts no liability or  responsibility for any use which may be made of any results, for the  results, nor for any reliance which may be placed on any such work or  results.



## System Requirements

- macOS or UNIX, the code was developed on macOS Mojave 10.14.
- R version  >= 3.6.1 (that s a wild guess based on Nature, mine is 4.0.5)

Package requirements are reported in the ``` covid19_brazil_hfr.yml``` file 

## Usage

To reproduce the main analyses from an UNIX environment, it is only necessary to move to the ``` inst``` directory and run:

```bash
$ Rscript start_local.R -outdir foo
```

specifying an output directory which will contain all the plots and Hamiltonian Monte Carlo posterior draws. 
The R script produces a bash script containing instructions to preprocess the data, fit the model and then analyse the results. We typically run 4 HMC chains in parallel, though the number of chains and the number of iterations can be modified in `inst/script/HFR.fit.cmdstan.R`  by modifying the arguments in:

```R
m_fit <- m$sample( 
                data=stan_data, seed=42,
                refresh=1e2, iter_warmup=5e2, iter_sampling=2e3, chains=4,
                parallel_chains=4, threads_per_chain = 1, save_warmup=TRUE,
                init= list(stan_init,stan_init,stan_init,stan_init)
        ) 
```



To run sensitivity analyses with other genomic data sources, it is first necessary to run the chains with the specified data by adding one of the flags ``` -random``` or ``` -fiocruz``` to the above command, or including the flag in the directory name.

```bash
$ # The following are equivalent:
$ Rscript start_local.R -outdir foo -random
$ Rscript start_local.R -outdir foo-random
```

In the first case, the `'-random'` flag will be appended to the output directory name, so that in both cases, the results will be stored in `foo-random` .

After the above analyses are completed, it is possible to reproduce the plots in the supplementary materials by moving to the ``` inst/scripts``` directory and running:

```bash
$ # if -random sensitivity:
$ Rscript sensitivity_analysis_ControlledSequenceDataCollection.R -dirs foo bar

$ # if -fiocruz sensitivity:
$ Rscript sensitivity_analysis_Fiocruz.R -dirs foo bar
```

where ```foo``` and ```bar``` are the directories containing the results from the standard and from the sensitivity analysis respectively.

We also include the `start.HPC.R` scripts that allows to run the models for the different cities in parallel.
This was designed to run on the Imperial College London Research Computing Service, but the script can be adapted for different set ups.
The script writes and queues bash scripts to run the preprocessing analyses, and then run the state capital models in parallel. 

## Acknowledgements

We thank all contributors to GISAID for making SARS‐CoV‐2 sequence data information publicly available as listed in`acknowledgments_GISAID`; all contributors to Rede Genomica Fiocruz for making SARS‐CoV‐2 variant frequency data publicly available; all members of the CADDE network for their comments throughout the project and earlier versions of the manuscript; Oliver G. Pybus, Andrew Rambaut and JT. McCrone for their
insightful comments on SARS‐CoV‐2 phylogenetic analyses; and the Imperial College Research Computing Service, DOI: 10.14469/hpc/2232, for providing the computational resources to perform this study

## Funding

This study was supported by the Medical Research Council-São Paulo Research Foundation (FAPESP) CADDE partnership award (MR/S0195/1 and FAPESP18/14389-0) (https://caddecentre.org), by the EPSRC through the EPSRC Centre for Doctoral Training in Modern Statistics and Statistical Machine Learning at Imperial and Oxford. RSA from the Rede Coronaômica BR MCTI/FINEP afﬁliated to RedeVírus/MCTI (FINEP 01.20.0029.000462/20, CNPq 404096/2020-4), from CNPq (312688/2017-2 and 439119/2018-9), MEC/CAPES (14/2020 - 23072.211119/2020-10), FINEP (0494/2001.20.0026.00); LSB acknowledges support from Inova Fiocruz (48401485034116); DSC is funded by the Clarendon Fund, University of Oxford Department of Zoology and Merton College; NF from Wellcome Trust and Royal Society (N.R.F.: Sir Henry Dale Fellowship. 204311/Z/16/Z) and from the Bill & Melinda Gates Foundation (INV-034540) ; CAP was supported by FAPESP (2019/21858-0), Fundação Faculdade de Medicina and Coordenação de Aperfeiçoamento de Pessoal de Nível Superior Brasil (CAPES) Finance Code 001 OTR from the Instituto de Salud Carlos III (Sara Borrell fellowship, CD19/00110); OR from the Bill & Melinda Gates Foundation (OPP1175094); ES from Bill & Melinda Gates Foundation (INV-034652); RPS from the Rede Coronaômica BR MCTI/FINEP afﬁliated to RedeVírus/MCTI (FINEP 01.20.0029.000462/20, CNPq 404096/2020-4), from CNPq (310627/2018-4), MEC/CAPES (14/2020 - 23072.211119/2020-10), FINEP (0494/20 01.20.0026.00), FAPEMIG (APQ-00475-20); WMS from FAPESP (2017/13981-0, 2019/24251-9) and the NIH (AI12094); CW acknowledges an MRC Doctoral Training partnership studentship
