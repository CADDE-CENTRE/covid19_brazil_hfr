cat(" \n -------------------------------- \n \n HFR.fit.cmdstan.R \n \n -------------------------------- \n")

library(data.table)
library(ggplot2)
library(ggsci)
library(viridis)
library(gtools)
library(MGLM)
library(fitdistrplus)
library(actuar)
library(cmdstanr)
library(posterior)
library(bayesplot)

data.table::setDTthreads(4)

if(0) # Andrea s 
{
	pkg.dir <- '~/git/covid19_brazil_hfr/'	 # Not really needed
	out.dir2 <- '~/Documents/P1Brazil/hfr-fit-210719h2-resources8-ex0-trvaluecs-hd1/'
	file.prefix2 <- basename(out.dir2)	
	loc_label <- 'sao paulo'
	fr.exclude.P1 <- 0
	use.only.hosp.deaths <- 1
}


# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-outdir')	
	stopifnot(args_line[[3]]=='-loc_label')
	stopifnot(args_line[[5]]=='-file_prefix')
	stopifnot(args_line[[7]]=='-pkgdir')
	stopifnot(args_line[[9]]=='-fr.exclude.P1')
	stopifnot(args_line[[11]]=='-use.only.hosp.deaths')	
	out.dir2 <- args_line[[2]]
	loc_label <- gsub('-|_',' ',args_line[[4]])
	file.prefix2 <- args_line[[6]]
	pkg.dir <- args_line[[8]]
	fr.exclude.P1 <- as.integer(args_line[[10]])
	use.only.hosp.deaths <- as.integer(args_line[[12]])
} 

# stop if output exists already
tmp <- paste0(file.path(out.dir2,file.prefix2),'_',gsub(' ','',loc_label),'_fit.rds')
cat('\nSearching for file ', tmp)
if( file.exists( tmp ) )
{
	cat('\noutput file exists')
	quit(save='no')
}

# load workspace
tmp <- paste0( file.path(out.dir2,file.prefix2),'_allloc_workspace.rda')
cat('\nAttempting to load ', tmp)
load( tmp )
plot.dir <- out.dir <- file.prefix <- NULL
gc()

# add to args the local arguments and update to new out.dir
args$out.dir <- out.dir2
args$out.base <- file.path(out.dir2,file.prefix2)
args$loc_label <- gsub('_', ' ', loc_label)
cat('\nSelected location ', args$loc_label, '\n')
stopifnot(args$loc_label %in% args$selected.locs)
args$fr.exclude.P1 <- fr.exclude.P1
args$use.only.hosp.deaths <- use.only.hosp.deaths

# this model does not support the previous alternative specification
stopifnot(args$use.only.hosp.deaths==1)
stopifnot(args$fr.exclude.P1==0)

####################################################################
# STAN MODEL
####################################################################

# subseting data table to the selected loc 
dd <- subset(dd, loc_label == args$loc_label)
dp <- subset(dp, loc_label == args$loc_label)
dsd <- subset(dsd, loc_label == args$loc_label)
drpop <- subset(drpop, loc_label == args$loc_label)
dh <- subset(dh, loc_label == args$loc_label)
dtd <- subset(dtd, loc_label == args$loc_label)
dhc <- subset(dhc, loc_label == args$loc_label)
dhf <- subset(dhf, loc_label == args$loc_label)
dha <- subset(dha, loc_label == args$loc_label)
dht <- subset(dht, loc_label == args$loc_label)
dicu <- subset(dicu, loc_label == args$loc_label)
dpop <- subset(dpop, loc_label == args$loc_label)
dr <- subset(dr, loc_label == args$loc_label)
dt0s <- subset(dt0s, loc_label == args$loc_label)

stan_data <- list()

#	add hospital admissions by age
dh[, week.idx := week-min(week)+1L]
dh[, age.idx := as.integer(age.label)]
stan_data$Th <- dh[, max(week.idx)]
stan_data$A <- dh[, max(age.idx)]
tmp <- dcast.data.table(dh, week.idx~age.idx, value.var='hosps')	
stan_data$age_hosp_admissions <- round(unname(as.matrix(tmp))[,-1])

# proportion of p1
stan_data$N <- nrow(dp)
stan_data$seq_timepoints <- dp[, date.idx]
stan_data$seq_sampled <- dp$n_sampled
stan_data$seq_P1 <- dp$n_positive
tmp <- subset(dp, n_sampled>=5)
tmp <- tmp[which.min(abs(tmp[n_sampled>=5,prop_P1]-0.5)),]
stan_data$propp1_logistic_midpoint_priormeansd <- c(tmp[, date.idx], 20)

# last date when proportion of p1 zero
tmp2 <- ddates[loc_label==args$loc_label, p1em.median]
tmp <- dw[Date == tmp2, week]
stan_data$seq_max_noP1_weekidx <- unique(dh[loc_label == args$loc_label & week == tmp, week.idx])

#	add size of susceptible population 3 weeks later at time of hosp admission
drpop[, week.idx := week-min(week)+1L]
drpop[, age.idx := as.integer(age.label)]
tmp <- dcast.data.table(drpop, week.idx~age.idx, value.var='rpop.alive.notvac.push3')
stan_data$age_pop_mat_underlying_hosps <- unname(as.matrix(tmp))[,-1]

#	add hospital collapse covariates
stan_data$P <- length(args$fr.predictors)
stopifnot( all(args$fr.predictors %in% colnames(dhc)) )
tmp <- subset(dhc, select=args$fr.predictors)
stopifnot( nrow(tmp)==stan_data$Th )
stan_data$fr_predictors <- lapply(1:stan_data$A, function(a) unname(as.matrix(tmp)) )

#	add hospital collapse data on fatalities corresponding to hospital adm in a week, during wildtype period
tmp <- dcast.data.table(dh, week.idx+week~age.idx, value.var='hosps')
tmp <- merge(subset(dhc, select=c(week,week.cat)), tmp, by='week')
tmp2 <- NA_integer_
if(nrow(tmp[week.cat=='initial']) > 0)
{
	tmp2 <- tmp[week.cat=='initial', max(week.idx)]
}
tmp <- unname(as.matrix(tmp[,-c(1,2,3)]))
if(!is.na(tmp2))
{
	tmp[1:tmp2,] <- 0
}
tmp <- as.numeric(t(tmp))
stan_data$Th_wildtype_fr_obs <- length(which(tmp>0))
stan_data$fr_wildtype_obsidx_row_major_order <- which(tmp>0)
stan_data$fr_wildtype_hosps_row_major_order <- tmp[tmp>0]
tmp <- dcast.data.table(dh, week.idx+week~age.idx, value.var='deaths_of_hospsi')
tmp <- as.numeric(t(unname(as.matrix(tmp))[,-c(1,2)]))
stan_data$fr_wildtype_deaths_row_major_order <- tmp[stan_data$fr_wildtype_obsidx_row_major_order]
stan_data$fr_shrinkage_df <- 1

#	set inits
stan_init <- list()
stan_init$log_age_hosp_rate_wildtype <- rep(0, stan_data$A)
stan_init$log_age_hosp_rate_P1 <- rep(0, stan_data$A)
stan_init$age_hosps_v_inflation <- 0.1
stan_init$logit_hosp_fatality_ratio_wildtype <- dhf$logit_hfr
stan_init$logit_hosp_fatality_ratio_P1_rnde <- rep(0, stan_data$A)
stan_init$logit_hosp_fatality_ratio_P1_rnde_sd <- 1
stan_init$propp1_logistic_growthrate <- .25
stan_init$propp1_logistic_midpoint <- stan_data$propp1_logistic_midpoint_priormeansd[1]
stan_init$propp1_overdisp_inv <- 0.05
stan_init$fr_regcoeff_overall <- rep(0.04, stan_data$P)
stan_init$fr_shrinkage <- rep(1, stan_data$P)
stan_init$fr_scale <- 1
stan_init$fr_overdisp_inv <- 0.01

# save workspace
save.image(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_workspace.rda'))

# copy Stan model file to working dir and compile model
file.copy(args$file.stan.model, file.path(args$out.dir, basename(args$file.stan.model)))
m <- cmdstanr::cmdstan_model(stan_file=file.path(args$out.dir, basename(args$file.stan.model)), force_recompile=TRUE)

# run Stan
if(0)
{
	m_fit <- m$sample( 
			data=stan_data, seed=42,
			refresh=10, iter_warmup=50, iter_sampling=20, chains=1,
			parallel_chains=1, threads_per_chain = 1, save_warmup=TRUE,
			init= list(stan_init)
		)		
}else{

m_fit <- m$sample( 
		data=stan_data, seed=42,
		refresh=1e2, iter_warmup=5e2, iter_sampling=2e3, chains=4,
		parallel_chains=4, threads_per_chain = 1, save_warmup=TRUE,
		init= list(stan_init,stan_init,stan_init,stan_init)
	)	
}
# save to be on safe side 
m_fit$save_object(file = paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_fit.rds'))
