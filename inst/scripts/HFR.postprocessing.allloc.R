cat(" \n -------------------------------- \n \n HFR.postprocessing.allloc.R \n \n -------------------------------- \n")

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

if(1) # Andrea s 
{
  out.dir2 <- '~/Documents/P1Brazil/hfr-fit-210719h2-resources8-ex0-trvaluecs-hd1/'
  file.prefix2 <- basename(out.dir2)
  threshold_check <-  0
  use.only.hosp.deaths <- 0
}

# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-outdir')	
  stopifnot(args_line[[3]]=='-file_prefix')
  stopifnot(args_line[[5]]=='-threshold_check')
  stopifnot(args_line[[7]]=='-use.only.hosp.deaths')
  out.dir2 <- args_line[[2]]
  file.prefix2 <- args_line[[4]]
  threshold_check <- as.integer(args_line[[6]])
  use.only.hosp.deaths <- as.integer(args_line[[8]])
} 

####################################################################
# Check that at least threshold number of chains are done running
####################################################################

tmp <- list.files(out.dir2, pattern='_fit.rds')
cat('\ntry to bugfix: start\n')
print(out.dir2)
print(tmp)
cat('\ntry to bugfix: end\n')
if(length(tmp) < threshold_check)
{
  stop("Not enough runs have completed yet!")
}

####################################################################
# search for available outputs
#################################################################### 

do <- data.table(FW=list.files(out.dir2, pattern='_workspace.rda'))
do[, PREFIX := do[, sapply(strsplit(FW, split='_'),'[[',1)]]
do[, loc_label2 := do[, sapply(strsplit(FW, split='_'),'[[',2)]]
tmp <- data.table(FF=list.files(out.dir2, pattern='_fit.rds'))
tmp[, PREFIX := tmp[, sapply(strsplit(FF, split='_'),'[[',1)]]
tmp[, loc_label2 := tmp[, sapply(strsplit(FF, split='_'),'[[',2)]]
do <- merge(do,tmp,by=c('PREFIX','loc_label2'))
tmp <- data.table(loc_label=c('belo horizonte','curitiba', "florianopolis", "goiania", "joao pessoa", "macapa", 'manaus', "natal", 'porto alegre', "porto velho", 'rio de janeiro', 'salvador', 'sao paulo', "sao luis"))
tmp[, loc_label2 := gsub(' ','',loc_label)]
do <- merge(do, tmp, by='loc_label2')

####################################################################
# load stan_data
####################################################################
stan_data <- list()
for(i in 1:nrow(do))
{	
	tmp <- file.path(out.dir2,do$FW[i])
	cat('\nreading ',tmp)
	load(tmp,  tmp_env <- new.env())
	tmp_list <- as.list.environment(tmp_env)	
	stan_data[[i]] <- tmp_list$stan_data
	rm(tmp_list)
}
names(stan_data) <- do[, loc_label]

####################################################################
# load workspace
####################################################################

tmp <- paste0( file.path(out.dir2,file.prefix2),'_allloc_workspace.rda')
cat('\nAttempting to load ', tmp)
load( tmp )
plot.dir <- out.dir <- file.prefix <- NULL
gc()

# add to args the local arguments and update to new out.dir
args$out.dir <- out.dir2
args$out.base <- file.path(out.dir2,file.prefix2)
args$city.palette <- grDevices::colorRampPalette(ggsci::pal_npg("nrc")(10))
args$use.only.hosp.deaths <- use.only.hosp.deaths

# process dh
tmp <- dh[, list(week=unique(week), week.idx= unique(week)-min(week)+1L), by='loc_label']
dh <- merge(dh, tmp, by=c('loc_label','week'))
dh[, age.idx := as.integer(age.label)]


####################################################################
# load fit
####################################################################
pos <- list()
for(i in 1:nrow(do))
{	
	tmp <- file.path(out.dir2,do$FF[i])
	cat('\nreading ',tmp)
	m_fit <- readRDS(tmp)
	
	# get HMC samples in wide array format
	pd <- m_fit$draws(inc_warmup = FALSE)	

	# bring samples into long array format
	select.chains <- seq_along(dimnames(pd)[['chain']])
	iters <- 1:(m_fit$metadata()[['iter_sampling']])
	age.idx.n <- stan_data[[i]][['A']]	
	#iters <- 1000:(m_fit$metadata()[['iter_sampling']])
	#pd <- pd[iters,,]
	
	po <- list()
	tmp <- pd[,,which(grepl('log_age_hosp_rate_wildtype',dimnames(pd)[[3]]))]
	po$log_age_hosp_rate_wildtype <- unname(apply(tmp[,select.chains,], 3, rbind))
	tmp <- pd[,,which(grepl('log_age_hosp_rate_P1',dimnames(pd)[[3]]))]
	po$log_age_hosp_rate_P1 <- unname(apply(tmp[,select.chains,], 3, rbind))
	tmp <- pd[,,which(grepl('age_prop_hosps_wildtype',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$age_prop_hosps_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('age_prop_hosps_P1',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$age_prop_hosps_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('age_prop_deaths_wildtype',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$age_prop_deaths_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('age_prop_deaths_P1',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$age_prop_deaths_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('^prop_P1',dimnames(pd)[[3]]) & !grepl('logit_prop_P1',dimnames(pd)[[3]]) & !grepl('age_prop_P1',dimnames(pd)[[3]])  & !grepl('prop_P1_fit',dimnames(pd)[[3]]))]
	po$prop_P1 <- unname(apply(tmp[,select.chains,], 3, rbind))
	tmp <- pd[,,which(grepl('logit_hosp_fatality_ratio_wildtype\\[',dimnames(pd)[[3]]))]
	po$logit_hosp_fatality_ratio_wildtype <- unname(apply(tmp[,select.chains,], 3, rbind))
	tmp <- pd[,,which(grepl('logit_hosp_fatality_ratio_P1_rnde\\[',dimnames(pd)[[3]]))]
	po$logit_hosp_fatality_ratio_P1_rnde <- unname(apply(tmp[,select.chains,], 3, rbind))
	tmp <- pd[,,which(grepl('age_hosps_v_inflation',dimnames(pd)[[3]]))]
	po$age_hosps_v_inflation <- as.vector(tmp[,select.chains,])
	tmp <- pd[,,which(grepl('exp_deaths_wildtype_in_hosp',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$exp_deaths_wildtype_in_hosp <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('exp_deaths_P1_in_hosp',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$exp_deaths_P1_in_hosp <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
 	tmp <- pd[,,which(grepl('exp_hosp_adm_wildtype',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$exp_hosp_adm_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('exp_hosp_adm_P1',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$exp_hosp_adm_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('hfr_overall',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$hfr_overall <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('hfr_wildtype',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$hfr_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('hfr_P1',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$hfr_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	tmp <- pd[,,which(grepl('fr_regcoeff_overall',dimnames(pd)[[3]]))]
	po$fr_regcoeff_overall <- unname(apply(tmp[,select.chains,], 3, rbind))
	tmp <- pd[,,which(grepl('fr_overdisp_inv',dimnames(pd)[[3]]))]
	po$fr_overdisp_inv <- unname(apply(tmp[,select.chains,], 3, rbind))
	tmp <- pd[,,which(grepl('^fr_multiplier_by_week\\[',dimnames(pd)[[3]]))]
	tmp <- apply(tmp[,select.chains,], 3, rbind)
	po$fr_multiplier_by_week <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
	
	#new to get \tau and \kappa_{l,i}
	tmp <- pd[,,which(grepl('fr_scale',dimnames(pd)[[3]]))]
	po$fr_scale <- apply(tmp[,select.chains,], 3, rbind)
	tmp <- pd[,,which(grepl('fr_shrinkage',dimnames(pd)[[3]]))]
	po$fr_shrinkage <- apply(tmp[,select.chains,], 3, rbind)
	
	rm(pd, m_fit)
	pos[[do$loc_label[i]]] <- po
	gc()
	rm(po)
}

####################################################################
# scatter plot overall deaths
####################################################################

pads <- list()
for(i in 1:length(pos))
{
	cat('\nprocessing',i,' ',names(pos)[i])
	pad <- as.data.table(reshape2::melt(pos[[i]][['exp_deaths_wildtype_in_hosp']]))
	setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.deaths.wildtype'))
	tmp <- as.data.table(reshape2::melt(pos[[i]][['exp_deaths_P1_in_hosp']]))
	setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.deaths.P1'))
	pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
	pad[, exp.deaths := exp.deaths.wildtype + exp.deaths.P1]
	tmp <- unique(subset(pad, is.na(exp.deaths), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pad <- subset(pad, !iterations%in%tmp)
	}
	pad <- pad[, list(value=quantile(exp.deaths, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx','age.idx')]
	pad[, loc_label := names(pos)[i]]
	pads[[i]] <- pad
}
pads <- do.call('rbind',pads)
tmp <- subset(dh, select=c(loc_label, age.label, week.start, week.idx, age.idx, deaths_of_hospsi))
pads <- merge(tmp, pads, by=c('loc_label','week.idx','age.idx'))
pads <- dcast.data.table(pads, loc_label+week.idx+week.start+age.idx+age.label+deaths_of_hospsi~stat, value.var='value')

p <- ggplot(pads, aes(x=deaths_of_hospsi)) +	
		geom_abline(intercept=0, slope=1) +
		geom_errorbar(aes(ymin=CL, ymax=CU), colour='grey50') +
		geom_point(aes(y=M, colour=age.label, shape=age.label)) +
		viridis::scale_colour_viridis(option='B', end=0.8, discrete = TRUE) +
		scale_shape_manual(values=seq_along(levels(pads$age.label))) +
		scale_x_continuous(expand=c(0,0)) +
		scale_y_continuous(expand=c(0,0)) +		
		facet_wrap(~change_city_label(loc_label), ncol=5, scales='free') +
		theme_bw() +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		theme(strip.background = element_blank(),
				legend.position='bottom') +
		labs(colour='age band', pch='age band', y='fitted weekly COVID-19 deaths', x='observed weekly COVID-19 deaths')
ggsave(file=paste0(args$out.base,'_allloc_deaths_ndeaths_fit_scatter.pdf'),p,w=10,h=10)
ggsave(file=paste0(args$out.base,'_allloc_deaths_ndeaths_fit_scatter.png'),p,w=10,h=10)


####################################################################
# scatter plot overall hospitalisations
####################################################################

pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pad <- as.data.table(reshape2::melt(pos[[i]][['exp_hosp_adm_wildtype']]))
	setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.hosp.adm.wildtype'))
	tmp <- as.data.table(reshape2::melt(pos[[i]][['exp_hosp_adm_P1']]))
	setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.hosp.adm.P1'))
	pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
	pad[, exp.hosp.adm := exp.hosp.adm.wildtype + exp.hosp.adm.P1]
	tmp <- unique(subset(pad, is.na(exp.hosp.adm), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pad <- subset(pad, !iterations%in%tmp)
	}	
	pad <- pad[, list(value=quantile(exp.hosp.adm, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx','age.idx')]
	pad[, loc_label := names(pos)[i]]
	pads[[i]] <- pad
}
pads <- do.call('rbind',pads)
tmp <- unique(subset(dh, select=c(loc_label, week.idx,age.idx,week,week.start,age.label,hosps)))
pads <- merge(tmp, pads, by=c('loc_label','week.idx','age.idx'))
pads <- dcast.data.table(pads, loc_label+week.idx+week+week.start+age.idx+age.label+hosps~stat, value.var='value')

p <- ggplot(pads, aes(x=hosps)) +	
		geom_abline(intercept=0, slope=1) +
		geom_errorbar(aes(ymin=CL, ymax=CU), colour='grey50') +
		geom_point(aes(y=M, colour=age.label, shape=age.label)) +
		viridis::scale_colour_viridis(option='D', end=0.8, discrete = TRUE) +
		scale_shape_manual(values=seq_along(levels(pads$age.label))) +
		scale_x_continuous(expand=c(0,0)) +
		scale_y_continuous(expand=c(0,0)) +
		facet_wrap(~change_city_label(loc_label), ncol=5, scales='free') +
		theme_bw() +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		theme(strip.background = element_blank(),
				legend.position='bottom') +
		labs(colour='age band', shape= 'age band', y='fitted weekly COVID-19 hospital admissions', x='actual COVID-19 hospital admissions')
ggsave(file=paste0(args$out.base,'_allloc_hospn_fit_scatter.pdf'),p,w=10,h=10)
ggsave(file=paste0(args$out.base,'_allloc_hospn_fit_scatter.png'),p,w=10,h=10)

####################################################################
# compare age distribution of hospital admissions P1 vs non-P1
####################################################################

pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]][['age_prop_hosps_wildtype']]))
	setnames(pa, 1:4, c('iterations','week.idx','age.idx','prop.age.in.hosps.wildtype'))
	tmp <- as.data.table(reshape2::melt(pos[[i]][['age_prop_hosps_P1']]))
	setnames(tmp, 1:4, c('iterations','week.idx','age.idx','prop.age.in.hosps.P1'))
	pa <- merge(pa,tmp,by=c('iterations','week.idx','age.idx'))
	pa[, ratio.P1.wildtype := prop.age.in.hosps.P1/prop.age.in.hosps.wildtype]
	pa <- melt(pa, id.vars=c('iterations','week.idx','age.idx'))
	setkey(pa, variable, iterations, week.idx, age.idx)
	tmp <- pa[variable!='ratio.P1.wildtype', list(age.idx=age.idx, value=cumsum(value)), by=c('variable','iterations','week.idx')]
	set(tmp, NULL, 'variable', tmp[,paste0('c',variable)])
	pa <- rbind(pa,tmp)
	tmp <- unique(subset(pa, is.na(value), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}		
	pa <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','week.idx','age.idx')]
	pa[, loc_label := names(pos)[i]]
	pads[[i]] <- pa
	
}
pads <- do.call('rbind',pads)
tmp <- unique(subset(dh, select=c(loc_label,age.idx, week.idx, week, week.start, age.label)))
pads <- merge(pads, tmp, by=c('loc_label','week.idx','age.idx'))
pads <- dcast.data.table(pads, variable+loc_label+age.idx+age.label+week+week.idx+week.start~stat, value.var='value')
pads[, variant := factor(grepl('P1', variable), levels=c(FALSE,TRUE), labels=c('non-P1','P1'))]
pads[, plot.facet := gsub('^([a-z]+)\\..*','\\1',variable)]

tmp <- subset(pads, week.idx==1 & plot.facet=='ratio')
tmp <- tmp[CU-CL<=3,]

p <- ggplot(tmp, aes(x=age.label)) +
		geom_hline(yintercept=1, colour='grey50') +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU, fill=change_city_label(loc_label)), stat='identity') +
		scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +					
		labs(x='',y='ratio of Gamma vs non-Gamma\nage composition of COVID-19 hospital admissions',fill='')
ggsave(file=paste0(args$out.base,'_allloc_hosps_ratio_age_composition.pdf'),p,w=10,h=5)
ggsave(file=paste0(args$out.base,'_allloc_hosps_ratio_age_composition.png'),p,w=10,h=5)


####################################################################
# compare age distribution of deaths P1 vs non-P1
####################################################################

pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]][['age_prop_deaths_wildtype']]))
	setnames(pa, 1:4, c('iterations','week.idx','age.idx','prop.age.in.deaths.wildtype'))
	tmp <- as.data.table(reshape2::melt(pos[[i]][['age_prop_deaths_P1']]))
	setnames(tmp, 1:4, c('iterations','week.idx','age.idx','prop.age.in.deaths.P1'))
	pa <- merge(pa,tmp,by=c('iterations','week.idx','age.idx'))
	pa[, ratio.P1.wildtype := prop.age.in.deaths.P1/prop.age.in.deaths.wildtype]
	pa <- melt(pa, id.vars=c('iterations','week.idx','age.idx'))
	setkey(pa, variable, iterations, week.idx, age.idx)
	tmp <- pa[variable!='ratio.P1.wildtype', list(age.idx=age.idx, value=cumsum(value)), by=c('variable','iterations','week.idx')]
	set(tmp, NULL, 'variable', tmp[,paste0('c',variable)])
	pa <- rbind(pa,tmp)
	tmp <- unique(subset(pa, is.na(value), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}		
	pa <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','week.idx','age.idx')]
	pa[, loc_label := names(pos)[i]]
	pads[[i]] <- pa
	
}
pads <- do.call('rbind',pads)
tmp <- unique(subset(dh, select=c(loc_label,age.idx, week.idx, week, week.start, age.label)))
pads <- merge(pads, tmp, by=c('loc_label','week.idx','age.idx'))
pads <- dcast.data.table(pads, variable+loc_label+age.idx+age.label+week+week.idx+week.start~stat, value.var='value')
pads[, variant := factor(grepl('P1', variable), levels=c(FALSE,TRUE), labels=c('non-P1','P1'))]
pads[, plot.facet := gsub('^([a-z]+)\\..*','\\1',variable)]

tmp <- subset(pads, week.idx==1 & plot.facet=='ratio')
tmp <- tmp[CU-CL<=3,]

p <- ggplot(tmp, aes(x=age.label)) +
		geom_hline(yintercept=1, colour='grey50') +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU, fill=change_city_label(loc_label)), stat='identity') +
		scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +					
		labs(x='',y='ratio of Gamma vs non-Gamma\nage composition of COVID-19 deaths',fill='')
ggsave(file=paste0(args$out.base,'_allloc_deaths_ratio_age_composition.pdf'),p,w=10,h=5)
ggsave(file=paste0(args$out.base,'_allloc_deaths_ratio_age_composition.png'),p,w=10,h=5)


####################################################################
# fitted HFR by variant in ref.week
####################################################################

pads <- list()
for(i in 1:length(pos))
{		
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]][['logit_hosp_fatality_ratio_wildtype']]))
	setnames(pa, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_wildtype'))
	tmp <- as.data.table(reshape2::melt(pos[[i]][['logit_hosp_fatality_ratio_P1_rnde']]))
	setnames(tmp, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_P1_rnde'))
	pa <- merge(pa,tmp,by=c('iterations','age.idx'))
	pa[, hosp_fatality_ratio_P1 := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype + logit_hosp_fatality_ratio_P1_rnde)]
	pa[, hosp_fatality_ratio_wildtype := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype)]
	pa[, HFR_ratio := hosp_fatality_ratio_P1 / hosp_fatality_ratio_wildtype]	
	pa <- melt(pa, id.vars=c('iterations','age.idx'), measure.vars=c('hosp_fatality_ratio_P1','hosp_fatality_ratio_wildtype','HFR_ratio'))
	tmp <- unique(subset(pa, is.na(value), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}			
	pa <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','age.idx')]
	pa[, loc_label := names(pos)[i]]	
	pads[[i]] <- pa	
}
pads <- do.call('rbind',pads)
tmp <- unique(subset(dh, select=c(age.idx, age.label)))
pads <- merge(pads, tmp, by='age.idx')
pads <- dcast.data.table(pads, loc_label+age.idx+age.label+variable~stat, value.var='value')
tmp <- dhs[, list(select.by.sample.size=prod(select.by.sample.size)), by=c('loc_label','age.label')]
pads <- merge(pads, tmp, by=c('loc_label','age.label'))

tmp <- subset(pads, variable=='HFR_ratio' & select.by.sample.size==1)
p.hfr.ratio <- ggplot(tmp, aes(x=age.label)) +
		geom_hline(yintercept=1, lwd=1, colour='grey50') +		
		geom_linerange(aes(ymin=CL, ymax=CU, colour=change_city_label(loc_label)), position=position_dodge(0.8)) +
		geom_point(aes(y=M, colour=change_city_label(loc_label)), position=position_dodge(0.8)) +
		scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
		scale_y_continuous(expand=c(0,0), breaks=1:100) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				legend.position='bottom') +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='',y='ratio of Gamma vs non-Gamma\nin-hospital fatality rates',colour='') 
#p.hfr.ratio <- ggplot(tmp, aes(x=age.label)) +
#		geom_hline(yintercept=1, lwd=1, colour='grey50') +
#		geom_boxplot(aes(ymin=CL, lower=IL, middle=M, upper=IU, ymax=CU, fill=change_city_label(loc_label)), stat='identity', width=0.7) +
#		scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
#		scale_y_continuous(expand=c(0,0), breaks=1:100) +
#		theme_bw() +
#		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +		
#		labs(x='',y='ratio of Gamma vs non-P1\nin-hospital fatality rates',fill='') 
ggsave(file=paste0(args$out.base,'_allloc_hfr_P1_vs_wildtype_ratio.pdf'),p.hfr.ratio,w=10,h=5)
ggsave(file=paste0(args$out.base,'_allloc_hfr_P1_vs_wildtype_ratio.png'),p.hfr.ratio,w=10,h=5)

tmp <- subset(pads, variable!='HFR_ratio' & select.by.sample.size==1)
tmp[, variant := factor(grepl('P1', variable), levels=c(FALSE,TRUE), labels=c('non-P1','P1'))]
tmp2 <- subset(pads, variable=='hosp_fatality_ratio_wildtype')
tmp2 <- tmp2[, list(M = sum(M)), by='loc_label']
tmp2 <- tmp2[order(M)]
tmp2[, hfr.cat := ceiling((1:nrow(tmp2))/4.7)]
set(tmp2, NULL, 'hfr.cat', tmp2[, factor(hfr.cat, levels=1:3, labels=c('low in-hospital fatality rate','mid in-hospital fatality rate','high in-hospital fatality rate'))])
set(tmp2, NULL, 'M', NULL)
tmp <- merge(tmp, tmp2, by='loc_label')

p.hfr.by.variants <- ggplot(tmp, aes(x=age.label)) +		
		geom_linerange(aes(ymin=CL, ymax=CU, colour=loc_label, lty=variant), stat='identity', position=position_dodge(.9)) +
		geom_point(aes(y=M, colour=loc_label, pch=variant), position=position_dodge(.9)) +
		scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
		scale_linetype_manual(values=c('non-P1'='solid', 'P1'='dotted')) +
		scale_y_continuous(expand=c(0,0),labels=scales::percent) +
		facet_grid(~hfr.cat) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				strip.background = element_blank()) +				
		labs(x='',y='in-hospital fatality rate',colour='',pch='',lty='') 
ggsave(file=paste0(args$out.base,'_allloc_hfr_P1_vs_wildtype.pdf'),p.hfr.by.variants,w=12,h=5)
ggsave(file=paste0(args$out.base,'_allloc_hfr_P1_vs_wildtype.png'),p.hfr.by.variants,w=12,h=5)



####################################################################
# ratio in fitted HFR by location in ref.week
####################################################################

pads <- list()
for(i in 1:length(pos))
{		
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]][['logit_hosp_fatality_ratio_wildtype']]))
	setnames(pa, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_wildtype'))
	pa[, hosp_fatality_ratio_wildtype := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype)]
	set(pa, NULL, 'logit_hosp_fatality_ratio_wildtype', NULL)
	tmp <- unique(subset(pa, is.na(hosp_fatality_ratio_wildtype), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}			
	pa[, loc_label := names(pos)[i]]	
	pads[[i]] <- pa	
}
pads <- do.call('rbind',pads)

#	prepare output for sep analysis
tmp <- pads[, list(value=quantile(hosp_fatality_ratio_wildtype, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label','age.idx')] 
tmp <- merge(tmp, unique(subset(dh, select=c(age.idx, age.label))), by='age.idx')
saveRDS(tmp, file=paste0(args$out.base,'_forsepanalysis_hfr_by_location_age_in_ref_week.rds'))
saveRDS(dhr, file=paste0(args$out.base,'_forsepanalysis_hfr_ref_week.rds'))

#	continue with analysis here
tmp <- pads[, list(M=median(hosp_fatality_ratio_wildtype)), by=c('loc_label','age.idx')]
tmp <- tmp[, list(M=mean(M)), by=c('loc_label')]
loc_label_min_hfr <- tmp[, loc_label[which.min(M)]]
cat('\nLocation with min in-hospital fatality rate in ref.week is ',loc_label_min_hfr)
tmp <- subset(pads, loc_label==loc_label_min_hfr, select=-loc_label)
setnames(tmp, 'hosp_fatality_ratio_wildtype', 'ref_hosp_fatality_ratio_wildtype')
pads <- merge(pads, tmp, by=c('iterations','age.idx'))
pads[, hosp_fatality_ratio := hosp_fatality_ratio_wildtype/ref_hosp_fatality_ratio_wildtype]
# For sexpr:
# tmp <- pads[, list(value=quantile(hosp_fatality_ratio, p=c(.025,.25,.5,.75,.975)), stat=c('CL', 'IL', 'M', 'IU', 'CU'))]
# tmp1 <- dpop[,unique(tpop), by = loc_label][,list(loc_label=loc_label, prop=V1/sum(V1))] NO: BH not counted
# tmp <- merge(pads, tmp1, by = 'loc_label')[,sum(ref_hosp_fatality_ratio_wildtype*prop) ,by = c('iterations', 'age.idx')]
# tmp <- tmp[, list(value=quantile(V1, p=c(.025,.25,.5,.75,.975)), stat=c('CL', 'IL', 'M', 'IU', 'CU'))]
pads <- pads[, list(value=quantile(hosp_fatality_ratio, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label','age.idx')]
pads <- dcast.data.table(pads, loc_label+age.idx~stat, value.var='value')
tmp <- unique(subset(dh, select=c(age.idx, age.label)))
pads <- merge(pads, tmp, by='age.idx')
tmp <- subset(dhs, variant=='non-Gamma', select=c(loc_label,age.label,select.by.sample.size))
pads <- merge(pads, tmp, by=c('loc_label','age.label'))

tmp <- subset(pads, select.by.sample.size==1 & loc_label!=loc_label_min_hfr)
p <- ggplot(tmp, aes(x=age.label)) +
		geom_hline(yintercept=1, lwd=1, colour='grey50') +
		#geom_boxplot(aes(ymin=CL, lower=IL, middle=M, upper=IU, ymax=CU, fill=change_city_label(loc_label)), stat='identity', width=0.7) +
		geom_linerange(aes(ymin=CL, ymax=CU, colour=change_city_label(loc_label)), position=position_dodge(0.8)) +
		geom_point(aes(y=M, colour=change_city_label(loc_label)), position=position_dodge(0.8)) +
		scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
		scale_y_continuous(expand=c(0,0), breaks=1:100) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				legend.position='bottom') +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='',y=paste0('ratio of in-hospital fatality rates in baseline week\nrelative to reference location ',loc_label_min_hfr),colour='') 
ggsave(file=paste0(args$out.base,'_allloc_hfr_location_ratio.pdf'),p,w=10,h=5)
ggsave(file=paste0(args$out.base,'_allloc_hfr_location_ratio.png'),p,w=10,h=5)


####################################################################
# expected deaths by age P1 and wildtype
####################################################################

pads <- list()
for(i in 1:length(pos))
{
	cat('\nprocessing',i,' ',names(pos)[i])
	pad <- as.data.table(reshape2::melt(pos[[i]][['exp_deaths_wildtype_in_hosp']]))
	setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.deaths.wildtype'))
	tmp <- as.data.table(reshape2::melt(pos[[i]][['exp_deaths_P1_in_hosp']]))
	setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.deaths.P1'))
	pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
	pad <- melt(pad, id.vars=c('iterations','week.idx','age.idx'), measure.vars=c('exp.deaths.wildtype','exp.deaths.P1'))
	tmp <- unique(subset(pad, is.na(value), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pad <- subset(pad, !iterations%in%tmp)
	}				
	pad <- pad[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('variable','week.idx','age.idx')]
	pad[, loc_label := names(pos)[i]]
	pads[[i]] <- pad
}
pads <- do.call('rbind',pads)
pads <- dcast.data.table(pads, loc_label+variable+week.idx+age.idx~stat, value.var='value')
pads[, variant := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\3',variable)]
set(pads, NULL, 'variant', pads[, factor(variant, levels=c('wildtype', 'P1'), labels=c('non-Gamma','Gamma'))])
pads[, variable := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\1.\\2',variable)]
tmp <- unique(subset(dh, select=c(loc_label,week.idx,age.idx,week,week.start,age.label))) # changed dd with dh, alternatively could add week,idx to dd at the start of the script
pads <- merge(tmp, pads, by=c('loc_label','week.idx','age.idx'))
setkey(pads, loc_label, variable, week.idx, variant, age.idx)
tmp <- pads[, list(variant=variant, age.idx=age.idx, sM=cumsum(M)), by=c('loc_label','variable','week.idx')]
pads <- merge(pads, tmp, by=c('loc_label','variable','week.idx','variant','age.idx'))
tmp <- dh[,list(tdeaths_of_hospsi=sum(deaths_of_hospsi)), by=c('loc_label','week.start')]
p <- ggplot(pads, aes(x=week.start)) +
		geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
		geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=age.label, alpha=variant)) +
		geom_point(data=tmp, aes(y=tdeaths_of_hospsi), pch=18) +
		viridis::scale_fill_viridis(option='B', end=0.8, discrete = TRUE) +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		scale_alpha_manual(values=c('non-Gamma'=0.7, 'Gamma'=1.)) +
		facet_wrap(~change_city_label(loc_label), scales='free', ncol=4) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
		labs(fill='age', y='estimated in-hospital COVID19 deaths by variant\n(posterior median)', x='week of hospital admission', alpha='')
ggsave(file=paste0(args$out.base,'_allloc_deaths_by_variant.pdf'),p, w=10,h=12)
ggsave(file=paste0(args$out.base,'_allloc_deaths_by_variant.png'),p, w=10,h=12)


p.deaths.by.variant <- ggplot(subset(pads, loc_label%in%c('belo horizonte','manaus','sao paulo')), aes(x=week.start)) +
		geom_vline(data=subset(ddates, loc_label%in%c('belo horizonte','manaus','sao paulo')), aes(xintercept=w2.start), colour='grey50') +
		geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=age.label, alpha=variant)) +
		geom_point(data=subset(tmp, loc_label%in%c('belo horizonte','manaus','sao paulo')), aes(y=tdeaths_of_hospsi), pch=18) +
		viridis::scale_fill_viridis(option='B', end=0.8, discrete = TRUE) +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		scale_alpha_manual(values=c('non-Gamma'=0.7, 'Gamma'=1.)) +
		facet_wrap(~change_city_label(loc_label), scales='free', nrow=1) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
			strip.background = element_blank(),
			legend.position='bottom') +
		guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
		labs(fill='age', y='estimated COVID19 deaths by variant\n(posterior median)', x='', alpha='')
ggsave(file=paste0(args$out.base,'_allloc_deaths_by_variant_v2.pdf'),p.deaths.by.variant,w=9,h=5)
ggsave(file=paste0(args$out.base,'_allloc_deaths_by_variant_v2.png'),p.deaths.by.variant,w=9,h=5)

####################################################################
# expected hospital admissions by age P1 and wildtype
####################################################################

pads <- list()
for(i in 1:length(pos))
{
	cat('\nprocessing',i,' ',names(pos)[i])
	pad <- as.data.table(reshape2::melt(pos[[i]][['exp_hosp_adm_wildtype']]))
	setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.hosps.wildtype'))
	tmp <- as.data.table(reshape2::melt(pos[[i]][['exp_hosp_adm_P1']]))
	setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.hosps.P1'))
	pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
	pad <- melt(pad, id.vars=c('iterations','week.idx','age.idx'), measure.vars=c('exp.hosps.wildtype','exp.hosps.P1'))
	tmp <- unique(subset(pad, is.na(value), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pad <- subset(pad, !iterations%in%tmp)
	}				
	pad <- pad[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('variable','week.idx','age.idx')]
	pad[, loc_label := names(pos)[i]]
	pads[[i]] <- pad
}
pads <- do.call('rbind',pads)
pads <- dcast.data.table(pads, loc_label+variable+week.idx+age.idx~stat, value.var='value')
pads[, variant := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\3',variable)]
set(pads, NULL, 'variant', pads[, factor(variant, levels=c('wildtype', 'P1'), labels=c('non-Gamma','Gamma'))])
pads[, variable := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\1.\\2',variable)]
tmp <- unique(subset(dh, select=c(loc_label,week.idx,age.idx,week,week.start,age.label))) # changed dd with dh, alternatively could add week,idx to dd at the start of the script
pads <- merge(tmp, pads, by=c('loc_label','week.idx','age.idx'))
setkey(pads, loc_label, variable, week.idx, variant, age.idx)
tmp <- pads[, list(variant=variant, age.idx=age.idx, sM=cumsum(M)), by=c('loc_label','variable','week.idx')]
pads <- merge(pads, tmp, by=c('loc_label','variable','week.idx','variant','age.idx'))
tmp <- unique(subset(dh, select=c(loc_label,week.start,thosps)))
p <- ggplot(pads, aes(x=week.start)) +
		geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
		geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=age.label, alpha=variant)) +
		geom_point(data=tmp, aes(y=thosps), pch=18) +
		viridis::scale_fill_viridis(option='D', end=0.8, discrete = TRUE) +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		scale_alpha_manual(values=c('non-Gamma'=0.7, 'Gamma'=1.)) +
		facet_wrap(~change_city_label(loc_label), scales='free', ncol=4) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
		labs(fill='age', y='estimated COVID-19 hospital admissions by variant\n(posterior median)', x='week of hospital admission', alpha='')
ggsave(file=paste0(args$out.base,'_allloc_hosps_by_variant.pdf'),p, w=10,h=12)
ggsave(file=paste0(args$out.base,'_allloc_hosps_by_variant.png'),p, w=10,h=12)


p.deaths.by.variant <- ggplot(subset(pads, loc_label%in%c('belo horizonte','manaus','sao paulo')), aes(x=week.start)) +
		geom_vline(data=subset(ddates, loc_label%in%c('belo horizonte','manaus','sao paulo')), aes(xintercept=w2.start), colour='grey50') +
		geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=age.label, alpha=variant)) +
		geom_point(data=subset(tmp, loc_label%in%c('belo horizonte','manaus','sao paulo')), aes(y=thosps), pch=18) +
		viridis::scale_fill_viridis(option='D', end=0.8, discrete = TRUE) +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		scale_alpha_manual(values=c('non-Gamma'=0.7, 'Gamma'=1.)) +
		facet_wrap(~change_city_label(loc_label), scales='free', nrow=1) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
		labs(fill='age', y='estimated COVID-19 hospital admissions by variant\n(posterior median)', x='week of hospital admission', alpha='')
ggsave(file=paste0(args$out.base,'_allloc_hosps_by_variant_v2.pdf'),p.deaths.by.variant,w=9,h=5)
ggsave(file=paste0(args$out.base,'_allloc_hosps_by_variant_v2.png'),p.deaths.by.variant,w=9,h=5)


####################################################################
# fitted prop P1 vs empirical
####################################################################

pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pad <- as.data.table(reshape2::melt(pos[[i]][['prop_P1']]))
	setnames(pad, 1:3, c('iterations','week.idx','value'))
	tmp <- unique(subset(pad, is.na(value), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pad <- subset(pad, !iterations%in%tmp)
	}	
	pad <- pad[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx')]
	pad <- dcast.data.table(pad, week.idx~stat, value.var='value')
	pad[, loc_label := names(pos)[i]]
	pads[[i]] <- pad
}
pads <- do.call('rbind',pads)
tmp <- unique(subset(dh, select=c(loc_label,week.idx, week, week.start)))
pads <- merge(pads, tmp, by=c('loc_label','week.idx'))
tmp <- subset(dp, select=c(loc_label, week, n_positive, n_sampled))
pads <- merge(pads, tmp, by=c('loc_label','week'), all.x=TRUE)

p <- ggplot(pads, aes(x=week.start)) +
		geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey50') +
		geom_line(aes(y=M)) +
		geom_point(data=subset(pads, !is.na(n_sampled)), aes(y=n_positive/n_sampled, size=n_sampled, fill=loc_label), pch=21) +
		scale_fill_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		facet_wrap(~change_city_label(loc_label), scales='free_y', ncol=5) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(fill='none') +
		labs(y='Gamma genotype frequency from GISAID metadata', x='week sampled in infected individuals', colour='', size='number of SARS-CoV-2 sequences')
ggsave(file=paste0(args$out.base,'_allloc_seq_propP1_over_time.pdf'),p,w=10,h=8)
ggsave(file=paste0(args$out.base,'_allloc_seq_propP1_over_time.png'),p,w=10,h=8)

tmp <- subset(pads, loc_label%in%c('belo horizonte','manaus','sao paulo'))
if(nrow(tmp))
{
  p.propP1 <- ggplot(tmp, aes(x=week.start)) +
    geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey50') +
    geom_line(aes(y=M)) +
    geom_point(data=subset(tmp, !is.na(n_sampled)), aes(y=n_positive/n_sampled, size=n_sampled, fill=loc_label), pch=21) +
    scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), scales='free_y', nrow=1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    guides(fill='none') +
    labs(y='Gamma genotype frequency from GISAID metadata', x='week sampled in infected individuals', colour='', size='number of SARS-CoV-2 sequences')
  ggsave(file=paste0(args$out.base,'_allloc_seq_propP1_over_time_v2.pdf'),p.propP1,w=9,h=5)
  ggsave(file=paste0(args$out.base,'_allloc_seq_propP1_over_time_v2.png'),p.propP1,w=9,h=5)
}

###############################################################################
# fitted HFR by variant and over time
##############################################################################

pads <- list()
for(i in 1:length(pos))
{	
  cat('\nprocessing',i,' ',names(pos)[i])
  pa <- as.data.table(reshape2::melt(pos[[i]]$hfr_wildtype))
  setnames(pa, 1:4, c('iterations','week.idx','age.idx','hfr.wildtype'))
  tmp <- as.data.table(reshape2::melt(pos[[i]]$hfr_overall))
  setnames(tmp, 1:4, c('iterations','week.idx','age.idx','hfr.overall'))
  pa <- merge(pa,tmp,by=c('iterations','week.idx','age.idx'))
  pa <- melt(pa, id.vars=c('iterations','week.idx','age.idx'), measure.vars=c('hfr.wildtype','hfr.overall'))
  tmp <- unique(subset(pa, is.na(value), select=iterations))[,iterations]
  if(length(tmp))
  {
	  cat('\nFound iterations with NA, removing:', tmp)
	  pa <- subset(pa, !iterations%in%tmp)
  }
  pas <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','age.idx','week.idx')]
  pas <- dcast.data.table(pas, variable + age.idx + week.idx  ~ stat, value.var='value')
  tmp <- unique(subset(dh, loc_label == names(pos)[[i]], select=c(age.idx, week.idx, week, week.start, age.label, hosps, deaths_of_hosps)))
  pas <- merge(pas, tmp, by=c('age.idx','week.idx'))
  pas[, loc_label := names(pos)[i]]
  pads[[i]] <- pas
  rm(pas,pa, tmp)
}
pads <- do.call('rbind',pads)
tmp <- subset(pads, variable=='hfr.overall' & age.label=='40-49')
tmp2 <- subset(pads, variable=='hfr.wildtype' & age.label=='40-49')

p <- ggplot(tmp, aes(x=week.start)) +
  geom_ribbon(aes(ymin=CL, ymax=CU, fill = loc_label),alpha=.3) +
  geom_line(data=tmp2, aes(y=M),colour='black',linetype='11') +
  geom_line(aes(y=M,colour = loc_label)) +
  scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
  scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  scale_x_date(breaks='4 months',expand=c(0,0)) +    
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.background = element_blank(),
        legend.position='bottom') +
  facet_wrap(~change_city_label(loc_label), ncol=5) +
  labs(x='',y='weekly in-hospital fatality rate in age group 40-49', color = '') + 
  guides(fill = FALSE)
ggsave(file=paste0(args$out.base,'_allloc_hfr_over_time_by_variant.png'),p,w=12,h=10)
ggsave(file=paste0(args$out.base,'_allloc_hfr_over_time_by_variant.pdf'),p,w=12,h=10)


###############################################################################
# min and max in age standardised HFR over time
##############################################################################

dage <- dpop[,list(pop=sum(pop)),by='age.label']
dage[, tpop:=sum(pop)]
dage[, ppop.all.cities:=pop/tpop]
setkey(dage,age.label)
dage[,age.idx:=1:nrow(dage)]

pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]]$hfr_overall))
	setnames(pa, 1:4, c('iterations','week.idx','age.idx','hfr.overall'))
	tmp <- as.data.table(reshape2::melt(pos[[i]]$hfr_wildtype))
	setnames(tmp, 1:4, c('iterations','week.idx','age.idx','hfr.wildtype'))
	pa <- merge(pa, tmp, by = c('iterations','week.idx','age.idx'))
	tmp <- unique(subset(pa, is.na(hfr.overall) | is.na(hfr.wildtype), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}
	if(1) # Up maximum and minimum considered up to 
	{ # could maybe use new column from ddates, but let s not risk to f up.
	  tmp <- dh[,.(week.idx,week.start,loc_label)][loc_label == names(pos)[i], unique(.SD)]
	  tmp2 <- ddates[, .(loc_label, w1.start, w2.start)]
	  tmp <- merge(tmp, tmp2, by='loc_label')
	}
	
	pa <- merge(pa,subset(dage,select=c(age.idx,ppop.all.cities)),by='age.idx')
	pa <- pa[,list(hfr.overall=sum(hfr.overall*ppop.all.cities), hfr.wildtype=sum(hfr.wildtype*ppop.all.cities)),by=c('week.idx','iterations')]
  
	pa$beforeP1 <- FALSE
	pa$before8months <- FALSE
	pa[week.idx %in% tmp[week.start <= w2.start, week.idx], beforeP1 := TRUE]
	pa[week.idx %in% tmp[week.start <= w1.start+35*7, week.idx], before8months := TRUE]
	
	# before P1 or before X months???
	tmp <- c('iterations', 'before8months')
	pa <- pa[,list(hfr.overall.min=min(hfr.overall),hfr.overall.max=max(hfr.overall),
	               hfr.wildtype.min=min(hfr.wildtype),hfr.wildtype.max=max(hfr.wildtype)),by=tmp]
	pa <- pa[, list(hfr.overall.min = hfr.overall.min[1], hfr.overall.max=max(hfr.overall.max),
	           hfr.wildtype.min=min(hfr.wildtype.min[1]),hfr.wildtype.max=max(hfr.wildtype.max)),by=c('iterations')]
	pa[, `:=` (hfr.overall.diff = hfr.overall.max-hfr.overall.min,
	           hfr.wildtype.diff = hfr.wildtype.max-hfr.wildtype.min )]
	pa <- melt(pa, id.vars='iterations')
	pa <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable')]	
	pa <- dcast.data.table(pa, variable ~ stat, value.var='value')
	pa[, loc_label := names(pos)[i]]
	pads[[i]] <- pa
	rm(pa, tmp)
}
pads <- do.call('rbind',pads)

pads[, LABEL:= paste0(format(round(M*100,1),nsmall=1),'% (',format(round(CL*100,1),nsmall=1),'%-',format(round(CU*100,1),nsmall=1),'%)')  ]
pas.hfr.range <- copy(pads)
pas.hfr.range[, age.idx:=0L]
pas.hfr.range[, age.label:='overall']
PAS.HFR.RANGE <- copy(pas.hfr.range)

###############################################################################
# age standardised HFR over time
##############################################################################

pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]]$hfr_overall))
	setnames(pa, 1:4, c('iterations','week.idx','age.idx','hfr.overall'))
	tmp <- unique(subset(pa, is.na(hfr.overall), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}
	pa <- merge(pa,subset(dage,select=c(age.idx,ppop.all.cities)),by='age.idx')
	pa <- pa[,list(hfr.overall=sum(hfr.overall*ppop.all.cities)),by=c('week.idx','iterations')]	
	pa <- pa[, list(value=quantile(hfr.overall, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('week.idx')]	
	pa <- dcast.data.table(pa, week.idx ~ stat, value.var='value')
	pa[, loc_label := names(pos)[i]]
	pads[[i]] <- pa
	rm(pa, tmp)
}
pads <- do.call('rbind',pads)
tmp <- unique(subset(dh, select=c(loc_label, week.idx, week, week.start)))
pads <- merge(pads,tmp,by=c('loc_label','week.idx'))

# Here no need to add time constraint on Gamma s emergence
# This is because ref week has been constrained in f2! FALSE!

tmp <- subset(dhr, select=c(loc_label,ref.week))
pads <- merge(pads,tmp,by='loc_label')
tmp <- ddates[, .(loc_label, w1.start)][, before8months:=w1.start+35*7]
pads[, months_threshold:=FALSE]
pads <- merge(pads,tmp,by='loc_label')[week.start<=before8months, months_threshold:=TRUE]

tmp <- pads[months_threshold==TRUE, list(min.week=week[which.min(M)], min.hfr=min(M)), by='loc_label']
pads <- merge(pads,tmp,by='loc_label')
dhr <- merge(dhr,tmp,by='loc_label')
setkey(dhr, min.hfr)
dhr[, loc_label_group := ceiling((1:nrow(dhr))/5)]
set(dhr,NULL,'loc_label_group',dhr[,factor(loc_label_group,levels=unique(loc_label_group),labels=c('lower tertile','mid tertile','higher tertile'))])
pads <- merge(pads,subset(dhr,select=c(loc_label,loc_label_group)),by='loc_label')
pads <- merge(pads,subset(dstates,select=c(loc_label,region_label)),by='loc_label')
set(pads,pads[,which(region_label%in%c('Central-West','North'))],'region_label','North + Central-West')

p <- ggplot(pads, aes(x=week.start)) +
		geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
		geom_vline(data=subset(pads,week==min.week), aes(xintercept=week.start), colour='black', lty='11') +
		geom_hline(data=subset(pads,week==min.week), aes(yintercept=M), colour='black', lty='11') +
		geom_ribbon(aes(ymin=CL, ymax=CU, fill=loc_label),alpha=.3) +		
		geom_line(aes(y=M,colour=loc_label)) +
		scale_colour_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_fill_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +    
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		facet_wrap(~change_city_label(loc_label), ncol=5) +
		labs(x='',y='weekly age-standardised in-hospital fatality rate', color='', fill='') + 
		guides(fill='none',colour='none')
ggsave(file=paste0(args$out.base,'_allloc_hfr_agestd_over_time.png'),p,w=12,h=10)
ggsave(file=paste0(args$out.base,'_allloc_hfr_agestd_over_time.pdf'),p,w=12,h=10)

p <- ggplot(pads, aes(x=week.start)) +
		geom_hline(data=subset(pads,week==min.week), aes(yintercept=M,colour=change_city_label(loc_label)), lty='11') +
		geom_ribbon(aes(ymin=CL, ymax=CU, fill=change_city_label(loc_label)),alpha=.3) +		
		geom_line(aes(y=M,colour=change_city_label(loc_label))) +
		scale_colour_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_fill_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +    
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		facet_grid(~region_label) +
		labs(x='',y='weekly age-standardised in-hospital fatality rate', color='', fill='') + 
		guides(fill='none',colour=guide_legend(nrow=2,byrow=TRUE))
ggsave(file=paste0(args$out.base,'_allloc_hfr_agestd_over_time_grouped.png'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_allloc_hfr_agestd_over_time_grouped.pdf'),p,w=10,h=6)


###############################################################################
# deaths that could have been averted with large resources 
##############################################################################

pads <- list()
pads2 <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]]$hfr_overall))
	setnames(pa, 1:4, c('iterations','week.idx','age.idx','hfr.overall'))
	pa[, loc_label := names(pos)[i]]
	tmp <- subset(dhr,select=c(loc_label,min.week))
	tmp2 <- unique(subset(dh,select=c(loc_label,week.idx,week)))
	setnames(tmp2,c('week.idx','week'),c('min.week.idx','min.week'))
	tmp <- merge(tmp,tmp2, by=c('min.week','loc_label'))	
	pa <- merge(pa,tmp,by='loc_label')
	pa <- subset(pa,week.idx==min.week.idx)
	set(pa,NULL,c('week.idx','min.week','min.week.idx'),NULL)	
	tmp <- unique(subset(pa, is.na(hfr.overall), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}
	tmp <- subset(dh,select=c(loc_label,week.idx,age.idx,hosps,deaths_of_hosps))
	pa <- merge(pa,tmp,by=c('loc_label','age.idx'),allow.cartesian=TRUE)
	pa[, hyp_deaths_of_hosps := hosps*hfr.overall]
	pa <- pa[,list(deaths_of_hosps=sum(deaths_of_hosps), hyp_deaths_of_hosps=sum(hyp_deaths_of_hosps)),by=c('iterations','loc_label','age.idx')]
	tmp <- pa[,list(deaths_of_hosps=sum(deaths_of_hosps), hyp_deaths_of_hosps=sum(hyp_deaths_of_hosps)),by=c('iterations','loc_label')]
	tmp[,age.idx:=0L]
	pa <- rbind(pa,tmp)	
	pa[, hyp_deaths_of_hosps_diff := deaths_of_hosps-hyp_deaths_of_hosps]
	pa[, hyp_deaths_of_hosps_pcred := 1-hyp_deaths_of_hosps/deaths_of_hosps]
	pa2 <- pa[age.idx==0,,]
	pa <- melt(pa,id.vars=c('iterations','loc_label','age.idx'))
	pa <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','age.idx','loc_label')]
	pa <- dcast.data.table(pa, loc_label + variable + age.idx  ~ stat, value.var='value')
	pads[[i]] <- pa
	pads2[[i]] <- pa2
	rm(pa, pa2, tmp)
}
pads <- do.call('rbind',pads)

# Compute quantiles for distribution of hyp deaths in the ensemble of cities
pads2 <- do.call('rbind', pads2)
pads2 <- pads2[,list(deaths_of_hosps=sum(deaths_of_hosps),
                     hyp_deaths_of_hosps=sum(hyp_deaths_of_hosps)),by = 'iterations']
pads2[, hyp_deaths_of_hosps_diff := deaths_of_hosps-hyp_deaths_of_hosps]
pads2[, hyp_deaths_of_hosps_pcred := 1-hyp_deaths_of_hosps/deaths_of_hosps]
pads2 <- melt(pads2,id.vars=c('iterations'))
pads2 <- pads2[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable')]
pads2 <- dcast.data.table(pads2, variable ~ stat, value.var='value')
pads2[, `:=` (loc_label="Total", age.idx=0, age.label='overall')]

# Back to normal
tmp <- unique(subset(dh, select=c(loc_label, age.idx, age.label)))
pads <- merge(pads, tmp, by=c('loc_label','age.idx'),all.x=TRUE)
set(pads,pads[,which(is.na(age.label))],'age.label','overall')


tmp <- subset(pads,variable=='hyp_deaths_of_hosps_diff')
tmp <- tmp[age.label!= 'overall']
p <- ggplot(tmp, aes(x=change_city_label(loc_label), fill=change_city_label(loc_label))) +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
		scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
		facet_wrap(~age.label,scales='free_y',ncol=4) +
		theme_bw() +
		theme(strip.background = element_blank(),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				panel.spacing = unit(1, "lines"),
				legend.position='bottom') +
		labs(x='', y='excess COVID-19 attributable deaths relative to minimum fatality rate', fill='') +
		guides(fill='none')
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_excessdeaths.pdf'),p,w=12,h=12)
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_excessdeaths.png'),p,w=12,h=12)

tmp <- subset(pads,variable=='hyp_deaths_of_hosps_pcred')
tmp <- tmp[age.label!= 'overall']
p <- ggplot(tmp, aes(x=change_city_label(loc_label), fill=change_city_label(loc_label))) +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
		scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
		facet_wrap(~age.label,scales='free_y',ncol=4) +
		scale_y_continuous(labels=scales::percent) +
		theme_bw() +
		theme(strip.background = element_blank(),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				panel.spacing = unit(1, "lines"),
				legend.position='bottom') +
		labs(x='', y='projected avoidable COVID-19 attributable deaths in hospitals\n in the absence of pandemic healthcare pressure', fill='') +
		guides(fill='none')
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_pcred.pdf'),p,w=12,h=12)
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_pcred.png'),p,w=12,h=12)

pas.hfr.range <- rbind(pas.hfr.range, pads, pads2, fill=TRUE)
tmp <- pas.hfr.range[,which(variable=='hyp_deaths_of_hosps_pcred')]
set(pas.hfr.range,tmp,'LABEL', pas.hfr.range[tmp,paste0(format(round(M*100,1),nsmall=1),'% (',format(round(CL*100,1),nsmall=1),'%-',format(round(CU*100,1),nsmall=1),'%)')])
tmp <- pas.hfr.range[,which(variable%in%c('hyp_deaths_of_hosps_diff','hyp_deaths_of_hosps'))]
set(pas.hfr.range,tmp,'LABEL', pas.hfr.range[tmp,paste0(format(round(M,0),nsmall=0),' (',format(round(CL,0),nsmall=0),'-',format(round(CU,0),nsmall=0),')')])
tmp <- pas.hfr.range[,which(variable=='deaths_of_hosps')]
set(pas.hfr.range,tmp,'LABEL', pas.hfr.range[tmp,paste0(format(round(M,0),nsmall=0))])

tmp <- subset(pas.hfr.range, age.idx==0 & variable%in%c('hfr.overall.min','hfr.overall.max','hfr.overall.diff','hyp_deaths_of_hosps_diff','hyp_deaths_of_hosps_pcred'))
tmp <- dcast.data.table(tmp, loc_label~variable, value.var='LABEL')
tmp[, loc_label:=change_city_label(loc_label)]
saveRDS(tmp, file=paste0(args$out.base,'_allloc_hfrmin_allagetable.rds'))

tmp <- subset(pas.hfr.range, variable%in%c('hyp_deaths_of_hosps_pcred'))
tmp <- dcast.data.table(tmp, loc_label~age.label, value.var='LABEL')
tmp[, loc_label:=change_city_label(loc_label)]
saveRDS(tmp, file=paste0(args$out.base,'_allloc_hfrmin_byage_pcred_table.rds'))

tmp <- subset(pas.hfr.range, variable%in%c('hyp_deaths_of_hosps_diff'))
tmp <- dcast.data.table(tmp, loc_label~age.label, value.var='LABEL')
tmp[, loc_label:=change_city_label(loc_label)]
saveRDS(tmp, file=paste0(args$out.base,'_allloc_hfrmin_byage_diff_table.rds'))

####################################################################
## Deaths that could have been averted if hfr was uniform minimum across locs
####################################################################
# reset pas.hfr.range
pas.hfr.range <- copy(PAS.HFR.RANGE)

i <- which(names(pos)==loc_label_min_hfr)
if(1)
{
  cat('Uniform minimum across locations')
  pa <- as.data.table(reshape2::melt(pos[[i]]$hfr_overall))
  setnames(pa, 1:4, c('iterations','week.idx','age.idx','hfr.overall'))
  pa[, loc_label := names(pos)[i]]
  tmp <- subset(dhr,select=c(loc_label,min.week))
  tmp2 <- unique(subset(dh,select=c(loc_label,week.idx,week)))
  setnames(tmp2,c('week.idx','week'),c('min.week.idx','min.week'))
  tmp <- merge(tmp,tmp2, by=c('min.week','loc_label'))	
  pa <- merge(pa,tmp,by='loc_label')
  pa <- subset(pa,week.idx==min.week.idx)
  set(pa,NULL,c('week.idx','min.week','min.week.idx'),NULL)	
  tmp <- unique(subset(pa, is.na(hfr.overall), select=iterations))[,iterations]
  if(length(tmp))
  {
    cat('\nFound iterations with NA, removing:', tmp)
    pa <- subset(pa, !iterations%in%tmp)
  }
  tmp <- subset(dh,select=c(loc_label,week.idx,age.idx,hosps,deaths_of_hosps))
  pa <- pa[, -"loc_label"]
  pa <- merge(tmp, pa,by=c('age.idx'), allow.cartesian = TRUE)
  pa[, hyp_deaths_of_hosps := hosps*hfr.overall]
  pa <- pa[,list(deaths_of_hosps=sum(deaths_of_hosps), hyp_deaths_of_hosps=sum(hyp_deaths_of_hosps)),by=c('iterations','loc_label','age.idx')]
  tmp <- pa[,list(deaths_of_hosps=sum(deaths_of_hosps), hyp_deaths_of_hosps=sum(hyp_deaths_of_hosps)),by=c('iterations','loc_label')]
  tmp[,age.idx:=0L]
  pa <- rbind(pa,tmp)	
  pa[, hyp_deaths_of_hosps_diff := deaths_of_hosps-hyp_deaths_of_hosps]
  pa[, hyp_deaths_of_hosps_pcred := 1-hyp_deaths_of_hosps/deaths_of_hosps]
  pa2 <- pa[age.idx==0,,]
  pa <- melt(pa,id.vars=c('iterations','loc_label','age.idx'))
  pa <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','age.idx','loc_label')]
  pa <- dcast.data.table(pa, loc_label + variable + age.idx  ~ stat, value.var='value')
}

# tmp <- pa[variable %in% c('hyp_deaths_of_hosps_diff', 'hyp_deaths_of_hosps_pcred')]
# tmp[age.idx == 0]
tmp1 <- pa[variable %in% c('hyp_deaths_of_hosps_diff','hyp_deaths_of_hosps_pcred') & age.idx == 0]
tmp1$value = '0'
tmp1[variable ==  'hyp_deaths_of_hosps_diff', value:= paste0(as.integer(M), ' (', as.integer(CL),'-',as.integer(CU),')')]
tmp1[variable ==  'hyp_deaths_of_hosps_pcred', value:= paste0(round(M*100,1), '% (', round(CL*100,1),'%-',round(CU*100,1),'%)')]
tmp1 <- dcast.data.table(tmp1, loc_label ~ variable, value.var = 'value')

tmp <- pa2[, list(hyp_deaths_of_hosps = sum(hyp_deaths_of_hosps),
                  deaths_of_hosps = sum(deaths_of_hosps))
           ,by = 'iterations']
tmp[, hyp_deaths_of_hosps_diff := deaths_of_hosps-hyp_deaths_of_hosps]
tmp[, hyp_deaths_of_hosps_pcred := 1-hyp_deaths_of_hosps/deaths_of_hosps]
tmp <- melt(tmp,id.vars=c('iterations'))
tmp <- tmp[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by='variable']
tmp <- tmp[variable %in% c('hyp_deaths_of_hosps_diff', 'hyp_deaths_of_hosps_pcred')]
tmp[variable == "hyp_deaths_of_hosps_diff", value:= as.integer(value)]
tmp[variable == "hyp_deaths_of_hosps_pcred", value:= round(value*100,2)]
tmp <- dcast.data.table(tmp,  variable  ~ stat, value.var='value')
tmp <- tmp[, lapply(.SD, as.character), .SDcols = c('variable', 'M', 'CL', 'CU')]
tmp[, lapply(.SD, function(x){gsub('^([0-9][0-9])([0-9][0-9][0-9])$','\\1,\\2',x)}), .SDcols = c('variable', 'M', 'CL', 'CU')]
tmp[variable == 'hyp_deaths_of_hosps_pcred', `:=` (M = paste0(M,'%') , CU = paste0(CU,'%'))]
tmp2 <- tmp[, paste0(M, ' (', CL, '-', CU, ')')]
names(tmp2) <- c('hyp_deaths_of_hosps_diff',  'hyp_deaths_of_hosps_pcred')

tmp2 <- list(diff = tmp2['hyp_deaths_of_hosps_diff'], pcred = tmp2['hyp_deaths_of_hosps_pcred'])
tmp[variable == 'hyp_deaths_of_hosps_pcred', CL := paste0(CL,'%')]
tmp[, value:= paste0(M, ' (', CL, '-', CU, ')')]
tmp[, loc_label := 'Total' ]
tmp <- dcast.data.table(tmp, loc_label ~ variable, value.var = 'value')
tmp <- rbind(tmp, tmp1)
tmp2$table <- tmp

saveRDS(tmp2, paste0(args$out.base,'_counterfactual_hfr_minimum_across_locations.rds') )

# Connect age-label and age.idx
tmp2 <-  unique(dh[, .(age.idx, Age_label)])
setnames(tmp2, 'Age_label', 'age.label')
tmp2 <- rbind(tmp2, data.table(age.idx=0, age.label = 'overall'))

tmp <- subset(pa,variable=='hyp_deaths_of_hosps_diff')
tmp <- merge(tmp, tmp2, by = 'age.idx')
tmp <- tmp[age.label!= 'overall']
p <- ggplot(tmp, aes(x=change_city_label(loc_label), fill=change_city_label(loc_label))) +
  geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
  scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
  facet_wrap(~age.label,scales='free_y',ncol=4) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.spacing = unit(1, "lines"),
        legend.position='bottom') +
  labs(x='', y='excess COVID-19 attributable deaths relative to minimum fatality rate', fill='') +
  guides(fill='none')
ggsave(file=paste0(args$out.base,'_allloc_hfrmin2_excessdeaths.pdf'),p,w=12,h=12)
ggsave(file=paste0(args$out.base,'_allloc_hfrmin2_excessdeaths.png'),p,w=12,h=12)

tmp <- subset(pa,variable=='hyp_deaths_of_hosps_pcred')
tmp <- merge(tmp, tmp2, by = 'age.idx')
tmp <- tmp[age.label!= 'overall']
p <- ggplot(tmp, aes(x=change_city_label(loc_label), fill=change_city_label(loc_label))) +
  geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
  scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
  facet_wrap(~age.label,scales='free_y',ncol=4) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.spacing = unit(1, "lines"),
        legend.position='bottom') +
  labs(x='', y='projected avoidable COVID-19 attributable deaths in hospitals\n in the absence of location inequities and pandemic healthcare pressure', fill='') +
  guides(fill='none')
ggsave(file=paste0(args$out.base,'_allloc_hfrmin2_pcred.pdf'),p,w=12,h=12)
ggsave(file=paste0(args$out.base,'_allloc_hfrmin2_pcred.png'),p,w=12,h=12)


####################################################################
# association of health care demand predictors with HFR
####################################################################

pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]]$fr_regcoeff_overall))
	setnames(pa, 1:3, c('iterations','p.idx','value'))
	tmp <- unique(subset(pa, is.na(value), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}	
	pas <- pa[, list(value=quantile(exp(value), p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('p.idx')]
	pas <- dcast.data.table(pas, p.idx ~ stat, value.var='value')
	tmp <- data.table(p.idx=seq_along(args$fr.predictors), p.names=args$fr.predictors)
	pas <- merge(pas, tmp, by='p.idx')
	pas[, loc_label := names(pos)[i]]
	pads[[i]] <- pas
	rm(pas,pa, tmp)
}
pads <- do.call('rbind',pads)

p <- ggplot(pads, aes(x=change_city_label(loc_label), fill=change_city_label(loc_label))) +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
		scale_fill_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		facet_grid(~p.names) +
		coord_flip() +
		theme_bw() +
		theme(strip.background = element_blank(),
				panel.spacing = unit(1, "lines"),
				legend.position='bottom') +
		labs(x='', y='regression coefficient', fill='') +
		guides(fill='none')
ggsave(file=paste0(args$out.base,'_allloc_frmult_regrcoeff.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_allloc_frmult_regrcoeff.png'),p,w=10,h=6)

p <- ggplot(pads, aes(x=p.names, fill=p.names)) +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
		scale_fill_manual(values = args$city.palette(length(unique(pads$p.names)))) +	
		facet_wrap(~change_city_label(loc_label), ncol=5) +
		coord_flip() +
		theme_bw() +
		theme(strip.background = element_blank(),
				panel.spacing = unit(1, "lines"),
				legend.position='bottom') +
		labs(x='', y='regression coefficient of standardised predictors', fill='') +
		guides(fill='none')
ggsave(file=paste0(args$out.base,'_allloc_frmult_regrcoeff_v2.pdf'),p,w=10,h=8)
ggsave(file=paste0(args$out.base,'_allloc_frmult_regrcoeff_v2.png'),p,w=10,h=8)

###############################################################################
## Significance of predictors
###############################################################################


# # Here 
# pads <- list()
# for (i in 1:length(pos))
# {
#   cat('\nprocessing',i,' ',names(pos)[i])
#   pa <- as.data.table(reshape2::melt(pos[[i]]$fr_shrinkage))
#   setnames(pa, 1:3, c('iterations','predictor.idx','value'))
#   pa[, gsub('*([[:digit:]]+)*','\\1',predictor.idx)]
#   
# }
# 


################################################################################
## relationships between health care demand predictors and population age-standardised hfr in reference week
###############################################################################

pads <- list()
for(i in 1:length(pos))
{		
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]][['logit_hosp_fatality_ratio_wildtype']]))
	setnames(pa, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_wildtype'))
	pa[, hosp_fatality_ratio_wildtype := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype)]
	set(pa, NULL, 'logit_hosp_fatality_ratio_wildtype', NULL)
	tmp <- unique(subset(pa, is.na(hosp_fatality_ratio_wildtype), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}	
	pa <- pa[, list(hosp_fatality_ratio_wildtype=mean(hosp_fatality_ratio_wildtype)), by='iterations']
	pas <- pa[, list(value=quantile(hosp_fatality_ratio_wildtype, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') )]
	pas[, loc_label := names(pos)[i]]	
	pads[[i]] <- pas
}
pads <- do.call('rbind',pads)
pads <- dcast.data.table(pads, loc_label~stat, value.var='value')
tmp <- subset(dhca, variable%in%c('p.nonres.hosps.00','r.nhosp.adm.00','r.nicu.adm.00'), select= c(variable,loc_label,week,week.start,ref.week,value))
set(tmp,NULL,'variable',tmp[,gsub('.00','',variable)])
tmp <- dcast.data.table(tmp, loc_label+week+week.start+ref.week~variable, value.var='value')
tmp <- merge(tmp,subset(dhcb, select=c(loc_label,week, r.physicians.all, r.physicians.specialist, r.icubeds.all, r.saribeds.all, r.beds.all, r.ventilator,  r.nurse, r.tech.nurse, r.physiotherapist)),by=c('loc_label','week'))
tmp <- melt(tmp, id.vars=c('loc_label','week','week.start','ref.week'))
pads <- merge(pads, tmp, by='loc_label',allow.cartesian=TRUE)
pad <- subset(pads, week==ref.week)
pad[, r.squared := NA_real_]
pad[, p.val := NA_real_]
z <- sort(unique(pad$variable))
for(i in seq_along(z))
{
	tmp <- subset(pad, variable==z[i]) 
	tmp <- summary(lm(M~value,data=tmp))	
	set(pad, pad[,which(variable==z[i])],'r.squared',tmp$r.squared)
	set(pad, pad[,which(variable==z[i])],'p.val',tmp$coefficients[2,4] )	
}
tmp <- data.table(
		variable=c('p.res.hosps','r.nhosp.adm','r.nicu.adm','r.physicians.all','r.physicians.specialist','r.icubeds.all','r.ventilator','r.nurse','r.tech.nurse','r.physiotherapist'),
		variable.label=c('prop residents in\nSARI hospital admissions','SARI hospital admissions','ICU admissions','physicians','critical care specialists','ICU beds','Ventilators','Nurses','Nurse assistants','Physiotherapist')
	)
pad <- merge(pad,tmp,by='variable')	
p <- ggplot(pad, aes(x=value)) +
		geom_smooth(aes(x=value,y=M), method='lm',formula='y~x',colour='black',lwd=0.5) +
		geom_linerange(aes(x=value,ymin=CL,ymax=CU,colour=change_city_label(loc_label))) +
		geom_point(aes(x=value,y=M,colour=change_city_label(loc_label))) +		
		geom_text(data=subset(pad, loc_label=='sao paulo'), aes(x=-Inf, y=Inf, label=paste0('p-val=',round(p.val, d=2))), colour="black", inherit.aes=FALSE, hjust=-0.2, vjust=1.5) +		
		#geom_text(aes(x= min(pad$value)*1.02, y= max(pad$value)*0.98, label= round(p.val, d=2)), hjust=0.5, vjust=0.5, colour='black') +
		scale_colour_manual(values = args$city.palette(length(unique(pad$loc_label)))) +			
		scale_y_continuous(expand=c(0,0), labels=scales::percent) +
		facet_wrap(~variable.label,scales='free',ncol=4) +
		theme_bw() +
		theme(legend.position='bottom',strip.background = element_blank()) +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='number per 100,000 population',y=paste0('in-hospital fatality rate in reference week'),colour='') 
ggsave(file=paste0(args$out.base,'_allloc_hfrrefweek_vs_predictors.pdf'),p,w=10,h=8)
ggsave(file=paste0(args$out.base,'_allloc_hfrrefweek_vs_predictors.png'),p,w=10,h=8)


####################################################################
# predicted baseline HFR at same low demand in each city 
####################################################################

#	start building data.table of all predictors as with fixed demand 
tmp.args <- copy(args)
tmp.args$make.preprocessing.plots <- 0
tmp.dicu <- copy(dicu)
tmp.dicu <- merge(tmp.dicu,unique(subset(dpop,select=c(loc_label,tpop))),by='loc_label')
set(tmp.dicu,NULL,'nicu.adm',tmp.dicu[,5*1e-5*tpop])
set(tmp.dicu,NULL,'tpop',NULL)
tmp.dhcadm <- copy(dhcadm)
tmp.dhcadm <- merge(tmp.dhcadm,unique(subset(dpop,select=c(loc_label,tpop))),by='loc_label')
set(tmp.dhcadm,NULL,'nhosp.adm',tmp.dhcadm[,10*1e-5*tpop])
set(tmp.dhcadm,NULL,'tpop',NULL)

dhr[ ,`:=` (hfr_loess.min=NULL, hfr_loess.max=NULL, hfr_loess.diff=NULL)]
tmp <- make_health_care_predictors(tmp.dhcadm, tmp.dicu, dooh, dhcnres, dhcb, dw, dr, dpop, ddates, dhr, tmp.args)
dhca2 <- copy(tmp$dhca)

# standardise across location 
tmp <- unique(subset(dhca,select=c(loc_label,variable,value_ref,value_log_ref,value_sd,value_log_sd)))
dhca2 <- merge(dhca2,tmp,by=c('loc_label','variable'))
tmp <- dhca2[, 
		list(
				loc_label = loc_label, 
				week = week,
				value_c = value-value_ref,
				value_cs = (value-value_ref)/value_sd,
				value_log = log(value+1),
				value_log_c = log(value+1)-value_log_ref,
				value_log_cs = (log(value+1)-value_log_ref)/value_log_sd
		), 
		by='variable']
dhca2 <- merge(dhca2, tmp, by=c('variable','loc_label','week'))
dhc2 <- subset(dhca2, variable%in%args$fr.predictors )
dhc2 <- dcast.data.table(dhc2, loc_label+week+week.start+week.cat~variable, value.var=args$fr.predictors.transformation)
dhc2 <- melt(dhc2, id.vars=c('loc_label','week'),measure.vars=args$fr.predictors)
setkey(dhc2,loc_label,week)
tmp <- unique(subset(dh, select=c(loc_label,week,week.idx)))
dhc2 <- merge(dhc2,tmp,by=c('loc_label','week')) 
tmp <- data.table(p.idx=seq_along(args$fr.predictors), variable=args$fr.predictors)
dhc2 <- merge(dhc2,tmp,by=c('variable'))
dhc2 <- subset(dhc2,week.idx==1)

# calculate baseline HFR, using inferred regression coefficients 
pas.ref.hfr <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]]$fr_regcoeff_overall))
	setnames(pa, 1:3, c('iterations','p.idx','coef'))	
	tmp <- unique(subset(pa, is.na(coef), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}	
	pa <- merge(pa, subset(dhc2, loc_label==names(pos)[i]), by='p.idx',allow.cartesian=TRUE)
	pa <- pa[, list(logit_fr_multiplier_by_week_hyp= sum(coef*value)),by=c('iterations','week.idx')]		
	tmp <- as.data.table(reshape2::melt(pos[[i]][['logit_hosp_fatality_ratio_wildtype']]))
	setnames(tmp, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_wildtype'))
	pa <- merge(pa,tmp,by=c('iterations'),allow.cartesian=TRUE)	
	pa[, baseline_hfr := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype+logit_fr_multiplier_by_week_hyp)]
	set(pa, NULL, c('logit_hosp_fatality_ratio_wildtype','logit_fr_multiplier_by_week_hyp'), NULL)		
	pa[, loc_label := names(pos)[i]]
	pas.ref.hfr[[i]] <- pa
	rm(pa, tmp)
}
pas.ref.hfr <- do.call('rbind',pas.ref.hfr)
# determine location with min baseline hfr 

tmp2 <- dpop[,list(pop=sum(pop)),by='age.label']
tmp2[, tpop:=sum(pop)]
tmp2[, ppop.all.cities:=pop/tpop]
setkey(tmp2,age.label)
tmp2[,age.idx:=1:nrow(tmp2)]
pas.ref.hfr <- merge(pas.ref.hfr,tmp2,by='age.idx')
tmp <- pas.ref.hfr[,list(baseline_hfr=sum(ppop.all.cities*baseline_hfr)),by=c('loc_label','iterations')]
tmp <- tmp[,list(baseline_hfr=median(baseline_hfr)),by=c('loc_label')]
loc_label.min.ref.hfr <- tmp[, loc_label[which.min(baseline_hfr)]]
# calculate ratio in baseline hfr
tmp <- subset(pas.ref.hfr, loc_label==loc_label.min.ref.hfr, select=c(iterations,age.idx,baseline_hfr))
setnames(tmp,'baseline_hfr','best_baseline_hfr')
pas.ref.hfr <- merge(pas.ref.hfr,tmp,by=c('iterations','age.idx'))
pas.ref.hfr[, baseline_hfr_ratio := baseline_hfr/best_baseline_hfr]
pads <- melt(pas.ref.hfr, id.vars=c('loc_label','week.idx','age.idx','iterations'), measure.vars=c('baseline_hfr','baseline_hfr_ratio'))
pads <- pads[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label','variable','week.idx','age.idx')]
pads <- dcast.data.table(pads, variable+loc_label+week.idx+age.idx ~ stat, value.var='value')
tmp <- subset(dh, select=c(loc_label,age.idx,age.label,week.idx,week,week.start))
pads <- merge(pads, tmp, by=c('loc_label','week.idx','age.idx'))


tmp <- subset(pads, week.idx==1 & variable=='baseline_hfr')
p <- ggplot(tmp, aes(x=age.label)) +		
		geom_linerange(aes(ymin=CL,ymax=CU,colour=change_city_label(loc_label)), alpha=0.3, position=position_dodge(.8)) +
		geom_point(aes(y=M,colour=change_city_label(loc_label)), position=position_dodge(.8)) +		
		scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +	
		#scale_fill_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +	
		scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
		scale_x_discrete(expand=c(0,0)) +		
		theme_bw() +
		theme(legend.position='bottom') +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='',y='baseline in-hospital fatality rate',colour='',fill='')
ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline.pdf'),p,w=10,h=8)
ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline.png'),p,w=10,h=8)


tmp <- subset(pads, week.idx==1 & variable=='baseline_hfr_ratio' & loc_label!=loc_label.min.ref.hfr & CU-CL<2)
p <- ggplot(tmp, aes(x=age.label)) +
		geom_hline(yintercept=1, colour='black') +
		geom_linerange(aes(ymin=CL,ymax=CU,colour=change_city_label(loc_label)), position=position_dodge(0.7), alpha=0.3) +
		geom_point(aes(y=M,colour=change_city_label(loc_label)), position=position_dodge(0.7)) +		
		scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +				
		scale_y_continuous(expand=c(0,0)) +				
		theme_bw() +
		theme(legend.position='bottom') +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='',y=paste0('ratio of baseline in-hospital fatality rates\nrelative to ',change_city_label(loc_label.min.ref.hfr)),colour='',fill='') 
ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline_ratio.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline_ratio.png'),p,w=10,h=6)


################################################################################
## relationships between health care demand predictors and age-specific hfr in baseline week
###############################################################################


# get resources in first observed week
tmp2 <- subset(dhca, variable%in%c('p.nonres.hosps.00','r.nhosp.adm.00','r.nicu.adm.00'), select= c(variable,loc_label,week,week.start,ref.week,value))
set(tmp2,NULL,'variable',tmp2[,gsub('.00','',variable)])
tmp2 <- dcast.data.table(tmp2, loc_label+week+week.start+ref.week~variable, value.var='value')
tmp2 <- merge(tmp2,subset(dhcb, select=c(loc_label,week, r.physicians.all, r.physicians.specialist, r.icubeds.all, r.saribeds.all, r.beds.all, r.ventilator,  r.nurse, r.tech.nurse, r.physiotherapist)),by=c('loc_label','week'))
tmp2 <- suppressWarnings(melt(tmp2, id.vars=c('loc_label','week','week.start')))
tmp2 <- merge(tmp2, tmp2[,list(week=min(week)),by='loc_label'],by=c('loc_label','week'))
tmp2 <- subset(tmp2, select=c(loc_label,variable,value))

# age specific fatality rates per baseline demand
pads <- pas.ref.hfr[, list(value=quantile(baseline_hfr, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label','age.idx')]
pads <- dcast.data.table(pads, loc_label+age.idx~stat, value.var='value')
tmp3 <- unique(subset(dh, select=c(age.idx, age.label)))
pads <- merge(pads, tmp3, by='age.idx')
pads <- merge(pads, tmp2, by='loc_label',allow.cartesian=TRUE)

z <- c('r.nhosp.adm','r.nicu.adm','r.physicians.all','r.physicians.specialist','r.icubeds.all','r.ventilator','r.nurse','r.tech.nurse','r.physiotherapist')
names(z) <- c('SARI hospital admissions','ICU admissions','physicians','critical care specialists','ICU beds','Ventilators','Nurses','Nurse assistants','Physiotherapist')
for(i in seq_along(z))
{
	tmp <- subset(pads, variable==z[i])	
	tmp[, r.squared := NA_real_]
	tmp[, p.val := NA_real_]
	zz <- sort(unique(pads$age.label))
	for(j in seq_along(zz))
	{
		tmp3 <- summary(lm(M~value,data=subset(tmp, age.label==zz[j])))	
		set(tmp, tmp[,which(age.label==zz[j])],'r.squared',tmp3$r.squared)
		set(tmp, tmp[,which(age.label==zz[j])],'p.val',tmp3$coefficients[2,4] )	
	}		
	p <- ggplot(tmp, aes(x=value)) +
			geom_smooth(aes(x=value,y=M), method='lm',formula='y~x',colour='black',lwd=0.5) +
			geom_linerange(aes(x=value,ymin=CL,ymax=CU,colour=change_city_label(loc_label))) +
			geom_point(aes(x=value,y=M,colour=change_city_label(loc_label))) +
			geom_text(data=subset(tmp, loc_label=='sao paulo'), aes(x=-Inf, y=Inf, label=paste0('p-val=',round(p.val, d=2))), colour="black", inherit.aes=FALSE, hjust=-0.2, vjust=1.5) +
			scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +			
			scale_y_continuous(expand=c(0,0), labels=scales::percent) +
			facet_wrap(~age.label,scales='free',ncol=4) +
			theme_bw() +
			theme(legend.position='bottom',strip.background = element_blank()) +
			guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
			labs(x=names(z)[i],y='baseline in-hospital fatality rate',colour='') 
	ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline_vs_predictor_',z[i],'.pdf'),p,w=10,h=8)
	ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline_vs_predictor_',z[i],'.png'),p,w=10,h=8)	
}

# all age fatality rates per baseline demand, standardising  fatality rates
tmp <- pas.ref.hfr[, list(baseline_hfr= sum(ppop.all.cities*baseline_hfr)), by=c('loc_label','iterations')]
pads <- tmp[, list(value=quantile(baseline_hfr, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label')]
pads <- dcast.data.table(pads, loc_label~stat, value.var='value')
pads <- merge(pads, tmp2, by='loc_label',allow.cartesian=TRUE)
pads <- merge(pads,data.table(variable=unname(z),p.names=gsub('per 100,000 ','\n',names(z))),by='variable')
pads[, r.squared := NA_real_]
pads[, p.val := NA_real_]
for(i in seq_along(z))
{
	tmp <- subset(pads, variable==z[i]) 
	tmp <- summary(lm(M~value,data=tmp))	
	set(pads, pads[,which(variable==z[i])],'r.squared',tmp$r.squared)
	set(pads, pads[,which(variable==z[i])],'p.val',tmp$coefficients[2,4] )	
}

p <- ggplot(pads, aes(x=value)) +
		geom_smooth(aes(x=value,y=M), method='lm',formula='y~x',colour='black',lwd=0.5) +
		geom_linerange(aes(x=value,ymin=CL,ymax=CU,colour=change_city_label(loc_label))) +
		geom_point(aes(x=value,y=M,colour=change_city_label(loc_label))) +		
		geom_text(data=subset(pads, loc_label=='sao paulo'), aes(x=-Inf, y=Inf, label=paste0('r2=',round(r.squared, d=2),' p-val=',round(p.val, d=3))), size=3, colour="black", inherit.aes=FALSE, hjust=-0.2, vjust=1.5) +		
		scale_colour_manual(values = args$city.palette(length(unique(pad$loc_label)))) +			
		scale_y_continuous(expand=c(0,0), labels=scales::percent) +
		facet_wrap(~p.names,scales='free',ncol=5) +
		theme_bw() +
		theme(legend.position='bottom',strip.background = element_blank()) +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='number per 100,000 population',y=paste0('baseline age-standardised in-hospital fatality rate'),colour='') 
ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline_vs_predictors.pdf'),p,w=10,h=8)
ggsave(file=paste0(args$out.base,'_allloc_hfrbaseline_vs_predictors.png'),p,w=10,h=8)



####################################################################
# health care demand attributable fatality rate multipliers across locations relative to min week
####################################################################


tmp2 <- dpop[,list(pop=sum(pop)),by='age.label']
tmp2[, tpop:=sum(pop)]
tmp2[, ppop.all.cities:=pop/tpop]
setkey(tmp2,age.label)
tmp2[,age.idx:=1:nrow(tmp2)]
pads <- list()
for(i in 1:length(pos))
{	
	cat('\nprocessing',i,' ',names(pos)[i])
	pa <- as.data.table(reshape2::melt(pos[[i]]$hfr_wildtype))
	setnames(pa, 1:4, c('iterations','week.idx','age.idx','hfr.wildtype'))	
	tmp <- unique(subset(pa, is.na(hfr.wildtype), select=iterations))[,iterations]
	if(length(tmp))
	{
		cat('\nFound iterations with NA, removing:', tmp)
		pa <- subset(pa, !iterations%in%tmp)
	}	
	tmp <- unique(subset(dh,loc_label==names(pos)[i],select=c(week.idx,week)))
	pa <- merge(pa,tmp,by='week.idx')	
	tmp <- subset(dhr,loc_label==names(pos)[i])[, min.week]
	tmp <- subset(pa, week==tmp, select=-c(week,week.idx))
	setnames(tmp,'hfr.wildtype','hfr.wildtype.min')
	pa <- merge(pa,tmp,by=c('age.idx','iterations'))
	pa[, hfr.multiplier := hfr.wildtype/hfr.wildtype.min]
	pa <- merge(pa,subset(tmp2,select=c(age.idx,ppop.all.cities)),by='age.idx')
	pa <- pa[,list(hfr.multiplier=sum(hfr.multiplier*ppop.all.cities)),by=c('iterations','week.idx')]	
	pa <- pa[, list(value=quantile(hfr.multiplier, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('week.idx')]
	pa <- dcast.data.table(pa, week.idx ~ stat, value.var='value')
	pa[, loc_label := names(pos)[i]]
	pads[[i]] <- pa
	rm(pa, tmp)
}
pads <- do.call('rbind',pads)
tmp <- unique(subset(dh, select=c(loc_label,week.idx,week,week.start)))
pads <- merge(pads, tmp, by=c('loc_label','week.idx'))

p <- ggplot(pads, aes(x=week.start)) +
		geom_hline(yintercept=1,colour='black') +
		geom_ribbon(aes(ymin=CL, ymax=CU,fill=change_city_label(loc_label)),alpha=0.5) +
		geom_line(aes(y=M,colour=change_city_label(loc_label))) +
		scale_colour_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_fill_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_x_date(breaks="2 months") +
		scale_y_continuous(expand=c(0,0)) +		
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='', y='multiplier to minimum in-hospital fatality rate\nbased on changing health care demand per available resources',colour='',fill='')
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_multiplier.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_multiplier.png'),p,w=10,h=6)

p <- ggplot(pads, aes(x=week.start)) +
		geom_hline(yintercept=1,colour='black') +
		geom_ribbon(aes(ymin=CL, ymax=CU,fill=change_city_label(loc_label)),alpha=0.5) +
		geom_line(aes(y=M,colour=change_city_label(loc_label))) +
		scale_colour_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_fill_manual(values = args$city.palette(length(unique(pads$loc_label)))) +
		scale_x_date(breaks="2 months") +
		scale_y_continuous(expand=c(0,0)) +		
		facet_wrap(~change_city_label(loc_label),ncol=5) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(fill='none',colour='none') +
		labs(x='', y='multiplier to minimum in-hospital fatality rate\nbased on changing health care demand per available resources',colour='',fill='')
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_multiplier2.pdf'),p,w=12,h=10)
ggsave(file=paste0(args$out.base,'_allloc_hfrmin_multiplier2.png'),p,w=12,h=10)


#########################################################################
# get smallest bulk ESS and plot corresponding trace
#########################################################################
cat('\nRead MCMC summaries to extract smallest bulk ESS \n')

tmp1 <- list.files(args$out.dir , pattern='_diagnostics.csv')

sexpr <- list(min_bulk_ESS = 10^10, max_Rhat = 0, min_bulk_ESS_loc = 'nowhere', max_Rhat_loc = 'nowhere')

for(file in tmp1)
{	
  tmp <- file.path(args$out.dir,file)
  cat('\nreading ',tmp)
  
  tmp <- as.data.table(read.csv(tmp), row.names = FALSE)
  tmp1 <- tmp[variable == 'min_ess_bulk', 'X1.']
  if (tmp1 < sexpr$min_bulk_ESS)
  {
    sexpr$min_bulk_ESS <- tmp1
    sexpr$min_bulk_ESS_loc <- stringr::str_match(file, "_\\s*(.*?)\\s*_diagnostics.csv")[,2]
  }

  tmp1 <- tmp[variable == 'max_rhat', 'X1.']
  if (tmp1 > sexpr$max_Rhat){
    sexpr$max_Rhat <- tmp1
    sexpr$max_Rhat_loc <- stringr::str_match(file, "_\\s*(.*?)\\s*_diagnostics.csv")[,2]}
}

cat(paste0('\nThe smallest bulk ESS was ', sexpr$min_bulk_ESS, ' and was obtained for ', sexpr$min_bulk_ESS_loc))
cat(paste0('\nThe maximum Rhat was ', sexpr$max_Rhat, ' and was obtained for ', sexpr$max_Rhat_loc))



# Now plot
tmp <- do[loc_label2 == sexpr$min_bulk_ESS_loc, FF]
m_fit <- readRDS(file.path(args$out.dir,tmp))
pd <- m_fit$draws(inc_warmup = FALSE)		
tmp <- as.data.table(summarise_draws(
  m_fit$draws(
    variables=c('log_age_hosp_rate_wildtype','log_age_hosp_rate_P1',                
                'logit_hosp_fatality_ratio_wildtype','logit_hosp_fatality_ratio_P1_rnde','logit_hosp_fatality_ratio_P1_rnde_sd',
                'age_hosps_v_inflation',
                'propp1_logistic_growthrate','propp1_logistic_midpoint','propp1_overdisp_inv',
                'fr_regcoeff_overall','fr_shrinkage','fr_scale','fr_overdisp_inv',
                'lp__'),
    inc_warmup = FALSE)
))
tmp <- tmp[ess_bulk == min(ess_bulk), variable]
pd <- m_fit$draws(inc_warmup = TRUE)

sexpr[['min_bulk_ESS']] <- as.character(as.integer(sexpr[['min_bulk_ESS']]) )
sexpr[['min_bulk_ESS_loc']] <- do[loc_label2 == sexpr[['min_bulk_ESS_loc']], loc_label]
sexpr[['min_bulk_ESS_loc']] <- change_city_label(as.character(sexpr[['min_bulk_ESS_loc']]))
sexpr[['max_Rhat']] <- as.character(round(sexpr[['max_Rhat']],2))
sexpr[['max_Rhat_loc']] <- do[loc_label2 == sexpr[['max_Rhat_loc']], loc_label]
sexpr[['max_Rhat_loc']] <- change_city_label(as.character(sexpr[['max_Rhat_loc']] ))
saveRDS(sexpr,paste0(args$out.base,'_substitute_expressions_HMC.rds'))

# can 
bayesplot:::color_scheme_set("mix-blue-pink")
p <- bayesplot:::mcmc_trace(pd,  
                            pars = tmp, 
                            n_warmup = 5e2)
ggsave(file=paste0(args$out.base,'_allloc_lowestESSbulk_traces.pdf'),p,w=8,h=4)



