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
require(scales)

if(1) # Andrea
{
  pkg.dir <- '~/git/covid19_brazil_hfr/'	
  out.dir <- '~/Documents/P1Brazil/hfr-fit-210719h4-resources8-ex0-trvaluecs-hd1'
  file.prefix <- basename(out.dir)
  out.base2 <- file.path(out.dir, file.prefix)
  plot.dir <- out.dir
  loc_label <- 'florianopolis'
}


# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-outdir')	
  stopifnot(args_line[[3]]=='-loc_label')
  stopifnot(args_line[[5]]=='-file_prefix')
  stopifnot(args_line[[7]]=='-pkgdir')
  plot.dir <- args_line[[2]]
  out.dir <- args_line[[2]]
  loc_label <- gsub('-|_',' ',args_line[[4]])
  file.prefix <- args_line[[6]]
  pkg.dir <- args_line[[8]]
} 
cat(loc_label); cat('\n')
####################################################################
# source functions
####################################################################

tmp <- list.files(file.path(pkg.dir,'inst/R'), full.names=TRUE)
cat('sourcing code from files\n', tmp)
invisible(lapply(tmp, source))

####################################################################
# source workspace
####################################################################

tmp <- paste0(file.path(out.dir,file.prefix), '_', gsub('_|-| ', '', loc_label),'_workspace.rda')
cat('\nAttempting to load ', tmp)
load(file = tmp)
if(dir.exists(dirname(out.base2)))
{
  args$out.base <- out.base2
}

########################################################################################
### Load fit and postprocess everything
########################################################################################

tmp <- paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_fit.rds')
cat('\nAttempting to load ', tmp)
m_fit <- readRDS( tmp)

# save HMC samples in array format
pd <- m_fit$draws(inc_warmup = FALSE)		

# save summary of parameters
su <- as.data.table(summarise_draws(
  m_fit$draws(
    variables=c('log_age_hosp_rate_wildtype','log_age_hosp_rate_P1',                
                'logit_hosp_fatality_ratio_wildtype','logit_hosp_fatality_ratio_P1_rnde','logit_hosp_fatality_ratio_P1_rnde_sd',
                'age_hosps_v_inflation',
                'propp1_logistic_growthrate','propp1_logistic_midpoint','propp1_overdisp_inv',
                'fr_regcoeff_overall','fr_shrinkage','fr_scale','fr_overdisp_inv',
                'lp__'),
    inc_warmup = FALSE)
))
write.csv(su, file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_convergence.csv'), row.names=TRUE)
su[,min(ess_bulk)]
su[,max(rhat)]

# save summary of diagnostics
diagnostics <- summarise_draws(posterior::as_draws_df(m_fit$sampler_diagnostics()), ~quantile(.x, probs = c(0.01, 0.5, 0.99)))
diagnostics <- rbind(diagnostics, c('min_ess_bulk', su[,min(ess_bulk)],0,0))
diagnostics <- rbind(diagnostics, c('max_rhat', su[,max(rhat)],0,0))
write.csv(diagnostics, file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_diagnostics.csv'), row.names=TRUE)

# Looks like this could be killing some jobs on the HPC? 
# null device
# 1
# /var/spool/PBS/mom_priv/jobs/4340131[4].pbs.SC: line 11: 2207921 Killed
if(0) {
  # traces
  tmp <- m_fit$draws(inc_warmup = TRUE)
  bayesplot:::color_scheme_set("mix-blue-pink")
  p <- bayesplot:::mcmc_trace(tmp,  
                              regex_pars = "^log_age_hosp_rate_wildtype|^log_age_hosp_rate_P1|^age_hosps_v_inflation$|^lp__$", 
                              n_warmup = 5e2,
                              facet_args = list(ncol = 1, labeller = label_parsed)
  )
  pdf(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_traces_agecomp.pdf'), w=10, h=70)
  print(p)
  dev.off()
  
  p <- bayesplot:::mcmc_trace(tmp,  
                              regex_pars = "^logit_hosp_fatality_ratio_wildtype|^logit_hosp_fatality_ratio_P1_rnde|^lp__$", 
                              n_warmup = 5e2,
                              facet_args = list(ncol = 1, labeller = label_parsed)
  )
  pdf(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_traces_hfr.pdf'), w=10, h=40)
  print(p)
  dev.off()	
  
  tmp <- m_fit$draws(variables=c("fr_regcoeff_overall",'fr_shrinkage','fr_scale',"fr_overdisp_inv","lp__"),inc_warmup = TRUE)
  dimnames(tmp)[[3]] <- gsub('\\[|,|\\]','_',dimnames(tmp)[[3]])
  p <- bayesplot:::mcmc_trace(tmp,  		 
                              n_warmup = 5e2,
                              facet_args = list(ncol = 1, labeller = label_parsed)
  )
  pdf(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_traces_frmultiplier.pdf'), w=10, h=50)
  print(p)
  dev.off()	
  
  
  
  p <- bayesplot:::mcmc_trace(tmp,  
                              regex_pars = "propp1_logistic_growthrate|propp1_logistic_midpoint|propp1_overdisp_inv|lp__$", 
                              n_warmup = 5e2,
                              facet_args = list(ncol = 1, labeller = label_parsed)
  )
  pdf(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_traces_seropossurvival_other.pdf'), w=10, h=25)
  print(p)
  dev.off()	
}
####################################################################
# get MCMC samples for particular chains into rstan::extract format 
####################################################################

select.chains <- seq_along(dimnames(pd)[['chain']])
iters <- 1:(m_fit$metadata()[['iter_sampling']])
#iters <- 1000:(m_fit$metadata()[['iter_sampling']])
#pd <- pd[iters,,]

po <- list()
tmp <- pd[,,which(grepl('log_age_hosp_rate_wildtype',dimnames(pd)[[3]]))]
po$log_age_hosp_rate_wildtype <- unname(apply(tmp[,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('log_age_hosp_rate_P1',dimnames(pd)[[3]]))]
po$log_age_hosp_rate_P1 <- unname(apply(tmp[,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('age_prop_hosps_wildtype',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$age_prop_hosps_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('age_prop_hosps_P1',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$age_prop_hosps_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('age_prop_deaths_wildtype',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$age_prop_deaths_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('age_prop_deaths_P1',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$age_prop_deaths_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
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
po$exp_deaths_wildtype_in_hosp <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('exp_deaths_P1_in_hosp',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$exp_deaths_P1_in_hosp <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('exp_hosp_adm_wildtype',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$exp_hosp_adm_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('exp_hosp_adm_P1',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$exp_hosp_adm_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('hfr_overall',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$hfr_overall <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('hfr_wildtype',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$hfr_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('hfr_P1',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$hfr_P1 <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))
tmp <- pd[,,which(grepl('fr_regcoeff_overall',dimnames(pd)[[3]]))]
po$fr_regcoeff_overall <- unname(apply(tmp[,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('fr_overdisp_inv',dimnames(pd)[[3]]))]
po$fr_overdisp_inv <- unname(apply(tmp[,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('^fr_multiplier_by_week\\[',dimnames(pd)[[3]]))]
tmp <- apply(tmp[,select.chains,], 3, rbind)
po$fr_multiplier_by_week <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/stan_data$A, stan_data$A))

pd <- NULL
gc()

if(1){
####################################################################
# fitted prop P1 vs empirical
####################################################################

pad <- as.data.table(reshape2::melt(po$prop_P1))
setnames(pad, 1:3, c('iterations','week.idx','value'))
pads <- pad[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx')]
pads <- dcast.data.table(pads, week.idx~stat, value.var='value')
tmp <- unique(subset(dh, select=c(week.idx, week, week.start)))
pads <- merge(pads, tmp, by='week.idx')
tmp <- subset(dp, select=c(week, n_positive, n_sampled))
pads <- merge(pads, tmp, by='week', all.x=TRUE)

p <- ggplot(pads, aes(x=week.start)) +
		geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey50') +
		geom_line(aes(y=M)) +
		geom_point(aes(y=n_positive/n_sampled, size=n_sampled), colour='orangered3') +		
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		scale_x_date(breaks='months',expand=c(0,0)) +		
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		labs(y='proportion P1\n(posterior median)', x='')
cat('\nplotting _seq_propP1_over_time.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_seq_propP1_over_time.pdf'),p,w=6,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_seq_propP1_over_time.png'),p,w=6,h=6)


####################################################################
# extract and compare age distribution of deaths
####################################################################

pa <- as.data.table(reshape2::melt(po$age_prop_deaths_wildtype))
setnames(pa, 1:4, c('iterations','week.idx','age.idx','prop.age.in.deaths.wildtype'))
tmp <- as.data.table(reshape2::melt(po$age_prop_deaths_P1))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','prop.age.in.deaths.P1'))
pa <- merge(pa,tmp,by=c('iterations','week.idx','age.idx'))
pa[, ratio.P1.wildtype := prop.age.in.deaths.P1/prop.age.in.deaths.wildtype]
pa <- melt(pa, id.vars=c('iterations','week.idx','age.idx'))
setkey(pa, variable, iterations, week.idx, age.idx)
tmp <- pa[variable!='ratio.P1.wildtype', list(age.idx=age.idx, value=cumsum(value)), by=c('variable','iterations','week.idx')]
set(tmp, NULL, 'variable', tmp[,paste0('c',variable)])
pa <- rbind(pa,tmp)

pas <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','week.idx','age.idx')]
tmp <- unique(subset(dh, select=c(age.idx, week.idx, week, week.start, age.label)))
pas <- merge(pas, tmp, by=c('week.idx','age.idx'))
pas <- dcast.data.table(pas, variable+age.idx+age.label+week+week.start~stat, value.var='value')
pas[, variant := factor(grepl('P1', variable), levels=c(FALSE,TRUE), labels=c('non-P1','P1'))]
pas[, plot.facet := gsub('^([a-z]+)\\..*','\\1',variable)]

p <- ggplot(subset(pas, week==min(week) & plot.facet=='prop'), aes(x=age.label, fill=variant)) +
  geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(labels=scales::percent) +
  ggsci::scale_fill_npg() +	
  labs(x='',y='share of age group among COVID-19 deaths in hospitals',fill='variant')
cat('\nplotting _deaths_prop_among_deaths_P1_vs_wildtype.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_prop_among_deaths_P1_vs_wildtype.pdf'),p,w=8,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_prop_among_deaths_P1_vs_wildtype.png'),p,w=8,h=6)

p <- ggplot(subset(pas, week==min(week) & plot.facet=='ratio'), aes(x=age.label)) +
  geom_hline(yintercept=1, colour='grey50') +
  geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +					
  labs(x='',y='prevalance ratio of age group among COVID-19 deaths\nP1 vs wildtype',fill='variant')
cat('\nplotting _deaths_prevratio_among_deaths_P1_vs_wildtype.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_prevratio_among_deaths_P1_vs_wildtype.pdf'),p,w=8,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_prevratio_among_deaths_P1_vs_wildtype.png'),p,w=8,h=6)

tmp <- subset(pas, week==min(week) & plot.facet=='cprop')
set(tmp, NULL, 'age.label', tmp[,gsub('\\[[0-9]+','\\[0',age.label)])
p <- ggplot(subset(tmp, age.label!='[0,150)'), aes(x=age.label, fill=variant)) +
  geom_bar(aes(y=M), stat='identity', position=position_dodge(0.8)) +
  geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.8), width=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  ggsci::scale_fill_npg() +	
  labs(x='',y='cumulated share of age group among COVID-19 deaths',fill='variant')
cat('\nplotting _deaths_cprop_among_deaths_P1_vs_wildtype.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_cprop_among_deaths_P1_vs_wildtype.pdf'),p,w=8,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_cprop_among_deaths_P1_vs_wildtype.png'),p,w=8,h=6)


####################################################################
# fitted share of age groups among deaths vs empirical
####################################################################

pa <- as.data.table(reshape2::melt(po$exp_deaths_wildtype_in_hosp))
setnames(pa, 1:4, c('iterations','week.idx','age.idx','exp_deaths_wildtype_in_hosp'))
tmp <- as.data.table(reshape2::melt(po$exp_deaths_P1_in_hosp))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp_deaths_P1_in_hosp'))
pa <- merge(pa, tmp, by=c('iterations','week.idx','age.idx'))
pa[,exp_deaths_in_hosp := exp_deaths_wildtype_in_hosp + exp_deaths_P1_in_hosp]
set(pa,NULL,c('exp_deaths_wildtype_in_hosp','exp_deaths_P1_in_hosp'),NULL)
tmp <- pa[,list(texp_deaths_in_hosp=sum(exp_deaths_in_hosp)),by=c('iterations','week.idx')]
pa <- merge(pa,tmp,by=c('iterations','week.idx'))
pa[, pexp_deaths_in_hosp:=exp_deaths_in_hosp/texp_deaths_in_hosp]
# need to NA remove something! (try easy fix atm)
pads <- pa[, list(value=quantile(pexp_deaths_in_hosp, p=c(0.025,0.25,.5,0.75,.975), na.rm=T),stat=c('CL','IL','M','IU','CU')),
           by=c('week.idx','age.idx')]
pads <- dcast.data.table(pads, week.idx+age.idx~stat, value.var='value')
tmp <- subset(dh, select=c(loc_label, age.label, week.start, week.idx, age.idx, deaths_of_hosps))
tmp2 <- tmp[,list(tdeaths_of_hosps=sum(deaths_of_hosps)),by=c('week.idx')]
tmp <- merge(tmp,tmp2,by='week.idx')
tmp[,pdeaths_of_hosps:=deaths_of_hosps/tdeaths_of_hosps]
pads <- merge(tmp, pads, by=c('week.idx','age.idx'))

# plot empirical proportion deaths (facets)
p <- ggplot(pads, aes(x=week.start)) +
  geom_ribbon(aes(ymin=CL, ymax=IL), fill='grey80') +
  geom_ribbon(aes(ymin=IU, ymax=CU), fill='grey80') +
  geom_ribbon(aes(ymin=IL, ymax=IU), fill='grey50') +
  geom_line(aes(y=M), colour='black') +		
  geom_point(aes(y=pdeaths_of_hosps, colour=age.label)) +		
  viridis::scale_colour_viridis(option='B', discrete = TRUE, end=0.8 ) +
  scale_x_date(breaks='2 months',expand=c(0,0)) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  labs(colour='age', y='share of age groups among COVID-19 deaths', x='') +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +
  facet_wrap(~age.label, ncol=6)
cat('\nplotting _deaths_pdeaths_vs_fitted.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_pdeaths_vs_fitted.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_pdeaths_vs_fitted.png'),p,w=10,h=6)


####################################################################
# expected deaths by age P1 and wildtype
####################################################################

pad <- as.data.table(reshape2::melt(po$exp_deaths_wildtype_in_hosp))
setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.deaths.wildtype'))
tmp <- as.data.table(reshape2::melt(po$exp_deaths_P1_in_hosp))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.deaths.P1'))
pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
pad <- melt(pad, id.vars=c('iterations','week.idx','age.idx'), measure.vars=c('exp.deaths.wildtype','exp.deaths.P1'))
setkey(pad, variable,iterations,age.idx,week.idx)
tmp <- pad[, list(week.idx=week.idx, value=cumsum(value)), by=c('variable','iterations','age.idx')]
set(tmp, NULL, 'variable', tmp[,gsub('deaths','cdeaths',variable)]) # might want to act on this tmp
cd <- tmp[, list(value = sum(value)), by = c("iterations", "age.idx", "week.idx")] #
cd <- cd[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx','age.idx')] #
pad <- rbind(pad,tmp)
pads <- pad[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('variable','week.idx','age.idx')]
pads[, variant := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\3',variable)]
pads[, variable := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\1.\\2',variable)]
tmp <- unique(subset(dh, select=c(week.idx,age.idx,week,week.start,age.label)))
pads <- merge(tmp, pads, by=c('week.idx','age.idx'))
pads <- dcast.data.table(pads, variable+variant+week.idx+week+week.start+age.idx+age.label~stat, value.var='value')
cd <- merge(tmp, cd, by=c('week.idx','age.idx')) #
cd <- dcast.data.table(cd, week.idx+week+week.start+age.idx+age.label~stat, value.var='value')  #
tmp <- pads[, list(age.idx=age.idx, sM=cumsum(M)), by=c('variable','variant','week.idx')]
pads <- merge(pads, tmp, by=c('variable','variant','week.idx','age.idx'))
set(pads, NULL, 'variant', pads[, factor(variant, levels=c('wildtype', 'P1'), labels=c('non-P1','P1'))])

pae <- subset(pads, variable=='exp.cdeaths', select=c(week.idx,age.idx,variant,M))
setkey(pae, week.idx,age.idx,variant)
tmp <- pae[, list(variant=variant, sM=cumsum(M)), by=c('week.idx','age.idx')]
pae <- merge(pae, tmp, by=c('week.idx','age.idx','variant'))
tmp <- subset(dh, select=c(loc_label, age.label, week.start, week.idx, age.idx, deaths_of_hospsi))
setkey(tmp, loc_label, age.idx, age.label, week.idx, week.start)
tmp <- tmp[, list(week.idx=week.idx, week.start=week.start, cdeaths_of_hospsi=cumsum(deaths_of_hospsi)),by=c('loc_label','age.idx','age.label')]
pae <- merge(tmp, pae, by=c('week.idx','age.idx'))

p <- ggplot(pae, aes(x=week.start)) +
  geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=variant)) +
  geom_linerange(data = subset(cd, week.idx%%4==1), aes(ymin = CL, ymax = CU), color = 'grey50') +
  geom_point(data = subset(cd, week.idx%%4==1), aes(y = M), color = 'grey50') +
  geom_line(data=tmp, aes(y=cdeaths_of_hospsi)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_date(breaks='3 months',expand=c(0,0)) +
  scale_fill_brewer() +
  facet_wrap(~age.label, ncol=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'bottom') +
  labs(fill='variant', y='expected COVID-19 in-hospital deaths by variant\n(posterior median)', x='date of hospital admission')
cat('\nplotting _deaths_cdeaths_by_variant_with_empirical.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_cdeaths_by_variant_with_empirical.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_cdeaths_by_variant_with_empirical.png'),p,w=10,h=6)


####################################################################
# scatter plot overall deaths
####################################################################


pad <- as.data.table(reshape2::melt(po$exp_deaths_wildtype_in_hosp))
setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.deaths.wildtype'))
tmp <- as.data.table(reshape2::melt(po$exp_deaths_P1_in_hosp))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.deaths.P1'))
pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
pad[, exp.deaths := exp.deaths.wildtype + exp.deaths.P1]
pads <- pad[, list(value=quantile(exp.deaths, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx','age.idx')]
tmp <- subset(dh, select=c(loc_label, age.label, week.start, week.idx, age.idx, deaths_of_hospsi))
pads <- merge(tmp, pads, by=c('week.idx','age.idx'))
pads <- dcast.data.table(pads, week.idx+week.start+age.idx+age.label+deaths_of_hospsi~stat, value.var='value')

p <- ggplot(pads, aes(x=deaths_of_hospsi)) +	
  geom_abline(intercept=0, slope=1) +
  geom_errorbar(aes(ymin=CL, ymax=CU), colour='grey50') +
  geom_point(aes(y=M, colour=age.label)) +
  viridis::scale_colour_viridis(option='B', discrete = TRUE) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw() +
  labs(colour='age', y='expected COVID-19 in-hospital deaths', x='actual COVID-19 in-hospital deaths')
cat('\nplotting _deaths_ndeaths_fit_scatter.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_ndeaths_fit_scatter.pdf'),p,w=7,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_deaths_ndeaths_fit_scatter.png'),p,w=7,h=6)


####################################################################
# expected hospitalisations by age P1 and wildtype
####################################################################


pad <- as.data.table(reshape2::melt(po$exp_hosp_adm_wildtype))
setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.hosp.adm.wildtype'))
tmp <- as.data.table(reshape2::melt(po$exp_hosp_adm_P1))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.hosp.adm.P1'))
pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
pad <- melt(pad, id.vars=c('iterations','week.idx','age.idx'), measure.vars=c('exp.hosp.adm.wildtype','exp.hosp.adm.P1'))
setkey(pad, variable,iterations,age.idx,week.idx)
tmp <- pad[, list(week.idx=week.idx, value=cumsum(value)), by=c('variable','iterations','age.idx')]
set(tmp, NULL, 'variable', tmp[,gsub('hosp','chosp',variable)]) #
cd <- tmp[, list(value = sum(value)), by = c("iterations", "age.idx", "week.idx")] #
cd <- cd[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx','age.idx')]
pad <- rbind(pad,tmp)
pads <- pad[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('variable','week.idx','age.idx')]
pads[, variant := gsub('^([a-z]+)\\.([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\4',variable)]
pads[, variable := gsub('^([a-z]+)\\.([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\1.\\2.\\3',variable)]
tmp <- unique(subset(dh, select=c(week.idx,age.idx,week,week.start,age.label)))
pads <- merge(tmp, pads, by=c('week.idx','age.idx'))
pads <- dcast.data.table(pads, variable+variant+week.idx+week+week.start+age.idx+age.label~stat, value.var='value')
cd <- merge(tmp, cd, by=c('week.idx','age.idx')) #
cd <- dcast.data.table(cd, week.idx+week+week.start+age.idx+age.label~stat, value.var='value')  #
tmp <- pads[, list(age.idx=age.idx, sM=cumsum(M)), by=c('variable','variant','week.idx')]
pads <- merge(pads, tmp, by=c('variable','variant','week.idx','age.idx')) # There doesn t seem to be a week.idx...
set(pads, NULL, 'variant', pads[, factor(variant, levels=c('wildtype','P1'), labels=c('non-P1','P1'))])

pae <- subset(pads, variable=='exp.chosp.adm', select=c(week.idx,age.idx,variant,M))
setkey(pae, week.idx,age.idx,variant)
tmp <- pae[, list(variant=variant, sM=cumsum(M)), by=c('week.idx','age.idx')]
pae <- merge(pae, tmp, by=c('week.idx','age.idx','variant'))
tmp <- unique(subset(dh, select=c(week.idx,age.idx,week,week.start,age.label, chosps)))
pae <- merge(tmp, pae, by=c('week.idx','age.idx'))

p <- ggplot(pae, aes(x=week.start)) +
  geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=variant)) +
  geom_linerange(data = subset(cd, week.idx%%4==1), aes(ymin = CL, ymax = CU), color = 'grey50') +
  geom_point(data = subset(cd, week.idx%%4==1), aes(y = M), color = 'grey50') +
  geom_line(data=dh, aes(y=chosps)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_date(breaks='6 months',expand=c(0,0)) +
  scale_fill_brewer() +
  facet_wrap(~age.label, ncol=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = "bottom") +
  labs(fill='variant', y='expected COVID-19 hospital admissions by variant\n(posterior median)', x='')
cat('\nplotting _hospc_by_variant_with_empirical.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hospc_by_variant_with_empirical.pdf'),p,w=10,h= 6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hospc_by_variant_with_empirical.png'),p,w=10,h= 6)


####################################################################
# scatter plot overall hospitalisations
####################################################################


pad <- as.data.table(reshape2::melt(po$exp_hosp_adm_wildtype))
setnames(pad, 1:4, c('iterations','week.idx','age.idx','exp.hosp.adm.wildtype'))
tmp <- as.data.table(reshape2::melt(po$exp_hosp_adm_P1))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','exp.hosp.adm.P1'))
pad <- merge(pad, tmp, by=c('iterations','week.idx','age.idx'))
pad[, exp.hosp.adm := exp.hosp.adm.wildtype + exp.hosp.adm.P1]

pads <- pad[, list(value=quantile(exp.hosp.adm, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx','age.idx')]
tmp <- unique(subset(dh, select=c(week.idx,age.idx,week,week.start,age.label,hosps)))
pads <- merge(tmp, pads, by=c('week.idx','age.idx'))
pads <- dcast.data.table(pads, week.idx+week+week.start+age.idx+age.label+hosps~stat, value.var='value')

p <- ggplot(pads, aes(x=hosps)) +	
  geom_abline(intercept=0, slope=1) +
  geom_errorbar(aes(ymin=CL, ymax=CU), colour='grey50') +
  geom_point(aes(y=M, colour=age.label, shape=age.label)) +
  viridis::scale_colour_viridis(option='D', discrete = TRUE) +
  scale_shape_manual(values=seq_along(levels(pads$age.label))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw() +
  labs(colour='age', shape= 'age', y='expected COVID-19 hospital admissions', x='actual COVID-19 hospital admissions')
cat('\nplotting _hospn_fit_scatter.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hospn_fit_scatter.pdf'),p,w=7,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hospn_fit_scatter.png'),p,w=7,h=6)

}

####################################################################
# fitted share of age groups among hospital admissions vs empirical
####################################################################
if(1)
{
pa <- as.data.table(reshape2::melt(po$age_prop_hosps_wildtype))
setnames(pa, 1:4, c('iterations','week.idx','age.idx','age.prop.hosps.wildtype'))
tmp <- as.data.table(reshape2::melt(po$age_prop_hosps_P1))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','age.prop.hosps.P1'))
pa <- merge(pa, tmp, by=c('iterations','week.idx','age.idx'))
tmp <- as.data.table(reshape2::melt(po$prop_P1))
setnames(tmp, 1:3, c('iterations','week.idx','pP1'))
pa <- merge(pa, tmp, by=c('iterations','week.idx'))
pa[, prop.age.hosps := (1-pP1)*age.prop.hosps.wildtype + pP1*age.prop.hosps.P1]
pads <- pa[, list(value=quantile(prop.age.hosps, p=c(0.025,0.25,.5,0.75,.975)),stat=c('CL','IL','M','IU','CU')), by=c('week.idx','age.idx')]
pads <- dcast.data.table(pads, week.idx+age.idx~stat, value.var='value')
tmp <- subset(dh, select=c(loc_label, age.label, week.start, week.idx, age.idx, phosps))
pads <- merge(tmp, pads, by=c('week.idx','age.idx'))

# plot empirical proportion hospital admissions (facets)
p <- ggplot(pads, aes(x=week.start)) +
    geom_point(aes(y=phosps, colour=age.label)) +	
		geom_ribbon(aes(ymin=CL, ymax=IL), fill='grey80') +
		geom_ribbon(aes(ymin=IU, ymax=CU), fill='grey80') +
		geom_ribbon(aes(ymin=IL, ymax=IU), fill='grey50') +
		geom_line(aes(y=M), colour='black') +		
		viridis::scale_colour_viridis(option='B', discrete = TRUE, end=0.8 ) +
		scale_x_date(breaks='2 months',expand=c(0,0), labels=date_format("%b %y") ) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		labs(colour='age', y='share of age groups\namong COVID-19 attributable hospital admissions', x='') +
		guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +
		facet_wrap(~age.label, ncol=6)
cat('\nplotting _hosps_phosps_vs_fitted.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_phosps_vs_fitted.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_phosps_vs_fitted.png'),p,w=10,h=6)
}else{
####################################################################
# extract and compare age distribution in hospital admissions
####################################################################

pa <- as.data.table(reshape2::melt(po$age_prop_hosps_wildtype))
setnames(pa, 1:4, c('iterations','week.idx','age.idx','prop.age.in.hosps.wildtype'))
tmp <- as.data.table(reshape2::melt(po$age_prop_hosps_P1))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','prop.age.in.hosps.P1'))
pa <- merge(pa,tmp,by=c('iterations','week.idx','age.idx'))
pa[, ratio.P1.wildtype := prop.age.in.hosps.P1/prop.age.in.hosps.wildtype]
pa <- melt(pa, id.vars=c('iterations','week.idx','age.idx'))
setkey(pa, variable, iterations, week.idx, age.idx)
tmp <- pa[variable!='ratio.P1.wildtype', list(age.idx=age.idx, value=cumsum(value)), by=c('variable','iterations','week.idx')]
set(tmp, NULL, 'variable', tmp[,paste0('c',variable)])
pa <- rbind(pa,tmp)
pas <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','week.idx','age.idx')]
tmp <- unique(subset(dh, select=c(age.idx, week.idx, week, week.start, age.label)))
pas <- merge(pas, tmp, by=c('week.idx','age.idx'))
pas <- dcast.data.table(pas, variable+age.idx+age.label+week+week.start~stat, value.var='value')
pas[, variant := factor(grepl('P1', variable), levels=c(FALSE,TRUE), labels=c('non-P1','P1'))]
pas[, plot.facet := gsub('^([a-z]+)\\..*','\\1',variable)]

p <- ggplot(subset(pas, week==min(week) & plot.facet=='prop'), aes(x=age.label, fill=variant)) +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		scale_y_continuous(labels=scales::percent) +
		ggsci::scale_fill_npg() +	
		labs(x='',y='share of age group among COVID-19 attributable hospital admissions',fill='variant')
cat('\nplotting _hosps_prop_among_hosps_P1_vs_wildtype.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_prop_among_hosps_P1_vs_wildtype.pdf'),p,w=8,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_prop_among_hosps_P1_vs_wildtype.png'),p,w=8,h=6)

p <- ggplot(subset(pas, week==min(week) & plot.facet=='ratio'), aes(x=age.label)) +
		geom_hline(yintercept=1, colour='grey50') +
		geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +					
		labs(x='',y='ratio of share of age group among COVID-19 attributable hospital admissions\nGamma vs non-Gamma',fill='variant')
cat('\nplotting _hosps_prevratio_among_hosps_P1_vs_wildtype.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_prevratio_among_hosps_P1_vs_wildtype.pdf'),p,w=8,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_prevratio_among_hosps_P1_vs_wildtype.png'),p,w=8,h=6)

tmp <- subset(pas, week==min(week) & plot.facet=='cprop')
set(tmp, NULL, 'age.label', tmp[,gsub('\\[[0-9]+','\\[0',age.label)])
p <- ggplot(subset(tmp, age.label!='[0,150)'), aes(x=age.label, fill=variant)) +
		geom_bar(aes(y=M), stat='identity', position=position_dodge(0.8)) +
		geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.8), width=0.3) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		ggsci::scale_fill_npg() +	
		labs(x='',y='cumulated share of age group among COVID-19 attributable hospital admissions',fill='variant')
cat('\nplotting _hosps_cprop_among_hosps_P1_vs_wildtype.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_cprop_among_hosps_P1_vs_wildtype.pdf'),p,w=8,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hosps_cprop_among_hosps_P1_vs_wildtype.png'),p,w=8,h=6)



####################################################################
# fitted HFR by variant
####################################################################

pa <- as.data.table(reshape2::melt(po$logit_hosp_fatality_ratio_wildtype))
setnames(pa, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_wildtype'))
tmp <- as.data.table(reshape2::melt(po$logit_hosp_fatality_ratio_P1_rnde))
setnames(tmp, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_P1_rnde'))
pa <- merge(pa,tmp,by=c('iterations','age.idx'))
pa[, hosp_fatality_ratio_P1 := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype + logit_hosp_fatality_ratio_P1_rnde)]
pa[, hosp_fatality_ratio_wildtype := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype)]
pa <- melt(pa, id.vars=c('iterations','age.idx'), measure.vars=c('hosp_fatality_ratio_P1','hosp_fatality_ratio_wildtype'))
setkey(pa, variable, iterations, age.idx)
tmp <- pa[, list(age.idx=age.idx, value=cumsum(value)), by=c('variable','iterations')]
set(tmp, NULL, 'variable', tmp[,paste0('c',variable)])
pa <- rbind(pa,tmp)
pas <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','age.idx')]
tmp <- unique(subset(dh, select=c(age.idx, age.label)))
pas <- merge(pas, tmp, by='age.idx')
pas <- dcast.data.table(pas, variable+age.idx+age.label~stat, value.var='value')
pas[, variant := factor(grepl('P1', variable), levels=c(FALSE,TRUE), labels=c('non-P1','P1'))]
pas[, plot.facet := gsub('_wildtype|_P1','',variable)]

p <- ggplot(subset(pas, plot.facet=='hosp_fatality_ratio'), aes(x=age.label, fill=variant)) +
  geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(labels=scales::percent) +
  ggsci::scale_fill_npg() +	
  labs(x='',y='hospital fatality ratio',fill='variant')
cat('\nplotting _hfr_P1_vs_wildtype.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_P1_vs_wildtype.pdf'),p,w=8,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_P1_vs_wildtype.png'),p,w=8,h=6)

pas <- dcast.data.table(pa, iterations+age.idx~variable, value.var='value')
pas[, HFR_ratio := hosp_fatality_ratio_P1 / hosp_fatality_ratio_wildtype]
pas <- pas[, list(value=quantile(HFR_ratio, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('age.idx')]
tmp <- unique(subset(dh, select=c(age.idx, age.label)))
pas <- merge(pas, tmp, by='age.idx')
pas <- dcast.data.table(pas, age.idx+age.label~stat, value.var='value')

p <- ggplot(pas, aes(x=age.label)) +
  geom_hline(yintercept=1, lwd=1, colour='grey50') +
  geom_boxplot(aes(ymin=CL, lower=IL, middle=M, upper=IU, ymax=CU), stat='identity') +
  scale_y_continuous(expand=c(0,0), breaks=1:100) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +		
  labs(x='',y='HFR ratio P1:wildtype') 
cat('\nplotting _hfr_P1_vs_wildtype_ratio.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_P1_vs_wildtype_ratio.pdf'),p,w=7,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_P1_vs_wildtype_ratio.png'),p,w=7,h=6)
}
####################################################################
# fitted HFR by variant and over time
####################################################################
if(1)
{

pa <- as.data.table(reshape2::melt(po$hfr_wildtype))
setnames(pa, 1:4, c('iterations','week.idx','age.idx','hfr.wildtype'))
tmp <- as.data.table(reshape2::melt(po$hfr_overall))
setnames(tmp, 1:4, c('iterations','week.idx','age.idx','hfr.overall'))
pa <- merge(pa,tmp,by=c('iterations','week.idx','age.idx'))
pa <- melt(pa, id.vars=c('iterations','week.idx','age.idx'), measure.vars=c('hfr.wildtype','hfr.overall'))
pas <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('variable','age.idx','week.idx')]
pas <- dcast.data.table(pas, variable + age.idx + week.idx  ~ stat, value.var='value')
tmp <- unique(subset(dh, select=c(age.idx, week.idx, week, week.start, age.label, hosps, deaths_of_hosps)))
pas <- merge(pas, tmp, by=c('age.idx','week.idx'))

tmp <- subset(pas, variable=='hfr.overall')
tmp2 <- subset(pas, variable=='hfr.wildtype')

p <- ggplot(tmp, aes(x=week.start)) +
  geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
  geom_line(data=tmp2, aes(y=M), linetype='11') +
  geom_line(aes(y=M)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.background = element_blank(),
        legend.position='bottom') +
  scale_x_date(breaks="3 months", expand=c(0,0),  labels=date_format("%b %y")) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  #viridis::scale_colour_viridis(option='B', end=0.85, labels = scales::percent, limits=c(0,1)) +
  facet_wrap(~age.label, ncol=6) +
  labs(x='',y='weekly in-hospital fatality rate', colour='Gamma attributable\nhospital admissions')
cat('\nplotting _hfr_over_time_by_variant.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_over_time_by_variant.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_over_time_by_variant.png'),p,w=10,h=6)


p <- ggplot(pas, aes(x=week.start)) +	
  geom_vline(xintercept=subset(ddates, loc_label==args$loc_label)[, w2.start], colour='grey50') +
  geom_point(aes(y=deaths_of_hosps/hosps, colour=age.label)) +
  geom_ribbon(data=tmp, aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
  geom_line(data=tmp2, aes(y=M), linetype='11') +
  geom_line(data=tmp, aes(y=M)) +	
  viridis::scale_colour_viridis(discrete=TRUE, option='B', end = .8) +
  facet_wrap(~age.label,ncol=6) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'bottom') +
  scale_x_date(breaks="3 months", expand=c(0,0), labels=date_format("%b %y")) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +			
  labs(x='', y="weekly in-hospital fatality rate", colour='age')
cat('\nplotting _hfr_over_time_by_variant_vs_empirical_noref.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_over_time_by_variant_vs_empirical_noref.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_over_time_by_variant_vs_empirical_noref.png'),p,w=10,h=6)


tmp3 <- subset(pas, variable=='hfr.wildtype' & week==dhr[loc_label==args$loc_label, ref.week])
p <- ggplot(pas, aes(x=week.start)) +	
	geom_vline(xintercept=subset(ddates, loc_label==args$loc_label)[, w2.start], colour='grey50') +
	geom_vline(xintercept=subset(dhr, loc_label==args$loc_label)[, ref.week.start], colour='black', lty='dotted') +
	geom_hline(data=tmp3, aes(yintercept=M), colour='black', lty='dotted') +
	geom_point(aes(y=deaths_of_hosps/hosps, colour=age.label)) +
	geom_ribbon(data=tmp, aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
	geom_line(data=tmp2, aes(y=M), linetype='11') +
	geom_line(data=tmp, aes(y=M)) +	
	viridis::scale_colour_viridis(discrete=TRUE, option='B', end = .8) +
	facet_wrap(~age.label,ncol=6) +
	guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'none') +
	scale_x_date(breaks="3 months", expand=c(0,0), labels=date_format("%b %y")) +
	scale_y_continuous(labels=scales::percent, expand=c(0,0)) +			
	labs(x='', y="weekly in-hospital fatality rate", colour='age')
cat('\nplotting _hfr_over_time_by_variant_vs_empirical.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_over_time_by_variant_vs_empirical.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_hfr_over_time_by_variant_vs_empirical.png'),p,w=10,h=6)
}

####################################################################
# loess fit over time
####################################################################
if(1)
{
tmp <- unique(subset(dh, select=c(loc_label, age.idx, week.idx, week, week.start, age.label, hosps, deaths_of_hosps)))
tmp[, hfr_empirical := deaths_of_hosps/hosps]

tmp2 <- data.frame()
for (age in unique(da$age.label)){
  tmp3 <-  tmp[age.label == age]
  tmp3_noNA <- tmp3[!is.na(hfr_empirical)]
  
  # If want some type of prediction intervals, work with this and expand
  # However assumption on normality of errors seems far fetched
  fit <- loess(hfr_empirical~week.idx, tmp3_noNA, span = 0.3)
  fit <- predict(fit, se=TRUE, tmp3$week.idx)
  fit <- data.table(week.idx = tmp3$week.idx, M = fit$fit, se = fit$se)
  fit <- fit[,`:=` (CL=M-1.96*se, CU=M+1.96*se, age.label=age)]
  
  tmp2 <- rbind(tmp2,fit)
  rm(tmp3_noNA, fit, tmp3)
}
tmp <- merge(tmp, tmp2, by = c("age.label", 'week.idx'))

p <- ggplot(tmp, aes(x=week.start)) +
  geom_ribbon(data=tmp, aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
  geom_line(aes(y=M)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.background = element_blank(),
        legend.position='bottom') +
  scale_x_date(breaks="3 months", expand=c(0,0), labels=date_format("%b %y") ) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  facet_wrap(~age.label, ncol=6) +
  labs(x='',y='smoothed weekly in-hospital fatality rate')
cat('\nplotting _loess_hfr_over_time_by_variant.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_loess_hfr_over_time_by_variant.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_loess_hfr_over_time_by_variant.png'),p,w=10,h=6)


p <- ggplot(tmp, aes(x=week.start)) +	
  geom_vline(xintercept=subset(ddates, loc_label==args$loc_label)[, w2.start], colour='grey50') +
  geom_point(aes(y=deaths_of_hosps/hosps, colour=age.label)) +
  geom_ribbon(data=tmp, aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
  geom_line(aes(y=M)) +	
  viridis::scale_colour_viridis(discrete=TRUE, option='B', end = .8) +
  facet_wrap(~age.label,ncol=6) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'bottom') +
  scale_x_date(breaks="3 months", expand=c(0,0), labels=date_format("%b %y")) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +			
  labs(x='', y="smoothed weekly in-hospital fatality rate", colour='age')
cat('\nplotting _loess_hfr_over_time_by_variant_vs_empirical_noref.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_loess_hfr_over_time_vs_empirical.pdf'),p,w=10,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_loess_hfr_over_time_vs_empirical.png'),p,w=10,h=6)

}else{

####################################################################
# extract and compare effect of predictors on logit linear predictor
####################################################################

pa <- as.data.table(reshape2::melt(po$fr_regcoeff_overall))
setnames(pa, 1:3, c('iterations','p.idx','value'))
set(pa, NULL, 'value', pa[,value])
pas <- pa[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('p.idx')]
tmp <- data.table(p.idx=seq_along(args$fr.predictors), p.label=names(args$fr.predictors))
pas <- merge(pas, tmp, by=c('p.idx'))
pas <- dcast.data.table(pas, p.idx+p.label~stat, value.var='value')
p <- ggplot(pas, aes(x=p.label, fill=p.label)) +
  geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU), stat='identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = args$city.palette(length(unique(pas$p.label)))) +
  guides(fill='none') +
  labs(x='', y='regression coeff', fill='predictor')
cat('\nplotting _frmult_regrcoeff.pdf ...')
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_frmult_regrcoeff.pdf'),p,w=6,h=6)
ggsave(file=paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_frmult_regrcoeff.png'),p,w=6,h=6)
}