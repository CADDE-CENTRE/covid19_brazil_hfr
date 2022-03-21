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
library(ggpubr) 
require(scales)

# source workspace
out.dir2 <- '/rds/general/project/ratmann_covid19/live/BrazilP1/hfr-fit-210719h4-resources8-ex0-trvaluecs-hd1'
out.dir2 <- '~/Documents/P1Brazil/hfr-fit-210719h4-resources8-ex0-trvaluecs-hd1/'
file.prefix <- basename(out.dir2)
out.base2 <- file.path(out.dir2, file.prefix)

tmp <- list.files(out.dir2, pattern='allloc_workspace.rda$', full.names=T )
cat('\nAttempting to load ', tmp)
load(file = tmp)
if(dir.exists(dirname(out.base2)))
{
  args$out.base <- out.base2
}


###############################################################################
# Helper Functions
###############################################################################

reqs <- theme(axis.text = element_text(size=5, family='sans'), text=element_text(size=7,family='sans'))
reqs2 <- theme(axis.text = element_text(size=2.5, family='sans'), text=element_text(size=3.5,family='sans'))
update_geom_defaults("point", list(size = 0.2))

load.env <- function(loc_label)
{
  tmp <- paste0(out.base2, '_', gsub('_|-| ', '', loc_label),'_workspace.rda')
  cat('\nAttempting to load ', tmp)
  load(file = tmp)
  if(dir.exists(dirname(out.base2)))
  {
    args$out.base <- out.base2
  }
  return(list(dh=dh, da=da, dhr=dhr))
}

load.chains <- function(loc_label)
{
  
  tmp <- paste0(out.base2, '_', gsub('_|-| ', '', loc_label),'_workspace.rda')
  cat('\nAttempting to load ', tmp)
  load(file = tmp)
  if(dir.exists(dirname(out.base2)))
  {
    args$out.base <- out.base2
  }
  
  ### Load fit and postprocess everything
  tmp <- paste0(args$out.base,'_',gsub(' ','',args$loc_label),'_fit.rds')
  cat('\nAttempting to load ', tmp)
  m_fit <- readRDS( tmp)
  pd <- m_fit$draws(inc_warmup = FALSE)		
  
  # get MCMC samples for particular chains into rstan::extract format 
  select.chains <- seq_along(dimnames(pd)[['chain']])
  iters <- 1:(m_fit$metadata()[['iter_sampling']])
  
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
  
  return(list(po=po, dh=dh, da=da, dhr=dhr))
}

plot.loess.over.time <- function(loc_label, out)
{
  po <- out$po; dh <- out$dh; da <- out$da; dhr <- out$dhr
  # loess fit over time
  
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
    geom_vline(xintercept=subset(ddates, loc_label==args$loc_label)[, w2.start], colour='grey50') +
    geom_ribbon(data=tmp, aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
    geom_point(aes(y=deaths_of_hosps/hosps, colour=age.label)) +
    geom_line(aes(y=M)) +	
    viridis::scale_colour_viridis(discrete=TRUE, option='B', end = .8) +
    facet_wrap(~age.label,ncol=6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'none') +
    scale_x_date(breaks="3 months", expand=c(0,0), labels=date_format("%b %y")) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +			
    coord_cartesian(ylim=c(0,1)) + 
    labs(x='', y=paste0("smoothed weekly in-hospital fatality rate in ", change_city_label(loc_label)))
  
  cat('\nplotting loess_hfr_over_time_vs_empirical_noref.pdf ...')
  txt <- '_loess_hfr_over_time_vs_empirical_noref'
  
  return(p)
}

plot.fittedhfr.byvar <- function(loc_label,out)
{
  po <- out$po; dh <- out$dh; da <- out$da; dhr <- out$dhr  # fitted HFR by variant and over time
  
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
  
  z <- which(dhr$loc_label == loc_label); ref.week <- dhr$ref.week
  tmp3 <- subset(pas, variable=='hfr.wildtype' & week==ref.week )
  
  pas$w2.start <- subset(ddates, loc_label==loc_label)[, unique(w2.start)[1]]
  pas$ref.week.start <- subset(dhr, loc_label==loc_label)[, unique(ref.week.start)[1]]
  
  p <- ggplot(pas, aes(x=week.start)) +	
    # geom_vline(xintercept=subset(ddates, loc_label==loc_label)[, w2.start], colour='grey50') +
    # geom_vline(xintercept=subset(dhr, loc_label==loc_label)[, ref.week.start], colour='black', lty='dotted') +
    geom_vline(aes(xintercept=w2.start), colour='grey50') +
    # geom_vline(aes(xintercept=ref.week.start), colour='black', lty='dotted') +
    # geom_hline(data=tmp3, aes(yintercept=M), colour='black', lty='dotted') +
    geom_point(aes(y=deaths_of_hosps/hosps, colour=age.label)) +
    geom_ribbon(data=tmp, aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
    geom_line(data=tmp2, aes(y=M), linetype='11') +
    geom_line(data=tmp, aes(y=M)) +	
    viridis::scale_colour_viridis(discrete=TRUE, option='B', end = .8) +
    facet_wrap(~age.label,ncol=6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'none') +
    scale_x_date(breaks="3 months", expand=c(0,0), labels=date_format("%b %y")) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +			
    labs(x='', y=paste0("weekly in-hospital fatality rate in ", change_city_label(loc_label)))
  
  # pas[, loc_label := LOC]
  # DFS[[LOC]] <- pas
  cat('\nplotting _hfr_over_time_by_variant_vs_empirical.pdf ...')
  txt <- '_hfr_over_time_by_variant_vs_empirical'
  return(p)
}

plot.hosps.fittedshare <- function(loc_label,out)
{
  po <- out$po; dh <- out$dh; da <- out$da; dhr <- out$dhr
  # fitted share of age groups among hospital admissions vs empirical
  
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
    # geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey80') +
    # geom_ribbon(aes(ymin=IU, ymax=CU), fill='grey80') +
    # geom_ribbon(aes(ymin=IL, ymax=IU), fill='grey50') +
    geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
    geom_line(aes(y=M), colour='black') +		
    geom_point(aes(y=phosps, colour=age.label)) +		
    viridis::scale_colour_viridis(option='B', discrete = TRUE, end=0.8 ) +
    scale_x_date(breaks='3 months',expand=c(0,0), labels=date_format("%b %y")) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    labs(y=paste0('share of age groups among COVID-19\nattributable hospital admissions in ', change_city_label(loc_label)), x='') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='none', strip.background = element_blank()) +
    facet_wrap(~age.label, ncol=6)
  # pads[, loc_label := LOC]
  # DFS[[LOC]] <- pads
  
  cat('\nplotting _hosps_phosps_vs_fitted.pdf ...')
  txt <- '_hosps_phosps_vs_fitted'
  return(p)
}

################################################################################
# Plot EDF
################################################################################

# out <- load.env('manaus')
# FIG2C <- plot.loess.over.time('manaus', out)
# FIG2C <- FIG2C +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#              panel.background = element_blank(), axis.line = element_line(colour = "black"),
#              strip.background = element_blank(),
#              strip.text.x = element_blank()) + labs(y="") +
#   geom_vline(aes(xintercept = ddates[loc_label == 'manaus', w2.start])) + 
#   scale_x_date(breaks="3 months", expand=c(0,0), labels=date_format("%b %y")) 
#   
# 
# ggsave(file=file.path(out.dir2, 'Brizzietal_mainfig1C.pdf'),
#        FIG2C,w=18,h=10, units='cm', dpi=350)



# select locations and then plot
LOCS <- c('goiania','natal','rio de janeiro', 'sao paulo')

EDF_oneloc <- function(loc_label)
{
  out <- load.chains(loc_label)

  p1 <- plot.loess.over.time(loc_label, out=out) # IDX=1
  # p2 <- plot.fittedhfr.byvar(loc_label, out=out)
  # p3 <- plot.hosps.fittedshare(loc_label, out=out)
  
  return(list(p1=p1))
}

PLOTS <- lapply(LOCS, EDF_oneloc)
names(PLOTS) <- LOCS


# EDF 3: Loess over time
p1 <- lapply(PLOTS, FUN = `[[`, "p1")
p1 <- ggpubr::ggarrange( NULL, p1[[1]] + reqs, NULL, p1[[2]] + reqs,
                         NULL, p1[[3]] + reqs, NULL, p1[[4]] + reqs,
                         widths = c(0.02, 1, 0.02, 1),
                         labels = c("","a","", "b","", "c","", "d"),
                         label.x = -0.02,
                         font.label=list(size=8, family='sans'),
                         ncol=4,
                         nrow=2)
ggsave(file=file.path(out.dir2, 'Brizzietal_extdatafig_3.pdf'),
       p1 ,w=24,h=17, units='cm', dpi=350)


# select locations and then plot
LOCS <- c('goiania','manaus','rio de janeiro', 'sao paulo')

EDF_oneloc <- function(loc_label)
{
  out <- load.chains(loc_label)
  
  p2 <- plot.fittedhfr.byvar(loc_label, out=out)
  p3 <- plot.hosps.fittedshare(loc_label, out=out)
  
  return(list(p1=p1, p2=p2, p3=p3))
}

PLOTS <- lapply(LOCS, EDF_oneloc)
names(PLOTS) <- LOCS


# EDF 9: Fitted HFR by variant
p2 <- lapply(PLOTS, FUN = `[[`, "p2")
p2 <- ggpubr::ggarrange( NULL, p2[[1]] + reqs, NULL, p2[[2]] + reqs,
                         NULL, p2[[3]] + reqs, NULL, p2[[4]] + reqs,
                         widths = c(0.02, 1, 0.02, 1),
                         labels = c("","a","", "b","", "c","", "d"),
                         label.x = -0.02,
                         font.label=list(size=8, family='sans'),
                         ncol=4,
                         nrow=2)
ggsave(file=file.path(out.dir2, 'Brizzietal_extdatafig_8.pdf'),
       p2 ,w=24,h=17, units='cm', dpi=350)


# EDF 10: Fitted share of hospital admissions by variant
p3 <- lapply(PLOTS, FUN = `[[`, "p3")
p3 <- ggpubr::ggarrange( NULL, p3[[1]] + reqs, NULL, p3[[2]] + reqs,
                         NULL, p3[[3]] + reqs, NULL, p3[[4]] + reqs,
                         widths = c(0.02, 1, 0.02, 1),
                         labels = c("","a","", "b","", "c","", "d"),
                         label.x = -0.02,
                         font.label=list(size=8, family='sans'),
                         ncol=4,
                         nrow=2)
ggsave(file=file.path(out.dir2, 'Brizzietal_extdatafig_9.pdf'),
       p3 ,w=24,h=17, units='cm', dpi=350)
