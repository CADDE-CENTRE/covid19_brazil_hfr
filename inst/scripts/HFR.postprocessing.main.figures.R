cat(" \n -------------------------------- \n \n HFR.postprocessing.main.figures.R \n \n -------------------------------- \n")

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
require(ggrepel)

if(1) # Andrea s 
{
  out.dir2 <- '~/Documents/P1Brazil/hfr-fit-210719h2-resources8-ex0-trvaluecs-hd1/'
  file.prefix2 <- basename(out.dir2)
  threshold_check <-  0
  pkg.dir2 <- '~/git/covid19_brazil_hfr/'
}


# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-outdir')	
  stopifnot(args_line[[3]]=='-file_prefix')
  stopifnot(args_line[[5]]=='-pkgdir')
  stopifnot(args_line[[7]]=='-threshold_check')
  out.dir2 <- args_line[[2]]
  file.prefix2 <- args_line[[4]]
  pkg.dir2 <- args_line[[6]]
  threshold_check <- as.integer(args_line[[8]])
} 

pkg.dir <- pkg.dir2
####################################################################
# Check that at least threshold number of chains are done running
####################################################################

tmp <- list.files(out.dir2, pattern='_fit.rds')
cat('\ntry to bugf<ix: start\n')
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

# process dd
tmp <- dd[, list(week=unique(week), week.idx= unique(week)-min(week)+1L), by='loc_label']
dd <- merge(dd, tmp, by=c('loc_label','week'))
dd[, age.idx := as.integer(age.label)]
dd[, adjdeathsi := as.integer(round(adjdeaths))]
dd[, deathshospadmi := as.integer(round(deaths_hospadm))]

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
		
	pos[[do$loc_label[i]]] <- po
	gc()
	rm(po, pd, m_fit)
}

####################################################################
# Figure 1
####################################################################

# Figure 1A (top left) -- map

# Figure 1B (top right) -- replacement dynamics

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
tmp <- unique(subset(dd, select=c(loc_label,week.idx, week, week.start)))
pads <- merge(pads, tmp, by=c('loc_label','week.idx'))
tmp <- subset(dp, select=c(loc_label, week, n_positive, n_sampled))
pads <- merge(pads, tmp, by=c('loc_label','week'), all.x=TRUE)

selectec_loc_label <- c('belo horizonte','manaus','sao paulo')
pads1B <- subset(pads, loc_label%in%selectec_loc_label)
tmp <- subset(pads1B, !is.na(n_sampled) & n_positive>0 & n_positive!=n_sampled)
p1B <- ggplot(pads1B, aes(x=week.start)) +
		geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey50') +
		geom_line(aes(y=M)) +
		geom_point(data=subset(pads1B, !is.na(n_sampled)), aes(y=n_positive/n_sampled, size=n_sampled, fill=loc_label), pch=21) +		
		scale_fill_manual(values = args$city.palette(length(unique(pads1B$loc_label)))) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		facet_wrap(~change_city_label(loc_label), scales='free_y', ncol=5) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(fill='none') +
		labs(y='Gamma variant frequency', x='', colour='', size='number of SARS-CoV-2 sequences')
ggsave(file=paste0(args$out.base,'_allloc_figure1Bv2.pdf'),p1B,w=8,h=4)
ggsave(file=paste0(args$out.base,'_allloc_figure1Bv2.png'),p1B,w=8,h=4)

p1B <- ggplot(pads1B, aes(x=week.start)) +
		geom_ribbon(aes(ymin=CL, ymax=CU), fill='grey50') +
		geom_line(aes(y=M)) +
		geom_point(data=subset(pads1B, !is.na(n_sampled)), aes(y=n_positive/n_sampled, size=n_sampled, fill=loc_label), pch=21) +
		ggrepel::geom_text_repel(data=tmp, aes(y=n_positive/n_sampled, label=n_sampled),max.overlaps=10) +
		scale_fill_manual(values = args$city.palette(length(unique(pads1B$loc_label)))) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		facet_wrap(~change_city_label(loc_label), scales='free_y', ncol=5) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				strip.background = element_blank(),
				legend.position='bottom') +
		guides(fill='none',size='none') +
		labs(y='Gamma genotype frequency', x='', colour='')
ggsave(file=paste0(args$out.base,'_allloc_figure1B.pdf'),p1B,w=8,h=4)
ggsave(file=paste0(args$out.base,'_allloc_figure1B.png'),p1B,w=8,h=4)


# Figure 1C (bottom) -- in hospital fatality rates BH

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
selectec_loc_label <- c('manaus')
pads1C <- subset(pads, loc_label%in%selectec_loc_label)

p1C <- ggplot(pads1C, aes(x=week.start)) +	
  geom_vline(xintercept=subset(ddates, loc_label%in%selectec_loc_label)[, w2.start], colour='black',linetype='11') +
  geom_point(data=subset(pads1C,variable=='hfr.overall'),aes(y=deaths_of_hosps/hosps, colour=age.label)) +
  geom_ribbon(data=subset(pads1C,variable=='hfr.overall'), aes(ymin=CL, ymax=CU), fill='grey50', alpha=.5) +
  geom_line(data=subset(pads1C,variable=='hfr.wildtype'), aes(y=M), linetype='11') +
  geom_line(data=subset(pads1C,variable=='hfr.overall'), aes(y=M)) +	
  viridis::scale_colour_viridis(discrete=TRUE, option='B', end = .8) +
  facet_wrap(~age.label,ncol=6) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'bottom') +
  scale_x_date(breaks="3 months", expand=c(0,0)) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +			
  labs(x='week of hospital admission', y=paste0("weekly in-hospital fatality rate\nin ",change_city_label(selectec_loc_label)), colour='age')

ggsave(file=paste0(args$out.base,'_allloc_figure1C.pdf'),p1C,w=10,h=6)
ggsave(file=paste0(args$out.base,'_allloc_figure1C.png'),p1C,w=10,h=6)

require(magick)
image1A <- image_read_pdf("~/Documents/P1Brazil/tmp/Figure1_Brazil_map.pdf")
image1B <- image_read_pdf(paste0(args$out.base,'_allloc_figure1B.pdf')) 
image1C <- image_read_pdf(paste0(args$out.base,'_allloc_figure1C.pdf')) 

tmp <- c(image1A, image1B)
g <- image_append(image_scale(tmp,'x1000'))
g <- image_append(image_scale(c(g,image1C), '1000'),stack = TRUE)
tmp <- list()
g <- image_annotate(g, 'A', location='+10+10', size = 18, weight = 700, font = 'mono')
g <- image_annotate(g, 'B', location='+350+10', size = 18, weight = 700, font = 'mono')
g <- image_annotate(g, 'C', location='+10+320', size = 18, weight = 700, font = 'mono')

image_write(g, path =paste0(args$out.base,'_allloc_figure1.png'), format = "png")

####################################################################
# Figure 4
####################################################################

# Figure 4A (top) -- weekly hosp admissions for BH, MAN, SP
if(all(c('belo horizonte', 'manaus', 'sao paulo')  %in% args$selected.locs ))
{
  
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
  set(pads, NULL, 'variant', pads[, factor(variant, levels=c('wildtype', 'P1'), labels=c('non-P1','P1'))])
  pads[, variable := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\1.\\2',variable)]
  tmp <- unique(subset(dh, select=c(loc_label,week.idx,age.idx,week,week.start,age.label))) 
  pads <- merge(tmp, pads, by=c('loc_label','week.idx','age.idx'))
  setkey(pads, loc_label, variable, week.idx, variant, age.idx)
  tmp <- pads[, list(variant=variant, age.idx=age.idx, sM=cumsum(M)), by=c('loc_label','variable','week.idx')]
  pads <- merge(pads, tmp, by=c('loc_label','variable','week.idx','variant','age.idx'))
  tmp <- unique(subset(dh, select=c(loc_label,week.start,thosps)))
  pads <- merge(pads, tmp, by=c('loc_label','week.start'))
  
  selectec_loc_label <- c('belo horizonte','manaus','sao paulo')
  pads4A <- subset(pads, loc_label%in%selectec_loc_label)
  p4A <- ggplot(pads4A, aes(x=week.start)) +
    geom_vline(data=subset(ddates, loc_label%in%selectec_loc_label), aes(xintercept=w2.start), colour='black', lty='11') +
    geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=age.label, alpha=variant)) +
    geom_point(data=unique(subset(pads4A,select=c(loc_label,week.start,thosps))), aes(y=thosps), pch=18) +
    viridis::scale_fill_viridis(option='D', end=0.8, discrete = TRUE) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_alpha_manual(values=c('non-P1'=0.7, 'P1'=1.)) +
    facet_wrap(~change_city_label(loc_label), scales='free', ncol=4) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE),alpha='none') +
    # guides(fill='none',alpha='none') +
    labs(fill='age', y='estimated COVID-19 hospital   \nadmissions by Gamma  \nand non-Gamma variant  ', x='', alpha='')
  
  
  
  # Figure 4B -- weekly deaths for BH, MAN, SP
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
  set(pads, NULL, 'variant', pads[, factor(variant, levels=c('wildtype', 'P1'), labels=c('non-P1','P1'))])
  pads[, variable := gsub('^([a-z]+)\\.([a-z]+)\\.([a-zA-Z0-9]+)$','\\1.\\2',variable)]
  tmp <- unique(subset(dh, select=c(loc_label,week.idx,age.idx,week,week.start,age.label))) # changed dd with dh, alternatively could add week,idx to dd at the start of the script
  pads <- merge(tmp, pads, by=c('loc_label','week.idx','age.idx'))
  setkey(pads, loc_label, variable, week.idx, variant, age.idx)
  tmp <- pads[, list(variant=variant, age.idx=age.idx, sM=cumsum(M)), by=c('loc_label','variable','week.idx')]
  pads <- merge(pads, tmp, by=c('loc_label','variable','week.idx','variant','age.idx'))
  tmp <- dh[,list(tdeaths_of_hospsi=sum(deaths_of_hospsi)), by=c('loc_label','week.start')]
  pads <- merge(pads, tmp, by=c('loc_label','week.start'))
  selectec_loc_label <- c('belo horizonte','manaus','sao paulo')
  pads4B <- subset(pads, loc_label%in%selectec_loc_label)
  
  p4B <- ggplot(pads4B, aes(x=week.start)) +
    geom_vline(data=subset(ddates, loc_label%in%selectec_loc_label), aes(xintercept=w2.start), colour='black', lty='11') +
    geom_ribbon(aes(ymin=sM-M, ymax=sM, fill=age.label, alpha=variant)) +
    geom_point(data=unique(subset(pads4B,select=c(loc_label,week.start,tdeaths_of_hospsi))), aes(y=tdeaths_of_hospsi), pch=18) +
    viridis::scale_fill_viridis(option='B', end=0.8, discrete = TRUE) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_alpha_manual(values=c('non-P1'=0.7, 'P1'=1.)) +
    facet_wrap(~change_city_label(loc_label), scales='free', nrow=1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE),alpha='none') +
    # guides(fill='none',alpha='none') +
    labs(fill='age', y='estimated COVID-19 deaths\nin hospitals by Gamma\nand non-Gamma variant', x='', alpha='')
  
  
  p4 <- ggpubr::ggarrange( p4A,
                           p4B,
                           legend="bottom",
                           common.legend=FALSE,
                           labels=c('A','B'),
                           font.label=list(size=12),
                           ncol=1,
                           heights=c(3,3),
                           widths=c(2,2))
  ggsave(file=paste0(args$out.base,'_allloc_figure4.pdf'),p4,w=10,h=7.5)
  ggsave(file=paste0(args$out.base,'_allloc_figure4.png'),p4,w=10,h=7.5)
}

####################################################################
# Baseline fold and Figure 5
####################################################################
# 
# # figure 3A -- baseline HFR versus resources
# tmp.args <- copy(args)
# tmp.args$make.preprocessing.plots <- 0
# tmp.dicu <- copy(dicu)
# tmp.dicu <- merge(tmp.dicu,unique(subset(dpop,select=c(loc_label,tpop))),by='loc_label')
# set(tmp.dicu,NULL,'nicu.adm',tmp.dicu[,5*1e-5*tpop])
# set(tmp.dicu,NULL,'tpop',NULL)
# tmp.dhcadm <- copy(dhcadm)
# tmp.dhcadm <- merge(tmp.dhcadm,unique(subset(dpop,select=c(loc_label,tpop))),by='loc_label')
# set(tmp.dhcadm,NULL,'nhosp.adm',tmp.dhcadm[,10*1e-5*tpop])
# set(tmp.dhcadm,NULL,'tpop',NULL)
# 
# tmp <- make_health_care_predictors(tmp.dhcadm, tmp.dicu, dooh, dhcnres, dhcb, dw, dr, dpop, ddates, dhr, tmp.args)
# 
# dhca2 <- copy(dhca)
# 
# # standardise across location 
# tmp <- unique(subset(dhca,select=c(loc_label,variable,value_ref,value_log_ref,value_sd,value_log_sd)))
# dhca2 <- merge(dhca2,tmp,by=c('loc_label','variable'))
# tmp <- dhca2[, 
# 		list(
# 				loc_label = loc_label, 
# 				week = week,
# 				value_c = value-value_ref,
# 				value_cs = (value-value_ref)/value_sd,
# 				value_log = log(value+1),
# 				value_log_c = log(value+1)-value_log_ref,
# 				value_log_cs = (log(value+1)-value_log_ref)/value_log_sd
# 		), 
# 		by='variable']
# dhca2 <- merge(dhca2, tmp, by=c('variable','loc_label','week'))
# dhc2 <- subset(dhca2, variable%in%args$fr.predictors )
# dhc2 <- dcast.data.table(dhc2, loc_label+week+week.start+week.cat~variable, value.var=args$fr.predictors.transformation)
# dhc2 <- melt(dhc2, id.vars=c('loc_label','week'),measure.vars=args$fr.predictors)
# setkey(dhc2,loc_label,week)
# tmp <- unique(subset(dh, select=c(loc_label,week,week.idx)))
# dhc2 <- merge(dhc2,tmp,by=c('loc_label','week')) 
# tmp <- data.table(p.idx=seq_along(args$fr.predictors), variable=args$fr.predictors)
# dhc2 <- merge(dhc2,tmp,by=c('variable'))
# dhc2 <- subset(dhc2,week.idx==1)
# 
# # calculate baseline HFR, using inferred regression coefficients 
# pas.ref.hfr <- list()
# for(i in 1:length(pos))
# {	
# 	cat('\nprocessing',i,' ',names(pos)[i])
# 	pa <- as.data.table(reshape2::melt(pos[[i]]$fr_regcoeff_overall))
# 	setnames(pa, 1:3, c('iterations','p.idx','coef'))	
# 	tmp <- unique(subset(pa, is.na(coef), select=iterations))[,iterations]
# 	if(length(tmp))
# 	{
# 		cat('\nFound iterations with NA, removing:', tmp)
# 		pa <- subset(pa, !iterations%in%tmp)
# 	}	
# 	pa <- merge(pa, subset(dhc2, loc_label==names(pos)[i]), by='p.idx',allow.cartesian=TRUE)
# 	pa <- pa[, list(logit_fr_multiplier_by_week_hyp= sum(coef*value)),by=c('iterations','week.idx')]		
# 	tmp <- as.data.table(reshape2::melt(pos[[i]][['logit_hosp_fatality_ratio_wildtype']]))
# 	setnames(tmp, 1:3, c('iterations','age.idx','logit_hosp_fatality_ratio_wildtype'))
# 	pa <- merge(pa,tmp,by=c('iterations'),allow.cartesian=TRUE)	
# 	pa[, baseline_hfr := gtools::inv.logit(logit_hosp_fatality_ratio_wildtype+logit_fr_multiplier_by_week_hyp)]
# 	set(pa, NULL, c('logit_hosp_fatality_ratio_wildtype','logit_fr_multiplier_by_week_hyp'), NULL)		
# 	pa[, loc_label := names(pos)[i]]
# 	pas.ref.hfr[[i]] <- pa
# 	rm(pa, tmp)
# }
# pas.ref.hfr <- do.call('rbind',pas.ref.hfr)
# # determine location with min baseline hfr 
# 
# tmp2 <- dpop[,list(pop=sum(pop)),by='age.label']
# tmp2[, tpop:=sum(pop)]
# tmp2[, ppop.all.cities:=pop/tpop]
# setkey(tmp2,age.label)
# tmp2[,age.idx:=1:nrow(tmp2)]
# pas.ref.hfr <- merge(pas.ref.hfr,tmp2,by='age.idx')
# tmp <- pas.ref.hfr[,list(baseline_hfr=sum(ppop.all.cities*baseline_hfr)),by=c('loc_label','iterations')]
# 
# if(tmp[, length(unique(loc_label))] == 14)
# {
#   # XX for folds in section 
#   folds <- list()
#   tmp1 <- dcast(tmp, iterations ~ loc_label, value.var = 'baseline_hfr')
#   tmp1 <- tmp1[, list( ref=`belo horizonte`/`belo horizonte`,`curitiba`/`belo horizonte`,`florianopolis`/`belo horizonte`, `goiania`/`belo horizonte`,`joao pessoa`/`belo horizonte`,`macapa`/`belo horizonte`,`manaus`/`belo horizonte`,`natal`/`belo horizonte`,
#                        `porto alegre`/`belo horizonte`,  `porto velho`/`belo horizonte`,`rio de janeiro`/`belo horizonte`,`salvador`/`belo horizonte`,`sao luis`/`belo horizonte`,      `sao paulo`/`belo horizonte`) ]
#   tmp1[, ref:= NULL]
#   tmp1 <- quantile(rowMeans(tmp1, na.rm = T),  p=c(.025,.5,.975), na.rm = T)
#   tmp1 <- round(tmp1, 2)
#   tmp1 <- paste0(tmp1[2],' (',tmp1[1],'-',tmp1[3],')')
#   folds[['baseline_ratio_overall']] <- tmp1 
#   
#   tmp1 <- tmp[loc_label != 'belo horizonte']
#   tmp1 <- tmp[, list(value=quantile(baseline_hfr, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label')]
#   tmp1 <- dcast(tmp1, loc_label ~ stat)
#   tmp1 <- merge(tmp1, dstates[,.(loc_label, region_label)], by = 'loc_label')
#   set(tmp1,tmp1[,which(region_label%in%c('Central-West','North'))],'region_label','North + Central-West')
#   
#   p_suppA <- ggplot(tmp1, aes(x = loc_label, color = change_city_label(loc_label) )) + 
#     geom_point(aes(y = M)) + 
#     geom_errorbar(aes(ymin = CL, ymax=CU)) +
#     scale_y_continuous(labels=scales::percent) +
#     facet_wrap(~region_label, scales = 'free_x') +
#     scale_colour_manual(values = args$city.palette(length(unique(tmp1$loc_label)))) +		
#     theme_bw() +
#     theme(legend.position='bottom',strip.background = element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) +
#     guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
#     labs(x='',y='baseline in-hospital fatality rate',colour='') 
# }
# 
# 
# tmp <- tmp[,list(baseline_hfr=median(baseline_hfr)),by=c('loc_label')]
# loc_label.min.ref.hfr <- tmp[, loc_label[which.min(baseline_hfr)]]
# # calculate ratio in baseline hfr
# tmp <- subset(pas.ref.hfr, loc_label==loc_label.min.ref.hfr, select=c(iterations,age.idx,baseline_hfr))
# setnames(tmp,'baseline_hfr','best_baseline_hfr')
# pas.ref.hfr <- merge(pas.ref.hfr,tmp,by=c('iterations','age.idx'))
# pas.ref.hfr[, baseline_hfr_ratio := baseline_hfr/best_baseline_hfr]
# pads <- melt(pas.ref.hfr, id.vars=c('loc_label','week.idx','age.idx','iterations'), measure.vars=c('baseline_hfr','baseline_hfr_ratio'))
# pads <- pads[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label','variable','week.idx','age.idx')]
# pads <- dcast.data.table(pads, variable+loc_label+week.idx+age.idx ~ stat, value.var='value')
# tmp <- subset(dh, select=c(loc_label,age.idx,age.label,week.idx,week,week.start))
# pads <- merge(pads, tmp, by=c('loc_label','week.idx','age.idx'))
# pads4B <- copy(pads)



if(0){# Load if sexpr was manually modified
  sexpr <- readRDS(paste0(args$out.base,'_substitute_expressions.rds'))
}
tmp <- sexpr$correlations
setcolorder(tmp, neworder =sort(names(tmp)[-1],))
tmp <- melt.data.table(tmp, id.vars='loc_label' )
# tmp[,variable:=gsub('admissions per','admissions per\n',variable)]
tmp[,variable:=gsub('2-wk ','',variable)]
setkey(tmp, variable)
# tmp[, round(as.numeric(value),2)]
tmp[,variable := factor(variable, levels = sort(unique(variable),decreasing=TRUE) )]
tmp[,variable := paste0("3-wk ", variable)]
tmp[,variable := gsub('intensivist', 'intensive care specialist', variable)]
tmp[,variable := gsub('critical bed', 'critical care bed', variable)]

tmp1 <- sort(unique(tmp$variable))

tmp[, variable := factor(variable, levels = tmp1  )]
tmp3B <- copy(tmp)


p3B <- ggplot(data=tmp, aes(x=loc_label,y=variable, fill = as.numeric(value))) +
  geom_tile() + scale_fill_gradient(low="blue", high="red") +
  labs(fill = 'Pearson\ncorrelation', y = '', x = '') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.background = element_blank(),
        legend.position='right')

# figure 3A -- age standardised HFR versus ICU adm per critical care specialist
tmp <- unique(subset(dh, select=c(loc_label, age.idx, week.idx, week, week.start, age.label, hosps, deaths_of_hosps)))
tmp[, hfr_empirical := deaths_of_hosps/hosps]

tmp2 <- data.frame()
for (loc in ddates$loc_label){
  for (age in unique(da$age.label)){
    
    tmp3 <-  tmp[loc_label == loc & age.label == age]
    tmp3_noNA <- tmp3[!is.na(hfr_empirical)]
    
    fit <- loess(hfr_empirical~week.idx, tmp3_noNA, span = 0.3)
    fit <- predict(fit, se=TRUE, tmp3$week.idx)
    fit <- data.table(week.idx = tmp3$week.idx, M = fit$fit, se = fit$se)
    fit <- fit[,`:=` (CL=M-1.96*se, CU=M+1.96*se, age.label=age, loc_label=loc)]
    
    tmp2 <- rbind(tmp2,fit)
    rm(tmp3_noNA, fit, tmp3)
    }
}
tmp <- merge(tmp, tmp2, by = c("loc_label", "age.label", 'week.idx'))

# Nearest-Neighbour attribution for NA hfr_empirical
setkey(tmp, loc_label, age.label, week.idx)
tmp2 <- tmp[,{
  a <- which(is.na(hfr_empirical)); b <- which(!is.na(hfr_empirical));
  l <- length(hfr_empirical);
  idx <- sapply(1:l, FUN = function(x){b[which.min(abs(b-x))]})
  list(hfr_empirical.adj=hfr_empirical[idx],week.idx=week.idx)
}, by = c('loc_label', 'age.label')]
tmp <- merge(tmp, tmp2, by=c('loc_label', 'age.label', 'week.idx') )

#Standardise by age groups
tmp2 <- dpop[,list(pop=sum(pop)),by='age.label']
tmp2[, tpop:=sum(pop)]
tmp2[, ppop.all.cities:=pop/tpop]
setkey(tmp2,age.label)
tmp2[,age.idx:=1:nrow(tmp2)]
tmp <- merge(tmp, tmp2[,.(age.idx, ppop.all.cities)], by = 'age.idx')
# se again computed assuming Gaussianity
tmp2 <- tmp[, list(age.idx=0,
                   age.label='overall',
                   hfr_empirical=sum(hfr_empirical*ppop.all.cities),
                   hfr_empirical.adj=sum(hfr_empirical.adj*ppop.all.cities),
                   M = sum(M*ppop.all.cities),
                   se = sqrt(sum(se^2*ppop.all.cities^2))
                  ),by = c('loc_label', 'week.idx', 'week.start')]
tmp2[, `:=` (CL=M-1.96*se, CU=M+1.96*se )]


# plot predictor against hfr-loess
z <- "icu.per.ventilator.02" 
names(z) <- '3-wk ICU admissions per ICU bed'
i<- 1
tmp <- subset(dhca, variable==z[i],select=-hfr)
tmp <- merge(tmp,tmp2,by=c('loc_label','week.start'))
setnames(tmp,'M','hfr_loess')
tmp2 <- tmp3B[variable == gsub('2-wk ','',names(z)), .(loc_label, value)] 
setkey(tmp2, loc_label)
tmp2[, loc_label :=sort(args$selected.locs)]
#tmp2 <- tmp[,{idx <- which(!is.na(hfr_loess));list(value=cor(value[idx], hfr_loess[idx]))},by = 'loc_label']
tmp3 <- tmp[, range(value)]
tmp4 <- tmp[, range(hfr_loess, na.rm=TRUE)*c(1,1.2)]
p3A <- ggplot(tmp, aes(x=week.start)) +
  geom_vline(data=ddates, aes(xintercept=w2.start), colour='black', lty='11') +
  geom_ribbon(data=subset(tmp, !is.na(hfr_loess)), aes(ymin=( CL - tmp4[1] ) / diff(tmp4) * diff(tmp3) + tmp3[1], ymax=( CU - tmp4[1] ) / diff(tmp4) * diff(tmp3) + tmp3[1]), fill = 'grey70', alpha=0.8) +
  geom_line(data=subset(tmp, !is.na(hfr_loess)), aes(y=( hfr_loess - tmp4[1] ) / diff(tmp4) * diff(tmp3) + tmp3[1]), color = 'black') +
  geom_step(aes(y=value, colour=change_city_label(loc_label)), lwd = .5) +
  # geom_point(aes(y=value, colour=change_city_label(loc_label))) +
  # geom_text(data=tmp2,hjust=0.5, vjust=0.5, colour='black',
  #           aes(x= min(tmp$week.start)+90, y= max(tmp$value)*.95,
  #               label= paste0('   r=',value))) + # round(value,d=2)
  # scale_x_date(breaks='3 months',expand=c(0,0)) +
  scale_y_continuous(expand=range(tmp$week.start)) +
  scale_y_continuous(sec.axis= sec_axis(~ (.-tmp3[1])/diff(tmp3)*diff(tmp4)+tmp4[1], name = "age standardised in-hospital fatality rate", labels = scales::percent), expand=c(0,0)) +
  scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
  facet_wrap(~change_city_label(loc_label), ncol=5) +
  labs(x='', y=ifelse(is.na(names(z)[i]), z[i], names(z)[i]), colour='') +
  guides(colour='none') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.background = element_blank(),
        legend.position='bottom',
        axis.line.y.right = element_line(color = 'black'), 
        axis.ticks.y.right = element_line(color = 'black'),
        axis.text.y.right = element_text(color = 'black'))


p3 <- ggpubr::ggarrange( p3A,
		p3B,
		# common.legend=TRUE,
		labels=c('A','B'),
		font.label=list(size=12),
		ncol=1,
		heights=c(4.5,3))

# Keep only p3B
ggsave(file=paste0(args$out.base,'_allloc_figure3.pdf'),p3,w=8,h=9)
ggsave(file=paste0(args$out.base,'_allloc_figure3.png'),p3,w=8,h=9)

		
####################################################################
# figure 5
####################################################################

# figure 5A age standardised HFR, grouped by region
tmp2 <- dpop[,list(pop=sum(pop)),by='age.label']
tmp2[, tpop:=sum(pop)]
tmp2[, ppop.all.cities:=pop/tpop]
setkey(tmp2,age.label)
tmp2[,age.idx:=1:nrow(tmp2)]

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
	pa <- merge(pa,subset(tmp2,select=c(age.idx,ppop.all.cities)),by='age.idx')
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

tmp <- subset(ddates, select = c(loc_label, w2.start))
# tmp <- subset(dhr, select=c(loc_label,ref.week))
pads <- merge(pads,tmp,by='loc_label')
tmp <- pads[week.start <= w2.start, list(min.week.idx=week.idx[which.min(M)], min.hfr=min(M)), by='loc_label']
pads <- merge(pads,tmp,by='loc_label')
dhr <- merge(dhr,tmp,by='loc_label')
setkey(dhr, min.hfr)
dhr[, loc_label_group := ceiling((1:nrow(dhr))/5)]
set(dhr,NULL,'loc_label_group',dhr[,factor(loc_label_group,levels=unique(loc_label_group),labels=c('lower tertile','mid tertile','higher tertile'))])
pads <- merge(pads,subset(dhr,select=c(loc_label,loc_label_group)),by='loc_label')
pads <- merge(pads,subset(dstates,select=c(loc_label,region_label)),by='loc_label')
set(pads,pads[,which(region_label%in%c('Central-West','North'))],'region_label','North + Central-West')
pads5A <- copy(pads)
p5A <- ggplot(pads5A, aes(x=week.start)) +
		geom_hline(data=subset(pads5A,week.idx==min.week.idx), aes(yintercept=M,colour=change_city_label(loc_label)), lty='11') +
		geom_ribbon(aes(ymin=CL, ymax=CU, fill=change_city_label(loc_label)),alpha=.3) +		
		geom_line(aes(y=M,colour=change_city_label(loc_label))) +
		scale_colour_manual(values = args$city.palette(length(unique(pads5A$loc_label)))) +
		scale_fill_manual(values = args$city.palette(length(unique(pads5A$loc_label)))) +
		scale_y_continuous(labels=scales::percent, expand=c(0,0), limits=c(0,NA)) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +    
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='bottom',
				strip.text.x = element_text(size = 12)) +
		facet_grid(~region_label) +
		labs(x='',y='weekly age-standardised in-hospital fatality rate', color='', fill='') + 
		guides(fill='none',colour=guide_legend(nrow=2,byrow=TRUE))


loc_label_min <- dhr[min.hfr == min(min.hfr), loc_label]
j <- which(names(pos) == loc_label_min)
cat('\nprocessing',j,' ',names(pos)[j])
pa_ref <- as.data.table(reshape2::melt(pos[[j]]$hfr_overall))
setnames(pa_ref, 1:4, c('iterations','week.idx','age.idx','hfr.overall_ref'))
tmp <- unique(subset(pa_ref, is.na(hfr.overall_ref), select=iterations))[,iterations]
if(length(tmp))
{
  cat('\nFound iterations with NA, removing:', tmp)
  pa_ref <- subset(pa_ref, !iterations%in%tmp)
}
pa_ref[,loc_label := loc_label_min]
pa_ref <- merge(pa_ref, dhr[,.(loc_label, min.week.idx)], by.x=c('loc_label','week.idx'), by.y=c('loc_label','min.week.idx'))
tmp <- merge(pa_ref, tmp2[, .(age.idx,ppop.all.cities)], by = 'age.idx')
tmp <- tmp[,list(loc_label =loc_label_min, age.idx =0,hfr.overall_ref =sum(ppop.all.cities*hfr.overall_ref)),by = 'iterations']
pa_ref <-  rbind(pa_ref[,-'week.idx'], tmp)

pads <- list()
pads2 <- list()
idx <- (1:length(pos))[-j]
for(i in idx)
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
  pa[,loc_label := names(pos)[i]]
  pa <- merge(pa, dhr[,.(loc_label, min.week.idx)], by.x=c('loc_label','week.idx'), by.y=c('loc_label','min.week.idx'))
  tmp <- merge(pa, tmp2[, .(age.idx,ppop.all.cities)], by = 'age.idx')
  tmp <- tmp[,list(loc_label =loc_label_min, age.idx =0,hfr.overall =sum(ppop.all.cities*hfr.overall)),by = 'iterations']
  pa <-  rbind(pa[,-'week.idx'], tmp)
  
  pa <- merge(pa, pa_ref[, .(iterations, age.idx, hfr.overall_ref)], by = c('iterations', 'age.idx'))
  pa[, baseline_ratio := hfr.overall/hfr.overall_ref]
  pads2[[i]] <- pa[age.idx == 0]
  
  pa <- pa[, list(baseline_ratio=quantile(baseline_ratio, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label','age.idx')]
  pa <- dcast(pa, loc_label + age.idx ~ stat, value.var = 'baseline_ratio')
  pa[, loc_label := names(pos)[i]]	
  pads[[i]] <- pa	
}
pads <- do.call('rbind',pads)
pads2 <- do.call('rbind',pads2)
pads <- merge(pads,  unique(dh[, .(age.idx, age.label)]), by = 'age.idx')
pads5B <- copy(pads)
# figure 5B ratio in fitted HFR by location in ref.week
#tmp <- subset(pads4B, variable=='baseline_hfr_ratio' & CU-CL<2 & loc_label!=loc_label.min.ref.hfr)

pads5B[, analysis:='location effect']
pads5B <- pads5B[CU-CL <= 3]
p5B <- ggplot(pads5B, aes(x=age.label)) +
		geom_hline(yintercept=1, lwd=1, colour='grey50') +		
		geom_linerange(aes(ymin=CL, ymax=CU, colour=change_city_label(loc_label)), position=position_dodge(0.5)) +
		geom_point(aes(y=M, colour=change_city_label(loc_label)), position=position_dodge(0.5)) +
		geom_boxplot(aes(y=M), width=0.5, outlier.shape = NA, fill='transparent') +
		scale_colour_manual(values = args$city.palette(14)[-1]) +
		scale_y_continuous(expand=c(0,0), breaks=1:100) +
		facet_wrap(~analysis,ncol=1) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				strip.background = element_blank(),
				legend.position='none',
				strip.text.x = element_text(size = 13)) +
		# guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='',y=paste0('ratio of baseline in-hospital fatality rates\nacross locations relative to ',change_city_label(loc_label_min)),colour='')

# XX for folds in section
folds <- list()
tmp1 <- pads2[, quantile(baseline_ratio, p=c(.025,.5,.975) )]
tmp1 <- round(tmp1, 2)
tmp1 <- paste0(tmp1[2],' (',tmp1[1],'-',tmp1[3],')')
folds[['baseline_ratio_overall']] <- tmp1


# figure 5C ratio in fitted HFR by variant
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
tmp <- unique(subset(dd, select=c(age.idx, age.label)))
pads <- merge(pads, tmp, by='age.idx')
pads <- dcast.data.table(pads, loc_label+age.idx+age.label+variable~stat, value.var='value')
tmp <- dhs[, list(select.by.sample.size=prod(select.by.sample.size)), by=c('loc_label','age.label')]
pads <- merge(pads, tmp, by=c('loc_label','age.label'))
pads5C <- subset(pads, variable=='HFR_ratio' & CU-CL<3)
pads5C[, analysis:='variant effect']
p5C <- ggplot(pads5C, aes(x=age.label)) +
		geom_hline(yintercept=1, lwd=1, colour='grey50') +		
		geom_linerange(aes(ymin=CL, ymax=CU, colour=change_city_label(loc_label)), position=position_dodge(0.5)) +
		geom_point(aes(y=M, colour=change_city_label(loc_label)), position=position_dodge(0.5)) +
		geom_boxplot(aes(y=M), width=0.5, outlier.shape = NA, fill='transparent') +
		scale_colour_manual(values = args$city.palette(length(unique(pads5C$loc_label)))) +
		scale_y_continuous(expand=c(0,0), breaks=1:100, lim=c(min(pads5B$CL),max(pads5B$CU))) +
		facet_wrap(~analysis,ncol=1) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				strip.background = element_blank(),
				legend.position='none',
				strip.text.x = element_text(size = 13)) +
		# guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='',y='ratio of Gamma in-hospital fatality rate\nrelative to non-Gamma in-hospital fatality rate',colour='')


#	figure 5D demand vs resource effect
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
	
	tmp <- subset(dhr,loc_label==names(pos)[i])[, min.week.idx]
	tmp <- subset(pa, week.idx==tmp, select=-c(week,week.idx))
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
pads5D <- copy(pads)
pads5D[, analysis:='healthcare pressure effect']
p5D <- ggplot(pads5D, aes(x=week.start)) +
		geom_hline(yintercept=1,colour='black') +
		geom_linerange(aes(ymin=CL, ymax=CU,colour=change_city_label(loc_label)),alpha=0.5, position=position_dodge(0.9)) +
		geom_point(aes(y=M,colour=change_city_label(loc_label)), position=position_dodge(0.9)) +
		scale_colour_manual(values = args$city.palette(length(unique(pads5D$loc_label)))) +
		scale_x_date(breaks="2 months",expand=c(0,0)) +
		scale_y_continuous(expand=c(0,0),lim=c(min(pads5B$CL),max(pads5B$CU))) +		
		facet_wrap(~analysis,ncol=1) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				strip.background = element_blank(),
				legend.position='none',
				strip.text.x = element_text(size = 13)) +
		# guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
		labs(x='', y='weekly multiplier to in-hospital fatality rates',colour='')


p5BD <- ggpubr::ggarrange( p5B,
		p5C,
		p5D,
		# legend="bottom",
		# common.legend=TRUE,
		labels=c('B','C','D'),
		font.label=list(size=12),
		ncol=3,
		align='h',
		heights=c(3,3,3),
		widths=c(3,3,3))

p5 <- ggpubr::ggarrange( p5A,
		p5BD,		
		legend="bottom",
		common.legend=TRUE,
		labels=c('A',''),
		font.label=list(size=12),
		ncol=1,
		heights=c(4,4),
		widths=c(3,3))
ggsave(file=paste0(args$out.base,'_allloc_figure5.pdf'),p5,w=12,h=11)
ggsave(file=paste0(args$out.base,'_allloc_figure5.png'),p5,w=12,h=11)




# median HFR ratio among all locs in 40-49(4) group
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
  # pa <- pa[age.idx == 4]
  pa[, loc_label := names(pos)[i]]
  pads[[i]] <- pa
}
pads <- do.call('rbind',pads)
tmp <- pads[, list(value = quantile(HFR_ratio, p=c(.025, .5, .975)), stat = c('2.5', '50', '97.5')), by = age.idx]
tmp[, value := round(value,2)]
tmp <- merge(tmp, unique(dd[,.(age.idx, age.label)]), 'age.idx')
tmp <- tmp[age.label == '40-49', paste0(value[2],' (',value[1],'-',value[3],')')]

folds[['hfr_ratio_4049']] <- tmp 


# peak times pandemic load
# get age-standardised hfr_multiplier for all locations
# quantile by week aggregating over all cities
# find week where median is highest?
#	figure 4D demand vs resource effect
tmp2 <- dpop[,list(pop=sum(pop)),by='age.label']
tmp2[, tpop:=sum(pop)]
tmp2[, ppop.all.cities:=pop/tpop]
setkey(tmp2,age.label)
tmp2[,age.idx:=1:nrow(tmp2)]

tmp3 <- pads5D[,{z <- which.max(M); list(M=M[z], week=week[z])},by = 'loc_label']

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
  tmp <- subset(dhr,loc_label==names(pos)[i])[, min.week.idx]
  tmp <- subset(pa, week.idx==tmp, select=-c(week,week.idx))
  setnames(tmp,'hfr.wildtype','hfr.wildtype.min')
  pa <- merge(pa,tmp,by=c('age.idx','iterations'))
  pa[, hfr.multiplier := hfr.wildtype/hfr.wildtype.min]
  pa <- merge(pa,subset(tmp2,select=c(age.idx,ppop.all.cities)),by='age.idx')
  pa <- pa[,list(hfr.multiplier=sum(hfr.multiplier*ppop.all.cities)),by=c('iterations','week')]	
  
  pa[, loc_label := names(pos)[i]]
  pads[[i]] <- pa
  rm(pa, tmp)
}
pads <- do.call('rbind',pads)
tmp <- pads5D[,{z <- which.max(M); list(week.idx = week.idx[z], M = max(M))},by = 'loc_label']
tmp <- merge(tmp, unique(dh[, .(week, week.idx, loc_label)]), by = c("loc_label", "week.idx"))
# mean(tmp$M) 
tmp <- merge(tmp, pads, by = c("loc_label", 'week') )
tmp <- tmp[, list(value = quantile(hfr.multiplier, p=c(.025,.5,.975)))]
tmp <- round(tmp$value, 2)

folds[['hfr_multiplier_overall']] <- sprintf("%2$1.2f (%1$1.2f-%3$1.2f)", tmp[1], tmp[2], tmp[3])
saveRDS(folds, paste0(args$out.base,'_allloc_XXfolds.rds'))
# readRDS(paste0(args$out.base,'_allloc_XXfolds.rds')) -> folds

