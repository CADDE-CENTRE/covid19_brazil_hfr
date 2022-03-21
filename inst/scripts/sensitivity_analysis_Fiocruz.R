######################################################################################
# Define the directories containing the results for Fiocruz and GISAID model running 
# -> Can be changed later when final outputs from the HPC have been run and are ready
######################################################################################
results_1_dir <- '~/Documents/P1Brazil/hfr-analyses-fiocruz'
results_2_dir <- '~/Documents/P1Brazil/hfr-analyses'
results_string = c('Rede Genomica Fiocruz', 'GISAID')
out.base2 <- file.path(results_2_dir, basename(results_2_dir)) 

# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) == 0)
{
  stop('Usage: Rscript sensitivity_analysis_Fiocruz.R -dirs foo bar
       Preferable that exactly one of foo and bar contain the string "fiocruz" ')
} 
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-dirs')	
  tmp <- c(args_line[[2]],args_line[[3]])
  if(sum(grepl('fiocruz', tmp)) != 1)
  {
    warning("No provided directories contain 'fiocruz' in their name\n
            Run without arguments to see instructions")
    results_1_dir <- tmp[1]; results_2_dir <- tmp[2]
  }else{
    results_1_dir <- grep('fiocruz', tmp, value = TRUE)
    results_2_dir <- tmp[!grepl('fiocruz', tmp)]
  }
  
} 

# Loading required libraries 
library(data.table); library(ggplot2); library(ggsci); library(viridis); library(gtools)
library(MGLM); library(fitdistrplus); library(actuar); library(cmdstanr); library(posterior)
library(bayesplot); library(ggpubr); require(ggrepel)

####################################################################
# search for available outputs
#################################################################### 

do <- data.table(FW=list.files(results_1_dir, pattern='_workspace.rda'))
do[, PREFIX := do[, sapply(strsplit(FW, split='_'),'[[',1)]]
do[, loc_label2 := do[, sapply(strsplit(FW, split='_'),'[[',2)]]
tmp <- data.table(FF=list.files(results_1_dir, pattern='_fit.rds'))
tmp[, PREFIX := tmp[, sapply(strsplit(FF, split='_'),'[[',1)]]
tmp[, loc_label2 := tmp[, sapply(strsplit(FF, split='_'),'[[',2)]]
do <- merge(do,tmp,by=c('PREFIX','loc_label2'))
tmp <- data.table(loc_label=c('belo horizonte','curitiba', "florianopolis", "goiania", "joao pessoa", "macapa", 'manaus', "natal", 'porto alegre', "porto velho", 'rio de janeiro', 'salvador', 'sao paulo', "sao luis"))
tmp[, loc_label2 := gsub(' ','',loc_label)]
do <- merge(do, tmp, by='loc_label2')
do1 <- do

do <- data.table(FW=list.files(results_2_dir, pattern='_workspace.rda'))
do[, PREFIX := do[, sapply(strsplit(FW, split='_'),'[[',1)]]
do[, loc_label2 := do[, sapply(strsplit(FW, split='_'),'[[',2)]]
tmp <- data.table(FF=list.files(results_2_dir, pattern='_fit.rds'))
tmp[, PREFIX := tmp[, sapply(strsplit(FF, split='_'),'[[',1)]]
tmp[, loc_label2 := tmp[, sapply(strsplit(FF, split='_'),'[[',2)]]
do <- merge(do,tmp,by=c('PREFIX','loc_label2'))
tmp <- data.table(loc_label=c('belo horizonte','curitiba', "florianopolis", "goiania", "joao pessoa", "macapa", 'manaus', "natal", 'porto alegre', "porto velho", 'rio de janeiro', 'salvador', 'sao paulo', "sao luis"))
tmp[, loc_label2 := gsub(' ','',loc_label)]
do <- merge(do, tmp, by='loc_label2')
do2 <- do

####################################################################
# load stan_data
####################################################################
do <- do1
stan_data <- list()
for(i in 1:nrow(do))
{	
  tmp <- file.path(results_1_dir,do$FW[i])
  cat('\nreading ',tmp)
  load(tmp,  tmp_env <- new.env())
  tmp_list <- as.list.environment(tmp_env)	
  stan_data[[i]] <- tmp_list$stan_data
  rm(tmp_list)
}
names(stan_data) <- do[, loc_label]
stan_data1 <- stan_data

do <- do2
stan_data <- list()
for(i in 1:nrow(do))
{	
  tmp <- file.path(results_2_dir,do$FW[i])
  cat('\nreading ',tmp)
  load(tmp,  tmp_env <- new.env())
  tmp_list <- as.list.environment(tmp_env)	
  stan_data[[i]] <- tmp_list$stan_data
  rm(tmp_list)
}
names(stan_data) <- do[, loc_label]
stan_data2 <- stan_data

####################################################################
# load workspace
####################################################################
# only one workspace needed
tmp <- paste0( out.base2,'_allloc_workspace.rda')
cat('\nAttempting to load ', tmp)
load( tmp )
# plot.dir <- out.dir <- file.prefix <- NULL
gc()

# add to args the local arguments and update to new out.dir

# args$city.palette <- grDevices::colorRampPalette(ggsci::pal_npg("nrc")(10))
# args$use.only.hosp.deaths <- use.only.hosp.deaths

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

load_samples <- function(do_tmp, dir2_tmp, stan_data_tmp){
  pos <- list()
  for(i in 1:nrow(do_tmp))
  {	
    tmp <- file.path(dir2_tmp,do_tmp$FF[i])
    cat('\nreading ',tmp)
    m_fit <- readRDS(tmp)
    
    # get HMC samples in wide array format
    pd <- m_fit$draws(inc_warmup = FALSE)	
    
    # bring samples into long array format
    select.chains <- seq_along(dimnames(pd)[['chain']])
    iters <- 1:(m_fit$metadata()[['iter_sampling']])
    age.idx.n <- stan_data_tmp[[i]][['A']]
    age.idx2.n <- stan_data_tmp[[i]][['A_SP']]

    po <- list()
    tmp <- pd[,,which(grepl('log_age_death_rate_wildtype',dimnames(pd)[[3]]))]
    po$log_age_death_rate_wildtype <- unname(apply(tmp[,select.chains,], 3, rbind))
    tmp <- pd[,,which(grepl('prop_P1',dimnames(pd)[[3]]) & !grepl('logit_prop_P1',dimnames(pd)[[3]]) & !grepl('age_prop_P1',dimnames(pd)[[3]])  & !grepl('prop_P1_fit',dimnames(pd)[[3]]))]
    po$prop_P1 <- unname(apply(tmp[,select.chains,], 3, rbind))
    tmp <- pd[,,which(grepl('logit_hosp_fatality_ratio_wildtype\\[',dimnames(pd)[[3]]))]
    po$logit_hosp_fatality_ratio_wildtype <- unname(apply(tmp[,select.chains,], 3, rbind))
    tmp <- pd[,,which(grepl('logit_hosp_fatality_ratio_P1_rnde\\[',dimnames(pd)[[3]]))]
    po$logit_hosp_fatality_ratio_P1_rnde <- unname(apply(tmp[,select.chains,], 3, rbind))
    tmp <- pd[,,which(grepl('^fr_multiplier_by_week\\[',dimnames(pd)[[3]]))]
    tmp <- apply(tmp[,select.chains,], 3, rbind)
    po$fr_multiplier_by_week <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
    tmp <- pd[,,which(grepl('hfr_wildtype',dimnames(pd)[[3]]))]
    tmp <- apply(tmp[,select.chains,], 3, rbind)
    po$hfr_wildtype <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
    tmp <- pd[,,which(grepl('hfr_overall',dimnames(pd)[[3]]))]
    tmp <- apply(tmp[,select.chains,], 3, rbind)
    po$hfr_overall <- array( tmp, dim=c(dim(tmp)[1], dim(tmp)[2]/age.idx.n, age.idx.n))
    tmp <- pd[,,which(grepl('fr_regcoeff_overall',dimnames(pd)[[3]]))]
    po$fr_regcoeff_overall <- unname(apply(tmp[,select.chains,], 3, rbind))
    
    pos[[do_tmp$loc_label[i]]] <- po
    gc()
    rm(po, pd, m_fit)
  }
  return(pos)
}

# memory.limit(size = 20000)
pos1 <- load_samples(do1, results_1_dir, stan_data1)
pos2 <- load_samples(do2, results_2_dir, stan_data2)

####################################################################
# figures
####################################################################
# counts of P1 vs non-P.1 genomes per week/month
tmp <- subset(dp, select=c(loc_label, n_positive, n_sampled, week, week.start))
tmp[, n_negative := n_sampled-n_positive]
tmp <- melt(tmp, id.vars=c('loc_label','week','week.start','n_sampled'))
set(tmp, NULL, 'variable', tmp[,factor(variable,levels=c('n_negative','n_positive'),labels=c('non-Gamma','Gamma'))])
p <- ggplot(tmp, aes(x=week.start)) +
  geom_bar(aes(y=value, fill=variable), stat='identity',position='stack') +
  ggsci::scale_fill_aaas() +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_date(breaks='1 month') +
  facet_wrap(~change_city_label(loc_label), ncol=5, scales='free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        legend.position = "bottom") +
  labs(x='', y='Number of weekly SARS-CoV-2 sequences', fill='') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title = "Lineage Assignment"))

ggsave(file=paste0(out.base2, '_fiocruz_genome_counts.pdf'),p,w=10,h=6)
ggsave(file=paste0(out.base2, '_fiocruz_genome_counts.png'),p,w=10,h=6)

# figure S.6.2.3 - P1 replacement dynamics
combine_pos <- function(p1, p2){
  p1$scenario <- results_string[1]
  p2$scenario <- results_string[2]
  p <- rbind(p1, p2)
  return(p)
}
pads_5D_processing <- function(pos){
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
  return(pads)
}

p5D1 <- pads_5D_processing(pos1)
p5D2 <- pads_5D_processing(pos2)
pads5D <- combine_pos(p5D1, p5D2)
pads5D$loc_label <- stringr::str_to_title(pads5D$loc_label)
tmp <- subset(pads5D, !is.na(n_sampled) & n_positive>0 & n_positive!=n_sampled)
p5D <- ggplot(pads5D, aes(x=week.start)) +
  geom_ribbon(aes(ymin=CL, ymax=CU, fill = scenario), alpha = 0.5) +
  geom_line(aes(y=M, colour = scenario), size = 1.25) +
  #geom_point(data=subset(pads5D, !is.na(n_sampled)), aes(y=n_positive/n_sampled, size=n_sampled, fill=scenario), pch=21) +		
  scale_fill_manual(values = c("#FFC759", "#9984D4")) +
  scale_colour_manual(values = c("#FFC759", "#9984D4")) +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  scale_x_date(breaks='2 months',expand=c(0,0)) +
  facet_wrap(. ~ loc_label, ncol=5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7), 
        strip.background = element_blank(),
        legend.position='bottom') +
  guides(colour=guide_legend(title="Source of SARS-CoV-2\nSequence Data"), fill = "none") +
  labs(y='SARS-CoV-2 Gamma Variant Frequency', x='', colour='', size='number of SARS-CoV-2 sequences')

ggsave(file=paste0(out.base2, '_gamma_frequency.pdf'),p5D,w=12,h=6)
ggsave(file=paste0(out.base2, '_gamma_frequency.png'),p5D,w=12,h=6)

# figure S.6.2.1 - ratio in fitted HFR by location in ref.week
# calculate baseline HFR, using inferred regression coefficients 
pads_5B_processing <- function(pos){
  
  # Compute age standardised HFRs and find location with min in ref,week (copied form make.main.figures)
  dhr_temp <- copy(dhr)
  
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
  pads <- merge(pads,tmp,by='loc_label')
  tmp <- pads[week.start <= w2.start, list(min.week.idx=week.idx[which.min(M)], min.hfr=min(M)), by='loc_label']
  pads <- merge(pads,tmp,by='loc_label')
  dhr_temp <- merge(dhr_temp,tmp,by='loc_label')
  setkey(dhr_temp, min.hfr)
  # dhr_temp[, loc_label_group := ceiling((1:nrow(dhr_temp))/5)]
  # set(dhr_temp,NULL,'loc_label_group',dhr_temp[,factor(loc_label_group,levels=unique(loc_label_group),labels=c('lower tertile','mid tertile','higher tertile'))])
  # pads <- merge(pads,subset(dhr_temp,select=c(loc_label,loc_label_group)),by='loc_label')
  # pads <- merge(pads,subset(dstates,select=c(loc_label,region_label)),by='loc_label')
  # set(pads,pads[,which(region_label%in%c('Central-West','North'))],'region_label','North + Central-West')
  # pads5A <- copy(pads)
  
  # Get posterior draws from ref week of ref location 
  loc_label_min <- dhr_temp[min.hfr == min(min.hfr), loc_label]
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
  pa_ref <- merge(pa_ref, dhr_temp[,.(loc_label, min.week.idx)], by.x=c('loc_label','week.idx'), by.y=c('loc_label','min.week.idx'))
  tmp <- merge(pa_ref, tmp2[, .(age.idx,ppop.all.cities)], by = 'age.idx')
  tmp <- tmp[,list(loc_label =loc_label_min, age.idx =0,hfr.overall_ref =sum(ppop.all.cities*hfr.overall_ref)),by = 'iterations']
  pa_ref <-  rbind(pa_ref[,-'week.idx'], tmp)
  
  # Go through the others (do not really need pads2 here, or do we?)
  pads <- list()
  # pads2 <- list()
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
    pa <- merge(pa, dhr_temp[,.(loc_label, min.week.idx)], by.x=c('loc_label','week.idx'), by.y=c('loc_label','min.week.idx'))
    tmp <- merge(pa, tmp2[, .(age.idx,ppop.all.cities)], by = 'age.idx')
    # tmp <- tmp[,list(loc_label =loc_label_min, age.idx =0,hfr.overall =sum(ppop.all.cities*hfr.overall)),by = 'iterations']
    # pa <-  rbind(pa[,-'week.idx'], tmp)
    
    pa <- merge(pa, pa_ref[, .(iterations, age.idx, hfr.overall_ref)], by = c('iterations', 'age.idx'))
    pa[, baseline_ratio := hfr.overall/hfr.overall_ref]
    # pads2[[i]] <- pa[age.idx == 0]
    
    pa <- pa[, list(baseline_ratio=quantile(baseline_ratio, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU') ), by=c('loc_label','age.idx')]
    pa <- dcast(pa, loc_label + age.idx ~ stat, value.var = 'baseline_ratio')
    pa[, loc_label := names(pos)[i]]	
    pads[[i]] <- pa	
  }
  pads <- do.call('rbind',pads)
  # pads2 <- do.call('rbind',pads2)
  pads <- merge(pads,  unique(dh[, .(age.idx, age.label)]), by = 'age.idx')
  pads5B <- copy(pads)
  #pads5B <- subset(pads5B, variable=='baseline_hfr_ratio' & loc_label!=loc_label.min.ref.hfr)
  pads5B[, analysis:='location effect']
  return(pads5B)
  
}

p5B1 <- pads_5B_processing(pos1)
p5B2 <- pads_5B_processing(pos2)
pads5B <- combine_pos(p5B1, p5B2)
pads5B$loc_label <- stringr::str_to_title(pads5B$loc_label)
p5B <- ggplot(pads5B, aes(x=age.label)) +
  geom_hline(yintercept=1, lwd=1, colour='grey50') +		
  geom_linerange(aes(ymin=CL, ymax=CU, colour=scenario), position=position_dodge(0.5)) +
  geom_point(aes(y=M, colour=scenario), position=position_dodge(0.5)) +
  scale_colour_manual(values = c("#FFC759", "#9984D4")) +
  facet_wrap(. ~ loc_label, ncol=7, scale = "free_y") +
  scale_y_continuous(expand=c(0,0), breaks=0:100, limits = c(0, NA)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7),
        strip.background = element_blank(),
        legend.position='bottom') +
  guides(colour=guide_legend(nrow=1,byrow=TRUE, title="Source of SARS-CoV-2\nSequence Data")) +
  labs(x='',y=paste0('Ratio of in-hospital fatality rates across locations\nrelative to Belo Horizonte'), colour='')
ggsave(file=paste0(out.base2, '_location_effect_sizes.pdf'),p5B,w=12,h=8)
ggsave(file=paste0(out.base2, '_location_effect_sizes.png'),p5B,w=12,h=8)

# figure S.6.2.2 - ratio in of P1 to non-P1 HFR by age-group and location
pads_5C_processing <- function(pos){
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
  pads4AC <- merge(pads, tmp, by=c('loc_label','age.label'))
  
  
  pads5C <- subset(pads4AC, variable=='HFR_ratio' & CU-CL<2)
  pads5C[, analysis:='variant effect']
  
  return(pads5C)
}

p5C1 <- pads_5C_processing(pos1)
p5C2 <- pads_5C_processing(pos2)
pads5C <- combine_pos(p5C1, p5C2)
pads5C$loc_label <- stringr::str_to_title(pads5C$loc_label)
p5C <- ggplot(pads5C, aes(x=age.label)) +
  geom_hline(yintercept=1, lwd=1, colour='grey50') +		
  geom_linerange(aes(ymin=CL, ymax=CU, colour=scenario), position=position_dodge(0.5)) +
  geom_point(aes(y=M, colour=scenario), position=position_dodge(0.5)) +
  scale_colour_manual(values = c("#FFC759", "#9984D4")) +
  scale_y_continuous(expand=c(0,0), breaks=c(0, 1, 2), lim=c(0,2.7)) +
  facet_wrap(. ~ loc_label, ncol=5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7),
        strip.background = element_blank(),
        legend.position='bottom') +
  guides(colour=guide_legend(nrow=1,byrow=TRUE, title="Source of SARS-CoV-2\nSequence Data")) +
  labs(x='',y='Ratio of Gamma in-hospital fatality rate\nrelative to non-Gamma in-hospital fatality rate',colour='')

ggsave(file=paste0(out.base2, '_P1_effect_sizes.pdf'),p5C,w=12,h=8)
ggsave(file=paste0(out.base2, '_P1_effect_sizes.png'),p5C,w=12,h=8)

# figure S.6.2.3 - health pressure effects
tmp2 <- dpop[,list(pop=sum(pop)),by='age.label']
tmp2[, tpop:=sum(pop)]
tmp2[, ppop.all.cities:=pop/tpop]
tmp2[,age.idx:=1:nrow(tmp2)]

# calculate weeks with minimum HFR for dataset 1
pads1 <- list()
pos <- pos1
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
  pads1[[i]] <- pa
  rm(pa, tmp)
}
pads1 <- do.call('rbind',pads1)
tmp <- unique(subset(dh, select=c(loc_label, week.idx, week, week.start)))
pads1 <- merge(pads1,tmp,by=c('loc_label','week.idx'))
tmp <- subset(dhr, select=c(loc_label,ref.week))
pads1 <- merge(pads1,tmp,by='loc_label')
tmp <- pads1[week <= ref.week, list(min.week=week[which.min(M)], min.hfr=min(M)), by='loc_label']
pads1 <- merge(pads1,tmp,by='loc_label')

dhr <- merge(dhr, tmp, by = "loc_label") 
dhr <- dhr %>%
  dplyr::rename(min_weeks_1 = min.week, min_hfr_1 = min.hfr)

# calculate weeks with minimum HFR for dataset 2
pads2 <- list()
pos <- pos2
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
  pads2[[i]] <- pa
  rm(pa, tmp)
}
pads2 <- do.call('rbind',pads2)
tmp <- unique(subset(dh, select=c(loc_label, week.idx, week, week.start)))
pads2 <- merge(pads2,tmp,by=c('loc_label','week.idx'))
tmp <- subset(dhr, select=c(loc_label,ref.week))
pads2 <- merge(pads2,tmp,by='loc_label')
tmp <- pads2[week <= ref.week, list(min.week=week[which.min(M)], min.hfr=min(M)), by='loc_label']
pads2 <- merge(pads2,tmp,by='loc_label')

dhr <- merge(dhr, tmp, by = "loc_label") 
dhr <- dhr %>%
  dplyr::rename(min_weeks_2 = min.week, min_hfr_2 = min.hfr)

# Plotting
pads_5A_processing <- function(pos, dataset) {
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
    if (dataset == 1) {
      tmp <- subset(dhr,loc_label==names(pos)[i])[, min_weeks_1]
    } else if (dataset == 2) {
      tmp <- subset(dhr,loc_label==names(pos)[i])[, min_weeks_2]
    }
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
  pads5A <- copy(pads)
  pads5A[, analysis:='demand vs. resource effect']
}

p5A1 <- pads_5A_processing(pos1, 1)
p5A2 <- pads_5A_processing(pos2, 2)
p5A1$genomic <- "Rede Genomica Fiocruz"
p5A2$genomic <- "GISAID"
p5A_overall <- rbind(p5A1, p5A2)
p5A_overall$loc_label <- stringr::str_to_title(p5A_overall$loc_label)
p5A_overall_plt <- ggplot(p5A_overall, aes(x=week.start)) +
  geom_hline(yintercept=1,colour='black') +
  geom_linerange(aes(ymin=CL, ymax=CU,colour=genomic), position=position_dodge(10)) +
  geom_point(aes(y=M,colour=genomic), position=position_dodge(10)) +
  scale_colour_manual(values = c("#FFC759", "#9984D4")) +
  scale_x_date(breaks="2 months",expand=c(0,0)) +
  lims(y = c(0, 4)) +
  #scale_y_continuous(expand=c(0,0),lim=c(min(pads4B$CL),max(pads4B$CU))) +		
  facet_wrap(~loc_label,ncol=5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.background = element_blank(),
        legend.position='bottom') +
  guides(colour=guide_legend(nrow=1,byrow=TRUE, title="Source of SARS-CoV-2\nSequence Data")) +
  labs(x='', y='Weekly multiplier to in-hospital fatality rates',colour='')

ggsave(file=paste0(out.base2, '_healthpressure_effect_sizes.pdf'),p5A_overall_plt,w=12,h=8)
ggsave(file=paste0(out.base2, '_healthpressure_effect_sizes.png'),p5A_overall_plt,w=12,h=8)


# p4 <- ggpubr::ggarrange( p4B,
#                          p4C,
#                          p4a_overall_plt,
#                          legend="bottom",
#                          common.legend=TRUE,
#                          labels=c('A','B','C'),
#                          font.label=list(size=12),
#                          ncol=1,
#                          align='v',
#                          heights=c(9,9,9),
#                          widths=c(3,3,3))
# 
# ggsave(file=paste0(out.base2, '_all_effect_sizes.pdf'),p4,w=9,h=13)
# ggsave(file=paste0(out.base2, '_all_effect_sizes.png'),p4,w=9,h=13)

