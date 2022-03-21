cat(" \n -------------------------------- \n \n HFR.postprocessing.tables.R \n \n -------------------------------- \n")

###########################################################################
# Get the RDS files and produce required tables for main text
###########################################################################

require(reshape2)
require(data.table)
require(xtable)
require(stringr)
require(tools)
require(lubridate)


if(1) # # Specify arguments "manual" run 
{
  out.dir3 <- '~/Documents/P1Brazil/hfr-fit-210719h4-resources8-ex0-trvaluecs-hd1/'
  file.prefix3 <- basename(out.dir3)
  threshold_check <-  0
  pkg.dir3 <-  '~/git/covid19_brazil_hfr/'
}

# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-outdir')	
  stopifnot(args_line[[3]]=='-file_prefix')
  stopifnot(args_line[[5]]=='-pkgdir')
  stopifnot(args_line[[7]]=='-threshold_check')
  
  out.dir3 <- args_line[[2]]
  file.prefix3 <- args_line[[4]]
  pkg.dir3 <- args_line[[6]]
  threshold_check <- as.integer(args_line[[8]])
} 

out.base3 <- file.path(out.dir3, file.prefix3)
############################################################
# Load workspace and helper functions
############################################################


tmp <- list.files(file.path(pkg.dir3, 'inst','R'), full.names=TRUE)
cat('sourcing code from files\n', tmp)
invisible(lapply(tmp, source))

write.date <- function(date){return( as.character(format(date, "%d %B %Y")))}

# read workspace
tmp <- list.files(out.dir3, pattern = 'allloc_workspace.rda', full.names = TRUE)[1]
load(tmp)
dha <- merge(dha,subset(ddates, select = c('loc_label','w1.start', 'w2.end')), by = 'loc_label')
dht <- merge(dht,subset(ddates, select = c('loc_label','w1.start', 'w2.end')), by = 'loc_label')
pkg.dir <- pkg.dir3
out.dir <- out.dir3
out.base <- out.base3


####################################################################
# Check that at least threshold number of chains are done running
####################################################################

tmp <- list.files(out.dir, pattern='_fit.rds')
cat('\ntry to bugfix: start\n')
print(out.dir)
print(tmp)
cat('\ntry to bugfix: end\n')
if(length(tmp) < threshold_check)
{
  stop("Not enough runs have completed yet!")
}


###############################################################################
# Settings
###############################################################################

# check if we count CLASSI_FIN = nA or not
tmp <- dht[ loc.res.cat=='resident' & DateAdmission >= w1.start & DateAdmission <= w2.end  & hosp.cat=='in hospital', sum(whosps)  ]
tmp_noNA <- dha[ DateAdmission >= w1.start & DateAdmission <= w2.end & hosp.cat=='in hospital' & loc.res.cat=='resident' &
                   adm.cat %in% c('COVID19 susp/conf')   , .N   ]
tmp_withNA <- dha[ DateAdmission >= w1.start & DateAdmission <= w2.end & hosp.cat=='in hospital' & loc.res.cat=='resident' &
                     adm.cat %in% c('COVID19 susp/conf', 'unknown')   , .N   ]

cat('check that Flowchart will be consistent')
stopifnot(tmp %in% c(tmp_noNA, tmp_withNA))
adm.classes <- c('COVID19 susp/conf')
if(tmp == tmp_withNA){adm.classes <- c('COVID19 susp/conf', 'unknown')}

# want P1 effect?
P1.eff <- TRUE
if(P1.eff)
{
  # look for rds files
  files <- list.files(out.dir, pattern = 'allloc(.*?)table.rds', full.names = TRUE)
  # tmp <- readRDS(files[1])
}else{
  # look for rds files
  files <- list.files(out.dir, pattern = 'wildtype.rds', full.names = TRUE)
}

################################################
##        Table 1 main text  
################################################
cat('\n\n==== Table 1 ====\n\n')

tmp1 <- readRDS(files[1])[, -c(5,6)]
tmp2 <- as.data.table(readRDS(files[2]))[,c(1,2)]
tmp3 <- as.data.table(readRDS(files[3]))[,c(1,2)]
tmp2 <- merge(tmp2, tmp3, by = 'loc_label')

tmp <- merge(tmp1, tmp2, by = 'loc_label')
tmp2 <- readRDS(paste0(out.base,'_minhfr_methods_table.rds'))
tmp2 <- tmp2[,.(loc_label,hfr_loess.min, hfr_loess.max, hfr_loess.fold)]
tmp2[,`:=`(loc_label=change_city_label(loc_label),
           hfr_loess.min=paste0(format(round(hfr_loess.min*100,1), nsmall=1),'%'),
           hfr_loess.max=paste0(format(round(hfr_loess.max*100,1), nsmall=1),'%'),
           hfr_loess.fold=as.character(format(round(hfr_loess.fold,2),nsmall=1),'%')
)]
tmp <- merge(tmp2, tmp, by='loc_label', all.y = TRUE)

# cumulative hosp admissions & CENSORING ADJUSTED deaths
# among Residents, chosen vax type, and COVID attr
tmp1 <- dht[Outcome2 %in% c('Death_Covid','Death_Covid_Estimated') & hosp.cat=='in hospital' & 
              loc.res.cat == 'resident' & grepl(args$keep.vac.type, vac_type) &
              DateAdmission >= w1.start & DateAdmission <= w2.end, ]

# Transfrom to integer later to avoid incongruences
tmp1 <- tmp1[,list(cens_adj_deaths = sum(wdeaths)),by = c('loc_label')]

tmp2 <- dha[ loc.res.cat == 'resident' & grepl(args$keep.vac.type, vac_type) &
               DateAdmission >= w1.start & DateAdmission <= w2.end &
               adm.cat %in% adm.classes & hosp.cat=='in hospital',  
             list(observed_admissions = .N),by = 'loc_label']

# instert observation period
tmp1 <- merge(tmp1, tmp2, by = 'loc_label')
tmp1[,loc_label:= change_city_label(loc_label)]
tmp1 <- rbind(tmp1, data.table(loc_label='Total', cens_adj_deaths=sum(tmp1$cens_adj_deaths),  observed_admissions=sum(tmp1$observed_admissions)))
tmp <- merge(tmp1, tmp, by = 'loc_label')

tmp2 <- ddates[,.(loc_label, w1.start, w2.end)]
tmp2 <- rbind(tmp2, data.table(loc_label = 'Total', w1.start=min(tmp2$w1.start), w2.end=max(tmp2$w2.end)))
tmp2[, obs_period := paste0(format(w1.start, '%d/%m/%y'),'-', format(w2.end, '%d/%m/%y'))]
tmp2[,`:=` (loc_label=change_city_label(loc_label),w1.start = NULL,w2.end = NULL ) ]
tmp <- merge(tmp2, tmp, by = 'loc_label', all.y = TRUE)

tmp <- rbind(tmp[loc_label=='Total'],tmp[loc_label!='Total'])
tmp[, cens_adj_deaths:= as.integer(cens_adj_deaths)]



# Load data from OTHER counterfactual
tmp2 <- readRDS(paste0(out.base,'_counterfactual_hfr_minimum_across_locations.rds') )
tmp2 <- as.data.table(tmp2$table)
tmp2[, loc_label := change_city_label(loc_label)]
tmp2[loc_label == 'Total', hyp_deaths_of_hosps_pcred := gsub('(\\.[0-9])[0-4]%', '\\1%',hyp_deaths_of_hosps_pcred)]
tmp3 <- tmp2[loc_label == 'Total', 
     { z <- as.numeric( gsub('^([0-9][0-9]\\.[0-9][0-9]).*$', '\\1',hyp_deaths_of_hosps_pcred));
       z <- paste0(as.character(round(z, 1)), '%') }]
tmp2[loc_label == 'Total', hyp_deaths_of_hosps_pcred :=  gsub('^([0-9][0-9]\\.[0-9][0-9])%', tmp3 ,hyp_deaths_of_hosps_pcred)]

tmp <- merge(tmp, tmp2, by = 'loc_label')

f <- function(s){
  if (is.na(s)){return(as.character(NA))}
  out <- strsplit(s, split = '\\(')[[1]]
  out[2] <- paste0('(', out[2])
  return(out)
}


if (P1.eff)
{
  tmp <- tmp[, list(obs_period = rep(obs_period, each=2),
                    observed_admissions = rep(observed_admissions, each=2),
                    cens_adj_deaths = rep(cens_adj_deaths, each=2),
                    hfr_loess.min=rep(hfr_loess.min, each=2),
                    hfr_loess.max=rep(hfr_loess.max, each=2),
                    hfr_loess.fold=rep(hfr_loess.fold, each=2),
                    hfr.overall.min = f(hfr.overall.min),
                    hfr.overall.max = f(hfr.overall.max),
                    hfr.overall.diff= f(hfr.overall.diff),
                    overall.x = f(overall.x),
                    overall.y = f(overall.y),
                    hyp_deaths_of_hosps_diff = f(hyp_deaths_of_hosps_diff),
                    hyp_deaths_of_hosps_pcred = f(hyp_deaths_of_hosps_pcred)), by = 'loc_label']
  tmp[, `:=`(hfr.overall.min=NULL, hfr.overall.max=NULL, hfr.overall.diff=NULL)]
}else{
  tmp <- tmp[, list(obs_period = rep(obs_period, each=2),
                    observed_admissions = rep(observed_admissions, each=2),
                    cens_adj_deaths = rep(cens_adj_deaths, each=2),
                    hfr_loess.min=rep(hfr_loess.min, each=2),
                    hfr_loess.max=rep(hfr_loess.max, each=2),
                    hfr_loess.fold=rep(hfr_loess.fold, each=2),
                    hfr.wildtype.min = f(hfr.wildtype.min),
                    hfr.wildtype.max = f(hfr.wildtype.max),
                    hfr.wildtype.diff= f(hfr.wildtype.diff),
                    overall.x = f(overall.x),
                    overall.y = f(overall.y)), by = 'loc_label']
  tmp[, `:=`(hfr.wildtype.min=NULL, hfr.wildtype.max=NULL, hfr.wildtype.diff=NULL)]
}

tmp1 <- dim(tmp)[1]/2
# tmp$obs_period[2*1:tmp1] <- ''
tmp$observed_admissions[2*1:tmp1] <- ''
tmp$cens_adj_deaths[2*1:tmp1] <- ''
# tmp$loc_label[2*1:tmp1] <- ''
tmp$hfr_loess.min[2*1:tmp1] <- ''
tmp$hfr_loess.max[2*1:tmp1] <- ''
tmp$hfr_loess.fold[2*1:tmp1] <- ''

# Want study period below location?
if(0)
{
  tmp1 <- tmp[, {z <- c(unique(loc_label), unique(obs_period)); list(loc_label2 = z)},by = 'loc_label']
  tmp <- cbind(tmp1[,.(loc_label2)], tmp)[, `:=` (loc_label = NULL, obs_period = NULL)] 
  setnames(tmp, 'loc_label2', 'loc_label')
}else{
  tmp$loc_label[2*1:tmp1] <- ''
  tmp$obs_period[2*1:tmp1] <- ''
}

## Latest adjustments:

# remove white spaces:
f <- function(x){gsub(' ','',x)}
tmp1 <- names(tmp)[! names(tmp) %in% c('loc_label','obs_period')]
tmp[,(tmp1) := lapply(.SD, FUN=f),.SDcols = tmp1]
# write commas
tmp1 <- c('observed_admissions','cens_adj_deaths' )
tmp[, (tmp1):=lapply(.SD, FUN=write_commas) ,.SDcols = tmp1]
f <- function(s){
  s <- str_match(s, '\\((.*?)-(.*?)\\)')[,c(2,3)]
  s <- write_commas(s)
  return(paste0('(', s[1], '-', s[2], ')'))
}
tmp[grepl('\\)', overall.x), overall.x:=sapply(overall.x,FUN = f)]
tmp[!grepl('\\)', overall.x), overall.x := sapply(overall.x, FUN=write_commas)]
tmp[grepl('\\)', hyp_deaths_of_hosps_diff), hyp_deaths_of_hosps_diff:=sapply(hyp_deaths_of_hosps_diff,FUN = f)]
tmp[!grepl('\\)', hyp_deaths_of_hosps_diff), hyp_deaths_of_hosps_diff := sapply(hyp_deaths_of_hosps_diff, FUN=write_commas)]

# change Total to All cities
tmp1 <- nrow(tmp)/2 - 1
tmp[loc_label == 'Total', loc_label:= sprintf('All %d cities', tmp1)]

# remove number columns
tmp[, `:=`(hyp_deaths_of_hosps_diff = NULL, overall.x=NULL)]

if(P1.eff)
{
  saveRDS(tmp, paste0(out.base, '_main_table1.rds'))
}else{
  saveRDS(tmp, paste0(out.base, '_main_table1_noP1.rds'))
}

if(0)
{
  tmp1 <- readRDS(paste0(out.base, '_main_table1.rds'))
  # tmp1[1:2, paste0(overall.x, collapse=' ')]
  # tmp1[1:2, paste0(gsub('%','\\\\%',overall.y), collapse=' ')]
  
  tmp1 <- print(xtable(tmp), include.rownames=FALSE,  include.colnames = FALSE)
  tmp1 <- gsub('\\\\begin\\{table\\}\\[ht\\]\\\n','',tmp1)
  tmp1 <- gsub('\\\\end\\{table\\}\\\n','',tmp1)
  tmp1 <- gsub('\\\\begin\\{tabular\\}\\{l(.*?)l\\}','',tmp1)
  tmp1 <- gsub('\\\\centering','',tmp1)
  tmp1 <- gsub('\\\\end\\{tabular\\}\\\n','',tmp1)
  tmp1 <- sub('\\\\hline\\\n','',tmp1) 
}


#################################################
##  Resource Table, supplementary table 2 
###############################################
cat('\n\n==== Supplementary Tables ====\n\n')

sexpr <- list.files(out.dir, pattern = "substitute_expressions.rds", full.names = TRUE)
sexpr <- readRDS(sexpr)

tmp <- sexpr[['Resources']]
tmp1 <- which(grepl('physicians.specialist', rownames(tmp)))
tmp <- tmp[-tmp1,]

colnames(tmp) <- change_city_label(colnames(tmp))
tmp <- as.data.table(tmp, keep.rownames = TRUE)
tmp1 <- names(tmp)[-1]
tmp[, (tmp1):=lapply(.SD, FUN=as.integer), .SDcols = tmp1]

# Compute percentages
tmp2 <- function(vec){
  l <- length(vec)/2
  out <- as.integer(vec)
  for (i in 1:l){
    out[2*i] <- paste0(vec[2*i], '(', round((vec[2*i] - vec[2*i-1])*100/vec[2*i-1],0),'%)')
  }
  return(out)
}
tmp[, (tmp1) := lapply(.SD, FUN = tmp2), .SDcols = tmp1]

tmp1 <- unique(str_match(tmp$rn, pattern = '_(.*?)$')[,2])
tmp1 <- paste(as.character(month(tmp1, label = TRUE, abbr=F)), year(tmp1))
tmp <- cbind(rep(tmp1, (nrow(tmp)/2)) , tmp)

tmp2 <- tmp$rn
tmp[, rn := NULL]
tmp <- as.matrix(tmp)
rownames(tmp) <- tmp2

rownames(tmp) <- str_match(rownames(tmp), pattern = 'r.(.*?)_')[,2]
rownames(tmp) <- gsub('\\.all', '', rownames(tmp))
rownames(tmp) <- gsub('(.*).specialist', 'specialist.\\1', rownames(tmp))
rownames(tmp) <- gsub('([^s])$', '\\1s', rownames(tmp))
rownames(tmp) <- gsub('\\.', ' ', rownames(tmp))
rownames(tmp) <- gsub('icubeds', 'ICU beds', rownames(tmp))
rownames(tmp) <- gsub('tech ', 'technical ', rownames(tmp))
rownames(tmp) <- toTitleCase(rownames(tmp))
rownames(tmp)[1:(nrow(tmp)/2)*2] <- ' '
rownames(tmp) <- gsub('Technical Nurses', 'Nurse assistants', rownames(tmp))
rownames(tmp) <- gsub('Intensivists', 'Intensive Care', rownames(tmp))
rownames(tmp) <- gsub('Specialist Physicians', 'Critical Care', rownames(tmp))
rownames(tmp)[which(grepl('Critical Care|Intensive Care',rownames(tmp))) + 1] <- 'Specialists'

rownames(tmp) <- gsub('Saribeds', 'CC Beds', rownames(tmp))
rownames(tmp)[which(grepl('Bedsvents',rownames(tmp))) + 1] <- 'ventilator'
rownames(tmp) <- gsub('Bedsvents', 'Beds with', rownames(tmp) )

f <- function(x){gsub('\\(NaN%\\)', '', x)} 
tmp <- apply(tmp, FUN=f, 2)


saveRDS(tmp, paste0(out.base, '_resource_table.rds'))



################################################
##         Supplement tables
################################################
change_age_label <- function(label){
  L <- length(label)
  out <- rep('O', L)
  for (i in 1:L)
  {
    if (label[i] %in% c('0-15', '16-29', '30-39', '40-49')){out[i] <- '0-49'}
    if (label[i] %in% c('50-59', '60-69', '70-74')){out[i] <-'50-74'}
    if (label[i] %in% c('75-79','80-84', '85-89', '90+')){out[i] <-'90+'}
  }
  if('0' %in% out){stop("Label does not correspond to classification")}
  return(out)
}
merge12 <- 1


# correlation table

# 1
tmp <- dht[Outcome2 %in% c('Death_Covid') & loc.res.cat == 'resident' &
      grepl(args$keep.vac.type, vac_type) & hosp.cat=='in hospital'  &
      DateAdmission >= w1.start & DateAdmission <= w2.end, ]

# tmp <- subset(dht, Outcome2 == 'Death_Covid' & loc.res.cat == 'resident' &
#                 DateAdmission >= w1.start & DateAdmission <= w2.end)

if(merge12){tmp <- tmp[, age.label := change_age_label(age.label)]}
tmp1 <- tmp[,list(observed_deaths = sum(wdeaths)),by = c('loc_label', 'age.label')]
tmp1 <- dcast.data.table(tmp1, loc_label~age.label,value.var = 'observed_deaths')
tmp1[,loc_label := change_city_label(loc_label)]
tmp1[, Total := rowSums(tmp1[,-1])]
setcolorder(tmp1, c(1,dim(tmp1)[2],2:(dim(tmp1)[2]-1)))
tmp1 <- tmp1[, lapply(.SD, as.integer), by=loc_label] # make all int
if(!merge12){
  saveRDS(tmp1, paste0(out.base, '_hosps_deaths_supp.rds'))
}

#2
tmp <- subset(dht, Outcome2 %in% c('Death_Covid','Death_Covid_Estimated') & hosp.cat=='in hospital'&
                loc.res.cat == 'resident' & grepl(args$keep.vac.type, vac_type) &
                DateAdmission >= w1.start & DateAdmission <= w2.end)

if(merge12){tmp <- tmp[, age.label := change_age_label(age.label)]}
tmp2 <- tmp[,list(deaths_censadj = as.integer(sum(wdeaths))),by = c('loc_label', 'age.label')]
tmp2 <- dcast(tmp2, loc_label~age.label,value.var = 'deaths_censadj')
tmp2[,loc_label := change_city_label(loc_label)]
tmp2[, Total := rowSums(tmp2[,-1])]
setcolorder(tmp2, c(1,dim(tmp2)[2],2:(dim(tmp1)[2]-1)))
tmp2 <- tmp2[, lapply(.SD, as.integer), by=loc_label]
if(!merge12){
saveRDS(tmp2, paste0(out.base, '_hosps_deaths_censadj_supp.rds'))
}

# 1,2
tmp <- merge(tmp1, tmp2, by = 'loc_label')
tmp <- tmp[, lapply(.SD, as.integer), by=loc_label]
saveRDS(tmp, paste0(out.base, '_hosps_deaths_and_censadj_supp.rds') )

# 
tmp <- readRDS(file=paste0(out.base,'_high_hfr_table.rds'))

tmp1 <- print(xtable(tmp), include.rownames=FALSE,  include.colnames = FALSE)
tmp1 <- gsub('^(.*?)begin\\{tabular\\}(.*?)\\}','',tmp1)
tmp1 <- gsub('\\\\end\\{tabular\\}(.*?)$','',tmp1)
tmp1 <- sub('\\\\hline','',tmp1) 


##################
# Additional XXs #
##################

#1: Introduction
intro <- list()

# Replacement dynamics in manaus
if('manaus' %in% args$selected.locs & 0)
{
  tmp <- paste0(out.base,'_manaus_fit.rds')
  cat('\nAttempting to load ', tmp)
  m_fit <- readRDS( tmp)
  pd <- m_fit$draws(inc_warmup = FALSE)		
  
  select.chains <- seq_along(dimnames(pd)[['chain']])
  iters <- 1:(m_fit$metadata()[['iter_sampling']])
  po <- list()
  tmp <- pd[,,which(grepl('^prop_P1',dimnames(pd)[[3]]) & !grepl('logit_prop_P1',dimnames(pd)[[3]]) & !grepl('age_prop_P1',dimnames(pd)[[3]])  & !grepl('prop_P1_fit',dimnames(pd)[[3]]))]
  po$prop_P1 <- unname(apply(tmp[,select.chains,], 3, rbind))
  pad <- as.data.table(melt(po$prop_P1))
  setnames(pad, 1:3, c('iterations','week.idx','value'))
  pad <- pad[, list(value=quantile(value, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('week.idx')]
  pad <- dcast.data.table(pad, week.idx~stat, value.var='value')
  tmp <- dh[loc_label == 'manaus', list(week=unique(week),week.start = unique(week.start), week.idx= unique(week)-min(week)+1L),]
  pad <- merge(pad, tmp, by = 'week.idx' )
  
  tmp <- pad[M > 0.5, min(week.start)] # want to translate this in terms of 
  tmp1 <- ddates[loc_label == 'manaus', w2.start]
  tmp - tmp1
}

# proportion of deaths
# dha[,range(DateSymptoms)]

tmp <- dha[ hosp.cat=='in hospital' & adm.cat %in% adm.classes &
              w1.start <= DateAdmission & DateAdmission <= w2.end]
intro[['SIVEP']][['tot_admissions_resnonres']] <- tmp[, write_commas(.N)]
tmp <- tmp[ loc.res.cat == 'resident']
intro[['SIVEP']][['tot_admissions_vaxnonvax']] <- tmp[, write_commas(.N)]
intro[['SIVEP']][['admissions_vax_diff']] <- tmp[! grepl(args$keep.vac.type, vac_type),write_commas(.N) ]
tmp <- tmp[ grepl(args$keep.vac.type, vac_type) ]
intro[['SIVEP']][['max_admission_date']] <- write.date(tmp[, max(DateAdmission)])
intro[['SIVEP']][['tot_admissions']] <-  tmp[, write_commas(.N)]
tmp1 <- tmp[, list(number=sum(Outcome2 == 'Death_Covid'), percentage = sum(Outcome2 == 'Death_Covid')/.N)]
intro[['SIVEP']][['N_deaths']] <- write_commas(tmp1$number)
intro[['SIVEP']][['p_deaths']] <- as.character(round(tmp1$percentage*100, 1))
tmp1 <- tmp[, list(number=sum(grepl('in_ICU|unknown',Outcome2)), percentage = sum(grepl('in_ICU|unknown',Outcome2))/.N)]
intro[['SIVEP']][['N_undetermined']] <- write_commas(tmp1$number)
intro[['SIVEP']][['p_undetermined']] <- as.character(round(tmp1$percentage*100, 1))
intro[['w2.end']] <- write.date(ddates[, unique(w2.end)])

tmp1 <- da[, paste0(unique(age.label), collapse = ', ')]
intro[['age_brackets']] <- gsub(",([^,]*)$"," and\\1",tmp1)


dha[ hosp.cat=='in hospital' & adm.cat %in% adm.classes &
       w1.start <= DateAdmission & DateAdmission <= w2.end &  loc.res.cat == 'resident' &
       grepl(args$keep.vac.type, vac_type), list(number=sum(Outcome2 == 'Death_Covid'), percentage = sum(Outcome2 == 'Death_Covid')/.N) ]
# 
health_care_indices <- list()
tmp <- dhcb[,list(ventilator=r.ventilator[1], week = week[1]),by = 'loc_label']
tmp <- tmp[ventilator %in% range(ventilator)][,loc_label := change_city_label(loc_label)]
health_care_indices[['min_ventilator']] <- write_commas(round(tmp$ventilator[1],1))
health_care_indices[['min_ventilator_loc']] <- tmp$loc_label[1]
health_care_indices[['max_ventilator']] <- write_commas(round(tmp$ventilator[2],1))
health_care_indices[['max_ventilator_loc']] <-tmp$loc_label[2]
tmp <- dhcb[,list(physicians=r.physicians.all[1], week = week[1]),by = 'loc_label']
tmp <- tmp[physicians %in% range(physicians)][,loc_label := change_city_label(loc_label)]
health_care_indices[['min_physicians']] <- write_commas(round(tmp$physicians[1],1))
health_care_indices[['min_physicians_loc']] <- tmp$loc_label[1]
health_care_indices[['max_physicians']] <- write_commas(round(tmp$physicians[2],1))
health_care_indices[['max_physicians_loc']] <-tmp$loc_label[2]



tmp <- c(intro, health_care_indices)

saveRDS(tmp, paste0(out.base, '_extra_XX.rds') )
# tmp <- readRDS(paste0(out.base, '_extra_XX.rds'))



###################################
# Flowchart Numbers
###################################

blocks_numbers <- list()

# missing the value across ALL BRAZIL
tmp <- file.path(pkg.dir, 'inst/data', basename(args$file.SIVEP))
SIVEP <- as.data.table(readRDS( tmp ))
tmp <- nrow(SIVEP[Admission == T &  DateAdmission <= unique(ddates$w2.end)])
blocks_numbers$Brazil_entries <- tmp
blocks_numbers$w2.end <- ddates$w2.end[1]
tmp <- dha[hosp.cat=='in hospital' & w1.start <= DateAdmission & DateAdmission <= w2.end]
blocks_numbers$SARI_admissions_14locs <- nrow(tmp)
tmp1 <- tmp[adm.cat %in% adm.classes] # Difference is in here I believe! Do we also use unknowns! or whattttt
blocks_numbers$COVID_admissions <- nrow(tmp1)
blocks_numbers$COVID_excluded <- nrow(tmp) - nrow(tmp1)
tmp <- tmp1[loc.res.cat == 'resident']
blocks_numbers$resident_admissions <- nrow(tmp)
blocks_numbers$resident_excluded <- nrow(tmp1)-nrow(tmp)
tmp1 <- tmp[grepl(args$keep.vac.type, vac_type)]
blocks_numbers$novax_admissions <- nrow(tmp1)
blocks_numbers$novax_excluded <- nrow(tmp) - nrow(tmp1)
tmp2 <- tmp1[Outcome2 %in% c('in_ICU', 'unknown')]
blocks_numbers$unknown_outcome <- nrow(tmp2)
blocks_numbers$observed_deaths <- tmp1[Outcome2 %in% c('Death_Covid', 'Death_Other'), .N]
blocks_numbers$observed_discharged <- tmp1[Outcome2 == 'discharged', .N]
tmp <- dht[w1.start <= DateAdmission & DateAdmission <= w2.end & hosp.cat=='in hospital' &
      loc.res.cat == 'resident' &  grepl(args$keep.vac.type, vac_type) ]
blocks_numbers$estimated_deaths <- tmp[, as.integer(sum(wdeaths))]
blocks_numbers$estimated_discharged <- tmp[, as.integer(sum(whosps)) - as.integer(sum(wdeaths)) ]

blocks_numbers$estimated_deaths <- blocks_numbers$estimated_deaths - blocks_numbers$observed_deaths
blocks_numbers$estimated_discharged <- blocks_numbers$estimated_discharged - blocks_numbers$observed_discharged

blocks_numbers$observed_deaths2 <- round(blocks_numbers$observed_deaths/blocks_numbers$novax_admissions*100,2)

blocks_numbers <- lapply(blocks_numbers, write_commas)
blocks_numbers$observed_deaths2 <-paste0(blocks_numbers$observed_deaths, ' (',blocks_numbers$observed_deaths2,'\\%)')

saveRDS(blocks_numbers, paste0(out.base, '_flowchart_numbers.rds'))


###########################
# Entirety of Brazil counterfactuals
###########################

dht <- SIVEP

setnames(dht, c('Municip','Municip_resid'), c('loc_label','loc.res.label'))
set(dht, NULL, 'DateSymptoms', dht[, as.Date(DateSymptoms, format='%Y-%m-%d')])
set(dht, NULL, 'DateAdmission', dht[, as.Date(DateAdmission, format='%Y-%m-%d')])
set(dht, NULL, 'DateICUstart', dht[, as.Date(DateICUstart, format='%Y-%m-%d')])
set(dht, NULL, 'DateDeath', dht[, as.Date(DateDeath, format='%Y-%m-%d')])
set(dht, NULL, 'loc_label', dht[, tolower(loc_label)])
set(dht, NULL, 'loc.res.label', dht[, tolower(loc.res.label)])
set(dht, dht[, which(is.na(loc.res.label))], 'loc.res.label', 'unknown')
dht[, loc.res.cat := factor(loc.res.label==loc_label,levels=c(TRUE,FALSE),labels=c('resident','not resident'))]
dht[, adm.cat := 'unknown']
set(dht, dht[, which(!is.na(CLASSI_FIN) & CLASSI_FIN%in%c(4,5))], 'adm.cat', 'COVID19 susp/conf')
set(dht, dht[, which(!is.na(CLASSI_FIN) & CLASSI_FIN%in%c(1,2,3))], 'adm.cat', 'other')
set(dht, NULL, 'adm.cat', dht[,factor(adm.cat, levels=c('COVID19 susp/conf', 'other', 'unknown'), labels=c('COVID19 susp/conf', 'other', 'unknown'))])
dht[, vac := 'unknown']
set(dht, dht[, which(!is.na(Vaccinated) & Vaccinated==TRUE)], 'vac', 'yes')
set(dht, dht[, which(!is.na(Vaccinated) & Vaccinated==FALSE)], 'vac', 'not vaccinated')
dht[, vac_type := vac]
set(dht, dht[, which(vac=='yes' & is.na(Date1stDose))], 'vac_type', 'vaccinated, date of dose unknown')
set(dht, dht[, which(vac=='yes' & !is.na(Date1stDose) & Date1stDose+14>DateSymptoms)], 'vac_type', 'vaccinated, 1 dose less than 2 weeks')
set(dht, dht[, which(vac=='yes' & !is.na(Date1stDose) & Date1stDose+14<=DateSymptoms)], 'vac_type', 'vaccinated, 1 dose at least 2 weeks')
set(dht, dht[, which(vac=='yes' & !is.na(Date2ndDose) & Date2ndDose+14<=DateSymptoms)], 'vac_type', 'vaccinated, 2 doses at least 2 weeks')
set(dht, NULL, 'vac_type', dht[, factor(vac_type, levels=c('not vaccinated','unknown','vaccinated, 1 dose less than 2 weeks','vaccinated, 1 dose at least 2 weeks','vaccinated, 2 doses at least 2 weeks','vaccinated, date of dose unknown'))])
set(dht, NULL, c('vac','Vaccinated'),NULL)
dht[, Outcome2 := Outcome]
set(dht, dht[,which(is.na(Outcome2))], 'Outcome2', 'unknown')
set(dht, dht[,which(Outcome2=='unknown' & !is.na(DateICUstart) & is.na(DateICUend))], 'Outcome2', 'in_ICU')
set(dht, NULL, 'Outcome2', dht[,factor(Outcome2, levels=c('Death_Covid','Death_Other','in_ICU','Cure','unknown'),labels=c('Death_Covid','Death_Other','in_ICU','discharged','unknown'))])
dht[, hosp.cat := 'in hospital']
tmp <- dht[, which( 
  ( is.na(DateAdmission) ) | 
    ( !is.na(DateDeath) & !is.na(DateAdmission) & DateAdmission==DateDeath )
)]
set(dht, tmp, 'hosp.cat', 'out of hospital')
set(dht, dht[, which(is.na(Hospital_type))], 'Hospital_type', 'Unknown')
set(dht, dht[,which(is.na(CLASSI_FIN) & Outcome2=='Death_Other')],'Outcome2','Death_Covid')
set(dht, dht[,which(is.na(CLASSI_FIN))],'adm.cat','COVID19 susp/conf')
dht <- subset(dht, adm.cat=='COVID19 susp/conf')
stopifnot(nrow(dht[Outcome2 == 'Death_Other']) == 0)

# dht[adm.cat=='COVID19 susp/conf' & DateAdmission <= unique(ddates$w2.end) & hosp.cat == "in hospital", table(Outcome2) ]
tmp <- dht[adm.cat=='COVID19 susp/conf' & DateDeath <= unique(ddates$w2.end)
           & hosp.cat == "in hospital" &  Outcome2 == 'Death_Covid', .N]


tmp1 <- as.data.table(readRDS(paste0(out.base, '_main_table1.rds')))
tmp1 <- tmp1[which(loc_label == 'All 14 cities') + 0:1, .(overall.y, hyp_deaths_of_hosps_pcred)]
tmp1 <- tmp1[, lapply(.SD, function(x){paste0(x, collapse = '')})]
tmp1 <- as.data.table(t(tmp1)); colnames(tmp1) <- 'value'

tmp1 <- tmp1[, list(M =as.numeric(gsub('^(.*?)\\%\\((.*?)\\%-(.*?)\\%\\)',"\\1",tmp1$value)), 
            CL=as.numeric(gsub('^(.*?)\\%\\((.*?)\\%-(.*?)\\%\\)',"\\2",tmp1$value)),
            CU=as.numeric(gsub('^(.*?)\\%\\((.*?)\\%-(.*?)\\%\\)',"\\3",tmp1$value)))]

# multiply by relevant percentage
counterfactual_all <- list(
  conservative = write_commas(as.integer(tmp * tmp1[1,] /100)),
  nonconservative = write_commas(as.integer(tmp * tmp1[2,] /100))
)
counterfactual_all <- lapply(counterfactual_all, function(vec){paste0(vec[1],' (', vec[2],'-', vec[3], ')')})
counterfactual_all$tot_hosp_deaths <- write_commas(tmp)

saveRDS(counterfactual_all, paste0(out.base, '_counterfactual_all.rds'))

