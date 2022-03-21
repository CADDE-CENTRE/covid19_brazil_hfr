cat(" \n -------------------------------- \n \n HFR.preprocessing.R \n \n -------------------------------- \n")

#### HFR.preprocessing.210719.R ####

# Includes data preprocessing described in supplementary text and outputs
# the workspace.rda file to be called prior to model fitting in HFR.fit.210719.cmdstan.R
# To perform "Unknown to discharged sensitivity analysis", the args$deaths.unknown2discharged
# parameter needs to set to 1 either manually or by including 'unkdis' to the output directory
####

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

data.table::setDTthreads(4)

if(1) # Andrea
{
  pkg.dir <- '~/git/covid19_brazil_hfr/'	
  out.dir <- '~/Documents/P1Brazil/hfr-fit-210719h4-resources8-ex0-trvaluecs-hd1'
  file.prefix <- basename(out.dir)
  plot.dir <- out.dir
  stan.model <- 'age_hfr_210719d.stan'
  genomic.data.type <- 'GISAID_darlan'
  fr.predictors <- "icu.per.ventilator.02,icu.per.beds.02,icu.per.nurse.02,icu.per.tech.nurse.02,icu.per.physiotherapist.02,icu.per.physicians.all.02,icu.per.intensivists.02,hosps.per.saribeds.02,hosps.per.saribedsvent.02,hosps.per.physicians.all.02,hosps.per.ventilator.02,hosps.per.nurse.02"
  fr.predictors.names <- gsub(' ','_',"2-wk ICU admissions per ventilator,2-wk ICU admissions per ICU bed,2-wk ICU admissions per nurse,2-wk ICU admissions per nurse assistant,2-wk ICU admissions per physiotherapist,2-wk ICU admissions per physician,2-wk ICU admissions per intensivist,2-wk SARI admissions per critical bed, 2-wk SARI admissions per critical bed with ventilator,2-wk SARI admissions per physician,2-wk SARI admissions per ventilator,2-wk SARI admissions per nurse")
  fr.predictors.transformation <- 'value_cs'
  keep.hospital.type <-'Private_Hospitals|Public_Hospitals|Unknown' # 'Private_Hospitals'
}

# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-outdir')	
  stopifnot(args_line[[3]]=='-file_prefix')
  stopifnot(args_line[[5]]=='-pkgdir')
  stopifnot(args_line[[7]]=='-stanmodel')
  stopifnot(args_line[[9]]=='-genomic_data_type')
  stopifnot(args_line[[11]]=='-fr_predictors')
  stopifnot(args_line[[13]]=='-fr_predictors_names')
  stopifnot(args_line[[15]]=='-fr_predictors_transformation')
  stopifnot(args_line[[17]]=='-keep.hospital.type')
  plot.dir <- args_line[[2]]
  out.dir <- args_line[[2]]
  file.prefix <- args_line[[4]]
  pkg.dir <- args_line[[6]]
  stan.model <- args_line[[8]]
  genomic.data.type <- gsub(',','|',args_line[[10]])
  fr.predictors <- args_line[[12]]
  fr.predictors.names <- args_line[[14]]
  fr.predictors.transformation <- args_line[[16]]
  keep.hospital.type <- gsub('_',' ',gsub(',','|',args_line[[18]]))
} 

stopifnot(fr.predictors.transformation %in% c('value_cs','value_log_cs'))
stopifnot(genomic.data.type %in% c('random_sample|GISAID_alfonso_lutz','genomica_fiocruz','GISAID_darlan'))

# set args 
args <- list()
args$out.dir <- out.dir
args$out.base <- file.path(out.dir,file.prefix)
args$file.brazil.states <- file.path(pkg.dir, 'inst/data', 'brazil_states.csv')
args$file.SIVEP <- file.path(pkg.dir, 'inst/data', 'SIVEP_hospital_31-01-2022-all.rds')

args$file.registry.death <- file.path(pkg.dir,'inst/data','registry_covid_detailed_09-08-2021.csv')
args$file.excess.burials <- file.path(pkg.dir, 'inst/data', 'private_public_burials_manaus_210706.csv')
args$file.p1.emergence <- file.path(pkg.dir, 'inst/data', '210627_phylo_P1_dates_emergence.csv')
args$file.demographics <- file.path(pkg.dir, 'inst/data', 'PNADc_populationpyramids_210617.csv')
args$file.vaccinations <- file.path(pkg.dir, 'inst/data', 'aggregated_vaccinations_210805.rds')
args$file.7_digit_IBGE_code <- file.path(pkg.dir, 'inst/data', 'municipios_brasil_dadoscompletos.csv')
args$file.hospital.demand.beds.physicians <- file.path(pkg.dir,'inst/data','IPEA_ICUbeds_physicians_210928.csv')
args$file.hospital.demand.ventilators <- file.path(pkg.dir, 'inst/data', 'MOH_CNES_ventilators_210628.csv')
args$file.hospital.collapse.times <- file.path(pkg.dir, 'inst/data', 'hospital_collapse_cities_list.csv')

# args$file.verity.ifr.estimates <- file.path(pkg.dir,'inst/data','Bob_latest_estimates_IFR_210430.csv')
args$file.genomic.data <- file.path(pkg.dir, 'inst/data','genomic_data_210702.csv')
# args$file.posterior.time.seroconversion.to.seroreversion.parameters <- file.path(pkg.dir,'inst/data','manaus_age_ifr_210428a_manaus_convergence.csv')
args$file.stan.model <- file.path(pkg.dir,'inst/stan-models',stan.model)
args$selected.locs <- c('belo horizonte','curitiba', "florianopolis", "goiania", "joao pessoa", "macapa", 'manaus', "natal", 'porto alegre', "porto velho", 'rio de janeiro', 'salvador', 'sao paulo', "sao luis")
# args$selected.locs <- c("aracaju","belem","belo horizonte" ,"boa vista","brasilia","campo grande","cuiaba","curitiba","florianopolis",
#                         "fortaleza","goiania","joao pessoa","macapa","maceio","manaus","natal","palmas","porto alegre", 
#                         "porto velho","recife","rio branco","rio de janeiro","salvador","sao luis","sao paulo","teresina","vitoria")
# args$selected.locs <- c('goiania', 'manaus', 'rio de janeiro', 'sao paulo')
tmp <- 74 + 9 # set this as the week.end for every state (before 74)?
args$select.deaths.week.end <- c("aracaju"=tmp, "belem"=tmp, "belo horizonte"=tmp ,"boa vista"=tmp, "brasilia"=tmp, "campo grande"=tmp, "cuiaba"=tmp, "curitiba"=tmp, "florianopolis"=tmp, 
                                 "fortaleza"=tmp, "goiania"=tmp, "joao pessoa"=tmp, "macapa"=tmp, "maceio"=tmp, "manaus"=tmp, "natal"=tmp, "palmas"=tmp, "porto alegre"=tmp,  
                                 "porto velho"=tmp, "recife"=tmp, "rio branco"=tmp, "rio de janeiro"=tmp, "salvador"=tmp, "sao luis"=tmp, "sao paulo"=tmp, "teresina"=tmp, "vitoria"=tmp)
args$deaths.weeks.increment <- 0L
args$deaths.adjust.for.civilregistrydeaths <- 'reg.deaths.all.exorcovid'
args$deaths.estimate.outcomes.for.those.unknown <- if(grepl('-nocens', args$out.base)){0}else{1}
args$deaths.unknown2discharged <- if(grepl('unkdis', args$out.base)){1}else{0}
args$select.vacc.increasing.vacc.rates <- 0
args$select.i2sc.max <- 45
args$select.h2d.max <- 65
args$select.s2d.max <- 65
args$select.only.SIVEP.patients.that.are.resident.in.loc <- 1
args$excess.death.zero.among.0.to.16 <- 0

args$vaccine.population <- 'loc_res_label'
args$vaccine.max.coverage <- 0.99
args$vaccine.astra.ve.1dose <- .40 
args$vaccine.astra.ve.2doses <- .91 
args$vaccine.astra.protection.delay <- 3L
args$vaccine.sinovac.ve.1dose <- .32 
args$vaccine.sinovac.ve.2doses <- .71 
args$vaccine.sinovac.protection.delay <- 2L
args$genomic.data.type <- genomic.data.type


args$make.preprocessing.plots <- if(length(args_line) > 0){1}else{0}
args$make.sexpr <-  if(length(args_line) > 0){1}else{0}
args$fr.predictors.transformation <- fr.predictors.transformation 
args$fr.predictors <- unlist(strsplit(fr.predictors, split=','))
names(args$fr.predictors) <-unlist(strsplit( gsub('_',' ',fr.predictors.names), split=','))
args$fr.predictors.ref.week.lower.bound <- '2020-01-01'
args$keep.hospital.type <- gsub('_',' ',gsub(',','|',keep.hospital.type))
#args$keep.vac.type <- "^unknown$|^no$|^vaccinated"
args$keep.vac.type <- "^unknown$|^not vaccinated$"
#args$keep.vac.type <- "^unknown$|^no$|^vaccinated, 1 dose less than 2 weeks$"

# By unclassified we mean COVID19 cases with CLASSI_FIN == 4
args$remove.unclassified.cases <-  if(grepl('-onlyconfirmed', args$out.base)){1}else{0}   
args$keep.CLASSIFIN.na <- 0

if(genomic.data.type == 'random_sample|GISAID_alfonso_lutz' )
{
  args$selected.locs <-  c('belo-horizonte', 'manaus', 'sao-paulo') 
}

args$city.palette <- grDevices::colorRampPalette(ggsci::pal_npg("nrc")(10))
args$nrows.plots <- ceiling(length(args$selected.locs)/5)
args$ncols.plots <- ceiling(length(args$selected.locs)/5)
# args$selected.locs.part1 <- args$selected.locs 
# args$selected.locs.part2 <- c()
# args$selected.locs.part3 <- c()
args$partition.locs <- 0
if(length(args$selected.locs) >= 6){
  args$partition.locs <- 1

  main.cities <- c('goiania', 'manaus', 'rio de janeiro', 'sao paulo')
  args$selected.locs.part1 <- main.cities
  
  tmp <- function(x) split(x, cut(seq_along(x), 2, labels = FALSE))
  
  args$selected.locs.part1 <- main.cities
  args$selected.locs.part2 <- tmp(args$selected.locs[! args$selected.locs %in% main.cities])[[1]]
  args$selected.locs.part3 <- tmp(args$selected.locs[! args$selected.locs %in% main.cities])[[2]]
}

print(args)

####################################################################
# Submission requirements
####################################################################
reqs <- theme(text=element_text(size=7, family = 'sans'), axis.text=element_text(size=5, family = 'sans'))

####################################################################
# source functions
####################################################################

tmp <- list.files(file.path(pkg.dir,'inst/R'), full.names=TRUE)
cat('sourcing code from files\n', tmp)
invisible(lapply(tmp, source))

####################################################################
# start building info for paper
####################################################################

if(args$make.sexpr)
{
  sexpr <- list()
  sexpr[['args']] <- args
  write.date <- function(date){return( as.character(format(date, "%d %B %Y")))}
  
  tmp <- c('belo horizonte', 'curitiba','florianopolis','goiania','joao pessoa','macapa', "manaus",'Natal','porto alegre','porto velho', 'rio de janeiro', 'salvador','sao luis', 'sao paulo')
  sexpr[['cities']][['names']]<- change_city_label(tmp)
  sexpr[['cities']][['number']] <- as.character(length(tmp))
  
  tmp <- stringr::str_match(args$file.vaccinations, pattern='vaccinations_(.*?).rds')[,2]
  tmp <- sub("(.{4})(.*)", "\\1/\\2", tmp)
  tmp <- sub("(.{2})(.*)", "\\1/\\2", tmp)
  tmp <- write.date(as.Date(tmp))
  sexpr[['vaccines']][['date_of_retrieval']] <- tmp # CHECK format...
  
  tmp <-  basename(args$file.hospital.demand.beds.physicians)
  sexpr[['resources']][['file']] <- gsub('_', '\\\\_', tmp)
}

####################################################################
# make age brackets
####################################################################

age.cat <- c(0,16,30,40,50,60,70,75,80,85,90,150)
da <- data.table(Age= 0:120)
da[, age.label := cut(Age, breaks=age.cat,right = FALSE)]
tmp2 <- unique(subset(da, select=age.label))
tmp2[, age.from := age.cat[-length(age.cat)]]
tmp2[, age.to := age.cat[-1]]
tmp2[, age.mid := (age.from+age.to)/2]
set(tmp2, nrow(tmp2), 'age.mid', tmp2[nrow(tmp2),age.from] )
da <- merge(da, tmp2, by='age.label')
set(da, NULL, 'age.label', da[, factor(gsub('90-149','90+',paste0(age.from,'-',age.to-1)))])

####################################################################
# make week time table
####################################################################

dw <- data.table(Date= seq(as.Date('2019-01-01') , as.Date('2021-12-31'), 1 ) )
dw[, week := as.integer(strftime(Date, format='%V'))]	
dw[, month := as.integer(strftime(Date, format='%m'))]
dw[, year := as.integer(strftime(Date, format='%Y'))]
tmp <- dw[,which(Date<='2019-12-29')]
set(dw, tmp, 'week',dw[tmp,week-52L])
tmp <- dw[,which(Date<='2019-12-31')]
set(dw, tmp, 'month',dw[tmp,month-12L])
tmp <- dw[,which(Date>'2021-01-03')]
set(dw, tmp, 'week',dw[tmp,week+53L])
tmp <- dw[,which(Date>='2021-01-01')]
set(dw, tmp, 'month',dw[tmp,month+12L])
tmp <- dw[,list(week.start=min(Date)),by='week']
dw <- merge(dw, tmp, by='week')
tmp <- dw[,list(month.start=min(Date)),by='month']
dw <- merge(dw, tmp, by='month')

# select only full weeks
# first week is padded below, so just need to make sure last week is complete
if( subset(dw, week==max(dw$week))[, max(Date)] < subset(dw, week==max(dw$week))[1, week.start+7-1] )
{
	dw <- subset(dw, week<max(dw$week))
	stopifnot( max(dw$week) >= args$select.deaths.week.end )
}

####################################################################
# define dates
####################################################################

cat('\ndefine dates ...')
ddates <- read.csv(args$file.hospital.collapse.times)
ddates <- data.table(
	  loc_label = tolower(ddates$City),
	  w1.start = as.Date(ddates$W1Start, '%d/%m/%Y'),
	  w1.end = as.Date(ddates$W1End,'%d/%m/%Y'),
	  w2.start = as.Date(ddates$W2Start,'%d/%m/%Y'),
	  w2.end = as.Date('2021-07-26')  
	)

set(ddates, ddates[,which(loc_label=='belo horizonte')],'w1.start',as.Date('2020-04-06'))
ddates[loc_label=='belo horizonte', w2.start:=w1.end + 1]
set(ddates, ddates[,which(loc_label=='sao paulo')],'w1.start',as.Date('2020-01-20'))
ddates[loc_label=='sao paulo', w2.start:=w1.end + 1]
tmp <- unique(subset(dw, select=c(week,week.start)))
tmp2<- copy(tmp)
setnames(tmp, c('week','week.start'),c('hosps.start','w1.start'))
ddates <- merge(ddates, tmp, by='w1.start')
setnames(tmp2, c('week', 'week.start'), c('p1.detection','w2.start'))

ddates[, w2.start.week := as.Date(cut(w2.start, 'week'))]
ddates <- merge(ddates, tmp2, by.x = 'w2.start.week', by.y='w2.start')
ddates[, w2.start.week := NULL]

# set start time of the death weeks in our model
ddates[, deaths.start:= hosps.start+args$deaths.weeks.increment]
tmp <- unique(subset(dw, select=c(week,week.start)))
setnames(tmp,c('week','week.start'),c('deaths.start','deaths.week.start'))
ddates <- merge(ddates, tmp, by='deaths.start')
# set end time of the death weeks in our model
tmp <- data.table( loc_label=names(args$select.deaths.week.end), deaths.end=as.numeric(args$select.deaths.week.end))
ddates <- merge(ddates, tmp, by='loc_label')
tmp <- unique(subset(dw, select=c(week,week.start)))
setnames(tmp,c('week','week.start'),c('deaths.end','deaths.week.end'))
ddates <- merge(ddates, tmp, by='deaths.end')

# set P.1 emergence times in our model
dstates <- as.data.table(read.csv(args$file.brazil.states, stringsAsFactors = FALSE ))
dstates <- subset(dstates, select=c(state, state_label, loc_label, loc_name,region_label))
tmp <- as.data.table(read.csv(args$file.p1.emergence))
tmp[, state := gsub('AM_to_','',Summary.Statistic)]
tmp <- merge(tmp, dstates, by='state')
setnames(tmp, c('median','lower','upper'), c('p1em.median','p1em.CL','p1em.CU'))
tmp <- subset(tmp, select=c(loc_label, p1em.median, p1em.CL, p1em.CU))
tmp <- melt(tmp, id.vars='loc_label')
tmp[, year:= paste0(floor(value),'-01-01')]
tmp[, day:= floor((value%%1) * 365)]
tmp[, date := as.Date(day, origin = year)]
tmp <- dcast.data.table(tmp, loc_label~variable, value.var='date')
ddates <- merge(ddates, tmp, by='loc_label',all.x=TRUE)
# set P.1 emergence time for Manaus in our model
set(ddates, ddates[, which(loc_label=='manaus')], 'p1em.median', as.Date('2020-12-04'))
set(ddates, ddates[, which(loc_label=='manaus')], 'p1em.CL', as.Date('2020-12-01'))
set(ddates, ddates[, which(loc_label=='manaus')], 'p1em.CU', as.Date('2020-12-07'))

# define subsets for plotting
ddates <- subset(ddates,loc_label %in% args$selected.locs)
if(args$partition.locs){
  ddates1 <- subset(ddates,loc_label %in% args$selected.locs.part1)
  ddates2 <- subset(ddates,loc_label %in% args$selected.locs.part2)
  ddates3 <- subset(ddates,loc_label %in% args$selected.locs.part3)
}

if(args$make.sexpr){
  sexpr[['ddates']] <- ddates  
  tmp <- dstates[loc_label %in% args$selected.locs, .(loc_label, state_label)] # table S1
  tmp <- merge(tmp, subset(ddates, select = c('loc_label','w1.start', 'w2.end', 'w2.start', 'p1em.median')), by = 'loc_label')
  set(tmp, NULL, 'loc_label', change_city_label(tmp$loc_label))
  tmp2 <- function(x){format(x, '%d/%m/%Y')}
  # tmp[, `:=` (w1.start = write.date(w1.start), w2.end = write.date(w2.end), w2.start = write.date(w2.start), p1em.median = write.date(p1em.median)) ,]
  tmp[, `:=` (w1.start = tmp2(w1.start), w2.end = tmp2(w2.end), w2.start = tmp2(w2.start), p1em.median = tmp2(p1em.median)),]
  sexpr[['dates']] <- tmp
  rm(tmp, tmp2)
}


####################################################################
# Read vaccination data
###################################################################

dv <- readRDS(args$file.vaccinations)
dv <- dv[loc_label %in% args$selected.locs]

if(args$make.sexpr)
{
  # tmp <- dv[week.start == max(week.start)]
  tmp <- ddates[,unique(w2.end)]
  tmp <- dv[week.start <= tmp]
  tmp <- tmp[, list(doses_administered = sum(vac1 + vac2)),by = c('vaccine')]
  tmp[, p_doses_administered := doses_administered / sum(tmp$doses_administered)]
  tmp <- round(tmp[vaccine == 'Pfizer', p_doses_administered]*100,2)
  sexpr[['Pop_at_risk']][['p_pfizer_administrations']] <- as.character(tmp)
}


####################################################################
# READ POPULATION COMPOSITION
####################################################################

cat('\nread population composition ...')
dpop <- as.data.table(read.csv(args$file.demographics, stringsAsFactors=FALSE, encoding="latin1"))
setnames(dpop, c('city','population'), c('loc_label','pop'))
set(dpop, NULL, 'loc_label', dpop[, tolower(gsub('^\\s|\\s$','',gsub('é','e',gsub('í','i',gsub('ó','o',gsub('á|ã|â','a',gsub('\xe3','a',loc_label)))))))])
tmp <- subset(da, select=c(Age, age.label))
dpop <- merge(dpop, tmp, by.x='age',by.y='Age')
dpop <- subset(dpop, loc_label %in% args$selected.locs)

#	adjust population so that vaccine coverage is at most args$vaccine.max.coverage 
tmp <- dpop[, list(pop=sum(pop)), by=c('loc_label','age.label')]
tmp2 <- dv[, list(cvac = sum(cvac1+cvac2)), by=c('loc_label','age.label','week')]
tmp2 <- tmp2[, list(cvac=cvac[which.max(week)]), by=c('loc_label','age.label')]
tmp <- merge(tmp, tmp2, by=c('loc_label','age.label'), all.x=TRUE)
set(tmp, tmp[,which(is.na(cvac))],'cvac',0)

cat('\nlocation age combinations where vaccinated pop is larger than 2020 pop projection:\n')
print(subset(tmp, cvac>pop))

cat('\nadjusting 2020 pop projection so that max vaccine coverage is ',args$vaccine.max.coverage,' ...')
tmp[, vacpop := cvac/args$vaccine.max.coverage]
tmp[, popmultiplier := pmax(1.,vacpop/pop)]
dpop <- merge(dpop, subset(tmp, select=c(loc_label,age.label,popmultiplier)), by=c('loc_label','age.label'))
setnames(dpop, 'pop', 'pop.PNADc')
dpop[, pop := pop.PNADc * popmultiplier]
dpop[, pop.vacc.correction := pop-pop.PNADc]
set(dpop, NULL, 'popmultiplier', NULL)
if(0){
  tmp <- copy(dpop)
  tmp[, pop.vacc.correction := NULL]
  setnames(tmp, 'pop', 'pop_adj')
  write.csv(tmp, file.path(pkg.dir,'inst/data','PNDAc_and_vaxadj_population.csv'))
  rm(tmp)
}
set(dpop, NULL, 'pop.PNADc', dpop[, as.numeric(pop.PNADc)])

if(args$make.preprocessing.plots)
{	
	tmp <- melt(dpop, id.vars=c('loc_label','age.label','age','sex'), measure.vars=c('pop.PNADc','pop.vacc.correction'))
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('pop.vacc.correction','pop.PNADc'), labels=c('Adjusted','PNADc 2020 projection'))])
	set(tmp, tmp[, which(sex=='Men')], 'value', tmp[sex=='Men', -value] )
	p <- ggplot(tmp, aes(x = age, y = value, fill = sex, alpha=variable)) +
			geom_bar(data = subset(tmp, sex == "Women"), stat = "identity") +
			geom_bar(data = subset(tmp, sex == "Men"), stat = "identity") +
			facet_wrap(~change_city_label(loc_label), ncol=4, scales='free') +
			coord_flip() +
			scale_x_continuous(expand=c(0,0)) +
			scale_y_continuous(labels=abs) +
			scale_alpha_manual(values=c('Adjusted'=1.0,'PNADc 2020 projection'=0.5)) +			
			ggsci::scale_fill_npg() +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
				strip.background = element_blank(), legend.position='bottom') +
			labs(y= '2020 population size', x='1-year age band', fill='', alpha='')
	ggsave(file=paste0(args$out.base,'_pop_pyramid.pdf'),p,w=10,h=15)
	ggsave(file=paste0(args$out.base,'_pop_pyramid.png'),p,w=10,h=15)
}

# aggregate pop counts by age.label, cumulate and add total
dpop <- dpop[, list(pop=sum(pop)), by=c('loc_label','age.label')]
tmp <- dpop[, list(tpop = sum(pop)), by='loc_label']
dpop <- merge(dpop,tmp,by='loc_label')
setkey(dpop, loc_label, age.label)
tmp <- dpop[,list(age.label=age.label, qpop=cumsum(pop)/sum(pop)), by=c('loc_label')]
dpop <- merge(dpop, tmp, by=c('loc_label','age.label'))


####################################################################
# plot vaccine coverage in adjusted population
####################################################################

cat('\nplot vaccine coverage in adjusted population ...')
dv <- merge(dv, dpop, by = c("loc_label", "age.label"))
dv[, pcvac1 := cvac1/pop]
dv[, pcvac2 := cvac2/pop]

if(args$make.preprocessing.plots)
{
	tmp <- subset(dv, week.start > '2021-01-01')  
	tmp <- melt(tmp, 
			id.vars = c("loc_label", "age.label", 'week', 'week.start', 'vaccine'), 
			measure.vars = c('pcvac1','pcvac2'),
			variable.name = 'Dose' , 
			value.name = 'prop')  
	set(tmp, NULL, 'Dose', tmp[, factor(Dose, levels=c('pcvac1','pcvac2'), labels=c('Only one','Two'))])
	set(tmp, NULL, 'vaccine', factor(tmp$vaccine, levels = c("Covishield", "Sinovac", "Pfizer", "Janssen")))  
	
	if(args$partition.locs){
	  tmp1 <- tmp[loc_label %in% args$selected.locs.part1,]
	  p <- ggplot(tmp1, aes(x = week.start)) +
	    geom_area(aes(y = prop,  alpha = Dose, fill = vaccine), position = 'stack') + # switch fill and alpha position depending on what u want
	    geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') + 
	    facet_grid(age.label~ change_city_label(loc_label)) +
	    scale_x_date(breaks='1 month',expand=c(0,0)) + 
	    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
	    coord_cartesian(ylim=c(0,1)) +
	    ggsci::scale_fill_lancet() + 
	    scale_alpha_manual(values=c('Only one'=0.3,'Two'=1.)) +
	    labs(fill='vaccine', y='proportions of population\nhaving received one or two vaccine doses', x='', alpha = 'doses administred') +
	    theme_bw() +
	    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
	          legend.position = 'bottom',
	          strip.background = element_blank(),
	          panel.spacing.y = unit(1, "lines"))
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part1.pdf'),p + reqs, w=15,h=18, dpi=350, units='cm')
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part1.png'),p + reqs, w=15,h=18, dpi=350, units='cm')
	  ggsave(file=file.path(args$out.dir,'Brizzietal_extdatafig_6.pdf'),p + reqs, w=15,h=18, dpi=350, units='cm')
	  
	  tmp1 <- tmp[loc_label %in% args$selected.locs.part2,]
	  p <- ggplot(tmp1, aes(x = week.start)) +
	    geom_area(aes(y = prop,  alpha = Dose, fill = vaccine), position = 'stack') + # switch fill and alpha position depending on what u want
	    geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') + 
	    facet_grid(age.label~ change_city_label(loc_label)) +
	    scale_x_date(breaks='1 month',expand=c(0,0)) + 
	    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
	    coord_cartesian(ylim=c(0,1)) +
	    ggsci::scale_fill_lancet() + 
	    scale_alpha_manual(values=c('Only one'=0.3,'Two'=1.)) +
	    labs(fill='vaccine', y='proportions of population \nhaving received one or two vaccine doses', x='', alpha = 'doses administred') +
	    theme_bw() +
	    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
	          legend.position = 'bottom',
	          strip.background = element_blank(),
	          panel.spacing.y = unit(1, "lines"))
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part2.pdf'),p + reqs, w=15,h=18, dpi=350, units='cm')
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part2.png'),p + reqs, w=15,h=18, dpi=350, units='cm')
	  
	  tmp1 <- tmp[loc_label %in% args$selected.locs.part3,]
	  p <- ggplot(tmp1, aes(x = week.start)) +
	    geom_area(aes(y = prop,  alpha = Dose, fill = vaccine), position = 'stack') + # switch fill and alpha position depending on what u want
	    geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') + 
	    facet_grid(age.label~ change_city_label(loc_label)) +
	    scale_x_date(breaks='1 month',expand=c(0,0)) + 
	    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
	    coord_cartesian(ylim=c(0,1)) +
	    ggsci::scale_fill_lancet() + 
	    scale_alpha_manual(values=c('Only one'=0.3,'Two'=1.)) +
	    labs(fill='vaccine', y='proportions of population \nhaving received one or two vaccine doses', x='', alpha = 'doses administred') +
	    theme_bw() +
	    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
	          legend.position = 'bottom',
	          strip.background = element_blank(),
	          panel.spacing.y = unit(1, "lines"))
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part3.pdf'),p + reqs, w=15,h=18, dpi=350, units='cm')
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part3.png'),p + reqs, w=15,h=18, dpi=350, units='cm')
	}else{
	  require(scales)
	  tmp1 <- copy(tmp)
	  p <- ggplot(tmp1, aes(x = week.start)) +
	    geom_area(aes(y = prop,  alpha = Dose, fill = vaccine), position = 'stack') + # switch fill and alpha position depending on what u want
	    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') + 
	    facet_grid(age.label~ change_city_label(loc_label)) +
	    scale_x_date(breaks='1 month',expand=c(0,0),  labels=date_format("%b %y")) + 
	    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
	    coord_cartesian(ylim=c(0,1)) +
	    ggsci::scale_fill_lancet() + 
	    scale_alpha_manual(values=c('Only one'=0.3,'Two'=1.)) +
	    labs(fill='vaccine', y='proportions of population \nhaving received one or two vaccine doses', x='', alpha = 'doses administred') +
	    theme_bw() +
	    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
	          legend.position = 'bottom',
	          strip.background = element_blank(),
	          panel.spacing.y = unit(0.5, "lines"),
	          legend.key.size = unit(10, 'pt'),
	          legend.spacing.y = unit(5, 'pt'),
	          axis.text=element_text(size=6.5, family = 'sans'),
	          text=element_text(size=7, family = 'sans')) 
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part.pdf'),p,w=15,h=18, dpi=350, units='cm')
	  ggsave(file=paste0(args$out.base,'_pop_prop_vaxcomp_part.png'),p,w=15,h=18, dpi=350, units='cm')
	}
}

#####################################################
# read death data (registry)
####################################################################

cat('\nread death data (registry) ...')

dr <- as.data.table(read.csv(args$file.registry.death, stringsAsFactors=FALSE))
setnames(dr, c('date','city'), c('Date','loc_label'))
set(dr, NULL, 'loc_label', dr[, tolower(loc_label)])
dr <- subset(dr, loc_label %in% args$selected.locs)
tmp <- subset(dpop, age.label%in%c('70-74','75-79','80-84','85-89'))
tmp[, age_group := as.character(factor(grepl('^7', age.label), levels=c(TRUE,FALSE), labels=c('70-79','80-89')))]
setkey(tmp, loc_label, age_group, age.label)
tmp <- tmp[,list(ppop.04 = pop[1]/sum(pop)),by=c('loc_label','age_group')]
dr <- merge(dr, tmp, by=c('loc_label','age_group'), all.x=TRUE)

tmp <- subset(dr, age_group%in%c('70-79','80-89'))
set(tmp, NULL, 'deaths_non_covid', tmp[, deaths_non_covid*ppop.04])
set(tmp, NULL, 'deaths_covid', tmp[, deaths_covid*ppop.04])
set(tmp, NULL, 'deaths_total', tmp[, deaths_total*ppop.04])
set(tmp, tmp[, which(grepl('^7',age_group))], 'age_group', '70-74')	
set(tmp, tmp[, which(grepl('^8',age_group))], 'age_group', '80-84')
dr <- rbind(dr, tmp)
tmp <- dr[, which(age_group%in%c('70-79','80-89'))]
set(dr, tmp, 'deaths_non_covid', dr[tmp, deaths_non_covid*(1-ppop.04)])
set(dr, tmp, 'deaths_covid', dr[tmp, deaths_covid*(1-ppop.04)])
set(dr, tmp, 'deaths_total', dr[tmp, deaths_total*(1-ppop.04)])
set(dr, dr[, which(age_group%in%c('70-79'))], 'age_group', '75-79')
set(dr, dr[, which(age_group%in%c('80-89'))], 'age_group', '85-89')
set(dr, NULL, 'ppop.04', NULL)
set(dr, dr[, which(age_group=="100+")], 'age_group', '100-149')
set(dr, dr[, which(age_group=="9-")], 'age_group', '1-9')
set(dr, NULL, 'Date', dr[, as.Date(Date, format='%Y-%m-%d')])

# allocate deaths with missing age_group to age groups according the observed proportion of deaths by age group in the previous 2 weeks
tmp <- dr[age_group != 'unknown', list(deaths_total=sum(deaths_total)), by=c('loc_label','Date','age_group') ]
tmp2 <- tmp[, list(Date=seq(min(Date),max(Date),1)), by='loc_label']
tmp2[, DUMMY:=1L]
tmp2 <- merge(tmp2, data.table(DUMMY=1L, age_group=tmp[,unique(age_group)]), by='DUMMY', allow.cartesian=TRUE)
tmp <- merge(tmp2, tmp, by=c('loc_label','age_group','Date'), all.x=TRUE)
set(tmp, tmp[, which(is.na(deaths_total))], 'deaths_total', 0L)
set(tmp, NULL, 'DUMMY', NULL)
tmp <- tmp[, 
           list( 
             Date=Date,
             deaths_total2 = data.table::frollmean(deaths_total, 14, algo="exact", align="right")		
           ) , 
           by=c('loc_label','age_group')]
set(tmp, tmp[, which(is.na(deaths_total2))], 'deaths_total2', 0.)
tmp <- tmp[, list(age_group=age_group, deaths_total2p = deaths_total2/sum(deaths_total2)), by=c('loc_label','Date')]
tmp <- merge(tmp, subset(dr, age_group == 'unknown', -age_group), by=c('loc_label','Date'))
tmp <- subset(tmp, !is.na(deaths_total2p))
set(tmp, NULL, 'deaths_non_covid', tmp[, deaths_non_covid*deaths_total2p])
set(tmp, NULL, 'deaths_covid', tmp[, deaths_covid*deaths_total2p])
set(tmp, NULL, 'deaths_total', tmp[, deaths_total*deaths_total2p])
set(tmp, NULL, 'deaths_total2p', NULL)
dr <- rbind(subset(dr, age_group!='unknown'), tmp )

# sum deaths by our age bands
set(dr, NULL, 'age.from', dr[, as.integer(gsub('^([0-9]+)-([0-9]+)$','\\1',age_group))])
set(dr, NULL, 'age.to', dr[, as.integer(gsub('^([0-9]+)-([0-9]+)$','\\2',age_group))])
dr <- merge(dr, subset(da, select=c(Age, age.label)), by.x='age.from', by.y='Age')
dr <- dr[, 
         list(
           reg.deaths.ooh = sum(deaths_total[place!='hospital']),
           reg.coviddeaths.all = sum(deaths_covid),
           reg.deaths.all = sum(deaths_total)
         ),
         by=c('Date','loc_label','age.label')]
setkey(dr, Date)
tmp <- dr[, list(Date= seq( min(Date), max(Date), 1 )), by='loc_label']
tmp[, DUMMY:=1L]
tmp <- merge(tmp, data.table(DUMMY=1, age.label= dr[, unique(age.label)]), by='DUMMY', allow.cartesian=TRUE)
set(tmp, NULL, 'DUMMY', NULL)
dr <- merge(tmp, dr, by=c('loc_label','age.label','Date'), all.x=TRUE)
set(dr, dr[, which(is.na(reg.deaths.ooh))], 'reg.deaths.ooh', 0L)
set(dr, dr[, which(is.na(reg.deaths.all))], 'reg.deaths.all', 0L)
set(dr, dr[, which(is.na(reg.coviddeaths.all))], 'reg.coviddeaths.all', 0L)

# 7 day rolling average and excess deaths defined as daily deaths - 7 day rolling average
setkey(dr, loc_label, age.label, Date)
tmp <- dr[, 
          list( 
            Date=Date,
            reg.deaths.all.7avg = data.table::frollmean(reg.deaths.all, 7, algo="exact", align="left")		
          ) , 
          by=c('loc_label','age.label')]
tmp[, year := as.integer(strftime(Date, format='%Y'))]
tmp <- subset(tmp, year==2019)
tmp[, Date2 := tmp[, gsub('2019-','',as.character(Date))]]
tmp2 <- subset(tmp, Date2=='02-28')
set(tmp2, NULL, 'Date2', '02-29')
tmp <- rbind(tmp, tmp2)
dr[, Date2 := gsub('^[0-9]+-','',as.character(Date))]
dr <- merge(dr, subset(tmp, select=c(loc_label, age.label, Date2, reg.deaths.all.7avg)), by=c('loc_label','age.label','Date2') )
dr[, reg.deaths.all.excess := reg.deaths.all - reg.deaths.all.7avg]

# smooth excess deaths with 7 day rolling average
setkey(dr, loc_label, age.label, Date)
tmp <- dr[, 
          list( 
            Date=Date,
            reg.deaths.all.excess.7avg = data.table::frollmean(reg.deaths.all.excess, 7, algo="exact", align="left")		
          ) , 
          by=c('loc_label','age.label')]
dr <- merge(dr, tmp, by=c('loc_label','age.label','Date') )
dr <- merge(dr, subset(dpop, select=c(loc_label,age.label,pop)), by=c('loc_label','age.label'))

if(args$excess.death.zero.among.0.to.16)
{
  cat('\nsetting excess death to zero among 0 to 15')
  set(dr, dr[, which(age.label=='0-15')], c('reg.deaths.all.excess','reg.deaths.all.excess.7avg'), 0.)
}

dr[, reg.deaths.all.exorcovid := pmax(reg.coviddeaths.all,reg.deaths.all.excess)]
dr <- merge(dr, dw, by=c('Date'))

if(0 & args$make.preprocessing.plots)
{		
  p <- ggplot(dr, aes(x=Date)) +
    geom_line(aes(y=reg.deaths.all.excess, colour=age.label)) +
    scale_x_date(breaks='6 months', expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
    facet_grid(change_city_label(loc_label)~age.label) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = 'bottom',
          strip.background = element_blank()) +
    labs(x='', y='excess deaths', colour='age')
  ggsave(file=paste0(args$out.base,'_pop_excessdeaths.pdf'),p,w=12,h=14)
  ggsave(file=paste0(args$out.base,'_pop_excessdeaths.png'),p,w=12,h=14)
}
if(args$make.preprocessing.plots)
{	
  if(args$partition.locs){
    # set 1
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part1)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_line(aes(y=reg.deaths.all.excess.7avg, colour=age.label)) +
      geom_vline(data=subset(ddates,loc_label %in% args$selected.locs.part1), aes(xintercept=w2.start), colour='grey50') +
      scale_x_date(breaks='3 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      facet_grid(age.label~change_city_label(loc_label)) +	
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg_part1.png'),p,w=10,h=12)
    
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part2)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_line(aes(y=reg.deaths.all.excess.7avg, colour=age.label)) +
      geom_vline(data=subset(ddates,loc_label %in% args$selected.locs.part2), aes(xintercept=w2.start), colour='grey50') +
      scale_x_date(breaks='3 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      facet_grid(age.label~change_city_label(loc_label)) +	
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg_part2.png'),p,w=10,h=12)
    
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part3)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_line(aes(y=reg.deaths.all.excess.7avg, colour=age.label)) +
      geom_vline(data=subset(ddates,loc_label %in% args$selected.locs.part3), aes(xintercept=w2.start), colour='grey50') +
      scale_x_date(breaks='3 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      facet_grid(age.label~change_city_label(loc_label)) +	
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg_part3.png'),p,w=10,h=12)
    
    #2nd set
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part1)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part1), aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part1), aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.all.excess/pop*1e5, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths per 100,000', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5_part1.png'),p,w=10,h=12)
    
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part2)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part2), aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part2), aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.all.excess/pop*1e5, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths per 100,000', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5_part2.png'),p,w=10,h=12)
    
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part3)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part3), aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part3), aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.all.excess/pop*1e5, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths per 100,000', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5_part3.png'),p,w=10,h=12)
    
    # 3rd set
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part1)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part1), aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part1), aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.ooh, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily out of hospital deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part1.png'),p,w=10,h=12)

    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part2)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part2), aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part2), aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.ooh, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily out of hospital deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part2.png'),p,w=10,h=12)
    
    tmp <- subset(dr, Date >= "2020-01-01" & loc_label %in% args$selected.locs.part3)
    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part3), aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=subset(ddates, loc_label %in% args$selected.locs.part3), aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.ooh, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily out of hospital deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part3.png'),p,w=10,h=12)	
    
  }else{
    tmp <- subset(dr, Date >= "2020-01-01")
    p <- ggplot(tmp, aes(x=Date)) +
      geom_line(aes(y=reg.deaths.all.excess.7avg, colour=age.label)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      scale_x_date(breaks='3 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      facet_grid(age.label~change_city_label(loc_label)) +	
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_7avg.png'),p,w=10,h=12)

    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=ddates, aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.all.excess/pop*1e5, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily excess deaths per 100,000', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_excessdeaths_1e5.png'),p,w=10,h=12)

    p <- ggplot(tmp, aes(x=Date)) +
      geom_vline(data=ddates, aes(xintercept=w1.start), colour='grey50') +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +			
      geom_line(aes(y=reg.deaths.ooh, colour=age.label)) +
      scale_x_date(breaks='6 months', expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      facet_grid(age.label~change_city_label(loc_label), scales = 'free_y') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = 'bottom',
            strip.background = element_blank()) +
      labs(x='', y='daily out of hospital deaths', colour='age')
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_oohdeaths_part3.png'),p,w=10,h=12)	
  }
}

#####################################################
# aggregate death data (registry) across ages
#####################################################

cat('\naggregate death data (registry) across ages ...')

# get all-age deaths
drt <- dr[, list(reg.deaths.ooh=sum(reg.deaths.ooh), reg.coviddeaths.all=sum(reg.coviddeaths.all), reg.deaths.all=sum(reg.deaths.all)), by=c('loc_label','Date','Date2')]
setkey(drt, loc_label, Date)
tmp <- drt[, 
           list( 
             Date=Date,
             reg.deaths.all.7avg = data.table::frollmean(reg.deaths.all, 7, algo="exact", align="left")		
           ) , 
           by=c('loc_label')]
tmp[, year := as.integer(strftime(Date, format='%Y'))]
tmp <- subset(tmp, year==2019)
tmp[, Date2 := tmp[, gsub('2019-','',as.character(Date))]]
tmp2 <- subset(tmp, Date2=='02-28')
set(tmp2, NULL, 'Date2', '02-29')
tmp <- rbind(tmp, tmp2)
drt <- merge(drt, subset(tmp, select=c(loc_label, Date2, reg.deaths.all.7avg)), by=c('loc_label','Date2') )
drt[, reg.deaths.all.excess := reg.deaths.all - reg.deaths.all.7avg]
setkey(drt, loc_label, Date)
tmp <- drt[, 
           list( 
             Date=Date,
             reg.deaths.all.excess.7avg = data.table::frollmean(reg.deaths.all.excess, 7, algo="exact", align="left")		
           ) , 
           by=c('loc_label')]
drt <- merge(drt, tmp, by=c('loc_label','Date') )
drt <- merge(drt, unique(subset(dpop, select=c(loc_label,tpop))), by=c('loc_label'))
drt <- merge(drt, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
drt <- merge(drt, tmp, by=c('loc_label','week'))

# Creating Manaus specific registry excess death dataset and merging (by date) with excess burials from Bruce Walker 
manaus_registry_excess <- subset(drt, loc_label=='manaus')
saopaulo_registry <- subset(drt, loc_label == 'sao paulo')
burials <- as.data.table(read.csv(args$file.excess.burials, stringsAsFactors=FALSE))
# burials_sp <- as.data.table(read.csv('~/Downloads/Burials_SP.csv')) 
burials[, Date:= as.Date(Date)]
# burials_sp[, Date:= as.Date(Date)]
manaus_registry_excess <- merge(manaus_registry_excess, burials, by = "Date")
# saopaulo_registry <- merge(saopaulo_registry, burials_sp, by = "Date")

if(args$make.preprocessing.plots)
{	
  # Manaus
  tmp <- manaus_registry_excess[,.(Date, reg.deaths.all, burials)]
  tmp <- tmp[,  `:=` (reg.deaths.all = data.table::frollmean(reg.deaths.all, 7, algo="exact", align="left"),
                      burials = data.table::frollmean(burials, 7, algo="exact", align="left")) ]
  tmp <- subset(tmp, !is.na(reg.deaths.all))
  tmp <- subset(tmp, !is.na(burials))
  tmp <- melt(tmp, id.vars = 'Date', value.name = 'deaths')
  
  p <- ggplot(tmp, aes(x = Date, y = deaths, color = variable)) +
    geom_line() +
    scale_x_date(breaks='1 months', expand=c(0,0)) +#, limits = c(as.Date("2020-04-01"), as.Date("2021-05-01"))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = 'bottom',
          strip.background = element_blank()) +
    labs(x='', y='7-day running mean of deaths', color = '') + 
    ggsci::scale_colour_npg(labels = c('all-cause deaths reported by the Brazilian Civil Registry for Manaus',
                                    "cemetery burials reported to the Mayor's Office of Manaus")) +
    guides(colour = guide_legend(nrow=2,byrow=TRUE) )
  ggsave(file=paste0(args$out.base,'_pop_deaths_reg_vs_burial_comp.pdf'),p,w=8,h=6)
  ggsave(file=paste0(args$out.base,'_pop_deaths_reg_vs_burial_comp.png'),p,w=8,h=6)
}
if(0)
{
  prop_res_nonres <- .83 # Derived in other script from dht (which is later)
  # prop_res_nonres <- prop_res_nonres * 0.88 # This is assuming ...
  # Sao Paulo. 
  tmp <- saopaulo_registry[,.(Date, reg.deaths.all, Deaths,Burials)]
  tmp <- tmp[,  `:=` (reg.deaths.all = data.table::frollmean(reg.deaths.all, 7, algo="exact", align="left"),
                      Deaths = data.table::frollmean(Deaths, 7, algo="exact", align="left"),
                      Burials = data.table::frollmean(Burials, 7, algo="exact", align="left")) ]
  tmp <- subset(tmp, !is.na(reg.deaths.all) & !is.na(Deaths) & !is.na(Burials))
  tmp <- melt(tmp, id.vars = 'Date', value.name = 'deaths')
  # tmp1 <- tmp[, list(deaths = sum(deaths)), by = 'variable']
  # burial_rate <- tmp1[variable == 'Burials', deaths]/tmp1[variable == 'Deaths', deaths]
  tmp[!variable == 'Burials', deaths := deaths * burial_rate * prop_res_nonres]
  tmp <- tmp[!variable == 'Deaths', ]
  
  p <- ggplot(tmp, aes(x = Date, y = deaths, color = variable)) +
    geom_line() +
    scale_x_date(breaks='1 months', expand=c(0,0)) +#, limits = c(as.Date("2020-04-01"), as.Date("2021-05-01"))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = 'bottom',
          strip.background = element_blank()) +
    labs(x='', y='7-day running mean of deaths', color = '') + 
    ggsci::scale_colour_npg(labels = c('all-cause deaths reported by the Brazilian Civil Registry for Sao Paulo x prop_residents_deaths x burials_rate',
                                       # 'all-cause deaths reported by the Funeral Service of the Municipality of São Paulo x prop_residents_deaths',
                                       "burials reported by the Funeral Service of the Municipality of São Paulo")) +
    guides(colour = guide_legend(nrow=3,byrow=TRUE) )
  ggsave(file=paste0(args$out.base,'_pop_deaths_reg_vs_burial_comp_3.pdf'),p,w=8,h=6)
  ggsave(file=paste0(args$out.base,'_pop_deaths_reg_vs_burial_comp_3.png'),p,w=8,h=6)
  
  # curious: what s the ratio between total civil registry deaths & Funeral Service Deaths
  tmp1 <- tmp[, list(deaths = sum(deaths)), by = 'variable']
  tmp1[variable == 'Deaths', deaths]/tmp1[variable == 'reg.deaths.all', deaths]
  
  # As if only 88% of deaths reported to the registry were also reported to the Funeral Service
  # Is there a high proportion of cremations? People that can t afford funerals?
}

if(0 & args$make.preprocessing.plots)
{		
  p <- ggplot(drt, aes(x=Date)) +
    geom_line(aes(y=reg.deaths.all.excess, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='3 months', expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggsci::scale_colour_npg() +
    facet_wrap(~change_city_label(loc_label), ncol=2, scales='free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = 'bottom',
          strip.background = element_blank()) +
    labs(x='', y='daily excess deaths', colour='cities') 
  ggsave(file=paste0(args$out.base,'_pop_excessdeathsallage.pdf'),p,w=8,h=8)
  ggsave(file=paste0(args$out.base,'_pop_excessdeathsallage.png'),p,w=8,h=8)
}
if(args$make.preprocessing.plots)
{							
  p <- ggplot(drt, aes(x=Date)) +
    geom_vline(data=ddates, aes(xintercept=w1.start), colour='grey50') +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +						
    geom_line(aes(y=reg.deaths.all.excess.7avg, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='4 months', expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(drt$loc_label)))) +
    facet_wrap(~change_city_label(loc_label), ncol=5, scales='free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = 'bottom',
          strip.background = element_blank()) +
    labs(x='', y='daily excess deaths', colour='cities') +
    guides(color='none')
  ggsave(file=paste0(args$out.base,'_pop_excessdeathsallage_7avg.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_pop_excessdeathsallage_7avg.png'),p,w=10,h=10)
	
  p <- ggplot(drt, aes(x=Date)) +
    geom_vline(data=ddates, aes(xintercept=w1.start), colour='grey50') +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +									
    geom_line(aes(y=reg.deaths.all.excess/tpop*1e5, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='4 months', expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(drt$loc_label)))) +
    facet_wrap(~change_city_label(loc_label), ncol=5, scales='free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = 'bottom',
          strip.background = element_blank()) +
    labs(x='', y='daily excess deaths per 100,000', colour='cities') +
    guides(color='none')
  ggsave(file=paste0(args$out.base,'_pop_excessdeathsallage_1e5.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_pop_excessdeathsallage_1e5.png'),p,w=10,h=10)

  p <- ggplot(drt, aes(x=Date)) +
    geom_vline(data=ddates, aes(xintercept=w1.start), colour='grey50') +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +												
    geom_line(aes(y=reg.deaths.ooh, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='4 months', expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(drt$loc_label)))) +
    facet_wrap(~change_city_label(loc_label), ncol=5, scales='free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = 'bottom',
          strip.background = element_blank()) +
    labs(x='', y='daily out of hospital deaths', colour='cities') +
    guides(color='none')
  ggsave(file=paste0(args$out.base,'_pop_oohdeathsallage.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_pop_oohdeathsallage.png'),p,w=10,h=10)	
}

tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dr <- merge(dr, tmp, by=c('loc_label','week'))
drt <- merge(drt, tmp, by=c('loc_label','week'))

dr <- subset(dr, Date>='2020-01-01', select=-c(Date2, reg.deaths.all.7avg))
drt <- subset(drt, Date>='2020-01-01', select=-c(Date2, reg.deaths.all.7avg))

####################################################################
# READ SIVEP AGE AND HOSPITALISATION DATA
####################################################################

cat('\nread SIVEP age and hospitalisation data ...')

# read hospitalisation data with symptom onset date, admission date, death date
dht <- as.data.table(readRDS(args$file.SIVEP))

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

if(args$remove.unclassified.cases) # By unclassified we mean COVID19 cases with CLASSI_FIN == 4 
{
  set(dht, dht[, which(!is.na(CLASSI_FIN) & CLASSI_FIN%in%c(5))], 'adm.cat', 'COVID19 susp/conf')
  set(dht, dht[, which(!is.na(CLASSI_FIN) & CLASSI_FIN%in%c(1,2,3,4))], 'adm.cat', 'other')
}else{
  set(dht, dht[, which(!is.na(CLASSI_FIN) & CLASSI_FIN%in%c(4,5))], 'adm.cat', 'COVID19 susp/conf')
  set(dht, dht[, which(!is.na(CLASSI_FIN) & CLASSI_FIN%in%c(1,2,3))], 'adm.cat', 'other')
}
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
dht <- subset(dht, loc_label %in% args$selected.locs)

if(args$make.preprocessing.plots)
{
  
  tmp <- dht[hosp.cat=='in hospital' & !is.na(DateAdmission),]
  tmp[, CLASSI_FIN := as.character(CLASSI_FIN)]
  tmp[CLASSI_FIN %in% c('1','2','3'), CLASSI_FIN := '1-3']
  tmp[is.na(CLASSI_FIN), CLASSI_FIN := 'unknown']
  tmp <- merge(tmp,dw[,.(Date, week.start)], by.x='DateAdmission', by.y='Date')
  
  tmp <- tmp[, list(.N),by=c('loc_label','week.start','CLASSI_FIN')]
  
  p <- ggplot(tmp, aes(x= week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_bar(aes(y=N, fill=CLASSI_FIN ), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_x_date(breaks='2 months') +
    scale_y_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +
    labs(x='weeks', y='Weekly COVID-19 attributable hospital admissions',
         fill='SIVEP Classification')
  
  ggsave(file=paste0(args$out.base,'_hospn_by_classification.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_classification.png'),p,w=10,h=12, limitsize = F)
  
  
  tmp <- dht[hosp.cat=='in hospital' & !is.na(DateAdmission) & !is.na(DateDeath),]
  tmp <- merge(tmp, ddates[,.(loc_label, w1.start, w2.end)], by ='loc_label')
  tmp <- tmp[DateAdmission >= w1.start & DateAdmission <= w2.end & grepl(args$keep.vac.type,vac_type)]
  tmp[, DateDiff := as.numeric(DateDeath - DateAdmission) ]
  tmp[, lengthCat := 'More than 3 weeks']
  tmp[DateDiff <= 7, lengthCat := '1 week or less']
  tmp[DateDiff <= 14 & DateDiff > 7,  lengthCat := '1 to 2 weeks']
  tmp[DateDiff <= 21 & DateDiff > 14,  lengthCat := '2 to 3 weeks']
  tmp <- merge(tmp,dw[,.(Date, week.start)], by.x='DateAdmission', by.y='Date')
  tmp <- tmp[, list(.N),by=c('loc_label','week.start','lengthCat')]
  tmp[,tot:= sum(N),by=c('loc_label','week.start')][, P:=N/tot]
  tmp[, lengthCat:= factor(lengthCat, levels = c('1 week or less', '1 to 2 weeks', '2 to 3 weeks', 'More than 3 weeks' ))]
  
  tmp1 <- dht[hosp.cat=='in hospital' & !is.na(DateAdmission) & !is.na(DateDeath),]
  tmp1 <- merge(tmp1, ddates[,.(loc_label, w1.start, w2.end)], by ='loc_label')
  tmp1 <- tmp1[DateAdmission >= w1.start & DateAdmission <= w2.end & grepl(args$keep.vac.type,vac_type)]
  tmp1 <- merge(tmp1,dw[,.(Date, week.start)], by.x='DateAdmission', by.y='Date')
  tmp1 <- tmp1[, list(.N),by=c('loc_label','week.start')]
  tmp1[, P := N/max(N), by='loc_label']
  
  p <- ggplot(tmp, aes(x= week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_bar(aes(y=P, fill=lengthCat), alpha=0.8, stat='identity', position=position_stack()) +
    geom_line(data=tmp1, aes(y=P)) + 
    ggsci::scale_fill_npg() +
    scale_x_date(breaks='2 months') +
    scale_y_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +
    labs(x='weeks', y='Weekly COVID-19 attributable hospital admissions',
         fill='Difference from time of hospital admission to death')
  
  ggsave(file=paste0(args$out.base,'_hosps_h2d_overtime.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hosps_h2d_overtime.png'),p,w=10,h=12, limitsize = F)

  # tmp <- merge(dht, da[, .(age.label, Age)], by='Age')
  # tmp <- tmp[!is.na(DateICUstart), {z <- table(age.label); list(variable= names(z), value=as.numeric(z)) } ,by='loc_label']
  # tmp <- dcast(tmp, loc_label ~variable, value.var = 'value')
  
  tmp <- dht[hosp.cat=='in hospital']
  tmp <- tmp[loc.res.cat == 'resident' & grepl(args$keep.vac.type, vac_type)]
  tmp <- subset(tmp, loc_label %in% args$selected.locs, select = c('loc_label', 'DateAdmission', 'CLASSI_FIN'))
  tmp <- merge(tmp, ddates[, .(loc_label, w1.start, w2.end)], by='loc_label')
  tmp <- merge(tmp, dw[, .(Date, week.start)], by.x='DateAdmission', by.y='Date')
  tmp[, CLASSI_FIN:=as.character(CLASSI_FIN)]
  tmp[is.na(CLASSI_FIN), CLASSI_FIN := 'Unreported']
  tmp[CLASSI_FIN %in% c('1','2','3'),CLASSI_FIN := '1-3',]
  tmp[week.start >= w1.start & week.start <= w2.end]
  tmp <- tmp[,list(.N, w2.end=unique(w2.end)),by=c('loc_label','week.start','CLASSI_FIN')]
  tmp[, classi := 'unreported']
  tmp[CLASSI_FIN == '1-3', classi := 'other cause']
  tmp[CLASSI_FIN == '4', classi:= 'suspected COVID-19']
  tmp[CLASSI_FIN == '5', classi:= 'confirmed COVID-19']
  
  p <- ggplot(tmp[week.start <= w2.end], aes(x= week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    # geom_vline(data=ddates, aes(xintercept=w1.start), colour='grey50') +
    geom_bar(aes(y=N, fill=CLASSI_FIN ), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_x_date(breaks='2 months', labels=scales:::date_format("%b %y") ) +
    scale_y_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +
    labs(x='', y='Weekly  hospital admissions',
         fill='SARI classification')
  
  ggsave(file=paste0(args$out.base,'_hospn_by_SRAG_class.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_SRAG_class.png'),p,w=10,h=12, limitsize = F)

  
  # .f <- function(loc)
  # {
  #   tmp <- dht[loc_label == loc & adm.cat %in% c('COVID19 susp/conf', 'unknown') & hosp.cat == 'in hospital']
  #   setkey(tmp, DateAdmission)
  #   i <- dpop[loc_label == loc, as.integer(2.5*unique(tpop)/100000)]
  #   c(loc, as.character(tmp[i, DateAdmission]), as.character(ddates[loc_label == loc, unique(w1.start)]))
  #   # paste0(tmp[i, DateAdmission], ' ---- ', ddates[loc_label == loc, unique(w1.start)])
  # }
  # v <- sapply(args$selected.locs,.f); v <- t(v)
  
  tmp <- dht[hosp.cat=='in hospital', .(loc_label, DateAdmission, Outcome2)]
  tmp <- merge(tmp, ddates, by='loc_label')
  tmp <- tmp[DateAdmission <= w2.end & DateAdmission >= w1.start]
  tmp[, round(100*sum(Outcome2 == 'unknown')/.N, 2), by='loc_label']
  
  # V1: GO: 4.8%;  JP 23.49%; NAT 14.43
  # V2: GO: 2.7%;  JP 22.78%; NAT 13.3%
  
  
  
}

if(args$make.sexpr)
{
  # Make table for SIVEP section
  tmp <- stringr::str_sub(args$file.SIVEP,nchar(args$file.SIVEP)-3) # get the file type
  tmp <- as.Date(stringr::str_match(args$file.SIVEP, paste0("hospital_(.*?)", tmp))[2], '%d-%m-%Y') # get the date of processing
  sexpr[['SIVEP']][['release_date']] <- write.date(tmp)
  
  tmp <- dht[!is.na(DateAdmission) & hosp.cat!='out of hospital',range(DateAdmission),]
  sexpr[['SIVEP']][['first_admission_date']] <- write.date(tmp[1])
  sexpr[['SIVEP']][['last_admission_date']] <- write.date(tmp[2])
  
  tmp <- dht[hosp.cat=='in hospital']
  tmp <- subset(tmp, loc_label %in% args$selected.locs, select = c('loc_label', 'DateAdmission', 'CLASSI_FIN'))
  tmp <- merge(tmp, ddates[, .(loc_label, w1.start, w2.end)], by='loc_label')
  
  tmp[, CLASSI_FIN:=as.character(CLASSI_FIN)]
  tmp[is.na(CLASSI_FIN), CLASSI_FIN := 'Unreported']
  tmp[CLASSI_FIN %in% c('1','2','3'),CLASSI_FIN := '1-3',]
  # tmp2 <- tmp[ !is.na(DateAdmission), list(Start = write.date(min(DateAdmission)), End = write.date(max(DateAdmission))) , by = 'loc_label']
  tmp <- dcast(tmp, loc_label ~ CLASSI_FIN, fun.aggregate = length)
  tmp[, Total := as.integer(rowSums(tmp[,-1,])),]
  tmp[, `:=` (`1-3`=format(round(`1-3`*100/Total,2), nsmall=2),
             `4`=format(round(`4`*100/Total,2),nsmall=2),
             `5`=format(round(`5`*100/Total,2),nsmall=2),
             Unreported=format(round(Unreported*100/Total,2), nsmall=2))  ]
  tmp[, `:=`(`1-3`=paste0(`1-3`, '%'),`4`= paste0(`4`, '%'), `5`=paste0(`5`, '%'), Unreported=paste0(`Unreported`, '%'))  ]
  #tmp <-merge(tmp2,tmp, by = 'loc_label')
  setcolorder(tmp, c('loc_label', 'Total', '1-3','4','5', 'Unreported'))
  setnames(tmp, 'loc_label', 'Cities')
  setnames(tmp, c('1-3','4','5'), c('1-3:Other resp.pathogen', '4:Suspected COVID-19', '5:Confirmed COVID-19'))
  tmp[, Cities := change_city_label(Cities)]
  sexpr[['SIVEP']][['composition_classes_table']] <- tmp
}

if(args$deaths.unknown2discharged){
  cat('\nsetting unknown reported outcomes to discharged ...')
  dht[Outcome2 == 'unknown', Outcome2:= 'discharged']
}
dha <- copy(dht)

if(args$keep.CLASSIFIN.na) # THIS SHOULD BE 0
{
  cat('\n Keeping COVID19 susp/conf + CLASSI_FIN = NA in dht \n')
  dht <- subset(dht, adm.cat=='COVID19 susp/conf' | adm.cat=='unknown')
  dht <- subset(dht, Outcome2 != 'Death_Other') # 10 people died of other causes
  
}else{
  cat('\n Only Keeping COVID19 susp/conf in dht \n')
  dht <- subset(dht, adm.cat=='COVID19 susp/conf')
}

stopifnot(nrow(dht[Outcome2 == 'Death_Other']) == 0)

if(args$make.preprocessing.plots){

  # Plot COVID admissions by residents and non-residents
  # and by vaccination status
  tmp <- dht[hosp.cat == 'in hospital' & !is.na(DateAdmission), ]
  setnames(tmp, 'DateAdmission', 'Date')
  tmp[, vac_type := gsub('^vaccinated(.*?)$', 'vaccinated with at least one dose', vac_type)]
  tmp <- merge(tmp, dw[, .(Date, week.start)], by = 'Date')
  tmp <- tmp[, .N, by = c('loc_label','week.start', 'vac_type', 'loc.res.cat') ]
  tmp[, loc.res.cat := factor(loc.res.cat, levels = c('not resident', 'resident')) ]
  tmp[, vac_type := factor(vac_type, c('vaccinated with at least one dose', 'unknown', 'not vaccinated'))]
  tmp <- merge(tmp, ddates[,.(loc_label, w1.start, w2.end)], by = 'loc_label')
  tmp <- tmp[week.start >= w1.start]
  tmp <- tmp[week.start <= w2.end]
  tmp1 <- tmp[loc.res.cat == 'resident' & vac_type %in%  c('unknown', 'not vaccinated'), list(Nres = sum(N)), by = c('week.start', 'loc_label', 'vac_type')]
  tmp1 <- tmp1[, list(Nres=sum(Nres)), by = c('week.start', 'loc_label')]

  p <- ggplot(tmp[week.start <= ddates$w2.end[1], ], aes(x= week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_bar(aes(y=N, alpha = loc.res.cat, fill=vac_type ), stat='identity', position=position_stack()) +
    geom_step(data=tmp1,aes(x = week.start - 3.5,y = Nres)) +
    ggsci::scale_fill_npg() +
    scale_x_date(breaks='2 months') +
    scale_y_continuous(expand=c(0,0)) +
    scale_alpha_manual(values=c('not resident'=0.5, 'resident'=1.)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +
    labs(x='', y='Weekly COVID-19 attributable hospital admissions',
         fill='Vaccination status', alpha = 'Residence status')

  ggsave(file=paste0(args$out.base,'_hospn_by_residence_vaxstatus.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_residence_vaxstatus.png'),p,w=10,h=12, limitsize = F)
  
  # Plot COVID-19 admissions that are resident, 
  # and by outcome
  tmp <- dht[hosp.cat == 'in hospital' & !is.na(DateAdmission) 
             & loc.res.cat == 'resident' & grepl(args$keep.vac.type,vac_type), ]
  setnames(tmp, 'DateAdmission', 'Date')
  # tmp[, vac_type := gsub('^vaccinated(.*?)$', 'vaccinated', vac_type)]
  tmp <- merge(tmp, dw[, .(Date, week.start)], by = 'Date')
  tmp <- tmp[, .N, by = c('loc_label','week.start', 'Outcome2') ]
  tmp[, Outcome2 := gsub('_', ' ', Outcome2)]
  tmp[Outcome2 == 'unknown', Outcome2 := 'Unknown']
  tmp[Outcome2 == 'discharged', Outcome2 := 'Discharged']
  
  p <- ggplot(tmp, aes(x= week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_bar(aes(y=N, fill = Outcome2), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_x_date(breaks='2 months') +
    scale_y_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +	
    labs(x='weeks', y='Weekly COVID-19 attributable hospital admissions',
         fill='Outcome')
  
  ggsave(file=paste0(args$out.base,'_hospn_by_resoutcome.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_resoutcome.png'),p,w=10,h=12, limitsize = F)
  
}

####################################################################
# Compare registry deaths and SIVEP deaths (esp for RdJ)
####################################################################

if(0)
{
  # Do not take into account residents, as there is no diff in civil registry

  tmp2 <- dht[!is.na(DateDeath),.N, by = c('loc_label', 'DateDeath')]
  tmp2 <- merge(tmp2, dw[, .(Date, week.start)], by.x = 'DateDeath', by.y='Date')
  tmp2 <- tmp2[, list(SIVEPdeaths=sum(N)), by=c('loc_label', 'week.start')]
  
  tmp <- drt[,list(reg.deaths.ooh=sum(reg.deaths.ooh),
                   reg.coviddeaths.all=sum(reg.coviddeaths.all)),
             by=c('loc_label','week.start')]
  
  tmp <- merge(tmp, tmp2, by=c('loc_label', 'week.start'))
  tmp <- melt(tmp, id.vars = c('loc_label', 'week.start'))
  
  ggplot(tmp[variable %in% c('SIVEPdeaths', 'reg.coviddeaths.all')], aes(x=week.start, color = variable)) +
    geom_line(aes(y = value)) + 
    facet_wrap(~change_city_label(loc_label))
    
}

####################################################################
# Exploratory plots
####################################################################


# plot SARI hosp admissions by residence
if(args$make.preprocessing.plots)
{	
  tmp <- subset(dha, !is.na(DateAdmission) & hosp.cat!='out of hospital')
  tmp[, week := as.integer(strftime(DateAdmission, format='%V'))]	
  tmp2 <- tmp[,which(DateAdmission>'2021-01-03')]
  set(tmp, tmp2, 'week',tmp[tmp2,week+53L])	
  tmp <- tmp[, list(n=length(Age)), by=c('week','loc_label','loc.res.cat')]	
  p <- ggplot(tmp, aes(x= week)) +
    geom_bar(aes(y=n, fill=loc.res.cat), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol= args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +	
    labs(x='weeks', y='Hospital admissions', fill='patient residence')
  ggsave(file=paste0(args$out.base,'_hospn_by_residence.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_residence.png'),p,w=10,h=12, limitsize = F)
  
  # plot SARI hosp admissions by residence
  tmp <- subset(dha, !is.na(DateAdmission) & hosp.cat!='out of hospital')
  tmp[, week := as.integer(strftime(DateAdmission, format='%V'))]	
  tmp2 <- tmp[,which(DateAdmission>'2021-01-03')]
  set(tmp, tmp2, 'week',tmp[tmp2,week+53L])	
  tmp <- tmp[, list(n=length(Age)), by=c('week','loc_label','adm.cat')]
  set(tmp, which(tmp$adm.cat == 'other'), 'adm.cat', "SARI cases other than COVID-19")
  tmp[, adm.cat := factor(tmp$adm.cat, levels = c('COVID19 susp/conf',  "SARI cases other than COVID-19", "unknown"))]
  p <- ggplot(tmp, aes(x= week)) +
    geom_bar(aes(y=n, fill=adm.cat), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +	
    labs(x='weeks', y='Hospital admissions', fill='reason hospital admission')
  ggsave(file=paste0(args$out.base,'_hospn_by_admcat.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_admcat.png'),p,w=10,h=12, limitsize = F)
  
  # plot SARI hosp admissions by clinical outcome
  tmp <- subset(dha, !is.na(DateAdmission) & hosp.cat!='out of hospital')
  tmp[, week := as.integer(strftime(DateAdmission, format='%V'))]	
  tmp2 <- tmp[,which(DateAdmission>'2021-01-03')]
  set(tmp, tmp2, 'week',tmp[tmp2,week+53L])
  tmp <- tmp[, list(n=length(Age)), by=c('week','loc_label','Outcome2')]
  p <- ggplot(tmp, aes(x= week)) +
    geom_bar(aes(y=n, fill=Outcome2), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +	
    labs(x='weeks', y='Hospital admissions', fill='outcome\nas of database closure')
  ggsave(file=paste0(args$out.base,'_hospn_by_outcome.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_outcome.png'),p,w=10,h=12, limitsize = F)
  
  # plot SARI hosp admissions among residents by clinical outcome
  tmp <- subset(dha, !is.na(DateAdmission) & loc.res.cat=='resident' & hosp.cat!='out of hospital')
  tmp[, week := as.integer(strftime(DateAdmission, format='%V'))]	
  tmp2 <- tmp[,which(DateAdmission>'2021-01-03')]
  set(tmp, tmp2, 'week',tmp[tmp2,week+53L])
  tmp <- tmp[, list(n=length(Age)), by=c('week','loc_label','Outcome2')]
  tmp[Outcome2 == 'discharged', Outcome2:= 'Discharged alive']
  tmp[Outcome2 == 'Death_Covid', Outcome2:= 'Died']
  tmp[Outcome2 == 'in_ICU', Outcome2:= 'in ICU']
  tmp[Outcome2 == 'unknown', Outcome2:= 'unknown/undetermined']
  tmp[Outcome2 == 'Death_Other', Outcome2:= 'unknown/undetermined']  
  p <- ggplot(tmp, aes(x= week)) +
    geom_bar(aes(y=n, fill=Outcome2), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol =args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +	
    labs(x='weeks', y='Weekly COVID-19 attributable hospital admissions', fill='outcome\nas of database closure')
  ggsave(file=paste0(args$out.base,'_hospn_by_outcome_residents.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_outcome_residents.png'),p,w=10,h=12, limitsize = F)
  
  tmp <- subset(dht, Outcome2 == 'Death_Covid' & loc.res.cat == 'resident' & !is.na(DateDeath))
  tmp <- tmp[,list(loc_label, DateDeath, DateAdmission,same_death_admission = as.integer(DateAdmission == DateDeath), after_admission = as.integer(DateAdmission < DateDeath)),]
  tmp[is.na(DateAdmission), `:=` (same_death_admission = 0L, after_admission = 0L)]
  tmp <- tmp[, list(same_death_admission = sum(same_death_admission), no_admission = sum(is.na(DateAdmission)), after_admission = sum(after_admission)), by = c('DateDeath', 'loc_label')]
  tmp <- melt(tmp, id.vars = c('DateDeath', 'loc_label'))
  tmp[, week := as.integer(strftime(DateDeath, format='%V'))]
  tmp2 <- min(tmp[,DateDeath])
  tmp[, week.start := tmp2  + 7 * (week - 1 ),]
  tmp <- tmp[, list(value = sum(value)),by =c('loc_label', 'week.start', 'variable')]
  tmp <- merge(tmp, ddates[,.(loc_label, w1.start)], by = 'loc_label')
  tmp[, variable := relevel(tmp$variable, 'after_admission')]
  p <- ggplot(tmp[week.start>w1.start,,], aes(x= week.start)) +
    geom_bar(aes(y=value, fill=as.factor(variable)), stat='identity', position=position_stack()) +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots ,scales='free') +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) + 
    ggsci::scale_fill_npg(labels = c("previously admitted and died  later", "admitted and died in the same date", "no hospital admission date")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position='bottom') +	
    labs(x='', y='Weekly COVID-19 attributable deaths reported to SIVEP-Gripe', fill='')  
  ggsave(file=paste0(args$out.base,'_hdeaths_by_admission_residents.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hdeaths_by_admission_residents.png'),p,w=10,h=12, limitsize = F)  
}

# make plots of in hospital fatality rates by vaccine coverage
if(args$make.preprocessing.plots)
{	  	
  dhv <- subset(dht, loc.res.cat=='resident' & hosp.cat=='in hospital')
  setnames(dhv, 'DateAdmission','Date')
  dhv <- merge(dhv,dw,by='Date')
  dhv <- merge(dhv,da,by='Age')
  
  # plot hospital admissions by detailed vaccination status
  tmp <- dhv[, list(n=length(Date)), by=c('loc_label','week','week.start','vac_type')]
  tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))
  p <- ggplot(tmp, aes(x= week.start)) +
    geom_bar(aes(y=n, fill=vac_type), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol= args$ncols.plots, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +	
    labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
  ggsave(file=paste0(args$out.base,'_hospn_by_vac2.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_vac2.png'),p,w=10,h=12, limitsize = F)
  
  # plot hospital admissions by detailed vaccination status and age
  if(args$partition.locs){
    
    # set 1
    tmp <- dhv[loc_label %in% args$selected.locs.part1, list(n=length(Date)), by=c('loc_label','week','week.start','age.label','vac_type')]
    tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position=position_stack()) +
      ggsci::scale_fill_npg() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            panel.spacing.y = unit(1, "lines"),
            strip.background = element_blank(), 
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_part1.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_part1.png'),p,w=10,h=12, limitsize = F)
    
    tmp <- dhv[loc_label %in% args$selected.locs.part2, list(n=length(Date)), by=c('loc_label','week','week.start','age.label','vac_type')]
    tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position=position_stack()) +
      ggsci::scale_fill_npg() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            panel.spacing.y = unit(1, "lines"),
            strip.background = element_blank(), 
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_part2.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_part2.png'),p,w=10,h=12, limitsize = F)
    
    tmp <- dhv[loc_label %in% args$selected.locs.part3, list(n=length(Date)), by=c('loc_label','week','week.start','age.label','vac_type')]
    tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position=position_stack()) +
      ggsci::scale_fill_npg() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            panel.spacing.y = unit(1, "lines"),
            strip.background = element_blank(), 
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_part3.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_part3.png'),p,w=10,h=12, limitsize = F)
    
    # set 2
    tmp <- dhv[loc_label %in% args$selected.locs.part1, list(n=length(Date)), by=c('loc_label','week','week.start','age.label','vac_type')]
    tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))	
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position='fill') +
      ggsci::scale_fill_npg() +
      scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions_part1.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions_part1.png'),p,w=10,h=12, limitsize = F)
    
    tmp <- dhv[loc_label %in% args$selected.locs.part2, list(n=length(Date)), by=c('loc_label','week','week.start','age.label','vac_type')]
    tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))	
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position='fill') +
      ggsci::scale_fill_npg() +
      scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions_part2.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions_part2.png'),p,w=10,h=12, limitsize = F)
    
    tmp <- dhv[loc_label %in% args$selected.locs.part3, list(n=length(Date)), by=c('loc_label','week','week.start','age.label','vac_type')]
    tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))	
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position='fill') +
      ggsci::scale_fill_npg() +
      scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions_part3.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions_part3.png'),p,w=10,h=12, limitsize = F)
  }else{
    tmp <- dhv[, list(n=length(Date)), by=c('loc_label','week','week.start','age.label','vac_type')]
    tmp <- merge(tmp, ddates[,list(week=seq(hosps.start,deaths.end,1)),by='loc_label'], by=c('loc_label','week'))
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position=position_stack()) +
      ggsci::scale_fill_npg() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            panel.spacing.y = unit(1, "lines"),
            strip.background = element_blank(), 
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2.png'),p,w=10,h=12, limitsize = F)
    
    p <- ggplot(tmp, aes(x= week.start)) +
      geom_bar(aes(y=n, fill=vac_type), stat='identity', position='fill') +
      ggsci::scale_fill_npg() +
      scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      facet_grid(age.label~change_city_label(loc_label), scales='free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position =  'bottom') +	
      labs(x='', y='COVID-19 attributable hospital admissions among residents', fill='vaccination status')
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions.pdf'),p,w=10,h=12, limitsize = F)
    ggsave(file=paste0(args$out.base,'_hospn_by_age_vac2_proportions.png'),p,w=10,h=12, limitsize = F)
  }
}

if(args$make.preprocessing.plots)
{
  dhv <- subset(dht, loc.res.cat=='resident' & hosp.cat=='in hospital')
  setnames(dhv, 'DateAdmission','Date')
  dhv <- merge(dhv,dw,by='Date')
  dhv <- merge(dhv,da,by='Age')
  
  # check if HFR among vaccinated is different 4 weeks after P1 was detected
  tmp <- subset(dhv, !Outcome2%in%c('unknown','in_ICU'))
  tmp <- merge(tmp, subset(ddates, select=c(loc_label,w2.start)),by='loc_label')
  tmp <- subset(tmp, w2.start+6*7<Date)
  tmp <- tmp[, list(hosps=length(Outcome2), deaths_of_hosps=length(which(!is.na(DateDeath)))), by=c('age.label','vac_type')]	
  tmp <- tmp[, list(hosps=hosps, deaths_of_hosps=deaths_of_hosps, value=c(deaths_of_hosps/hosps,suppressWarnings(as.numeric(prop.test(deaths_of_hosps,hosps,conf.level=0.95, correct = FALSE)$conf.int))), stat=c('M','CL','CU') ), by=c('age.label','vac_type')]
  tmp <- dcast.data.table(tmp, age.label+vac_type+hosps+deaths_of_hosps~stat, value.var='value')
  p <- ggplot(tmp, aes(x= age.label)) +
    geom_linerange(aes(ymin=CL,ymax=CU,colour=vac_type),position=position_dodge(0.8)) +
    geom_point(aes(y=M, colour=vac_type,size=hosps),position=position_dodge(0.8)) +
    ggsci::scale_colour_npg() +
    scale_y_continuous(label=scales::percent,expand=c(0,0)) +
    scale_x_discrete() +		
    theme_bw() +
    theme(legend.position= 'bottom', legend.box = "vertical", legend.justification = c(0, 0)) +	
    labs(x='', y='in-hospital fatality rate among\nhospital admissions 6 weeks after Gamma detection', colour='vaccination\nstatus', size='hospital admissions') +
    guides(colour=guide_legend(nrow=4,byrow=TRUE))
  ggsave(file=paste0(args$out.base,'_hfr_by_vac2_6wksafterP1.pdf'),p,w=10,h=7, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hfr_by_vac2_6wksafterP1.png'),p,w=10,h=7, limitsize = F)
  
  # check if HFR among vaccinated is different in 2021
  tmp <- subset(dhv, !Outcome2%in%c('unknown','in_ICU'))
  tmp <- merge(tmp, subset(ddates, select=c(loc_label,w2.start)),by='loc_label')
  tmp <- subset(tmp, Date>'2021-01-01')
  tmp <- tmp[, list(hosps=length(Outcome2), deaths_of_hosps=length(which(!is.na(DateDeath)))), by=c('age.label','vac_type')]
  tmp <- tmp[, list(hosps=hosps, deaths_of_hosps=deaths_of_hosps, value=c(deaths_of_hosps/hosps,suppressWarnings(as.numeric(prop.test(deaths_of_hosps,hosps,conf.level=0.95, correct = FALSE)$conf.int))), stat=c('M','CL','CU') ), by=c('age.label','vac_type')]
  tmp <- dcast.data.table(tmp, age.label+vac_type+hosps+deaths_of_hosps~stat, value.var='value')
  p <- ggplot(tmp, aes(x= age.label)) +
    geom_linerange(aes(ymin=CL,ymax=CU,colour=vac_type),position=position_dodge(0.8)) +
    geom_point(aes(y=M, colour=vac_type,size=hosps),position=position_dodge(0.8)) +
    ggsci::scale_colour_npg() +
    scale_y_continuous(label=scales::percent,expand=c(0,0)) +
    scale_x_discrete() +		
    theme_bw() +
    theme(legend.position= 'bottom', legend.box = "vertical", legend.justification = c(0, 0)) +	
    labs(x='', y='in-hospital fatality rate among\nhospital admissions afer January 1, 2021', colour='vaccination\nstatus', size='hospital admissions') +
    guides(colour=guide_legend(nrow=4,byrow=TRUE))
  ggsave(file=paste0(args$out.base,'_hfr_by_vac2_2021.pdf'),p,w=10,h=7, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hfr_by_vac2_2021.png'),p,w=10,h=7, limitsize = F)
  
  
  # check if HFR among vaccinated is different in 2021, for each location
  tmp <- subset(dhv, !Outcome2%in%c('unknown','in_ICU'))
  tmp <- merge(tmp, subset(ddates, select=c(loc_label,w2.start)),by='loc_label')
  tmp <- subset(tmp, Date>'2021-01-01')
  tmp <- tmp[, list(hosps=length(Outcome2), deaths_of_hosps=length(which(!is.na(DateDeath)))), by=c('loc_label','age.label','vac_type')]
  tmp <- tmp[, list(hosps=hosps, deaths_of_hosps=deaths_of_hosps, value=c(deaths_of_hosps/hosps,suppressWarnings(as.numeric(prop.test(deaths_of_hosps,hosps,conf.level=0.95, correct = FALSE)$conf.int))), stat=c('M','CL','CU') ), by=c('loc_label','age.label','vac_type')]
  tmp <- dcast.data.table(tmp, loc_label+age.label+vac_type+hosps+deaths_of_hosps~stat, value.var='value')
  p <- ggplot(tmp, aes(x= age.label)) +
    geom_linerange(aes(ymin=CL,ymax=CU,colour=vac_type),position=position_dodge(0.8)) +
    geom_point(aes(y=M, colour=vac_type,size=hosps),position=position_dodge(0.8)) +
    ggsci::scale_colour_npg() +
    scale_y_continuous(label=scales::percent,expand=c(0,0)) +
    scale_x_discrete() +	
    facet_wrap(~change_city_label(loc_label), ncol= args$ncols.plots) +
    theme_bw() +
    theme(legend.position= 'bottom', 
          legend.box = "vertical", 
          legend.justification=c(0, 0),
          strip.background = element_blank(), 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +	
    labs(x='', y='in-hospital fatality rate among\nhospital admissions afer January 1, 2021', colour='vaccination\nstatus', size='hospital admissions') +
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  ggsave(file=paste0(args$out.base,'_hfr_by_vac2_2021_byloc.pdf'),p,w=12,h=14, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hfr_by_vac2_2021_byloc.png'),p,w=12,h=14, limitsize = F)
}



####################################################################
# calculate proportion of deaths in past 15 days among known outcomes by date of hospital admission
# we do this to resolve unknown outcomes in SIVEP records by date of hospital admission 
# and add to individual-level SIVEP records
####################################################################


cat('\ncalculate proportion of deaths among known outcomes ...')
dhk <- subset(dht, !grepl('unknown|in_ICU', Outcome2) & 
	hosp.cat!='out of hospital' &
	grepl(args$keep.hospital.type, Hospital_type) &
	grepl(args$keep.vac.type, vac_type))
setnames(dhk, 'DateAdmission', 'Date')
dhk <- merge(dhk, da, by='Age')
dhk <- dhk[, 
           list(
             deaths=sum(as.numeric(!is.na(DateDeath))),
             hosps=length(DateDeath)
           ), 
           by=c('Date','age.label','loc_label')]
tmp <- dht[!is.na(DateAdmission), list(Date= seq.Date(min(DateAdmission),max(DateAdmission),1)), by='loc_label']
tmp[, DUMMY:=1L]
tmp2 <- as.data.table(expand.grid(DUMMY=1L, age.label=unique(da$age.label)))
tmp <- merge(tmp, tmp2, by='DUMMY', allow.cartesian=TRUE)
set(tmp, NULL, 'DUMMY', NULL)
dhk <- merge(tmp, dhk, all.x=TRUE, by=c('loc_label','age.label','Date'))
set(dhk, dhk[, which(is.na(deaths))], 'deaths', 0)
set(dhk, dhk[, which(is.na(hosps))], 'hosps', 0)
setkey(dhk, loc_label, age.label, Date)
dhk <- dhk[, 
           list( 
             Date=Date,
             deaths2 = data.table::frollmean(deaths, 14, algo="exact", align="right"),
             hosps2 = data.table::frollmean(hosps, 14, algo="exact", align="right")
           ) , 
           by=c('loc_label','age.label')]
dhk[, pdeaths2 := deaths2/hosps2]
set(dhk, dhk[, which(hosps2==0)], 'pdeaths2', 0.)
set(dhk, NULL, c('deaths2','hosps2'), NULL)
dhk <- merge(dhk, subset(da, select=c(age.label,Age)), by='age.label', allow.cartesian=TRUE)
set(dhk, NULL, 'age.label', NULL)
setnames(dhk, 'Date', 'DateAdmission')
dht <- merge(dht, dhk, by=c('loc_label','Age','DateAdmission'), all.x=TRUE)


####################################################################
# TIME FROM HOSPITALISATION TO DEATH
# among all hospitalised patients with a fatal outcome, regardless if resident or not
####################################################################
cat('\ntime from hospitalisation to death ...')
dt <- subset(dht, 
		!is.na(DateDeath) & 		
		hosp.cat!='out of hospital' &
		grepl(args$keep.hospital.type, Hospital_type) &
		grepl(args$keep.vac.type, vac_type)
	)
dt[, h2d := as.numeric(DateDeath - DateAdmission) ]
setnames(dt, 'DateDeath', 'Date')
dt <- merge(dt, da, by='Age')
dt <- merge(dt, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label']
dt <- merge(dt, tmp, by=c('loc_label','week'))
dta <- dt[, 
          list(
            h2d.median=as.numeric(median(h2d)), 
            h2d.il=as.numeric(quantile(h2d, p=.25)), 
            h2d.iu=as.numeric(quantile(h2d, p=.75)), 
            h2d.mean=as.numeric(mean(h2d)), 
            h2d.sd=as.numeric(sd(h2d)) 
          ), 
          by=c('loc_label', 'month','month.start','age.label')]

if(args$make.preprocessing.plots)
{	
  if(args$partition.locs){
    # Set 1
    tmp <- dt[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=Date, y=h2d)) +	
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age') +
      coord_cartesian(ylim=c(0,75))
    ggsave(file=paste0(args$out.base,'_timeh2d_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_part1.png'),p,w=10,h=12)
    
    tmp <- dt[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=Date, y=h2d)) +	
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age') +
      coord_cartesian(ylim=c(0,75))
    ggsave(file=paste0(args$out.base,'_timeh2d_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_part2.png'),p,w=10,h=12)
    
    tmp <- dt[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=Date, y=h2d)) +	
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age') +
      coord_cartesian(ylim=c(0,75))
    ggsave(file=paste0(args$out.base,'_timeh2d_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_part3.png'),p,w=10,h=12)
    
    # set 2
    tmp <- dta[loc_label %in% args$selected.locs.part1,]
    p <- ggplot(tmp, aes(x=month.start, y=h2d.median)) +
      geom_errorbar(aes(ymin=h2d.il, ymax=h2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,30)) +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_timeh2d_summary_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_summary_part1.png'),p,w=10,h=12)
    
    tmp <- dta[loc_label %in% args$selected.locs.part2,]
    p <- ggplot(tmp, aes(x=month.start, y=h2d.median)) +
      geom_errorbar(aes(ymin=h2d.il, ymax=h2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,30)) +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_timeh2d_summary_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_summary_part2.png'),p,w=10,h=12)  
    
    tmp <- dta[loc_label %in% args$selected.locs.part3,]
    p <- ggplot(tmp, aes(x=month.start, y=h2d.median)) +
      geom_errorbar(aes(ymin=h2d.il, ymax=h2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,30)) +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_timeh2d_summary_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_summary_part3.png'),p,w=10,h=12)  
  }else{
    tmp <- copy(dt)
    p <- ggplot(tmp, aes(x=Date, y=h2d)) +	
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age') +
      coord_cartesian(ylim=c(0,75))
    ggsave(file=paste0(args$out.base,'_timeh2d.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d.png'),p,w=10,h=12)
    
    tmp <- copy(dta)
    p <- ggplot(tmp, aes(x=month.start, y=h2d.median)) +
      geom_errorbar(aes(ymin=h2d.il, ymax=h2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,30)) +
      labs(x='', y='time from hospitalisation to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_timeh2d_summary.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_summary.png'),p,w=10,h=12)  
  }
  
}

####################################################################
# ESTIMATE TIME FROM HOSPITALISATION TO DEATH DISTRIBUTION
# among all hospitalised patients with a fatal outcome, regardless if resident or not
####################################################################

cat('\nestimate time from hospitalisation to death distribution ...')
set.seed(12)
tmp <- ddates[, 
	{
    	z1 <- seq.Date(w1.start, w1.end, by='day')
        z2 <- seq.Date(w2.start, w2.end, by='day')
        list(Date=c(z1,z2))
	}, 
	by='loc_label']
tmp <- unique(tmp)
dt0 <- merge(dt, tmp, by = c('Date', 'loc_label'))
dt0s <- dt0[, list(value=quantile(h2d, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('loc_label', 'age.label')]
dt0s <- dcast.data.table(dt0s, loc_label+age.label~stat, value.var='value')
dt0[, h2dc := h2d + runif(nrow(dt0), 0, 1)]

if(0 & args$make.preprocessing.plots  )
{	
  # could also make this 3 different plots
  p <- ggplot(dt0s, aes(x=age.label)) +
    geom_boxplot(aes(ymin=CL, lower=IL, middle=M, upper=IU, ymax=CU, fill=loc_label), stat='identity') +
    scale_fill_manual(values = args$city.palette(length(unique(dt0s$loc_label))))+
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +    	
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'bottom') +
    labs(x='', y='days from hospital admission to \nCOVID-19 attributable death', fill='') +
    coord_cartesian(ylim=c(0,75)) 
  ggsave(file=paste0(args$out.base,'_timeh2d_age_boxplots.pdf'),p,w=12,h=10)
  ggsave(file=paste0(args$out.base,'_timeh2d_age_boxplots.png'),p,w=12,h=10)
  
  for(m in unique(dt0$loc_label)){
    # find distributions to fit to time from hospitalisation to death in each age band
    pdf(file=paste0(args$out.base,'_timeh2d_wildtype_fits_', stringr::str_replace_all(m, stringr::fixed(" "), "_"), '.pdf'),w=8,h=50)
    dt0_m = subset(dt0, loc_label == m)
    par(mfrow = c(dt0_m[, length(levels(age.label))], 1))
    for(x in dt0_m[, levels(age.label)])
    {	
      fg <- fitdistrplus::fitdist(dt0_m[age.label==x,h2dc], "gamma")
      fln <- fitdistrplus::fitdist(dt0_m[age.label==x,h2dc], "lnorm")
      fw <- fitdistrplus::fitdist(dt0_m[age.label==x,h2dc], "weibull")
      #fll <- fitdist(dt0[age.label==x,h2dc], "llogis", start = list(shape = 1, scale = 1))
      fp <- fitdist(dt0_m[age.label==x,h2dc], "pareto", start = list(shape = 1, scale = 1))	
      fitdistrplus::qqcomp(list(fw, fln, fg, fp), xlogscale = TRUE, ylogscale = TRUE, main=x)	
    }
    dev.off()
  }
  
  #	--> we choose gamma
}

#	get gamma parameters and prob that admission occurred on day 0,1,2,.. before date of death
#	get weibull parameters and prob that symptom onset occurred on day 0,1,2,.. before date of death
dtd <- as.data.table(expand.grid(age.label = dt0[, levels(age.label)], loc_label = dt0[, unique(loc_label)], shape = NA_real_, rate = NA_real_))
Age_label <- dt0[, levels(age.label)]
for(i in seq_along(Age_label))
{	
  for(m in dt0[, unique(loc_label)]){
    fg <- fitdistrplus::fitdist(dt0[age.label==Age_label[i] & loc_label == m,h2dc], "gamma")
    estimate =  unname(fg$estimate)
    dtd[age.label==Age_label[i] & loc_label == m,shape := estimate[1]] 
    dtd[age.label==Age_label[i] & loc_label == m,rate := estimate[2]] 
  }
  
}

dtd[, age.label := factor(age.label, levels = Age_label)]
dtd[, mean := shape/rate]
dtd[, upper.lim := qgamma(1-1e-2, shape=shape, rate=rate)]
tmp <- dtd[,{
  z <- pgamma(1:(args$select.h2d.max+1), shape=shape, rate=rate)
  z <- c(z[1], diff(z))
  z <- z/sum(z)
  list(day=0:args$select.h2d.max, p=z)
}, 	 
by= c('loc_label', 'age.label')]
dtd <- merge(dtd, tmp, by=c('loc_label', 'age.label'))

if(args$make.preprocessing.plots)
{		
  p <- ggplot(dtd, aes(x=day)) +
    geom_point(aes(y=p, colour=age.label)) +
    viridis::scale_colour_viridis(option='D', end=0.8, discrete = TRUE) +	
    theme_bw() +
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.background = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          legend.position='bottom') +		
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots, scales='free_y') +
    labs(y='probability mass function', x='time from hospitalisation to death (days)', colour='age')
  ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf.pdf'),p,w=10,h=12)
  ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf.png'),p,w=10,h=12)
  
  
  dtde <- dt0[, list(n=length(Date)), by=c('loc_label', 'age.label','h2d')]
  dtde <- merge(dtde, dtde[, list(total=sum(n)),by=c('loc_label', 'age.label')], by=c('loc_label', 'age.label'))
  dtde[,p:=n/total]
  dtde[,stat:='empirical']
  setnames(dtde, 'h2d','day')
  set(dtde, NULL, c('n','total'), NULL)
  tmp <- subset(dtd, select=c(loc_label,age.label,day,p))
  tmp[, stat:='Gamma']
  dtde <- rbind(dtde, tmp)
  
  if(args$partition.locs){
    tmp <- dtde[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='days from hospital admission to COVID-19 attributable death', fill='') +
      scale_fill_discrete(labels = c("empirical data", "estimated gamma distribution"))
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical_part1.png'),p,w=10,h=12)
    
    tmp <- dtde[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='days from hospital admission to COVID-19 attributable death', fill='') +
      scale_fill_discrete(labels = c("empirical data", "estimated gamma distribution"))
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical_part2.png'),p,w=10,h=12)
    
    tmp <- dtde[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='days from hospital admission to COVID-19 attributable death', fill='') +
      scale_fill_discrete(labels = c("empirical data", "estimated gamma distribution"))
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical_part3.png'),p,w=10,h=12)
  }else{
    tmp <- copy(dtde)
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='days from hospital admission to COVID-19 attributable death', fill='') +
      scale_fill_discrete(labels = c("empirical data", "estimated gamma distribution"))
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_timeh2d_age_discretized_pmf_vs_empirical.png'),p,w=10,h=12)
  }
}

####################################################################
# TIME FROM SYMPTOM ONSET TO DEATH
# among all symptomatic patients with a fatal outcome, regardless if resident or not
# note: this INCLUDES out of hosp admissions
####################################################################

cat('\ntime from symptom onset to death ...')
ds <- subset(dht, 
		!is.na(DateDeath) & 
		!is.na(DateSymptoms) &
		grepl(args$keep.hospital.type, Hospital_type) &
		grepl(args$keep.vac.type, vac_type)
	)

ds[, s2d := as.numeric(DateDeath-DateSymptoms) ]	  
setnames(ds, c('DateDeath'), c('Date'))
ds <- merge(ds, da, by='Age')
ds <- merge(ds, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label']
ds <- merge(ds, tmp, by=c('loc_label','week'))

dsa <- ds[, 
	list(
    	s2d.median=as.numeric(median(s2d)), 
        s2d.il=as.numeric(quantile(s2d, p=.25)), 
        s2d.iu=as.numeric(quantile(s2d, p=.75)), 
        s2d.mean=as.numeric(mean(s2d)), 
        s2d.sd=as.numeric(sd(s2d)) 
	), 
    by=c('loc_label', 'month','month.start','age.label')]

if(args$make.preprocessing.plots)
{	
  if(args$partition.locs){
    # set 1
    tmp <- ds[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=Date, y=s2d)) +	
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +    
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,100)) +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_part1.png'),p,w=10,h=12)
    
    tmp <- ds[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=Date, y=s2d)) +	
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +    
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,100)) +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_part2.png'),p,w=10,h=12)
    
    tmp <- ds[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=Date, y=s2d)) +	
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +    
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,100)) +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_part3.png'),p,w=10,h=12)
    
    # set 2
    tmp <- dsa[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=month.start, y=s2d.median)) +
      geom_errorbar(aes(ymin=s2d.il, ymax=s2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d_summary_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_summary_part1.png'),p,w=10,h=12)
    
    tmp <- dsa[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=month.start, y=s2d.median)) +
      geom_errorbar(aes(ymin=s2d.il, ymax=s2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d_summary_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_summary_part2.png'),p,w=10,h=12)
    
    tmp <- dsa[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=month.start, y=s2d.median)) +
      geom_errorbar(aes(ymin=s2d.il, ymax=s2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d_summary_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_summary_part3.png'),p,w=10,h=12)
  }else{
    tmp <- copy(ds)
    p <- ggplot(tmp, aes(x=Date, y=s2d)) +	
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(colour=age.label)) +
      geom_smooth(lwd=1, colour='black') +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +    
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      coord_cartesian(ylim=c(0,100)) +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d.png'),p,w=10,h=12)
    
    tmp <- copy(dsa)
    p <- ggplot(tmp, aes(x=month.start, y=s2d.median)) +
      geom_errorbar(aes(ymin=s2d.il, ymax=s2d.iu)) +
      geom_point(aes(colour=age.label)) +		
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='time from symptom onset to death\nmeasured backwards from date of death', colour='age')
    ggsave(file=paste0(args$out.base,'_times2d_summary.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_summary.png'),p,w=10,h=12)
  }
}

####################################################################
# ESTIMATE TIME FROM SYMPTOM ONSET TO DEATH DISTRIBUTION
# among all symptomatic patients with a fatal outcome, regardless if resident or not
####################################################################

cat('\nestimate time from symptom onset to death distribution ...')
set.seed(12)
tmp <- ddates[, 
	{
    	z1 <- seq.Date(w1.start, w1.end, by='day')
		z2 <- seq.Date(w2.start, w2.end, by='day')
		list(Date=c(z1,z2))
	}, 
	by='loc_label']
tmp <- unique(tmp) # Why did it work without this before??
ds0 <- merge(ds, tmp, by = c('Date', 'loc_label'))
ds0s <- ds0[, list(value=quantile(s2d, p=c(.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('loc_label', 'age.label')]
ds0s <- dcast.data.table(ds0s, loc_label+age.label~stat, value.var='value')
ds0[, s2dc := s2d + runif(nrow(ds0), 0, 1)]

if(args$make.preprocessing.plots)
{	
  p <- ggplot(ds0s, aes(x=age.label)) +
    geom_boxplot(aes(ymin=CL, lower=IL, middle=M, upper=IU, ymax=CU, fill=change_city_label(loc_label)), stat='identity', position=position_dodge()) +
    scale_colour_manual(values = args$city.palette(length(unique(ds0s$loc_label)))) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'bottom') +     
    labs(x='', y='time from symptom onset to death (days)', fill='cities') +
    coord_cartesian(ylim=c(0,100)) + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  ggsave(file=paste0(args$out.base,'_times2d_age_boxplots.pdf'),p,w=10,h=6)
  ggsave(file=paste0(args$out.base,'_times2d_age_boxplots.png'),p,w=10,h=6)
  
  # find distributions to fit to time from hospitalisation to death in each age band
  for(m in unique(ds0s$loc_label)){
    
    pdf(file=paste0(args$out.base,'_times2d_wildtype_fits_', stringr::str_replace_all(m, stringr::fixed(" "), "_"), '.pdf', collapse = '_'),w=8,h=50)
    ds0_m = subset(ds0, loc_label == m)
    par(mfrow = c(ds0_m[, length(levels(age.label))], 1))
    for(x in ds0_m[, levels(age.label)])
    {	
      fg <- fitdistrplus::fitdist(ds0_m[age.label==x,s2dc], "gamma")
      fln <- fitdistrplus::fitdist(ds0_m[age.label==x,s2dc], "lnorm")
      fw <- fitdistrplus::fitdist(ds0_m[age.label==x,s2dc], "weibull")
      fll <- fitdist(ds0_m[age.label==x,s2dc], "llogis", start = list(shape = 1, scale = 1))
      fp <- fitdist(ds0_m[age.label==x,s2dc], "pareto", start = list(shape = 1, scale = 1))	
      fitdistrplus::qqcomp(list(fw, fln, fll, fg, fp), xlogscale = TRUE, ylogscale = TRUE, main=x)	
    }
    dev.off()
  }
  #	--> we choose weibull
}

#	get weibull parameters and prob that symptom onset occurred on day 0,1,2,.. before date of death
dsd <- as.data.table(expand.grid(
  age.label = ds0[, levels(age.label)], 
  loc_label = ds0[, unique(loc_label)], 
  shape = NA_real_, 
  scale = NA_real_))
Age_label <- ds0[, levels(age.label)]

for(i in seq_along(Age_label))
{	
  for(m in ds0[, unique(loc_label)]){
    fg <- fitdistrplus::fitdist(ds0[age.label==Age_label[i] & loc_label == m,s2dc], "weibull")
    estimate =  unname(fg$estimate)
    dsd[age.label==Age_label[i] & loc_label == m,shape := estimate[1]] 
    dsd[age.label==Age_label[i] & loc_label == m,scale := estimate[2]] 
  }
}

dsd[, age.label := factor(age.label, levels = Age_label)]
dsd[, mean := scale*gamma(1+1/shape)]
dsd[, upper.lim := qweibull(1-1e-2, shape=shape, scale=scale)]
tmp <- dsd[,
           {
             z <- pweibull(1:(args$select.s2d.max+1), shape=shape, scale=scale)
             z <- c(z[1], diff(z))
             z <- z/sum(z)
             list(day=0:args$select.s2d.max, p=z)
           }, 	 
           by=c('loc_label', 'age.label')]
dsd <- merge(dsd, tmp, by=c('loc_label', 'age.label'))

if(args$make.preprocessing.plots)
{		
  p <- ggplot(dsd, aes(x=day)) +
    geom_point(aes(y=p, colour=age.label)) +
    viridis::scale_colour_viridis(option='E', end=0.8, discrete = TRUE) +	
    theme_bw() +
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.background = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          legend.position='bottom') +		
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots, scales='free_y') +
    labs(y='probability mass function', x='time from symptom onset to death (days)', colour='age')
  ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf.pdf'),p,w=10,h=12)
  ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf.png'),p,w=10,h=12)
  
  dsde <- ds0[, list(n=length(Date)), by=c('loc_label', 'age.label','s2d')]
  dsde <- merge(dsde, dsde[, list(total=sum(n)),by= c('loc_label', 'age.label')], by=c('loc_label', 'age.label'))
  dsde[,p:=n/total]
  dsde[,stat:='empirical']
  setnames(dsde, 's2d','day')
  set(dsde, NULL, c('n','total'), NULL)
  tmp <- subset(dsd, select=c(loc_label,age.label,day,p))
  tmp[, stat:='Weibull']
  dsde <- rbind(dsde, tmp) 
  
  if(args$partition.locs){
    tmp <- dsde[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='time from symptom onset to death (days)', fill='')
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical_part1.png'),p,w=10,h=12)
    
    tmp <- dsde[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='time from symptom onset to death (days)', fill='')
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical_part2.png'),p,w=10,h=12)
    
    tmp <- dsde[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='time from symptom onset to death (days)', fill='')
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical_part3.png'),p,w=10,h=12)
  }else{
    tmp <- copy(dsde)
    p <- ggplot(tmp, aes(x=day)) +
      geom_bar(aes(y=p, fill=stat), stat='identity', position=position_dodge()) +
      coord_cartesian(xlim=c(0,args$select.s2d.max)) +
      scale_y_continuous(expand=c(0,0)) +
      ggsci::scale_fill_aaas() +
      theme_bw() +	
      theme(strip.background = element_blank(),legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(y='probability mass function', x='time from symptom onset to death (days)', fill='')
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_times2d_age_discretized_pmf_vs_empirical.png'),p,w=10,h=12)
  }
}

####################################################################
# ESTIMATE TIME FROM INFECTION TO DEATH DISTRIBUTION
# using all symptomatic patients with a fatal outcome, regardless if resident or not
####################################################################

cat('\nestimate time from infection to death distribution ...')

# Evidence for longer incubation period among older individuals
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7484300/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7753659/
# https://www.dovepress.com/front_end/cr_data/cache/pdf/download_1616369757_6057d85d99e47/rmhp-257907-probable-longer-incubation-period-for-elderly-covid-19-cases.pdf
tmp <- data.table(
  age.from=c(0, 30, 60),
  age.to=c(30, 60, 150),
  incubation_shape= c(1.3, 1.6, 2.4),
  incubation_scale= c(5.3, 7.3, 8.9)
)
tmp2 <- unique(subset(da, select=c(age.label, age.from, age.to)))
tmp <- tmp2[, 
	{
  		z <- which(age.from>=tmp$age.from & age.to<=tmp$age.to)
  		list(incubation_shape=tmp$incubation_shape[z], incubation_scale=tmp$incubation_scale[z])
	}, 
	by='age.label']
dsd <- unique(subset(dsd, select=c(loc_label,age.label, shape, scale)))
setnames(dsd, c('shape','scale'), c('o2d_shape','o2d_scale'))
dsd <- merge(dsd, tmp, by= 'age.label')
dsd[, incubation_mean := incubation_scale*gamma(1+1/incubation_shape)]
dsd[, o2d_mean := o2d_scale*gamma(1+1/o2d_shape)]

set.seed(42)
tmp <- dsd[,
	{
		z1 <- rweibull(1e6, shape=incubation_shape, scale=incubation_scale)
		z2 <- rweibull(1e6, shape=o2d_shape, scale=o2d_scale)
		ecdf_i2d <- ecdf(z1+z2)			
		i2d <- sapply(1:(args$select.s2d.max+1), function(s) ecdf_i2d(s) - ecdf_i2d(s-1.0))
		list(day=0:args$select.s2d.max, i2d_p= i2d)			
	},
    by=c('loc_label', 'age.label')]
dsd <- merge(dsd,tmp,by=c('loc_label', 'age.label'))

####################################################################
# use time from hospital admission to death to 
# attribute deaths to hospitalised patients with unknown outcomes 
####################################################################


cat('\nattribute deaths to hospitalised patients with unknown outcomes ...')
dhunknown <- subset(dht, 
		grepl('unknown|in_ICU',Outcome2) & 
		hosp.cat!='out of hospital' &
		!is.na(pdeaths2)
	)
dhunknown <- merge(dhunknown, da, by='Age')
dhunknown <- dhunknown[, list(whosps=length(pdeaths2), pdeaths2=pdeaths2[1]), by=c('loc_label','loc.res.cat','DateAdmission','age.label','Hospital_type','vac_type')]
dhunknown <- dhunknown[, 
	{		
		z <- which(dtd[['loc_label']]==loc_label & dtd[['age.label']]==age.label)
		z <- dtd[['p']][z]
		list(
			DateDeath = DateAdmission + seq.int(0,args$select.h2d.max),
			Outcome2 = 'Death_Covid_Estimated',
			whosps = c(whosps, rep(0, args$select.h2d.max)),
			wdeaths = whosps * pdeaths2 * z 
		)			
	}, 
	by=c('loc_label','loc.res.cat','DateAdmission','age.label','Hospital_type','vac_type')]
dhunknown[, hosp.cat := 'in hospital']
if(args$deaths.estimate.outcomes.for.those.unknown==0)
{
  cat('\nNot adjusting SIVEP deaths for individuals with unknown outcomes ...')
  set(dhunknown, NULL, 'wdeaths', 0.)
}
if(args$deaths.estimate.outcomes.for.those.unknown==1)
{
  cat('\nAdjusting SIVEP deaths for individuals with unknown outcomes ...')
}

if(0) # dhunknown for curitiba
{
  tmp <- dhunknown[loc_label == 'curitiba']
  tmp[, cat(sum(wdeaths), 'out of ', sum(whosps), ' unknown hosps \n')]
  
  tmp  <- merge(tmp, ddates[, .(loc_label, w1.start, w2.end)], by='loc_label')
  tmp[DateAdmission >= w1.start & DateAdmission <= w2.end  & loc.res.cat == 'resident' & grepl(args$keep.vac.type, vac_type),
      cat(sum(wdeaths), 'out of ', sum(whosps), ' unknown hosps \n')]
  
}


####################################################################
# combine patients with known outcomes and patients with estimated resolved outcomes 
# still:
#	include out of hosp records
#	include all hosp types
####################################################################

dht <- subset(dht, 
		!grepl('unknown|in_ICU',Outcome2), 
		select=c(loc_label, loc.res.cat, Age, DateSymptoms, DateAdmission, DateICUstart, DateICUend, DateDeath, Outcome2, Hospital_type, vac_type, hosp.cat)
	)
dht[, wdeaths := as.numeric(Outcome2=='Death_Covid')]
dht[, whosps := as.numeric(!is.na(DateAdmission))]
dht <- merge(dht, subset(da, select=c(Age, age.label)), by='Age')
dht <- rbind(dht, dhunknown, fill=TRUE)

####################################################################
# prepare all other data sets based on the patients with known outcomes and estimated resolved outcomes  
####################################################################

cat('\nprepare all other data sets based on the patients with known outcomes and estimated resolved outcomes ...')

# make death data including non-residents
# 	regardless of hospital admission (COVID19)
#	regardless of hospital type
#	regardless of vaccine status
dda <- subset(dht, !is.na(DateDeath))
setnames(dda, 'DateDeath', 'Date')
dda <- dda[, list(Deaths=sum(wdeaths)),by=c('Date', 'age.label', 'loc_label','loc.res.cat')]
dda <- merge(dda, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dda <- merge(dda, tmp, by=c('loc_label','week'))

# keep only resident
if(args$select.only.SIVEP.patients.that.are.resident.in.loc)
{
  cat('\nselecting only SIVEP patients that are residents in the location')
  dht <- subset(dht, loc.res.cat == 'resident')
}

# make death data among residents
#	regardless of hospital admission (COVID19)
#	by hospital type
#	by in/out of hospital
dd <- subset(dht, !is.na(DateDeath))
setnames(dd, 'DateDeath', 'Date')
dd <- dd[, 
	list(
		deaths=sum(wdeaths),
		deaths_ntcns = sum(wdeaths[Outcome2=='Death_Covid']),
		deaths_hospadm= sum(wdeaths[hosp.cat!='out of hospital'])		
	),
	by=c('loc_label','age.label','Date','Hospital_type','vac_type','hosp.cat')]
setkey(dd, loc_label, age.label, Hospital_type, vac_type, hosp.cat, Date)
dd <- merge(dd, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dd <- merge(dd, tmp, by=c('loc_label','week'))

# make hosp admission data (COVID19)
# note: 
# here we adjust for censored outcomes in hospitals
# we do not adjust for under-reporting of COVID-19 deaths in the population through
# civil registry data, since the in-hospital fatality rate is not dependent on deaths
# occurring outside the hospital
dh <- dht[
	hosp.cat!='out of hospital', 
	list(
    	hosps= sum(whosps), 
    	deaths_of_hosps= sum(wdeaths) 
  	), 
  	by=c('loc_label','age.label','DateAdmission', 'Hospital_type','vac_type')]

tmp <- dht[
	hosp.cat!='out of hospital' & Outcome2!='Death_Covid_Estimated', 
  	list(
    	hosps_ntcns= sum(whosps), 
    	deaths_of_hosps_ntcns= sum(wdeaths) 
  	), 
  	by=c('loc_label','age.label','DateAdmission', 'Hospital_type','vac_type')]
dh <- merge(dh, tmp, by=c('loc_label','age.label','DateAdmission', 'Hospital_type','vac_type'))
set(dh, dh[, which(is.na(hosps_ntcns))],'hosps_ntcns',0)
set(dh, dh[, which(is.na(deaths_of_hosps_ntcns))],'deaths_of_hosps_ntcns',0)
setnames(dh, 'DateAdmission', 'Date')
dh <- merge(dh, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dh <- merge(dh, tmp, by=c('loc_label','week'))

####################################################################
# AGGREGATE COVID19 class 4+5 DEATHS (residents and non-residents) 
####################################################################

dda <- dda[, list(deaths=sum(Deaths)), by=c('week','age.label', 'loc_label','loc.res.cat')]
tmp <- ddates[,list(week=seq(hosps.start, deaths.end, 1)),by='loc_label']
tmp[, DUMMY:=1L]
tmp2 <- as.data.table(expand.grid(DUMMY=1L, age.label=unique(da$age.label), loc.res.cat=unique(dda$loc.res.cat) ))
tmp <- merge(tmp, tmp2, by='DUMMY', allow.cartesian=TRUE)
set(tmp, NULL, 'DUMMY', NULL)
tmp <- merge(tmp, unique(subset(dw, select=c(week,week.start))), by='week')
dda <- merge(tmp, dda, all.x=TRUE, by=c('week','age.label','loc_label','loc.res.cat'))
set(dda, dda[, which(is.na(deaths))], 'deaths', 0L)

if(args$make.preprocessing.plots)
{				
  if(args$partition.locs){
    tmp <- dda[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_bar(aes(y=deaths, fill=loc.res.cat), stat='identity', position='stack') +    
      ggsci::scale_fill_npg() +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(fill='', y='COVID19 deaths', x='') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_resnonres_part1.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_resident_part1.png'),p,w=10,h=12,limitsize = FALSE)
    
    tmp <- dda[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_bar(aes(y=deaths, fill=loc.res.cat), stat='identity', position='stack') +    
      ggsci::scale_fill_npg() +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(fill='', y='COVID19 deaths', x='') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_resident_part2.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_resident_part2.png'),p,w=10,h=12,limitsize = FALSE)
    
    tmp <- dda[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_bar(aes(y=deaths, fill=loc.res.cat), stat='identity', position='stack') +    
      ggsci::scale_fill_npg() +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(fill='', y='COVID19 deaths', x='') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_resident_part3.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_resident_part3.png'),p,w=10,h=12,limitsize = FALSE)
  }else{
    tmp <- copy(dda)
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_bar(aes(y=deaths, fill=loc.res.cat), stat='identity', position='stack') +    
      ggsci::scale_fill_npg() +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(fill='', y='COVID19 deaths', x='') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_resident.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_resident.png'),p,w=10,h=12,limitsize = FALSE)
  }
}

####################################################################
# AGGREGATE REGISTRY DEATHS (residents and non-residents)
####################################################################

cat('\naggregate registry deaths ...')

# aggregate
dr <- melt(dr, id.vars=c('week','age.label','loc_label','pop','Date'), measure.vars= c('reg.deaths.ooh','reg.coviddeaths.all', 'reg.deaths.all','reg.deaths.all.excess','reg.deaths.all.excess.7avg'))
dr <- dr[, list(value=sum(value)), by=c('loc_label','age.label','pop','week','variable')]
tmp <- ddates[,list(week=seq(hosps.start, deaths.end, 1)),by='loc_label']
tmp[, DUMMY:=1L]
tmp2 <- as.data.table(expand.grid(DUMMY=1L, age.label=unique(da$age.label), variable=unique(dr$variable)))
tmp <- merge(tmp, tmp2, by='DUMMY', allow.cartesian=TRUE)
set(tmp, NULL, 'DUMMY', NULL)
tmp <- merge(tmp, unique(subset(dw, select=c(week,week.start))), by='week')
dr <- merge(tmp, dr, all.x=TRUE, by=c('loc_label','age.label','week','variable'))
set(dr, dr[, which(is.na(value))], 'value', 0L)

# take proportion corresponding to residents seen in SIVEP data set
tmp <- dcast.data.table(dda, week+age.label+loc_label~loc.res.cat, value.var='deaths')
setnames(tmp, 'not resident', 'not_resident')
tmp[, president := resident / (resident+not_resident)]
set(tmp, tmp[, which(0==(resident+not_resident))], 'president', 1.)
set(tmp, NULL, c('resident','not_resident'), NULL)
dr <- merge(dr, tmp, by=c('loc_label','age.label','week'))
set(dr, NULL, 'value', dr[, value*president])

# bring into wide format
dr <- dcast.data.table(dr, loc_label+age.label+week+week.start+pop ~ variable, value.var='value')
dr[, reg.deaths.all.exorcovid:= pmax(reg.coviddeaths.all, reg.deaths.all.excess.7avg)]

####################################################################
# make death adjustment ratio, reflecting the ratio of excess deaths vs reported deaths among residents and non-residents
####################################################################

cat('\nmake death adjustment ratio...')

tmp <- subset(dd, select=c(loc_label,age.label,week,deaths))
tmp <- tmp[, list(deaths=sum(deaths)), by=c('loc_label','age.label','week')]
dr <- merge(dr, tmp, by=c('loc_label','age.label','week'), all.x=TRUE)
set(dr, dr[, which(is.na(deaths))], 'deaths', 0.)
dr[, adjdeaths := deaths]

if(!is.na(args$deaths.adjust.for.civilregistrydeaths))
{
	cat('\nadjusting total deaths to ',args$deaths.adjust.for.civilregistrydeaths)
	set(dr, NULL, 'adjdeaths', pmax( dr[, adjdeaths], dr[[args$deaths.adjust.for.civilregistrydeaths]]))		
}

if(args$make.preprocessing.plots){
  
  tmp <- c('loc_label', 'age.label','week.start','deaths', 'reg.coviddeaths.all','reg.deaths.all.excess')
  tmp <- dr[, ..tmp]
  tmp <- melt(tmp, id.vars=c('loc_label','age.label', 'week.start'))
  tmp <- tmp[, list(value=sum(value)), by=c('loc_label', 'week.start', 'variable')]
  tmp[, variable:= as.factor(variable)]
  levels(tmp$variable) <- c('SIVEP censoring-adjusted deaths', 'Registry COVID deaths', 'Registry Excess deaths')

  tmp2 <- dht[loc.res.cat == 'resident' & Outcome2 == 'Death_Covid' & !is.na(DateDeath)]
  tmp2 <- merge(tmp2, dw[,.(Date, week.start)], by.x = 'DateDeath', by.y='Date')
  tmp2 <- tmp2[,list( variable='SIVEP observed deaths',value=.N), by=c('loc_label', 'week.start')]
  tmp <- rbind(tmp2, tmp)
  
  p <- ggplot(tmp, aes(x=week.start, y=value, color=variable)) + 
    geom_step() + 
    facet_wrap(~change_city_label(loc_label), scale = 'free_y') + 
    ggsci::scale_colour_npg() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    scale_x_date(breaks="3 months", expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    #viridis::scale_colour_viridis(option='B', end=0.85, labels = scales::percent, limits=c(0,1)) +
    labs(x='',y='weekly in-hospital deaths', colour='source') 
  ggsave(file=paste0(args$out.base, '_deaths_censvsregistry_2.png'), p, w=10, h=12)
  ggsave(file=paste0(args$out.base, '_deaths_censvsregistry_2.pdf'), p, w=10, h=12)
}

# note:
# we can adjust the deaths only on the total across hospital types since we don t have registry 
# data by private / public hospital

# plots: compare censoring adjustment to underreporting+censoring adjustment
if (args$make.preprocessing.plots)
{  
  if(args$partition.locs){
    # set 1
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    tmp <- tmp1[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') +
      geom_step(aes(y=value, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y='weekly COVID-19 attributable deaths', colour='') 
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n_part1.png'),p,w=10,h=12)
    
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    tmp <- tmp1[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') +
      geom_step(aes(y=value, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y='weekly COVID-19 attributable deaths', colour='') 
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n_part2.png'),p,w=10,h=12)
    
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    tmp <- tmp1[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') +
      geom_step(aes(y=value, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y='weekly COVID-19 attributable deaths', colour='') 
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n_part3.png'),p,w=10,h=12)
    
    # Set 2 
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp <- tmp[, list(age.label=age.label, pvalue=value/sum(value)), by=c('variable','week','week.start','loc_label')]
    set(tmp, tmp[, which(is.nan(pvalue))], 'pvalue', 0.)
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    tmp <- tmp1[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=pvalue, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y="share of age group among weekly COVID-19 attributable deaths", colour='')  
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p_part1.png'),p,w=10,h=12)
    
    
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp <- tmp[, list(age.label=age.label, pvalue=value/sum(value)), by=c('variable','week','week.start','loc_label')]
    set(tmp, tmp[, which(is.nan(pvalue))], 'pvalue', 0.)
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    tmp <- tmp1[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=pvalue, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y="share of age group among weekly COVID-19 attributable deaths", colour='')  
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p_part2.png'),p,w=10,h=12)
    
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp <- tmp[, list(age.label=age.label, pvalue=value/sum(value)), by=c('variable','week','week.start','loc_label')]
    set(tmp, tmp[, which(is.nan(pvalue))], 'pvalue', 0.)
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    tmp <- tmp1[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=pvalue, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y="share of age group among weekly COVID-19 attributable deaths", colour='')  
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p_part3.png'),p,w=10,h=12)	
  }else{
    
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    p <- ggplot(tmp1[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      geom_step(aes(y=value, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y='weekly COVID-19 attributable deaths', colour='') 
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_n.png'),p,w=10,h=12)
    
    tmp <- melt(dr, id.vars=c('week','week.start','age.label','loc_label'), measure.vars= c('deaths','adjdeaths'))
    set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('adjdeaths','deaths'), labels=c( "underreporting and censoring adjusted","censoring adjusted"))])
    tmp <- tmp[, list(age.label=age.label, pvalue=value/sum(value)), by=c('variable','week','week.start','loc_label')]
    set(tmp, tmp[, which(is.nan(pvalue))], 'pvalue', 0.)
    tmp1 <- merge(tmp, ddates[,c("loc_label", 'w1.start')], by = 'loc_label')
    
    p <- ggplot(tmp1[week.start > w1.start], aes(x= week.start)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=pvalue, colour=variable)) +    
      ggsci::scale_colour_aaas() +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', strip.background = element_blank()) +	
      labs(x='', y="share of age group among weekly COVID-19 attributable deaths", colour='')  
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_deaths_adjcomp_p.png'),p,w=10,h=12)	
  }
}

####################################################################
# AGGREGATE COVID19 class 4+5 DEATHS (residents only)
# note:
#	deaths are regardless of in hospital or not
#	deaths_hospadm are only in hospital
####################################################################

cat('\naggregate COVID-19 susp/conf deaths (residents only) ...')
dd <- dd[, 
	list(
		deaths=sum(deaths), 
		deaths_ntcns=sum(deaths_ntcns), 
		deaths_hospadm=sum(deaths_hospadm[hosp.cat!='out of hospital'])
	), 
    by=c('week','age.label','loc_label','Hospital_type','vac_type')]
tmp <- ddates[,list(week=seq(hosps.start, deaths.end, 1)),by='loc_label']
tmp[, DUMMY:=1L]
tmp2 <- as.data.table(expand.grid(
		DUMMY=1L, 
		age.label=unique(da$age.label), 
		Hospital_type=unique(dd$Hospital_type),
		vac_type=unique(dd$vac_type)
	))
tmp <- merge(tmp, tmp2, by='DUMMY', allow.cartesian=TRUE)
set(tmp, NULL, 'DUMMY', NULL)
tmp <- merge(tmp, unique(subset(dw, select=c(week,week.start))), by='week')
dd <- merge(tmp, dd, all.x=TRUE, by=c('loc_label', 'Hospital_type','vac_type','age.label','week'))
set(dd, dd[, which(is.na(deaths))], 'deaths', 0)
set(dd, dd[, which(is.na(deaths_ntcns))], 'deaths_ntcns', 0L)
set(dd, dd[, which(is.na(deaths_hospadm))], 'deaths_hospadm', 0)
set(dd, NULL, 'deaths_ooh', dd[, deaths-deaths_hospadm])

####################################################################
# plot source of all age deaths
####################################################################

if(args$make.preprocessing.plots)
{	
	tmp <- dd[, 
		list(
			deaths=sum(deaths), 
			deaths_ntcns=sum(deaths_ntcns), 
			deaths_hospadm=sum(deaths_hospadm)
		), 
		by=c('loc_label','Hospital_type','age.label','week','week.start')]
	tmp <- copy(dd)
	tmp[, deaths_hospadm_cns := deaths_hospadm - (deaths_ntcns - deaths_ooh)]
	tmp[, deaths_hospadm_ntcns := deaths_ntcns - deaths_ooh]
	tmp <- melt(tmp, id.vars=c('loc_label','Hospital_type','age.label','week','week.start'), measure.vars=c('deaths_hospadm_ntcns','deaths_ooh','deaths_hospadm_cns'))
	set(tmp, tmp[, which(grepl('ooh', variable))], 'Hospital_type', 'out of hospital')
	set(tmp, tmp[, which(grepl('_cns', variable))], 'Hospital_type', 'censoring')
	tmp <- tmp[, list(value=sum(value)), by=c('loc_label','Hospital_type','week','week.start')]
	
	tmp2 <- copy(dr)
	tmp2[, underreporting := adjdeaths-deaths]
	tmp2 <- tmp2[, list(value=sum(underreporting)), by=c('loc_label','week','week.start')]
	tmp2[, Hospital_type := 'underreporting']
	tmp <- rbind(tmp, tmp2)
	set(tmp, NULL, 'Hospital_type', tmp[, factor(Hospital_type, 
		levels=c('underreporting','out of hospital','Private Hospitals','Public Hospitals','Unknown','censoring'), 
		labels=c('difference to excess deaths', 
				 'reported out-of-hospital deaths', 
				 'in private hospitals',
				 'in public hospitals',
				 'in private or public hospitals',
				 'expected deaths among COVID-19 hospitalised\npatients with as of yet unknown outcome'))])
	tmp <- merge(tmp, ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label' ], by=c('loc_label','week'))
	
	p <- ggplot(tmp) +    	
			# geom_rect(data=ddates, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
			# geom_vline(data=ddates, aes(xintercept=p1em.median), colour='black', linetype='solid') +
			geom_vline(data=ddates, aes(xintercept=w2.start), colour='black', linetype='dotted') +	
			geom_bar(aes(x=week.start, y=value, fill=Hospital_type), width=6, stat='identity', position='stack') +		
			scale_y_continuous(expand=c(0,0)) +
			scale_x_date(breaks='months',expand=c(0,0), labels=date_format("%b %y") ) +
			ggsci::scale_fill_nejm() +
			theme_bw() +
			facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots, scale = 'free_y') +
			guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
					strip.background = element_blank(),
					legend.position='bottom') +		
			labs(x='', y='observed and expected COVID-19 attributable deaths', fill='') + reqs
	ggsave(file=paste0(args$out.base,'_deathsn_over_time_byexpobs_2.pdf'),p,w=15,h=18, dpi=350, units='cm')
	ggsave(file=paste0(args$out.base,'_deathsn_over_time_byexpobs_2.png'),p,w=15,h=18, dpi=350, units='cm')
}

####################################################################
# after this plot, we only need: 
#	underreporting + censoring adjusted deaths
#	censoring adjusted deaths after hosp admission in hospital type
####################################################################

dd <- subset(dd, 
		grepl(args$keep.hospital.type, Hospital_type) &
		grepl(args$keep.vac.type, vac_type)
	)
dd <- dd[, list(deaths_hospadm=sum(deaths_hospadm)),by=c('week','age.label', 'loc_label')]
tmp <- subset(dr, select=c(loc_label,age.label,week,week.start,adjdeaths))
dd <- merge(tmp, dd, by=c('loc_label','age.label','week'))
dd <- dd[order(loc_label,age.label,week),]
tmp <- dd[, list(week=week, cadjdeaths=cumsum(adjdeaths), cdeathshospadm=cumsum(deaths_hospadm)), by=c('age.label', 'loc_label')]
dd <- merge(dd, tmp, by=c('week','age.label', 'loc_label'))
tmp <- dd[,list(age.label=age.label, scadjdeaths=cumsum(cadjdeaths), scdeathshospadm=cumsum(cdeathshospadm)), by=c('week', 'loc_label')]
dd <- merge(dd, tmp, by=c('week','age.label', 'loc_label'))
tmp <- dd[,list(tadjdeaths=sum(adjdeaths), tdeathshospadm=sum(deaths_hospadm)), by=c('week', 'loc_label')]
dd <- merge(dd, tmp, by=c('week', 'loc_label'))
dd[, padjdeaths:=adjdeaths/tadjdeaths]
dd[, pdeathshospadm:=deaths_hospadm/tdeathshospadm]
tmp <- dd[,list(age.label=age.label, spadjdeaths=cumsum(padjdeaths), spdeathshospadm=cumsum(pdeathshospadm) ), by=c('week', 'loc_label')]
dd <- merge(dd, tmp, by=c('week','age.label', 'loc_label'))
set(dd, dd[, which(is.nan(padjdeaths))], 'padjdeaths', 0.)
set(dd, dd[, which(is.nan(spadjdeaths))], 'spadjdeaths', 0.)
set(dd, dd[, which(is.nan(pdeathshospadm))], 'pdeathshospadm', 0.)
set(dd, dd[, which(is.nan(spdeathshospadm))], 'spdeathshospadm', 0.)


if(args$make.preprocessing.plots)
{	
  	# plot empirical cum deaths
  	p <- ggplot(dd) +
		geom_rect(data=ddates, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
		geom_vline(data=ddates, aes(xintercept=p1em.median), colour='black', linetype='solid') +
		geom_vline(data=ddates, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  		  
    	geom_ribbon(aes(x=week.start, ymin=scadjdeaths-cadjdeaths, ymax=scadjdeaths, fill=age.label)) +
    	viridis::scale_fill_viridis(option='B', end=0.85, discrete = TRUE) +
    	scale_x_date(breaks='2 months',expand=c(0,0)) +
    	scale_y_continuous(expand=c(0,0)) +
    	labs(fill='age', y='cumulated, weekly underreporting and censoring adjusted COVID-19 attributable deaths', x='') +
    	theme_bw() +
    	facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots , scale = 'free_y') +
    	guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') 
  	ggsave(file=paste0(args$out.base,'_deaths_adj_c.pdf'),p,w=10,h=12,limitsize = FALSE)
  	ggsave(file=paste0(args$out.base,'_deaths_adj_c.png'),p,w=10,h=12,limitsize = FALSE)
  
  	# plot empirical cum deaths
  	p <- ggplot(dd) +
		geom_rect(data=ddates, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
		geom_vline(data=ddates, aes(xintercept=p1em.median), colour='black', linetype='solid') +
		geom_vline(data=ddates, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  		  
		geom_ribbon(aes(x=week.start,ymin=scdeathshospadm-cdeathshospadm, ymax=scdeathshospadm, fill=age.label)) +
		viridis::scale_fill_viridis(option='B', end=0.85, discrete = TRUE) +
		scale_x_date(breaks='2 months',expand=c(0,0)) +
		scale_y_continuous(expand=c(0,0)) +
		labs(fill='age', y='cumulated, weekly censoring adjusted COVID-19 attributable deaths in hospitals', x='') +
		theme_bw() +
		facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots , scale = 'free_y') +
		guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				  strip.background = element_blank(),
				  legend.position='bottom') 
  	ggsave(file=paste0(args$out.base,'_deaths_inhosp_c.pdf'),p,w=10,h=12,limitsize = FALSE)
  	ggsave(file=paste0(args$out.base,'_deaths_inhosp_c.png'),p,w=10,h=12,limitsize = FALSE)
}
if(args$make.preprocessing.plots)
{	  

  if(args$partition.locs){
    tmp <- dd[loc_label %in% args$selected.locs.part1]
    tmp <- merge(tmp, ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label' ], by=c('loc_label','week'))
    tmp2 <- ddates[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp) +
      geom_rect(data=tmp2, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
      geom_vline(data=tmp2, aes(xintercept=p1em.median), colour='black', linetype='solid') +
      geom_vline(data=tmp2, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  
      geom_point(aes(x=week.start, y=padjdeaths, colour=age.label, pch=change_city_label(loc_label))) +    
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(colour='age', y='share of age groups among underreporting and censoring adjusted COVID-19 attributable deaths', x='', pch='cities') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      guides(colour='none') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet_part1.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet_part1.png'),p,w=10,h=12,limitsize = FALSE)
    
    tmp <- dd[loc_label %in% args$selected.locs.part2]
    tmp <- merge(tmp, ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label' ], by=c('loc_label','week'))
    tmp2 <- ddates[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp) +
      geom_rect(data=tmp2, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
      geom_vline(data=tmp2, aes(xintercept=p1em.median), colour='black', linetype='solid') +
      geom_vline(data=tmp2, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  
      geom_point(aes(x=week.start, y=padjdeaths, colour=age.label, pch=change_city_label(loc_label))) +    
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(colour='age', y='share of age groups among underreporting and censoring adjusted COVID-19 attributable deaths', x='', pch='cities') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      guides(colour='none') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet_part2.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet_part2.png'),p,w=10,h=12,limitsize = FALSE)
    
    tmp <- dd[loc_label %in% args$selected.locs.part3]
    tmp <- merge(tmp, ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label' ], by=c('loc_label','week'))
    tmp2 <- ddates[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp) +
      geom_rect(data=tmp2, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
      geom_vline(data=tmp2, aes(xintercept=p1em.median), colour='black', linetype='solid') +
      geom_vline(data=tmp2, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  
      geom_point(aes(x=week.start, y=padjdeaths, colour=age.label, pch=change_city_label(loc_label))) +    
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(colour='age', y='share of age groups among underreporting and censoring adjusted COVID-19 attributable deaths', x='', pch='cities') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      guides(colour='none') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet_part3.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet_part3.png'),p,w=10,h=12,limitsize = FALSE)
  }else{
    tmp <- copy(dd)
    tmp <- merge(tmp, ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label' ], by=c('loc_label','week'))
    p <- ggplot(tmp) +
      geom_rect(data=ddates, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
      geom_vline(data=ddates, aes(xintercept=p1em.median), colour='black', linetype='solid') +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  
      geom_point(aes(x=week.start, y=padjdeaths, colour=age.label, pch=change_city_label(loc_label))) +    
      viridis::scale_colour_viridis(option='B', end=0.85, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(colour='age', y='share of age groups among underreporting and censoring adjusted COVID-19 attributable deaths', x='', pch='cities') +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label), scale = 'free_y') +
      guides(colour='none') +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='top') 	
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet.pdf'),p,w=10,h=12,limitsize = FALSE)
    ggsave(file=paste0(args$out.base,'_deaths_adj_p_facet.png'),p,w=10,h=12,limitsize = FALSE)
  }
  
}
if(args$make.preprocessing.plots)
{	  	 
	tmp <- merge(dd, ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label' ], by=c('loc_label','week'))
	p <- ggplot(tmp) +
			geom_rect(data=ddates, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
			geom_vline(data=ddates, aes(xintercept=p1em.median), colour='black', linetype='solid') +
			geom_vline(data=ddates, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  
			geom_bar(aes(x=week.start, y=adjdeaths, fill=age.label), width=6, stat='identity', position='stack') +		
			scale_y_continuous(expand=c(0,0)) +
			scale_x_date(breaks='months',expand=c(0,0)) +
			viridis::scale_fill_viridis(option='B', end=0.85, discrete = TRUE) +
			theme_bw() +
			facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plot, scale = 'free_y') +
			guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
					strip.background = element_blank(),
					legend.position='bottom') +		
			labs(x='', y='weekly underreporting and censoring adjusted COVID-19 attributable deaths', fill='age')
	ggsave(file=paste0(args$out.base,'_deaths_adj_n_over_time.pdf'),p,w=10,h=12,limitsize = FALSE)
	ggsave(file=paste0(args$out.base,'_deaths_adj_n_over_time.png'),p,w=10,h=12,limitsize = FALSE)	
		
	p <- ggplot(tmp) +
			geom_rect(data=ddates, aes(xmin=p1em.CL, xmax=p1em.CU, ymin=-Inf, ymax=Inf), fill='grey80', show.legend = FALSE) +
			geom_vline(data=ddates, aes(xintercept=p1em.median), colour='black', linetype='solid') +
			geom_vline(data=ddates, aes(xintercept=w2.start), colour='black', linetype='dotted') +			  
			geom_bar(aes(x=week.start, y=deaths_hospadm, fill=age.label), width=6, stat='identity', position='stack') +		
			scale_y_continuous(expand=c(0,0)) +
			scale_x_date(breaks='months',expand=c(0,0)) +
			viridis::scale_fill_viridis(option='B', end=0.85, discrete = TRUE) +
			theme_bw() +
			facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plot, scale = 'free_y') +
			guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
					strip.background = element_blank(),
					legend.position='bottom') +		
			labs(x='', y='weekly censoring adjusted COVID-19 attributable deaths in hospitals', fill='age')
	ggsave(file=paste0(args$out.base,'_deaths_inhosp_n_over_time.pdf'),p,w=10,h=12,limitsize = FALSE)
	ggsave(file=paste0(args$out.base,'_deaths_inhosp_n_over_time.png'),p,w=10,h=12,limitsize = FALSE)	
}


####################################################################
# calculate population protected through vaccines
####################################################################

cat('\ncalculate population protected through vaccines...')
tmp <- subset(dv, select=c(week, loc_label, age.label, vaccine, pcvac1, pcvac2))
tmp[, prtct := 0.]
tmp2 <- which(tmp$vaccine %in% c("Covishield", "Pfizer", "Janssen"))
set(tmp, tmp2, 'week', tmp[tmp2, week+args$vaccine.astra.protection.delay])
set(tmp, tmp2, 'prtct', tmp[tmp2, pcvac1*args$vaccine.astra.ve.1dose + pcvac2*args$vaccine.astra.ve.2doses])
tmp2 <- which(tmp$vaccine %in% c("Sinovac"))
set(tmp, tmp2, 'week', tmp[tmp2, week+args$vaccine.sinovac.protection.delay])
set(tmp, tmp2, 'prtct', tmp[tmp2, pcvac1*args$vaccine.sinovac.ve.1dose + pcvac2*args$vaccine.sinovac.ve.2doses])
set(tmp, NULL, c('pcvac1','pcvac2'), NULL)

dv <- merge(dv, tmp, by=c('loc_label','week','age.label','vaccine'), all.x=TRUE)
set(dv, dv[,which(is.na(prtct))],'prtct',0.)
dv[, cprtct := prtct*pop]

if(args$make.preprocessing.plots)
{
  tmp <- dv[,.(week.start, loc_label, age.label, vaccine,prtct)]
  tmp <- tmp[,list(prtct=sum(prtct)), by = c('week.start', 'loc_label', 'age.label')]
  tmp <- tmp[prtct > .001,,]
  
  p <- ggplot(tmp, aes(x=week.start)) +
    # geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_line(aes(y = prtct, color = age.label)) +    
    scale_x_date(breaks='3 weeks',expand=c(0,0)) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    viridis::scale_color_viridis(discrete=TRUE, option='B', end=.8, expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    labs(color = "Age group", y = 'estimated proportion of population protected from severe COVID-19 requiring hospitalisation by vaccine roll-out', x = '')  
  ggsave(file=paste0(args$out.base,'_pop_prop_vaxprotected.pdf'),p,w=10,h= 12)
  ggsave(file=paste0(args$out.base,'_pop_prop_vaxprotected.png'),p,w=10,h= 12)
  
  p <- ggplot(tmp, aes(x=week.start)) +
    # geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_line(aes(y = 1-prtct, color = age.label)) +    
    scale_x_date(breaks='3 weeks',expand=c(0,0)) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    viridis::scale_color_viridis(discrete=TRUE, option='B', end=.8, expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    labs(color = "Age group", y = 'mortality adjusted population at risk of severe COVID-19 requiring hospitalisation', x = '')  
  ggsave(file=paste0(args$out.base,'_pop_prop_vaxprotected2.pdf'),p,w=10,h= 12)
  ggsave(file=paste0(args$out.base,'_pop_prop_vaxprotected2.png'),p,w=10,h= 12)
}


####################################################################
# calculate remaining population (rpop.alive)
####################################################################

cat('\ncalculate remaining population...')
drpop <- subset(dd, select=c(loc_label, age.label, week, cadjdeaths)) 
drpop <- merge(drpop, subset(dpop,select=c(loc_label,age.label, pop)), by=c('age.label', 'loc_label'))
drpop[, rpop.alive := pop - cadjdeaths ]
set(drpop, NULL, c('cadjdeaths'), NULL)

####################################################################
# calculate remaining unprotected population (rpop.alive.notvac)
####################################################################

tmp <- dv[,
          list(
            cvac1 = as.numeric(sum(cvac1+cvac2)), 
            cvac2 = as.numeric(sum(cvac2)),         
            cprtct = sum(cprtct)
          ),
          by = c('loc_label', 'age.label', 'week')]
tmp <- subset(tmp, select=c(loc_label, age.label, week, cvac1, cvac2, cprtct))

drpop <- merge(drpop, tmp, by=c('loc_label','age.label','week'), all.x=TRUE)
set(drpop, drpop[, which(is.na(cvac1))], 'cvac1', 0)
set(drpop, drpop[, which(is.na(cvac2))], 'cvac2', 0)
set(drpop, drpop[, which(is.na(cprtct)) ], 'cprtct', 0)

# make consistent: remaining pop size must be larger or equal than number having received 1st dose
tmp <- drpop[, list(week = max(week)), by='loc_label']
tmp <- merge(drpop, tmp, by=c('loc_label','week'))
tmp[, vac1.adj := pmin(1.,rpop.alive/cvac1)]
tmp[, vac1.rate := cvac1/rpop.alive*vac1.adj]
write.csv(tmp, file=paste0(args$out.base,'_pop_vaccinated.csv'))

# assume vaccination rate cannot exceed args$vaccine.max.coverage
tmp2 <- tmp[, which(vac1.rate>args$vaccine.max.coverage)]
set(tmp, tmp2, 'vac1.adj', tmp[tmp2, args$vaccine.max.coverage*vac1.adj])

# adjust remaining population for vaccination status
drpop <- merge(drpop, subset(tmp, select=c(loc_label,age.label,vac1.adj)), by=c('loc_label','age.label'))
drpop[, rpop.alive.notvac := rpop.alive - cprtct*vac1.adj]
drpop <- subset(drpop, select=c(loc_label, week, age.label, pop, rpop.alive, rpop.alive.notvac))
drpop[, p.alive := rpop.alive/pop]
drpop[, p.alive.notvac := rpop.alive.notvac/pop]
setkey(drpop, loc_label,week, age.label)
tmp <- drpop[, 
             list(
               age.label=age.label, 
               crpop.alive.notvac=cumsum(rpop.alive.notvac)
             ), by=c('loc_label', 'week') ]
drpop <- merge(drpop, tmp, by=c('loc_label', 'week','age.label'))
drpop <- merge(drpop, unique(subset(dw, select=c(week, week.start))), by=c('week'))

# to estimate the death composition, we need to push forward 
# the remaining, unvaccinated population by approximately 3 weeks 
# (1 week infection to symptom onset, 2 weeks symptom onset to death) 
tmp <- subset(drpop, select=c(loc_label, age.label, week, rpop.alive.notvac))
set(tmp, NULL, 'week', tmp[,week+3L])
tmp2 <- drpop[, 
              {
                z <- week==min(week)
                list(week= week[z]+0:2, pop=pop[z])
              }, by=c('loc_label','age.label')]
setnames(tmp, 'rpop.alive.notvac', 'rpop.alive.notvac.push3')
setnames(tmp2, 'pop', 'rpop.alive.notvac.push3')
tmp <- rbind(tmp, tmp2)
drpop <- merge(drpop, tmp, by=c('loc_label','age.label','week'))
drpop[, p.alive.notvac.push3 := rpop.alive.notvac.push3/pop]

if(args$make.preprocessing.plots)
{
  p <- ggplot(drpop, aes(x=week.start)) +
    geom_ribbon(aes(ymin=crpop.alive.notvac-rpop.alive.notvac, ymax=crpop.alive.notvac, fill=age.label)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    viridis::scale_fill_viridis(discrete=TRUE, option='E', end=0.8, expand=c(0,0)) +
    labs(fill='age', y='population\n(alive, not protected from a severe outcome through vaccination)', x='') +
    theme_bw() +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots, scale = 'free_y') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')
  ggsave(file=paste0(args$out.base,'_pop_over_time.pdf'),p,w=10,h=12)
  ggsave(file=paste0(args$out.base,'_pop_over_time.png'),p,w=10,h=12)
}

if(args$make.preprocessing.plots)
{	  
  if(args$partition.locs){
    # set 1
    tmp <- drpop[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac, fill=age.label)) +
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y='proportion of initial population alive\nand not protected from a severe outcome through vaccination', x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_part1.png'),p,w=10,h=12)
    
    tmp <- drpop[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac, fill=age.label)) +
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y='proportion of initial population alive\nand not protected from a severe outcome through vaccination', x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_part2.png'),p,w=10,h=12)
    
    tmp <- drpop[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac, fill=age.label)) +
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y='proportion of initial population alive\nand not protected from a severe outcome through vaccination', x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_part3.png'),p,w=10,h=12)
    
    # set 2 
    tmp <- drpop[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac.push3, fill=age.label)) +
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y=paste0('proportion of initial population alive\nand not protected from a severe outcome through vaccination\n', args$vaccine.protection.delay, ' weeks'), x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom') +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) 
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks_part1.png'),p,w=10,h=12)
    
    tmp <- drpop[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac.push3, fill=age.label)) +
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y=paste0('proportion of initial population alive\nand not protected from a severe outcome through vaccination\n', args$vaccine.protection.delay, ' weeks'), x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom') +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) 
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks_part2.png'),p,w=10,h=12)
    
    tmp <- drpop[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac.push3, fill=age.label)) +
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y=paste0('proportion of initial population alive\nand not protected from a severe outcome through vaccination\n', args$vaccine.protection.delay, ' weeks'), x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom') +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) 
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks_part3.png'),p,w=10,h=12)
  }else{
    tmp <- copy(drpop)
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac, fill=age.label)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y='proportion of initial population alive\nand not protected from a severe outcome through vaccination', x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')
    ggsave(file=paste0(args$out.base,'_pop_prop_alive.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive.png'),p,w=10,h=12)
    
    tmp <- copy(drpop)
    p <- ggplot(tmp, aes(x=week.start)) +	
      geom_ribbon(aes(ymin=0, ymax=p.alive.notvac.push3, fill=age.label)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='black') +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      viridis::scale_fill_viridis(discrete=TRUE, option='E', end=.8, expand=c(0,0)) +
      labs(fill='Age', y=paste0('proportion of initial population alive\nand not protected from a severe outcome through vaccination\n', args$vaccine.protection.delay, ' weeks'), x='', fill = "age") +
      theme_bw() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom') +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) 
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_pop_prop_alive_push3weeks.png'),p,w=10,h=12)
  }

  p <- ggplot(drpop, aes(x=week.start)) +
	geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_line(aes(y = 1 - p.alive, color = age.label)) +    
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    viridis::scale_color_viridis(discrete=TRUE, option='B', end=.8, expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    labs(color = "Age group", y = 'cumulative COVID-19 attributable mortality in percent of 2020 population', x = '')  
  ggsave(file=paste0(args$out.base,'_pop_prop_coviddeaths.pdf'),p,w=10,h= 12)
  ggsave(file=paste0(args$out.base,'_pop_prop_coviddeaths.png'),p,w=10,h= 12)
  
  ##### 
  
  p <- ggplot(drpop, aes(x=week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_line(aes(y = p.alive, color = age.label)) +    
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    viridis::scale_color_viridis(discrete=TRUE, option='B', end=.8, expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    labs(color = "Age group", y = 'population still alive in percent of 2020 population', x = '')  
  ggsave(file=paste0(args$out.base,'_pop_prop_coviddeaths2.pdf'),p,w=10,h= 12)
  ggsave(file=paste0(args$out.base,'_pop_prop_coviddeaths2.png'),p,w=10,h= 12)
  
}

if(args$make.sexpr & 'manaus' %in% args$selected.locs){
  tmp <- round(drpop[week == max(week) & loc_label == 'manaus' & age.label %in% c('85-89', '90+'),1-p.alive,]*100,2)
  sexpr[['Pop_at_risk']][['manaus_covid_deaths']][['85_89']] <- as.character(tmp[1])
  sexpr[['Pop_at_risk']][['manaus_covid_deaths']][['90+']] <- as.character(tmp[2])
  sexpr[['Pop_at_risk']][['last date']] <- write.date(max(drpop$week.start))
}

#	now we can subset the deaths to the time frame needed
tmp <- ddates[, list(week=seq(deaths.start, deaths.end, 1)), by='loc_label']
dd <- merge(dd, tmp, by=c('loc_label','week'))
set(dd, NULL, c('cadjdeaths','cdeathshospadm','scadjdeaths','scdeathshospadm'), NULL)



####################################################################
# investigate reference week & fit loess smoother
####################################################################

cat('\ninvestigate reference week ...')

# Define as useful even in postprocessing
dage <- dpop[,list(pop=sum(pop)),by='age.label']
dage[, tpop:=sum(pop)]
dage[, ppop.all.cities:=pop/tpop]
setkey(dage,age.label)
dage[,age.idx:=1:nrow(dage)]

# Left out comparison of UNVAXED and UNVAXED+UNKNOWN pops, but easy to do
# Compute HFR both in UNVACCINATED and in the UNVACCINATED|UNKNOWN populations
tmp <- subset(dh, grepl(args$keep.hospital.type, Hospital_type) & grepl(args$keep.vac.type, vac_type))
tmp <- tmp[, list(hosps=sum(hosps),deaths_of_hosps=sum(deaths_of_hosps)), by=c('loc_label','age.label','week')]
tmp1 <- tmp[,list(week=min(week):max(week)),by='loc_label'][,list(age.label=da[,unique(age.label)]),by=c('loc_label','week')]
tmp <- merge(tmp, tmp1, by=c('loc_label', 'week','age.label'), all.y=TRUE)
tmp[, hfr_empirical := deaths_of_hosps/hosps]
tmp[is.na(hosps), `:=`(hosps=0, deaths_of_hosps=0)]

tmp2 <- data.frame()
for (loc in ddates$loc_label){
  for (age in unique(da$age.label)){
    
    tmp3 <-  tmp[loc_label == loc & age.label == age]
    tmp3_noNA <- tmp3[!is.na(hfr_empirical)]
    tmp3 <- data.table(loc_label=loc, age.label=age, week=min(tmp$week):max(tmp$week),
                       hfr_loess30=predict(loess(hfr_empirical~week, tmp3_noNA, span = 0.3), min(tmp$week):max(tmp$week)),
                       hfr_loess40=predict(loess(hfr_empirical~week, tmp3_noNA, span = 0.4), min(tmp$week):max(tmp$week)),
                       hfr_loess50=predict(loess(hfr_empirical~week, tmp3_noNA, span = 0.5), min(tmp$week):max(tmp$week)),
                       hfr_loess75=predict(loess(hfr_empirical~week, tmp3_noNA, span = 0.75), min(tmp$week):max(tmp$week)))
    tmp2 <- rbind(tmp2, tmp3)
  }
}

# adjustment for weeks and age groups with 0 cases
tmp <- merge(tmp, tmp2, by = c("loc_label", "age.label", 'week'))
setkey(tmp, loc_label, age.label, week)
tmp2 <- tmp[,{
  a <- which(is.na(hfr_empirical)); b <- which(!is.na(hfr_empirical));
  l <- length(hfr_empirical);
  idx <- sapply(1:l, FUN = function(x){b[which.min(abs(b-x))]})
  list(week=week,
    hfr_empirical.adj=hfr_empirical[idx])
}, by = c('loc_label', 'age.label')]
tmp <- merge(tmp, tmp2, by=c('loc_label', 'age.label', 'week') )


# Write table with longest period of hfr > .5 per age&loc label
f <- function(x){
  y <- rle(as.numeric(x))
  len <- y$lengths[y$value == 1]
  if (length(len) > 0) max(len) else as.integer(0)
}
tmp1 <- tmp[,.(loc_label, age.label, week,hfr_loess30)]
tmp1 <- tmp1[is.na(hfr_loess30), hfr_loess30 := 0]
tmp1 <- tmp1[, list(V1 = f(hfr_loess30>.5)) , by = c('loc_label','age.label')]
tmp1 <- dcast(tmp1, loc_label ~ age.label, value.var = 'V1')
tmp1[,loc_label:= change_city_label(loc_label)]
tmp1[, (names(tmp1)):=lapply(.SD, as.character)]
saveRDS(tmp1, file=paste0(args$out.base,'_high_hfr_table.rds'))

tmp <- merge(tmp, unique(dw[,.(week.start,week)]), by='week')
tmp <- merge(tmp, dage[,.(age.label, ppop.all.cities)], by = 'age.label')

TMP <- dpop[,list(p=pop/tpop),by=c('loc_label','age.label')]
tmp <- merge(tmp, TMP, by=c('loc_label','age.label'))

tmp <- tmp[, list(hfr_loess30 = sum(hfr_loess30*ppop.all.cities),
                  hfr_loess30_nonadj = sum(hfr_loess30*p),
                  hfr_loess40 = sum(hfr_loess40*ppop.all.cities),
                  hfr_loess50 = sum(hfr_loess50*ppop.all.cities),
                  hfr_loess75 = sum(hfr_loess75*ppop.all.cities),
                  hfr_empirical=sum(hfr_empirical*ppop.all.cities),
                  hfr_empirical.adj=sum(hfr_empirical.adj*ppop.all.cities)),by = c('loc_label', 'week', 'week.start')]

rm(tmp2, tmp3, TMP)
dloess <- copy(tmp)

if(args$make.preprocessing.plots)
{
  tmp1 <- melt(tmp, measure.vars =  c('hfr_empirical', 'hfr_empirical.adj'))
  tmp2 <- tmp1[variable=='hfr_empirical' & !is.na(value), list(min=min(week),max=max(week)) ,by='loc_label']
  tmp1 <- merge(tmp1, tmp2, by='loc_label')[week > max | week < min, hfr_loess30 := NaN]

  
  
  p <- ggplot(data = tmp1, aes(x = week.start, color = loc_label)) + 
    geom_vline(data=ddates, aes(xintercept=w2.start), pch='11', colour='black', linetype='dashed') +
    geom_point(aes(y=value, pch = variable, size = variable)) + 
    geom_line(aes(y=hfr_loess30)) +# , linetype=shQuote('smoothed fit using the loess method'))) +
    # geom_line(data=pads, aes(y=M, linetype='posterior fit using health care indicators'))+
    facet_wrap(~change_city_label(loc_label)) + 
    theme_bw() +
    guides(colour=FALSE, size=FALSE, pch=guide_legend(nrow=2, byrow=TRUE), linetype=guide_legend(nrow=2, byrow=TRUE)) +
    scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    scale_size_manual(values=c(2,1)) + 
    scale_shape_manual(values=c(16,15),labels=c('empirical data','empirical data, no hospitalisations in at least one age group'))+
    scale_x_date(breaks="2 months") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.background = element_blank(),
          legend.position='bottom') +
    labs(x = '', y = 'weekly age-standardised COVID-19 attributable in-hospital fatality rates', shape='', linetype='')
  
  ggsave(file=paste0(args$out.base,'_allloc_loess_empirical.pdf'),p,w=10,h=12)
  ggsave(file=paste0(args$out.base,'_allloc_loess_empirical.png'),p,w=10,h=12)
}

# I need to store reference week as well as min hfr
tmp1 <- tmp[!is.na(hfr_empirical), list(min=min(week), max=max(week)) ,by=loc_label]
tmp1 <- merge(tmp, tmp1, by='loc_label')
tmp1 <- tmp1[week <= max & week >= min]
tmp1 <- merge(tmp1, ddates[,.(loc_label, w2.start)])
dhr <- tmp1[week.start<=w2.start, 
             {
              z <- which.min(hfr_loess30); 
              list(hfr_loess.min=min(hfr_loess30), ref.week=week[z], ref.week.start=week.start[z])
             }, by='loc_label']
if(0) # comparison between age-standardised relative to aggregate pop and city-specific city
{
  dhr2 <- tmp1[week.start<=w2.start, 
               {
                 z <- which.min(hfr_loess30_nonadj); 
                 list(hfr_loess.min_nonadj=min(hfr_loess30_nonadj), ref.week=week[z], ref.week.start=week.start[z])
               }, by='loc_label']
  dhr2 <- merge(dhr, dhr2, by='loc_label')
  dhr2 <- dhr2[loc_label == 'manaus', .(hfr_loess.min, hfr_loess.min_nonadj)]
}


tmp1 <- merge(tmp1[,list(hfr_loess.max=max(hfr_loess30)),by='loc_label'], dhr, by='loc_label')
tmp1[,hfr_loess.fold:=hfr_loess.max/hfr_loess.min,by='loc_label']
saveRDS(tmp1, file=paste0(args$out.base,'_minhfr_methods_table.rds'))


if(args$make.preprocessing.plots)
{
  
  p <- ggplot(tmp, aes(x=week.start)) + 
    geom_vline(data=ddates, aes(xintercept=w1.start), colour='grey50') +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_line(aes(y=hfr_loess30)) + 
    geom_hline(data=tmp1, aes(yintercept=hfr_loess.min)) +
    facet_wrap(~change_city_label(loc_label), ncol=4, scale='free') +
    scale_x_date(breaks='2 months') +
    expand_limits(y=0) +
    labs(x='', y='smoothed in-hospital fatality rate') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')
  ggsave(file=paste0(args$out.base,'_allloc_loessmins.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_allloc_loessmins.png'),p,w=10,h=10)
  
}

####################################################################
# make hospital collapse predictors
####################################################################

gc()

dhfrt <- subset(dh, 
		grepl(args$keep.hospital.type, Hospital_type) & 
		grepl(args$keep.vac.type, vac_type)
	)
dhfrt[, biweek := week%/%2]
tmp <- dhfrt[, 
	list(
		week.start= week.start[1],
		hfr= ifelse(sum(hosps)==0, 0., sum(deaths_of_hosps) / sum(hosps)) 
	), 
	by=c('loc_label','biweek', 'age.label')]
# age standardise biweekl in-hospital fatality rates
# at the moment we are not considering the problem that 0 deaths in age group loc underestimate as hfr
dage <- dpop[,list(pop=sum(pop)),by='age.label']
dage[, tpop:=sum(pop)]
dage[, ppop.all.cities:=pop/tpop][,.(age.label,ppop.all.cities)]
setkey(dage,age.label)
tmp <- merge(tmp, dage, by = 'age.label')
tmp <- tmp[, list(hfr = sum(hfr*ppop.all.cities)), by=c('loc_label', 'biweek')]
dhfrt <- dhfrt[, 
             list(
               week.start= week.start[1],
               hfr_notsd= ifelse(sum(hosps)==0, 0., sum(deaths_of_hosps) / sum(hosps)) 
             ), 
             by=c('loc_label','biweek')]
dhfrt <- merge(dhfrt,tmp,  by = c('loc_label', 'biweek'))

# icu admissions (all cause) in selected hospital type
cat('\nmake hospital collapse predictor ICU admissions ...')
dicu <- subset(dha, 
		hosp.cat!='out of hospital' &	
		!is.na(DateICUstart) &
		grepl(args$keep.hospital.type, Hospital_type)
	)
setnames(dicu, 'DateICUstart', 'Date')
dicu <- merge(dicu, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dicu <- merge(dicu, tmp, by=c('loc_label','week'))
dicu <- dicu[, list( nicu.adm = as.numeric(length(Date))), by=c('week','loc_label')]

# hosp admissions 75+ and hosp admissions (all cause) in selected hospital type
cat('\nmake hospital collapse predictor hospital admissions ...')
dhcadm <- subset(dha, 
    hosp.cat =='in hospital' &	
		# !is.na(DateAdmission) &
		grepl(args$keep.hospital.type, Hospital_type) # sens analysis doesn t make so much sense because we are changing demand but not supply in pandemic loads
	)
setnames(dhcadm, 'DateAdmission', 'Date')
dhcadm <- merge(dhcadm, dw, by=c('Date'))
dhcadm <- dhcadm[, list( nhosp.adm= as.numeric(length(Date)), nhosp.adm.75plus= as.numeric(length(Date[Age>=75])) ), by=c('week','loc_label')]
tmp <- ddates[,list(week=seq(hosps.start, deaths.end, 1)),by='loc_label']
dhcadm <- merge(tmp, dhcadm, by=c('week','loc_label'), all.x=TRUE)

# ICU beds, physicians, nurses, ventilators
cat('\nmake hospital collapse predictor ICU beds and physicians ...')
dhcb <- as.data.table(read.csv(args$file.hospital.demand.beds.physicians, stringsAsFactors=FALSE))

set(dhcb, NULL, c('X'), NULL)
setnames(dhcb, c('capital', 'competen', 'physician', 'physb', 'anestesio', 'nurse','tec_nurse','physicalterapist', 'vent', 'ICUBR', 'IU'),
         c('loc_label','month.start','physicians.all', 'physicians.specialist' ,'anesthesiologist', 'nurse','tech.nurse','physiotherapist', 'ventilator', 'ICUBr', 'intermed_beds'),
         skip_absent = TRUE) # suggest to comment this out while testing
setnames(dhcb, c('interm_beds', 'isola_beds', 'beds'), c('interm.beds.all', 'isola.beds.all', 'beds.all' ),
         skip_absent = TRUE) # not setting other names as we would have set them in d1 because of incoherences


# # Compute input ventilators as Luciana suggested.
# # To the 'ventilator' column, add the number of ventilators associated to the
# # minimum requirements for ICU beds II and III and intermediate
tmp <- function(x){ if(is.na(x)){0}else{x}}

if(!'ICUBr' %in% names(dhcb) )
{
  dhcb[, ICUBr := sapply(typeII_ICU,tmp) + sapply(typeIII_ICU,tmp) +sapply(typeII_covid,tmp) ]
}

# add ventilators reserved to ICU beds and intermediate beds to ventilators:
if(1)
{
  dhcb[, ventilator := ventilator+ ceiling(ICUBr/2) +ceiling(interm.beds.all/3) ]
}
# saribeds WITH or WITHOUT ventilator
dhcb[,saribeds:= typeII_covid + typeI_ICU + typeII_ICU + typeIII_ICU + interm.beds.all]

tmp <- c('typeI_ICU', 'typeII_ICU','typeIII_ICU','typeII_covid', 'isola.beds.all', 'interm.beds.all', 'intcrit_cbeds', 'intcrit_beds')
tmp <- c(tmp, 'LSVP', 'size_totalbeds')
tmp <- tmp[which(tmp %in% names(dhcb))]
dhcb[, (tmp):= NULL,]

# set icu/saribeds.all as the icu/sari beds with available ventilator
setnames(dhcb, c('saribeds','icubeds.with.ventilator','saribeds.with.ventilator'),c('saribeds.all','icubeds.all', 'saribedsvent.all'))

dhcb <- dhcb[loc_label %in% args$selected.locs]
set(dhcb, NULL, 'month.start', dhcb[, gsub('([0-9]{4})([0-9]{2})','\\1-\\2-01', month.start)] )
set(dhcb, NULL, 'month.start', dhcb[, as.Date(month.start)])
stopifnot( dhcb[, !any(is.na(physicians.all))] )
stopifnot( dhcb[, !any(is.na(physicians.specialist))] )
stopifnot( dhcb[, !any(is.na(intensivists))] )
stopifnot( dhcb[, !any(is.na(icubeds.all))] )
stopifnot( dhcb[, !any(is.na(ventilator))] )
stopifnot( dhcb[, !any(is.na(nurse))] )
stopifnot( dhcb[, !any(is.na(tech.nurse))] )
stopifnot( dhcb[, !any(is.na(physiotherapist))] )
tmp <- unique(subset(dpop, select=c(loc_label, tpop)))
dhcb <- merge(dhcb, tmp, by='loc_label')
dhcb[, r.physicians.all := physicians.all / tpop * 1e5]
dhcb[, r.physicians.specialist := physicians.specialist / tpop * 1e5]
dhcb[, r.intensivists := intensivists / tpop *1e5]
dhcb[, r.icubeds.all := icubeds.all / tpop * 1e5]
dhcb[, r.saribeds.all := saribeds.all / tpop * 1e5]
dhcb[, r.saribedsvent.all := saribedsvent.all / tpop * 1e5]
dhcb[, r.beds.all := beds.all / tpop * 1e5]
dhcb[, r.ventilator := ventilator / tpop * 1e5]
dhcb[, r.nurse := nurse / tpop * 1e5]
dhcb[, r.tech.nurse := tech.nurse / tpop * 1e5]
dhcb[, r.physiotherapist := physiotherapist / tpop * 1e5]

if(args$make.sexpr)
{
  tmp <- min(ddates$w2.end)
  tmp <- as.Date(gsub("^(.*?)-(.*?)-(.*?)$", "\\1-\\2-01", tmp))
  tmp <- dhcb[month.start == '2020-03-01' | month.start == tmp]
  # For each location, extract resources per 10^5 inhabitants in March 2020 and July 2021
  # But we have for April 2020 to May 2020
  
  # first personnel, (alphabetically), then equipment (alphabetically)
  tmp1 <- sort(grep(pattern = '^r.',names(tmp), value = T))
  tmp2 <- c("r.saribeds.all", "r.saribedsvent.all","r.icubeds.all", "r.ventilator")  # rm beds.all?
  tmp1 <- tmp1[!tmp1 %in% tmp2]
  tmp1 <- unique(tmp[, c('loc_label', 'month.start', ..tmp1)])
  tmp2 <- unique(tmp[, c('loc_label', 'month.start', ..tmp2)])
  tmp <- merge(tmp1, tmp2, by = c('loc_label', 'month.start'))
  tmp1 <- names(tmp)[!names(tmp) %in% c('loc_label', 'month.start')]
  tmp[,  (tmp1) := lapply(.SD, FUN=as.integer, MARGIN = 2) ,.SDcols=tmp1]
  tmp <- dcast(tmp,  loc_label ~ month.start, value.var=tmp1)
  tmp <- t(tmp)
  colnames(tmp) <- tmp['loc_label',] 
  tmp <- tmp[2:nrow(tmp),]
  
  sexpr[['Resources']] <- tmp
}

tmp <- unique(subset(dw, select=c(week, week.start, month.start)))
tmp <- tmp[, list(week.start=week.start[1], month.start=min(month.start)), by='week']
dhcb <- merge(dhcb,tmp,by='month.start',allow.cartesian=TRUE)
dhcb <- subset( dhcb, loc_label%in%args$selected.locs )
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dhcb <- merge(tmp, dhcb, by=c('loc_label','week'), all.x=TRUE)
dhcb <- merge(dhcb, dicu, by=c('loc_label','week'), all.x=TRUE)
dhcb <- merge(dhcb, dhcadm, by=c('loc_label','week'), all.x=TRUE)
set(dhcb, dhcb[, which(is.na(nicu.adm))], 'nicu.adm', 0)
set(dhcb, dhcb[, which(is.na(nhosp.adm))], 'nhosp.adm', 0)
set(dhcb, dhcb[, which(is.na(nhosp.adm.75plus))], 'nhosp.adm.75plus', 0)
dhcb[, icu.per.beds := nicu.adm / icubeds.all]
dhcb[, icu.per.physicians.specialist := nicu.adm / physicians.specialist]
dhcb[, icu.per.intensivists := nicu.adm / intensivists]
dhcb[, icu.per.physicians.all := nicu.adm / physicians.all]
dhcb[, icu.per.ventilator := nicu.adm / ventilator]
dhcb[, icu.per.nurse := nicu.adm / nurse]
dhcb[, icu.per.tech.nurse := nicu.adm / tech.nurse]
dhcb[, icu.per.physiotherapist := nicu.adm / physiotherapist]
dhcb[, hosps.per.physicians.all := nhosp.adm / physicians.all]
dhcb[, hosps.per.saribeds := nhosp.adm / saribeds.all ]
dhcb[, hosps.per.saribedsvent := nhosp.adm / saribedsvent.all ]
dhcb[, hosps.per.allbeds := nhosp.adm / beds.all ]
dhcb[, hosps.per.ventilator := nhosp.adm / ventilator]
dhcb[, hosps.per.nurse := nhosp.adm / nurse]
dhcb[, hosps.per.tech.nurse := nhosp.adm / tech.nurse]
dhcb[, hosps.per.physiotherapist := nhosp.adm / physiotherapist]

tmp <- unlist(lapply(dhcb, is.integer))
tmp <- names(tmp[tmp == T])
tmp <- tmp[!(tmp %in% c('week'))]
dhcb[, (tmp):=lapply(.SD, FUN = as.numeric), .SDcols = tmp]
tmp <- unlist(lapply(dhcb, is.numeric))
tmp <- names(tmp[tmp == T])
tmp <- tmp[!(tmp %in% c('week'))]
for(x in tmp)
{
  set(dhcb, which(is.nan(dhcb[[x]])), x, 0.)	
}

# proportion of non-residents in all-cause hospital admissions
cat('\nmake hospital collapse predictor proportion of non-residents ...')
dhcnres <- subset(dha, hosp.cat!='out of hospital' & grepl(args$keep.hospital.type, Hospital_type))
setnames(dhcnres, 'DateAdmission', 'Date')
dhcnres <- merge(dhcnres, dw, by=c('Date'))
dhcnres <- dhcnres[, list(hosps=as.numeric(length(Sex))), by=c('loc_label','week','loc.res.cat')]
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
tmp[, DUMMY:=1L]
tmp <- merge(tmp,data.table(DUMMY=1L,loc.res.cat=unique(dhcnres$loc.res.cat)),by='DUMMY',allow.cartesian=TRUE)
dhcnres <- merge(tmp, dhcnres, by=c('loc_label','week','loc.res.cat'), all.x=TRUE)
set(dhcnres,dhcnres[,which(is.na(hosps))],'hosps',0)
set(dhcnres,NULL,'DUMMY',NULL)
tmp <- dhcnres[, list(thosps=sum(hosps)), by=c('loc_label','week')]
dhcnres <- merge(dhcnres, tmp, by=c('loc_label','week'))
dhcnres[,phosps:=hosps/thosps]
set(dhcnres, dhcnres[, which(is.nan(phosps))],'phosps',0)
dhcnres <- subset(dhcnres, loc.res.cat=='not resident', select=-c(hosps,thosps,loc.res.cat))
setnames(dhcnres, 'phosps', 'p.nonres.hosps')
dhcnres[, p.res.hosps := 1-p.nonres.hosps]

# out of hospital deaths
cat('\nmake hospital collapse predictor out of hospital deaths ...')
dooh <- subset(dha, !is.na(DateDeath) & hosp.cat=='out of hospital')
setnames(dooh, 'DateDeath', 'Date')
dooh <- dooh[, list(sivep.deaths.ooh = length(Outcome)), by=c('loc_label','Date')]
dooh <- merge(dooh, dw, by=c('Date'))
dooh <- dooh[, list(sivep.deaths.ooh = sum(sivep.deaths.ooh)), by=c('loc_label','week')]
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dooh <- merge(tmp, dooh, by=c('loc_label','week'), all.x=TRUE)
set(dooh, dooh[, which(is.na(sivep.deaths.ooh))], 'sivep.deaths.ooh', 0)
tmp <- dr[, list( reg.deaths.ooh = sum(reg.deaths.ooh)), by=c('week', 'loc_label')]
dooh <- merge(dooh, tmp, by=c('week', 'loc_label'), all.x=TRUE )
set(dooh, dooh[, which(is.na(reg.deaths.ooh))], 'reg.deaths.ooh', 0)
dooh <- merge(dooh, unique(subset(dw, select=c(week,week.start))), by='week')
dooh[, deaths.ooh := pmax(reg.deaths.ooh,sivep.deaths.ooh)]
tmp <- unique(subset(dpop, select=c(loc_label, tpop)))
dooh<- merge(dooh, tmp, by='loc_label')
dooh[, r.deaths.ooh := deaths.ooh/tpop*1e5]

if(args$make.preprocessing.plots)
{		
  tmp <- merge(dhcb, dhfrt, by=c('loc_label','week.start'), all.x=TRUE)
  tmp <- merge(tmp,subset(dooh,select=c(loc_label,week,r.deaths.ooh)),by=c('loc_label','week'))
  tmp <- merge(tmp,subset(dhcnres,select=c(loc_label,week,p.res.hosps)),by=c('loc_label','week'))
  tmp <- subset(tmp, select=c(loc_label,week,week.start,month.start,r.physicians.all, r.intensivists, r.icubeds.all,r.saribeds.all, r.ventilator, r.nurse, r.tech.nurse, r.physiotherapist))
  tmp <- melt(tmp,id.vars=c('loc_label','week','week.start','month.start'))
  tmp2 <- data.table(variable=c('r.physicians.all','r.intensivists','r.icubeds.all','r.saribeds.all','r.ventilator','r.nurse','r.tech.nurse','r.physiotherapist'),
                     variable_label=factor(c('Physicians','Intensive Care\nspecialists','adult ICU\nbeds','Critical care\nbeds','Ventilators','Nurses','Nurse assistants','Physiotherapists'),
                                           levels = c('Physicians','Intensive Care\nspecialists','adult ICU\nbeds','Critical care\nbeds','Ventilators','Nurses','Nurse assistants','Physiotherapists')))
  tmp <- merge(tmp,tmp2,by='variable')
  tmp2 <- subset(dstates,select=c(loc_label,region_label))
  set(tmp2,tmp2[,which(region_label%in%c('Central-West','North'))],'region_label','North + Central-West')
  tmp <- merge(tmp,tmp2,by='loc_label')
  tmp <- unique(tmp,by=c('loc_label','month.start','variable'))
  p <- ggplot(tmp, aes(x=month.start)) +			
    geom_step(aes(y=value, colour=change_city_label(loc_label)),lwd=1) +
    #geom_point(aes(y=value, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous() +
    scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    facet_grid(variable_label~region_label, scales='free_y') +			
    guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    labs(x='',y='resource per 100,000',colour='') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          panel.spacing.x = unit(1, "lines"),
          legend.position='bottom')							
  ggsave(file=paste0(args$out.base,'_hfr_resources.pdf'),p,w=10,h=12)
  ggsave(file=paste0(args$out.base,'_hfr_resources.png'),p,w=10,h=12)					
}



###############################################################################
# Build table of all predictors
###############################################################################

#	start building data.table of all predictors
if('hfr_notsd' %in% names(dhfrt)){
  dhfrt_notsd <- copy(dhfrt)
  dhfrt[, hfr_notsd := NULL]
}

tmp <- make_health_care_predictors(dhcadm, dicu, dooh, dhcnres, dhcb, dw, dr, dpop, ddates, dhr, args)
dhca <- copy(tmp$dhca)
dhcc <- copy(tmp$dhcc)
dhccw <- copy(tmp$dhccw)



if(args$make.sexpr)
{
  
  tmp <- dhca[ variable %in% args$fr.predictors]
  tmp <- merge(tmp, dloess[,.(loc_label, week.start, hfr_empirical, hfr_empirical.adj)], by=c('loc_label', 'week.start'))
  
  tmp1 <- tmp[, list(cor.adj = cor(value ,hfr_empirical.adj),
                     cor     = cor(value, hfr_empirical, use="na.or.complete")), by = c('loc_label', 'variable') ]
  
  tmp <-  dcast(tmp1, loc_label ~ variable, value.var = 'cor')
  tmp1 <- colnames(tmp)[2:ncol(tmp)]
  for(i in 1:length(tmp1)){
    tmp2 <- names(which(args$fr.predictors == tmp1[i])); tmp1[i] <- tmp2
  }
  colnames(tmp)[2: length(colnames(tmp))] <- tmp1
  
  tmp[,loc_label:= change_city_label(loc_label)]
  colnames(tmp) <- gsub('^ ','', colnames(tmp) )
  tmp[,2:ncol(tmp)] <- round(tmp[,2:ncol(tmp)] , 2)
  setcolorder(tmp)
  
  sexpr[['correlations']] <- tmp
}

# get log predictors
tmp <- dhca[, list(loc_label = loc_label, week = week, value_log = log(value+1)), by='variable']
dhca <- merge(dhca, tmp, by=c('variable','loc_label','week'))

# center each variable relative to reference week with lowest HFR in each location
tmp <- dhca[, 
            {
              z <- which(week.cat=='before.P1')		
              value_ref <- value[week==ref.week]
              value_log_ref <- value_log[week==ref.week]
              stopifnot(!is.na(value_ref))
              list(
                week=week,
                value_ref = value_ref,
                value_c= value-value_ref,
                value_log_ref = value_log_ref,
                value_log_c= value_log-value_log_ref						
              )	
            }, 
            by=c('variable','loc_label')]
dhca <- merge(dhca, tmp, by=c('variable','loc_label','week'))

# standardise across location to range of -3,3 for each variable
tmp <- dhca[, 
            list(		 
              week = week,
              value_sd = sd(value_c),
              value_cs = value_c/sd(value_c),
              value_log_sd = sd(value_log_c),
              value_log_cs = value_log_c/sd(value_log_c)
            ), 
            by=c('variable','loc_label')]
dhca <- merge(dhca, tmp, by=c('variable','loc_label','week'))

if(args$make.preprocessing.plots)
{
  tmp <- subset(dhca, variable%in%args$fr.predictors)
  tmp <- merge(tmp,data.table(variable=args$fr.predictors, variable.label=names(args$fr.predictors)),by='variable')
  
  p <- ggplot(tmp, aes(x=week.start)) +
    geom_line(aes(y=value, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    labs(x='', y='health care demand predictor', colour='') +
    facet_wrap(~variable.label, ncol=3, scales='free') +	
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')
  ggsave(file=paste0(args$out.base,'_hfr_predictors.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_hfr_predictors.png'),p,w=10,h=10)
  
  p <- ggplot(tmp, aes(x=week.start)) +
    geom_line(aes(y=value_c, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    labs(x='', y='centered standardised health care demand predictor', colour='') +
    facet_wrap(~variable.label, ncol=3, scales='free') +	
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')
  ggsave(file=paste0(args$out.base,'_hfr_predictors_centered.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_hfr_predictors_centered.png'),p,w=10,h=10)
  
  p <- ggplot(tmp, aes(x=week.start)) +
    geom_line(aes(y=value_log_c, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    labs(x='', y='centered standardised log health care demand predictor', colour='') +
    facet_wrap(~variable.label, ncol=3, scales='free') +	
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')
  ggsave(file=paste0(args$out.base,'_hfr_predictors_centeredlog.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_hfr_predictors_centeredlog.png'),p,w=10,h=10)
  
  p <- ggplot(tmp, aes(x=week.start)) +
    geom_line(aes(y=value_cs, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    labs(x='', y='standardised health care demand predictor', colour='') +
    facet_wrap(~variable.label, ncol=3, scales='free') +	
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')
  ggsave(file=paste0(args$out.base,'_hfr_predictors_std.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_hfr_predictors_std.png'),p,w=10,h=10)
  
  p <- ggplot(tmp, aes(x=week.start)) +
    geom_line(aes(y=value_log_cs, colour=change_city_label(loc_label))) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
    labs(x='', y='standardised log health care demand predictor', colour='') +
    facet_wrap(~variable.label, ncol=3, scales='free') +	
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')
  ggsave(file=paste0(args$out.base,'_hfr_predictors_stdlog.pdf'),p,w=10,h=10)
  ggsave(file=paste0(args$out.base,'_hfr_predictors_stdlog.png'),p,w=10,h=10)	
}

dhc <- subset(dhca, variable%in%args$fr.predictors )
dhc <- dcast.data.table(dhc, loc_label+week+week.start+week.cat~variable, value.var=args$fr.predictors.transformation)
stopifnot( all(args$fr.predictors %in%colnames(dhc)) )

####################################################################
# PLOT HFRs
####################################################################

# time trends in biweekly HFR by type
# and boxplots in biweekly HFR
cat('\nplot empirical in-hospital fatality rates ...')

if(args$make.preprocessing.plots)
{      
  tmp <- dh[Hospital_type!='Unknown']
  tmp[, biweek := week%/%2]	  
  tmp <- tmp[,
             list(
               week.start = week.start[1],
               deaths_of_hosps_bw = sum(deaths_of_hosps),
               deaths_of_hosps_ntcs_bw = sum(deaths_of_hosps_ntcns),
               hosps_bw = sum(hosps),
               hosps_ntcns_bw = sum(hosps_ntcns)
             ), 
             by = c('loc_label', 'biweek', 'age.label', 'Hospital_type')]		  
  tmp[, hfr_bw := deaths_of_hosps_bw/hosps_bw]
  tmp[, hfr_ntcns_bw := deaths_of_hosps_ntcs_bw/hosps_bw]
  tmp <- merge(tmp, subset(ddates, select = c('loc_label','w1.start', 'w2.start')), by = 'loc_label')
  
  if(args$partition.locs){
    tmp1 <- tmp[week.start>w1.start & loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp1, aes(x = week.start)) +
      geom_vline(aes(xintercept = w2.start), col = 'grey50') +
      geom_point(aes(y=hfr_bw, color = Hospital_type)) +
      scale_x_date(breaks='2 months', expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      ggsci::scale_colour_npg() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')	+
      labs(x = '', y = 'biweekly in-hospital fatality rate', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time_part1.pdf'),p,w=6,h=8)
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time_part1.png'),p,w=6,h=8)
    
    tmp1 <- tmp[week.start>w1.start & loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp1, aes(x = week.start)) +
      geom_vline(aes(xintercept = w2.start), col = 'grey50') +
      geom_point(aes(y=hfr_bw, color = Hospital_type)) +
      scale_x_date(breaks='2 months', expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      ggsci::scale_colour_npg() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')	+
      labs(x = '', y = 'biweekly in-hospital fatality rate', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time_part2.pdf'),p,w=6,h=8)
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time_part2.png'),p,w=6,h=8)
    
    tmp1 <- tmp[week.start>w1.start & loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp1, aes(x = week.start)) +
      geom_vline(aes(xintercept = w2.start), col = 'grey50') +
      geom_point(aes(y=hfr_bw, color = Hospital_type)) +
      scale_x_date(breaks='2 months', expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      ggsci::scale_colour_npg() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')	+
      labs(x = '', y = 'biweekly in-hospital fatality rate', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time_part3.pdf'),p,w=6,h=8)
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time_part3.png'),p,w=6,h=8)
  }else{
    p <- ggplot(tmp, aes(x = week.start)) +
      geom_vline(aes(xintercept = w2.start), col = 'grey50') +
      geom_point(aes(y=hfr_bw, color = Hospital_type)) +
      scale_x_date(breaks='2 months', expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      ggsci::scale_colour_npg() +
      facet_grid(age.label~change_city_label(loc_label)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            strip.background = element_blank(),
            legend.position='bottom')	+
      labs(x = '', y = 'biweekly in-hospital fatality rate', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time.pdf'),p,w=6,h=8)
    ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type_over_time.png'),p,w=6,h=8)
  }
  
  tmp <- tmp[, list(value=quantile(hfr_ntcns_bw, p=c(0.025,.25,.5,.75,.975)), stat=c('CL','IL','M','IU','CU')), by=c('loc_label','age.label','Hospital_type')]
  tmp <- dcast.data.table(tmp, loc_label+age.label+Hospital_type~stat, value.var='value')
  p <- ggplot(tmp, aes(x=age.label)) +		  
    geom_boxplot(aes(min=CL, lower=IL, middle=M, upper=IU, max=CU, fill=Hospital_type), stat='identity') +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    ggsci::scale_fill_aaas() +
    facet_wrap(~change_city_label(loc_label), ncol=5, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')	+
    labs(x = '', y = 'biweekly in-hospital fatality rate', fill='')
  ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type.pdf'),p,w=12,h=10)
  ggsave(file=paste0(args$out.base,'_hfr_by_hospital_type.png'),p,w=12,h=10)      
}

####################################################################
# after the figures we can now keep only the selected hospital types
####################################################################

cat('\nAggregate hospital admissions ...')
dh <- subset(dh, 
             grepl(args$keep.hospital.type, Hospital_type) &
               grepl(args$keep.vac.type, vac_type)				
)
dh <- dh[, 
         list(
           hosps=sum(hosps), 
           deaths_of_hosps=sum(deaths_of_hosps),
           hosps_ntcns=sum(hosps_ntcns),
           deaths_of_hosps_ntcns=sum(deaths_of_hosps_ntcns)	
         ), 
         by=c('loc_label', 'week','age.label')]
tmp <- ddates[,list(week=seq(hosps.start, deaths.end, 1)),by='loc_label']
tmp[, DUMMY:=1L]
tmp2 <- as.data.table(expand.grid(DUMMY=1L, age.label=unique(da$age.label)))
tmp <- merge(tmp, tmp2, by='DUMMY', allow.cartesian=TRUE)
set(tmp, NULL, 'DUMMY', NULL)
tmp <- merge(tmp, unique(subset(dw, select=c(week,week.start))), by=c('week'))
dh <- merge(tmp, dh, all.x=TRUE, by=c('loc_label', 'week','age.label'))
set(dh, dh[, which(is.na(hosps))], 'hosps', 0L)
set(dh, dh[, which(is.na(deaths_of_hosps))], 'deaths_of_hosps', 0L)
set(dh, dh[, which(is.na(hosps_ntcns))], 'hosps_ntcns', 0L)
set(dh, dh[, which(is.na(deaths_of_hosps_ntcns))], 'deaths_of_hosps_ntcns', 0L)
dh <- dh[order(loc_label,age.label,week),]
tmp <- dh[, list(week=week, chosps=cumsum(hosps)), by=c('loc_label', 'age.label')]
dh <- merge(dh, tmp, by=c('loc_label', 'week','age.label'))
tmp <- dh[,list(age.label=age.label, schosps=cumsum(chosps)), by=c('loc_label', 'week')]
dh <- merge(dh, tmp, by=c('loc_label', 'week','age.label'))
tmp <- dh[,list(thosps=sum(hosps)), by=c('loc_label', 'week')]
dh <- merge(dh, tmp, by=c('loc_label', 'week'))
dh[, phosps:=hosps/thosps]
tmp <- dh[,list(age.label=age.label, sphosps=cumsum(phosps)), by=c('loc_label', 'week')]
dh <- merge(dh, tmp, by=c('loc_label', 'week','age.label'))
dh[, deaths_of_hospsi := as.integer(round(deaths_of_hosps))]

if(args$make.preprocessing.plots)
{
  # plot empirical cum deaths
  p <- ggplot(dh, aes(x=week.start)) +
    geom_ribbon(aes(ymin=schosps-chosps, ymax=schosps, fill=age.label)) +
    viridis::scale_fill_viridis(option='D', end=0.8, discrete = TRUE) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(fill='age', y='cumulated COVID-19 hospitalisations', x='') +
    theme_bw() +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots, scale = 'free_y') + 
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom')	
  ggsave(file=paste0(args$out.base,'_hospsc.pdf'),p,w=8,h=10)
  ggsave(file=paste0(args$out.base,'_hospsc.png'),p,w=8,h=10)
  
  if(args$partition.locs){
    
    # plot empirical proportion deaths (facets)
    tmp <- dh[loc_label %in% args$selected.locs.part1,]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=phosps, colour=age.label, pch=change_city_label(loc_label))) +
      viridis::scale_colour_viridis(option='D', end=0.8, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(y='share of age groups among COVID-19 hospitalisations', x='') +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y')
    ggsave(file=paste0(args$out.base,'_hospsp_facet_part1.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_hospsp_facet_part1.png'),p,w=10,h=12)
    
    tmp <- dh[loc_label %in% args$selected.locs.part2,]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=phosps, colour=age.label, pch=change_city_label(loc_label))) +
      viridis::scale_colour_viridis(option='D', end=0.8, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(y='share of age groups among COVID-19 hospitalisations', x='') +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y')
    ggsave(file=paste0(args$out.base,'_hospsp_facet_part2.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_hospsp_facet_part2.png'),p,w=10,h=12)
    
    tmp <- dh[loc_label %in% args$selected.locs.part3,]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=phosps, colour=age.label, pch=change_city_label(loc_label))) +
      viridis::scale_colour_viridis(option='D', end=0.8, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(y='share of age groups among COVID-19 hospitalisations', x='') +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y')
    ggsave(file=paste0(args$out.base,'_hospsp_facet_part3.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_hospsp_facet_part3.png'),p,w=10,h=12)
  }else{

    p <- ggplot(dh, aes(x=week.start)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=phosps, colour=age.label, pch=change_city_label(loc_label))) +
      viridis::scale_colour_viridis(option='D', end=0.8, discrete = TRUE) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      labs(y='share of age groups among COVID-19 hospitalisations', x='') +
      guides(colour= 'none', pch='none') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +
      facet_grid(age.label~change_city_label(loc_label), scales='free_y')
    ggsave(file=paste0(args$out.base,'_hospsp_facet.pdf'),p,w=10,h=12)
    ggsave(file=paste0(args$out.base,'_hospsp_facet.png'),p,w=10,h=12)
  }
 
  ### 
  tmp <- dh[, list(deaths= sum(deaths_of_hosps), discharged= sum(hosps)-sum(deaths_of_hosps) ) , by=c('loc_label','week','week.start')]
  tmp <- melt(tmp, id.vars=c('loc_label','week','week.start'), measure.vars=c('deaths','discharged'))
  set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('deaths','discharged'), labels=c('COVID19 death','discharged'))])
  p <- ggplot(tmp, aes(x= week.start)) +
    geom_bar(aes(y=value, fill=variable), stat='identity', position=position_stack()) +
    ggsci::scale_fill_npg() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols , scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position = 'bottom') +	
    labs(x='epiweek', y='Hospital admissions among residents', fill='expected outcome\nas of database closure')
  ggsave(file=paste0(args$out.base,'_hospn_by_expoutcome_residents.pdf'),p,w=8,h=10, limitsize = F)
  ggsave(file=paste0(args$out.base,'_hospn_by_expoutcome_residents.png'),p,w=8,h=10, limitsize = F)
  
  p <- ggplot(dh, aes(x=week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_bar(aes(y=hosps, fill=age.label), width=6, stat='identity', position='stack') +		
    scale_y_continuous(expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    viridis::scale_fill_viridis(option='D', end=.8, discrete = TRUE) +
    theme_bw() +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots, scale = 'free_y') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.background = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          legend.position='bottom') +			
    labs(x='', y='hospital admissions', fill='age')
  ggsave(file=paste0(args$out.base,'_hospsn_over_time.pdf'),p,w=8,h=10,limitsize = FALSE)
  ggsave(file=paste0(args$out.base,'_hospsn_over_time.png'),p,w=8,h=10,limitsize = FALSE)
  
  
  tmp <- subset(dh, select=c(week, age.label, loc_label, week.start, hosps, hosps_ntcns))
  tmp[, hosps_cns := hosps-hosps_ntcns]
  tmp <- melt(tmp, id.vars=c('week','age.label','loc_label','week.start'), measure.vars=c('hosps_cns','hosps_ntcns'))
  set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('hosps_ntcns','hosps_cns'), labels=c('hospital admissions with observed outcome','hospital admissions with unknown outcome'))])
  p <- ggplot(tmp, aes(x=week.start)) +
    geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
    geom_bar(aes(y=value, fill=variable), width=6, stat='identity', position='stack') +		
    scale_y_continuous(expand=c(0,0)) +
    #scale_x_date(breaks='month Porto Alegre (RS), Rio de Janeiro (RJ) and Salvador (BA).s',expand=c(0,0)) +
    ggsci::scale_fill_npg() +
    theme_bw() +
    facet_wrap(~change_city_label(loc_label), ncol = args$ncols.plots, scale = 'free_y') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          strip.background = element_blank(),
          legend.position='bottom') +		
    labs(x='', y='COVID-19 hospital admissions', fill='')
  ggsave(file=paste0(args$out.base,'_hospn_over_time_bycensoring.pdf'),p,w=8,h=10,limitsize = FALSE)
  ggsave(file=paste0(args$out.base,'_hospn_over_time_bycensoring.png'),p,w=8,h=10,limitsize = FALSE)
}

####################################################################
# EMPIRICAL HOSPITAL FATALITY RATIO OVER TIME
####################################################################

cat('\ncalculate empirical in-hospital fatality rate over time ...')
dh[, hfr := deaths_of_hosps/hosps]
dh[, hfr_ntcns := deaths_of_hosps_ntcns/hosps_ntcns]

# biweekly hfr
dh[, biweek := week %/%2]
tmp <- dh[, 
          list(
            week.start = week.start[1],
            deaths_of_hosps_bw = sum(deaths_of_hosps),
            deaths_of_hosps_ntcs_bw = sum(deaths_of_hosps_ntcns),
            hosps_bw = sum(hosps),
            hosps_ntcns_bw = sum(hosps_ntcns)
          ), 
          by = c('loc_label', 'age.label', 'biweek')]
tmp[,hfr_bw := deaths_of_hosps_bw/hosps_bw]
tmp[,hfr_ntcns_bw := deaths_of_hosps_ntcs_bw/hosps_bw]
set(dh, NULL, 'biweek', NULL)

if(0 & args$make.preprocessing.plots)
{
  # manaus HFR over time by age group
  tmp1 <- dh[loc_label == 'manaus' & deaths_of_hospsi!=0]
  p <- ggplot(tmp1, aes(x=week.start)) +
    geom_vline(data=ddates[loc_label == 'manaus'], aes(xintercept=w2.start), colour='grey50') +
    # geom_line(aes(y=hfr, colour=age.label)) +
    geom_point(aes(y=hfr, colour=age.label, size=deaths_of_hospsi)) +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    scale_x_date(breaks='2 months',expand=c(0,0)) +
    viridis::scale_colour_viridis(option='B', end=.85, discrete = TRUE) +
    coord_cartesian(ylim=c(0,1)) +
    theme_bw() +
    guides(colour=guide_legend(nrow=2,byrow=TRUE), size=guide_legend(color='white')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.background = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          legend.position='bottom') +		
    facet_wrap(~age.label, ncol = 6) +
    labs(x='', y='hospital fatality ratio', colour='age', size='weekly deaths')
  ggsave(file=paste0(args$out.base,'manaus_hfr_age.png'),p,w=12,h=8)
  ggsave(file=paste0(args$out.base,'manaus_hfr_age.pdf'),p,w=12,h=8)
  
  if (args$partition.locs){
    # Set 1
    tmp <- dh[loc_label %in% args$selected.locs.part1,]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=hfr, colour=age.label)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label)) +
      labs(x='', y='hospital fatality ratio', colour='age')
    ggsave(file=paste0(args$out.base,'_hfr_over_time_part1.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time_part1.png'),p,w=10,h=10)
    
    tmp <- dh[loc_label %in% args$selected.locs.part2,]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=hfr, colour=age.label)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label)) +
      labs(x='', y='hospital fatality ratio', colour='age')
    ggsave(file=paste0(args$out.base,'_hfr_over_time_part2.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time_part2.png'),p,w=10,h=10)
    
    tmp <- dh[loc_label %in% args$selected.locs.part3,]
    p <- ggplot(tmp, aes(x=week.start)) +
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=hfr, colour=age.label)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label)) +
      labs(x='', y='hospital fatality ratio', colour='age')
    ggsave(file=paste0(args$out.base,'_hfr_over_time_part3.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time_part3.png'),p,w=10,h=10)
    
    # Set 2
    tmp <- subset(dh, select=c(week, age.label, loc_label, week.start, hfr, hfr_ntcns))  
    tmp <- melt(tmp, id.vars=c('week','age.label','loc_label','week.start'), measure.vars=c('hfr','hfr_ntcns'))
    
    tmp1 <- tmp[loc_label %in% args$selected.locs.part1]
    p <- ggplot(tmp1, aes(x=week.start)) +
      geom_vline(data=ddates1, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=value, colour=variable)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      ggsci::scale_colour_npg() +		  
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='hospital fatality ratio', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring_part1.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring_part1.png'),p,w=10,h=10)
    
    tmp1 <- tmp[loc_label %in% args$selected.locs.part2]
    p <- ggplot(tmp1, aes(x=week.start)) +
      geom_vline(data=ddates2, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=value, colour=variable)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      ggsci::scale_colour_npg() +		  
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='hospital fatality ratio', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring_part2.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring_part2.png'),p,w=10,h=10) 
    
    tmp1 <- tmp[loc_label %in% args$selected.locs.part3]
    p <- ggplot(tmp1, aes(x=week.start)) +
      geom_vline(data=ddates3, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=value, colour=variable)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      ggsci::scale_colour_npg() +		  
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='hospital fatality ratio', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring_part3.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring_part3.png'),p,w=10,h=10) 
  }else{

    p <- ggplot(dh, aes(x=week.start)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=hfr, colour=age.label)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      viridis::scale_colour_viridis(option='D', end=.8, discrete = TRUE) +
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label)) +
      labs(x='', y='hospital fatality ratio', colour='age')
    ggsave(file=paste0(args$out.base,'_hfr_over_time.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time.png'),p,w=10,h=10)
    
    tmp <- subset(dh, select=c(week, age.label, loc_label, week.start, hfr, hfr_ntcns))  
    tmp <- melt(tmp, id.vars=c('week','age.label','loc_label','week.start'), measure.vars=c('hfr','hfr_ntcns'))

    p <- ggplot(tmp, aes(x=week.start)) +
      geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
      geom_point(aes(y=value, colour=variable)) +    
      scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
      scale_x_date(breaks='2 months',expand=c(0,0)) +
      ggsci::scale_colour_npg() +		  
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            strip.background = element_blank(),
            panel.spacing.y = unit(1, "lines"),
            legend.position='bottom') +		
      facet_grid(age.label~change_city_label(loc_label), scales='free_y') +
      labs(x='', y='hospital fatality ratio', colour='')
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring.pdf'),p,w=10,h=10)
    ggsave(file=paste0(args$out.base,'_hfr_over_time_bycensoring.png'),p,w=10,h=10)
  }
}


####################################################################
# SELECTION CRITERIA FOR HFR CALCULATIONS
####################################################################

dhs <- merge(dh, ddates, by='loc_label')
dhs[, variant := factor(week.start<w2.start, levels=c(TRUE,FALSE), labels=c('non-Gamma','Gamma'))]		
dhs <- dhs[, list(hosps=sum(hosps)), by=c('loc_label','age.label','variant')]
dhs <- merge(dhs, dpop, by=c('loc_label','age.label'))
dhs[, select.by.sample.size := as.integer(hosps>=100)]

if(args$make.preprocessing.plots)
{
  p <- ggplot(dhs, aes(x=loc_label, y=hosps, fill=variant)) +
    geom_bar(stat='identity', position=position_dodge()) +
    facet_wrap(~age.label, ncol=4, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.background = element_blank()) +
    scale_y_log10() +
    ggsci::scale_fill_aaas() +
    labs(x='', y='total COVID-19 attributable hospitalisations', fill='')
  ggsave(file=paste0(args$out.base,'_hfr_sample_size.pdf'),p,w=12,h=10)
  ggsave(file=paste0(args$out.base,'_hfr_sample_size.png'),p,w=12,h=10) 
}

####################################################################
# EMPIRICAL HOSPITAL FATALITY RATIO FOR FIRST WAVE
####################################################################

cat('\ncalculate empirical hospital fatality rate for first wave...')
tmp <- ddates[, 
              {
                z1 <- seq.Date(w1.start, w1.end, by='day')
                z2 <- seq.Date(w2.start, w2.end, by='day')
                list(DateAdmission=c(z1,z2))
              }, 
              by='loc_label']
dhf <- merge(dht, tmp, by = c('DateAdmission', 'loc_label'))
dhf <- subset(dht, 		
              hosp.cat!='out of hospital' &
                grepl(args$keep.hospital.type,Hospital_type)
)
dhf <- dhf[, list( logit_hfr= gtools::logit( sum(wdeaths)/sum(whosps) ) ), by= c('loc_label', 'age.label')]

####################################################################
# READ PROPORTION P1 AND MAKE PRIOR
####################################################################

cat('\nread genomic data...')
dp <- as.data.table( read.csv(args$file.genomic.data, stringsAsFactors = FALSE) )
dp <- subset(dp, grepl(args$genomic.data.type,data_type))
dp <- subset(dp, n_sampled > 0) # To avoid NAs!
dp[, Date := as.Date(week, format = '%d/%m/%Y')]
tmp <- data.table(  
  location = c('manaus', 'porto alegre', 'belo horizonte', 'sao paulo city', 'acre', 'amapa', 'amazonas', 'parana', 'roraima', 'alagoas', 'bahia', 'maranhao', 'paraiba', 'pernambuco', 'piaui', 'rio grande do norte', 'distrito federal', 'mato grosso do sul', 'mato grosso', 'espirito santo', 'minas gerais', 'rio de janeiro', 'sau paulo state', 'rio grande do sul', 
               'santa catarina', 'tocantins', 'ceara', 'goias', 'rondonia', 'sergipe', 'sao paulo'), # The last sao paulo here refers to the state
  loc_label = c('manaus', 'porto alegre', 'belo horizonte', 'sao paulo', 'rio branco', 'macapa', 'manaus', 'curitiba', 'boa vista', 'maceio', 'salvador', 'sao luis', 'joao pessoa', 'recife', 'teresina', 'natal', 'brasilia', 'campo grande', 'cuiaba', 'vitoria', 'belo horizonte', 'rio de janeiro', 'sao paulo', 'porto alegre', 
                'florianopolis', 'palmas', 'fortaleza', 'goiania', 'porto velho', 'aracaju', 'sao paulo')  		
)
dp <- merge(dp, tmp, by = 'location')
dp <- subset(dp, select = c(-week, -location))
dp <- merge(dp,dw, by= 'Date')
tmp <- subset(ddates, select=c(loc_label, hosps.start))
setnames(tmp, 'hosps.start', 'week')
tmp <- merge(tmp, unique(subset(dw, select=c(week, week.start))), by='week')
setnames(tmp, 'week.start', 'hosps.week.start')
dp <- merge(dp, subset(tmp, select=c(loc_label, hosps.week.start)), by='loc_label')
dp[, date.idx := as.numeric(Date - hosps.week.start + 1L)]
dp[, prop_P1 := n_positive/n_sampled]

if(args$make.preprocessing.plots)
{
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
          strip.background = element_blank()) +
    labs(x='', y='number of weekly SARS-CoV-2 sequences', fill='')
  ggsave(file=paste0(args$out.base,'_seq_sample_size.pdf'),p,w=12,h=10)
  ggsave(file=paste0(args$out.base,'_seq_sample_size.png'),p,w=12,h=10)
  
  tmp <- subset(dp, select = c('loc_label', 'week.start', 'n_positive', 'n_sampled', 'prop_P1'))
  tmp <- merge(tmp, subset(ddates, select = c('loc_label', 'w2.start')), by = 'loc_label')
  tmp <- merge(tmp, dstates, by = 'loc_label')
  tmp$state_label <- factor(tmp$state_label, levels=c('Minas Gerais','Paraná','Santa Catarina','Goiás',
                                                      'Paraíba', 'Amapá','Amazonas', 'Rio Grande do Norte',
                                                      'Rio Grande do Sul', 'Rondônia', 'Rio de Janeiro', 'Bahia',
                                                      'Maranhão','São Paulo' ))
  
  p <- ggplot(tmp, aes(x = week.start, color = loc_label)) + 
    geom_point(aes(y = prop_P1, size = n_sampled)) +
    facet_wrap(.~state_label, ncol = args$ncols.plots) + 
    # viridis::scale_color_viridis(option='B', end=0.85, discrete = TRUE) +
    scale_colour_manual(values = args$city.palette(14)) +
    scale_x_date(breaks='4 weeks',expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    guides(colour = 'none') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background = element_blank(), legend.position =  'bottom') +	
    labs(x='weeks', y='Proportion of sampled P1 sequences', size = 'Total weekly samples', color = '')
  
  ggsave(file=paste0(args$out.base,'_p1_prop.pdf'),p,w=10,h=12, limitsize = F)
  ggsave(file=paste0(args$out.base,'_p1_prop.png'),p,w=10,h=12, limitsize = F)
}

if(args$make.sexpr)
{
  saveRDS(sexpr,paste0(args$out.base,'_substitute_expressions.rds'))
  #sexpr <- readRDS(paste0(args$out.base,'_substitute_expressions.rds'))
}

cat('\nSaving to ', paste0(args$out.base,'_allloc_workspace.rda'))
save.image(file=paste0(args$out.base,'_allloc_workspace.rda'))

cat('\nDone HFR.preprocessing.210719.R\n')

# erased empirical manaus, but if needed get it from before

if(args$make.preprocessing.plots)
{
  tmp <- dht[hosp.cat == 'in hospital'  & Outcome2 != 'discharged',, ]
  tmp <- merge(tmp, subset(ddates, select = c('loc_label', 'w1.start', 'w1.end')), by = 'loc_label')
  tmp <- tmp[DateDeath >= w1.start & DateDeath < w1.end,
             list(wdeaths = sum(wdeaths)),
             by = c('loc_label', 'age.label')]
  tmp[, t_wdeaths := sum(wdeaths),by = 'loc_label']
  tmp <- tmp[, list(p_wdeaths = unname(as.numeric(Hmisc::binconf(wdeaths,t_wdeaths))), stat=c('M','CL','CU')), by=c('loc_label','age.label')]
  tmp <- dcast(tmp, loc_label + age.label ~ stat, value.var = "p_wdeaths")
  tmp[, type := 'in-hospital deaths']
  
  tmp1 <- dr[,c('loc_label', 'age.label', 'reg.deaths.ooh', 'week.start')]
  tmp1 <- tmp1[, list(reg.deaths.ooh = sum(reg.deaths.ooh)), by = c('loc_label','age.label', 'week.start')]
  tmp1 <- merge(tmp1, subset(ddates, select = c('loc_label', 'w1.start', 'w1.end')), by = 'loc_label')
  tmp1 <- tmp1[week.start >= w1.start & week.start < w1.end,]
  tmp1 <- tmp1[,list(ooh_deaths = sum(reg.deaths.ooh)),by = c('loc_label', 'age.label')]
  tmp1[, t_ooh_deaths := sum(ooh_deaths),by = 'loc_label']
  tmp1 <- tmp1[, list(p_ooh_deaths = unname(as.numeric(Hmisc::binconf(ooh_deaths,t_ooh_deaths))), stat=c('M','CL','CU')), by=c('loc_label','age.label')]
  tmp1 <- dcast(tmp1, loc_label + age.label ~ stat, value.var = "p_ooh_deaths")
  tmp1[, type := 'out of hospital deaths']
  
  tmp <- rbind(tmp, tmp1)
  
  p <- ggplot(tmp, aes(x = age.label, color = type)) + 
    geom_point(aes(y = M))+
    geom_errorbar(aes(ymin = CL, ymax = CU)) + 
    facet_wrap(~change_city_label(loc_label), ncol= args$ncols.plots) + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    #scale_x_date(breaks='2 months',expand=c(0,0)) +
    ggsci::scale_colour_npg() +		  
    #coord_cartesian(ylim=c(0,1)) +
    theme_bw() +
    guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.background = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          legend.position='bottom') +		
    labs(x='', y='proportion of COVID-19 attributable deaths by age group', colour='')
  
  ggsave(file=paste0(args$out.base,'_deathcompage_inouthosp.pdf'),p,w=8,h=10,limitsize = FALSE)
  ggsave(file=paste0(args$out.base,'_deathcompage_inouthosp.png'),p,w=8,h=10,limitsize = FALSE)
}


