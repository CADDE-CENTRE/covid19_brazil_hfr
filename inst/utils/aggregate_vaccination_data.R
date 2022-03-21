cat(" \n -------------------------------- \n \n aggregate_vaccination_data.R \n \n -------------------------------- \n")

require(data.table)

###############################################################################
###  Requirements
################################################################################
# In order to run this script, it is first neccessary to download the vaccination data
# from the Brazilian ministry of health, and then run preprocess_vaccination_data.R in inst/utils

if(1)
{
  pkg.dir <- '~/git/covid19.P1.IFR/'
  file.vaccinations <- file.path(pkg.dir, 'inst/data', 'saudegovbr_vaccinations_210805.rds')
}

# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-filename')
  file.vaccinations <- args_line[[2]]
  stopifnot(args_line[[3]]=='-pkgdir')
  pkg.dir <- args_line[[4]]
} 

###############
# set args
###############
args <- list()
args$selected.locs <- c('belo horizonte','curitiba', "florianopolis", "goiania", "joao pessoa", "macapa", 'manaus', 
                        "natal", 'porto alegre', "porto velho", 'rio de janeiro', 'salvador', 'sao paulo', "sao luis")
args$vaccine.population <- 'loc_res_label'
args$file.hospital.collapse.times <- file.path(pkg.dir, 'inst/data', 'hospital_collapse_cities_list.csv')
args$deaths.weeks.increment <- 0L
tmp <- 74 + 9 # set this as the week.end for every state (before 74)?
args$select.deaths.week.end <- c("aracaju"=tmp, "belem"=tmp, "belo horizonte"=tmp ,"boa vista"=tmp, "brasilia"=tmp, "campo grande"=tmp, "cuiaba"=tmp, "curitiba"=tmp, "florianopolis"=tmp, 
                                 "fortaleza"=tmp, "goiania"=tmp, "joao pessoa"=tmp, "macapa"=tmp, "maceio"=tmp, "manaus"=tmp, "natal"=tmp, "palmas"=tmp, "porto alegre"=tmp,  
                                 "porto velho"=tmp, "recife"=tmp, "rio branco"=tmp, "rio de janeiro"=tmp, "salvador"=tmp, "sao luis"=tmp, "sao paulo"=tmp, "teresina"=tmp, "vitoria"=tmp)


####################################################################
#  make age brackets and weeks datatable
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
if( subset(dw, week==max(dw$week))[, max(Date)] < subset(dw, week==max(dw$week))[1, week.start+7-1] )
{
  dw <- subset(dw, week<max(dw$week))
  stopifnot( max(dw$week) >= args$select.deaths.week.end )
}

ddates <- read.csv(args$file.hospital.collapse.times)
ddates <- data.table(
  loc_label = tolower(ddates$City),
  w1.start = as.Date(ddates$W1Start, '%d/%m/%Y'),
  w1.end = as.Date(ddates$W1End,'%d/%m/%Y'),
  w2.start = as.Date(ddates$W2Start,'%d/%m/%Y'),
  w2.end = as.Date('2021-07-26') 
)
# exception to BH, the zero deaths are just weird:
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

# set start and end time of the death weeks in our model
ddates[, deaths.start:= hosps.start+args$deaths.weeks.increment]
tmp <- unique(subset(dw, select=c(week,week.start)))
setnames(tmp,c('week','week.start'),c('deaths.start','deaths.week.start'))
ddates <- merge(ddates, tmp, by='deaths.start')
tmp <- data.table( loc_label=names(args$select.deaths.week.end), deaths.end=as.numeric(args$select.deaths.week.end))
ddates <- merge(ddates, tmp, by='loc_label')
tmp <- unique(subset(dw, select=c(week,week.start)))
setnames(tmp,c('week','week.start'),c('deaths.end','deaths.week.end'))
ddates <- merge(ddates, tmp, by='deaths.end')

#################################################### 
# Aggregate the individual level vaccination data
####################################################
cat('\nread vaccination data ...')

# dv <- readRDS(args$file.vaccinations)
if(file.exists(file.vaccinations))
{
  dv <- readRDS(file.vaccinations)
}else{
  tmp1 <- gsub('(.*?).rds$','\\1_1.rds', file.vaccinations)
  tmp2 <- gsub('(.*?).rds$','\\1_2.rds', file.vaccinations)
  tmp1 <- readRDS(tmp1)
  tmp2 <- readRDS(tmp2)
  dv <- rbind(tmp1, tmp2)
  rm(tmp1, tmp2)
}

set(dv, dv[, which(vaccine %in% c('Covishield', 'AstraZeneca'))], 'vaccine', 'Covishield')
dv[is.na(loc_label), loc_res_label:= loc_label ,]
dv <- subset(dv, age<=120)
dv <- merge(da, dv, by.x='Age', by.y='age')
if(args$vaccine.population=='loc_res_label and loc_label')
{
  dv <- subset(dv, loc_res_label%in%args$selected.locs & loc_res_label==loc_label)
  dv <- dv[, list(vac=length(patient)), by=c('loc_label','age.label','Date', 'vaccine')]
  dv <- dcast(dv, loc_res_label + age.label+ Date +  vaccine ~  dose, value.var = 'vac')
  setnames(dv, c('1','2'), c("vac1", "vac2"))
  set(dv, dv[,which(is.na(vac1))],'vac1',0L)
  set(dv, dv[,which(is.na(vac2))],'vac2',0L)
  # --> manaus  289134
}
if(args$vaccine.population=='loc_res_label')
{
  dv <- subset(dv, loc_res_label%in%args$selected.locs)
  dv <- dv[, list(vac=length(patient)), by=c('loc_res_label','age.label','Date', 'dose','vaccine')]
  dv <- dcast(dv, loc_res_label + age.label+ Date +  vaccine ~  dose, value.var = 'vac')
  setnames(dv, c('1','2'), c("vac1", "vac2"))
  set(dv, dv[,which(is.na(vac1))],'vac1',0L)
  set(dv, dv[,which(is.na(vac2))],'vac2',0L)
  setnames(dv, 'loc_res_label', 'loc_label')
  # --> manaus  312873
}
if(args$vaccine.population=='loc_label')
{
  dv <- subset(dv, loc_label%in%args$selected.locs)
  dv <- dv[, list(vac=length(patient)), by=c('loc_label','age.label','Date', 'vaccine')]
  # --> manaus  318098
  dv <- dcast(dv, loc_res_label + age.label+ Date +  vaccine ~  dose, value.var = 'vac')
  setnames(dv, c('1','2'), c("vac1", "vac2"))
  set(dv, dv[,which(is.na(vac1))],'vac1',0L)
  set(dv, dv[,which(is.na(vac2))],'vac2',0L)
}

# some birthdates were most likely reported as vaccination dates, remove them:
# (3 entries)
dv <- subset(dv, Date>="2020-01-01")
dv <- merge(dv, dw, by=c('Date'))
tmp <- ddates[, list(week=seq(hosps.start, deaths.end, 1)), by='loc_label']
dv <- merge(dv, tmp, by=c('loc_label','week'))

####################################################################
# aggregate vaccinations
####################################################################

cat('\naggregate vaccinations...')
dv <- dv[,list(vac1=sum(vac1), vac2 = sum(vac2)),by=c('loc_label','age.label','week', 'vaccine')]
tmp <- dv[, as.data.table(expand.grid(week=seq(min(week),max(week),1), age.label=unique(age.label), vaccine = unique(vaccine))), by='loc_label']
dv <- merge(tmp, dv, by=c('loc_label','age.label','week', 'vaccine'), all.x=TRUE)
set(dv, dv[, which(is.na(vac1))], 'vac1', 0L)
set(dv, dv[, which(is.na(vac2))], 'vac2', 0L)
setkey(dv, loc_label, age.label, week)
tmp <- dv[, list(week=week, cvac1=cumsum(vac1), cvac2=cumsum(vac2)),by=c('loc_label','age.label', 'vaccine')] 
dv <- merge(dv, tmp, by=c('loc_label','age.label', 'vaccine','week'))
tmp <- dv[,which(cvac1<cvac2)]
if(length(tmp))
{
  warning('unexpected cvac1<cvac2, cumulative second doses cannot be larger than cumulative first doses')
  print(dv[tmp])
  set(dv, tmp, 'cvac2', dv[tmp,cvac1])
}
set(dv, NULL, 'cvac1', dv[, cvac1-cvac2])
dv <- merge(dv, unique(dw[,.(week, week.start)]), by = 'week')

tmp1 <- gsub('saudegovbr_vaccinations_(.*?).rds', '\\1', basename(file.vaccinations))
tmp <- paste0('aggregated_vaccinations_',tmp1, '.rds')

saveRDS(dv, file.path(pkg.dir, 'inst/data', tmp))

