require(data.table)
require(ggplot2)
library("readxl")

pkg.dir <- "~/git/covid19_brazil_hfr/"

################################################################################
# 28 of September
###############################################################################

# Need to update previous data to include intensivists from Luciana s latest mail
tmp <- file.path(pkg.dir, 'inst/data/IPEA_ICUbeds_physicians_210906.csv')
tmp <- as.data.table(read.csv(tmp))
tmp1 <- file.path(pkg.dir, 'inst/data/physicians_intensivists_professionals.xlsx')
tmp1 <- as.data.table(read_xlsx(tmp1))
tmp1 <- melt(tmp1,id.vars=c('codibge','Capital'), variable.name = 'competen',value.name = 'intensivists')
tmp1[, capital:= tolower(gsub('^[0-9]+ ','', Capital)) ]
tmp1[, capital:= tolower(gsub('^\\s|\\s$','',gsub('é','e',gsub('í','i',gsub('ó','o',gsub('á|ã|â','a',gsub('\xe3','a',capital)))))))]
tmp1[, Capital:= NULL]
set(tmp1, NULL, 'competen', tmp1[, gsub('([0-9]{4})([0-9]{2})','\\1-\\2-01', competen)] )
tmp <- merge(tmp,tmp1, by=c('capital', 'competen'))

if(0)
{
  # Exploratory plot:
  dpop <- as.data.table(read.csv(file.path(pkg.dir, 'inst/data', 'PNADc_populationpyramids_210617.csv'),
                                 stringsAsFactors=FALSE, encoding="latin1"))
  setnames(dpop, c('city','population'), c('loc_label','pop'))
  set(dpop, NULL, 'loc_label', dpop[, tolower(gsub('^\\s|\\s$','',gsub('é','e',gsub('í','i',gsub('ó','o',gsub('á|ã|â','a',gsub('\xe3','a',loc_label)))))))])
  dpop <- dpop[,list(pop=sum(pop)),by='loc_label']
  
  tmp1 <- merge(tmp, dpop, by.x='capital', by.y='loc_label')
  tmp1[, r.intensivists:=intensivists/pop*1e5]
    
  p <- ggplot(tmp1, aes(x=as.Date(competen), y=r.intensivists)) +
    geom_point() + 
    facet_wrap(~capital) + 
    labs(x='', y='Intensivists per 100,000 inhabitants')
  
  ggsave('~/Documents/P1Brazil/tmp/intensivists_per_pop.png',p,  w=10, h=8)
}

write.csv(tmp, file.path(pkg.dir, 'inst/data/IPEA_ICUbeds_physicians_210928.csv'))


###############################################################################
# 6th of September stuff
###############################################################################
file <- file.path(pkg.dir,  "/instdata/Imperial_beds_vent_cnes_20210904.xlsx")

tmp <- excel_sheets(file)
dhc_all <- as.data.table(read_excel(file, sheet = "beds_vent_cnes_all"))
dhc_comp <- as.data.table(read_excel(file, sheet = "beds_vent_cnes_complement"))

names(dhc_all)
names(dhc_comp)
dhc_all[, table(unit_type1)]
dhc_comp[, table(unit_type)]
# Comp has no NA unit_types!
nrow(dhc_all[is.na(unit_type1)])
nrow(dhc_comp[is.na(unit_type)])
# LSVP = Leitos de suporte ventilatorio polmunar.

# Let s work with comp, I think it should be more appropriate
NAto0 <- function(x){if(is.na(x)){return(0)};return(x)}

set(dhc_comp, NULL, 'capital', dhc_comp[, tolower(gsub('^[0-9]+ ','',capital))])
set(dhc_comp, NULL, 'capital', dhc_comp[, tolower(gsub('^\\s|\\s$','',gsub('é','e',gsub('í','i',gsub('ó','o',gsub('á|ã|â','a',gsub('\xe3','a',capital)))))))])
set(dhc_all, NULL, 'capital', dhc_all[, tolower(gsub('^[0-9]+ ','',capital))])
set(dhc_all, NULL, 'capital', dhc_all[, tolower(gsub('^\\s|\\s$','',gsub('é','e',gsub('í','i',gsub('ó','o',gsub('á|ã|â','a',gsub('\xe3','a',capital)))))))])


# Take minimum after aggregation
tmp <- c("typeII_covid","typeI_ICU","typeII_ICU", "typeIII_ICU", "interm_beds", "vent_inuse")
tmp1 <- dhc_comp[, lapply(.SD, FUN=function(x){sum(x,na.rm=T)}) , .SDcols=tmp, ,by = c('capital', 'competen')]
tmp1[,saribeds:= typeII_covid + typeI_ICU + typeII_ICU + typeIII_ICU + interm_beds]
tmp1[,icubeds:= typeII_covid  + typeII_ICU + typeIII_ICU ]
tmp1 <- tmp1[, list(saribeds.with.ventilator = NAto0(ceiling(typeII_covid/2)) +  NAto0(ceiling(typeII_ICU/2)) +  NAto0(ceiling(typeIII_ICU/2)) + NAto0(ceiling(interm_beds/3)),
            saribeds.without.ventilator = NAto0(typeI_ICU) +  NAto0(floor(typeII_covid/2)) +  NAto0(floor(typeII_ICU/2)) +  NAto0(floor(typeIII_ICU/2)) + NAto0(floor(2*interm_beds/3)),
            icubeds.with.ventilator = NAto0(ceiling(typeII_covid/2)) +  NAto0(ceiling(typeII_ICU/2)) +  NAto0(ceiling(typeIII_ICU/2)),
            icubeds.without.ventilator = NAto0(floor(typeII_covid/2)) +  NAto0(floor(typeII_ICU/2)) +  NAto0(floor(typeIII_ICU/2)),
            vent_inuse = vent_inuse),
      by =c('capital', 'competen')]

sprintf("Number of inadequate icubeds: %i",tmp1[icubeds.without.ventilator>vent_inuse, sum(icubeds.without.ventilator - vent_inuse)] )
sprintf("Number of inadequate saribeds: %i",tmp1[saribeds.without.ventilator>vent_inuse, sum(saribeds.without.ventilator - vent_inuse)] )

tmp1[, `:=` (saribeds.without.ventilator = min(saribeds.without.ventilator, vent_inuse),
             icubeds.without.ventilator = min(icubeds.without.ventilator, vent_inuse)),
     by =  c('capital','competen') ]
tmp1 <- tmp1[, list(saribeds.with.ventilator = saribeds.with.ventilator + saribeds.without.ventilator,
             icubeds.with.ventilator = icubeds.with.ventilator + icubeds.without.ventilator)
     ,by =  c('capital','competen')]


# Take minimum before aggregation
tmp2 <- dhc_comp[, list( saribeds.with.ventilator = NAto0(ceiling(typeII_covid/2)) +  NAto0(ceiling(typeII_ICU/2)) +  NAto0(ceiling(typeIII_ICU/2)) + NAto0(ceiling(interm_beds/3)),
                         saribeds.without.ventilator = NAto0(typeI_ICU) +  NAto0(floor(typeII_covid/2)) +  NAto0(floor(typeII_ICU/2)) +  NAto0(floor(typeIII_ICU/2)) + NAto0(floor(2*interm_beds/3)),
                         icubeds.with.ventilator = NAto0(ceiling(typeII_covid/2)) +  NAto0(ceiling(typeII_ICU/2)) +  NAto0(ceiling(typeIII_ICU/2)),
                         icubeds.without.ventilator = NAto0(floor(typeII_covid/2)) +  NAto0(floor(typeII_ICU/2)) +  NAto0(floor(typeIII_ICU/2)),
                         vent_inuse = vent_inuse),
                 by = c('capital','cnes','competen') ]
sprintf("Number of inadequate icubeds: %i",tmp2[icubeds.without.ventilator>vent_inuse, sum(icubeds.without.ventilator - vent_inuse)] )
sprintf("Number of inadequate saribeds: %i",tmp2[saribeds.without.ventilator>vent_inuse, sum(saribeds.without.ventilator - vent_inuse)] )
# "Number of inadequate icubeds: 10472";  "Number of inadequate saribeds: 29763"

tmp2[, `:=` (saribeds.without.ventilator = min(saribeds.without.ventilator, vent_inuse),
             icubeds.without.ventilator = min(icubeds.without.ventilator, vent_inuse)),
     by =  c('capital','cnes','competen') ]
tmp2 <- tmp2[, list(saribeds.with.ventilator = saribeds.with.ventilator + saribeds.without.ventilator,
                    icubeds.with.ventilator = icubeds.with.ventilator + icubeds.without.ventilator),
             by =  c('capital','cnes','competen') ]
tmp2 <- tmp2[, list(saribeds.with.ventilator = sum(saribeds.with.ventilator, na.rm = T),
                    icubeds.with.ventilator = sum(icubeds.with.ventilator, na.rm = T)),
             by =  c('capital','competen') ]


set(tmp1, NULL, 'competen', tmp1[, gsub('([0-9]{4})([0-9]{2})','\\1-\\2-01', competen)] )
set(tmp2, NULL, 'competen', tmp2[, gsub('([0-9]{4})([0-9]{2})','\\1-\\2-01', competen)] )
tmp1[, min := as.factor('after')]
tmp2[, min := as.factor('before')]
rbind(tmp1, tmp2) -> tmp


ggplot(data=tmp, aes(x = competen, y = saribeds.with.ventilator, color = min)) + 
  geom_point() + 
  facet_wrap(~capital, scale = 'free_y') + 
  scale_y_continuous(limits = c(0,NA)) + 
  geom_point(aes(x = competen, y = saribeds.with.ventilator))

tmp1 <- dcast.data.table(tmp, capital + competen ~ min, value.var = 'saribeds.with.ventilator')
tmp1 <- tmp1[, ratio := after/before][, max(ratio), by = capital]
setkey(tmp1, V1)
setnames(tmp1,'V1','sariratio')

tmp2 <- dcast.data.table(tmp, capital + competen ~ min, value.var = 'icubeds.with.ventilator')
tmp2 <- tmp2[, ratio := after/before][, max(ratio), by = capital]
setkey(tmp2, V1)
setnames(tmp2,'V1','icuratio')

tmp1 <- merge(tmp1, tmp2, by = 'capital')
if(0){
# most unequal by income according to Geographic Access:
tmp3 <- data.table(  capital = c( 'brasilia', 'Sao Paulo', 'belo horizonte', 'Fortaleza', 'Salvador', 'rio de janeiro', 'Goiania', 'Curitiba', 'Porto Alegre', 'Sao Luis', 'Natal', 'Manaus'))
tmp3[, geo_rank := 1:.N]
tmp3[, capital := tolower(capital)]

tmp1 <- merge(tmp3, tmp1, by = 'capital')
tmp1[, sarirank := order(sariratio, decreasing = T)]
tmp1[, icurank := order(icuratio, decreasing = T)]

ggplot(tmp1, aes(x=icurank, y=geo_rank, col=capital)) + 
  geom_point()+ 
  scale_y_continuous(limits=c(1,12)) + 
  scale_x_continuous(limits=c(1,12)) 
}


# aggregate micro data and join new columns
tmp <- tmp[min == 'before', -'min']

set(dhc_all, NULL, 'competen', dhc_all[, gsub('([0-9]{4})([0-9]{2})','\\1-\\2-01', competen)] )
set(dhc_comp, NULL, 'competen', dhc_comp[, gsub('([0-9]{4})([0-9]{2})','\\1-\\2-01', competen)] )

tmp_names <- names(dhc_all)[!names(dhc_all) %in% c("capital", "codibge", "cnes", "competen", "size_c", "size_d", "unit_type1", "unit_typ2")]
tmp1 <- dhc_all[, lapply(.SD, FUN=function(x){sum(x, na.rm = T)}), .SDcols=tmp_names, c("capital", "competen")]
tmp2 <- dhc_comp[, lapply(.SD, FUN=function(x){sum(x, na.rm = T)}), .SDcols=tmp_names, c("capital", "competen")]

# In summary, the difference btw all and comp will be exclusively in the # of vents in use 
# (and I guess health units sizes, which anyways we don t have here)
# In ptcular, all will have more vent_inuse
# I would use all, as this should lead to eq # of vents as before
tmp <- merge(tmp1, tmp, by = c("capital", "competen"))
tmp1 <- as.data.table(read.csv(paste0(pkg.dir,'inst/data/IPEA_ICUbeds_physicians_210824.csv') ))
set(tmp1, NULL, 'competen', tmp1[, gsub('([0-9]{4})([0-9]{2})','\\1-\\2-01', competen)] )
set(tmp1, NULL, 'capital', tmp1[, tolower(gsub('^[0-9]+ ','',capital))])
set(tmp1, NULL, 'capital', tmp1[, tolower(gsub('^\\s|\\s$','',gsub('é','e',gsub('í','i',gsub('ó','o',gsub('á|ã|â','a',gsub('\xe3','a',capital)))))))])
merge( tmp1[,.(capital,competen,vent)], tmp[,.(capital,competen,vent_inuse)]  ,by = c('capital','competen'))
# Indeed they are coherent!

# merge staff data:
tmp1 <- tmp1[,.(capital, competen, beds,physician,physicians.specialist,nurse,tec_nurse,physiotherapist)]
tmp <- merge(tmp, tmp1)

# Save
setnames(tmp, 'vent_inuse', 'vent')
write.csv(tmp, paste0(pkg.dir, 'inst/data/IPEA_ICUbeds_physicians_210906.csv'))

################################################################################
# 24 of August files
################################################################################

# files we use in prepro. Where did we get those from?
dprep1 <- as.data.table(read.csv(paste0(pkg.dir,'/inst/data/IPEA_ICUbeds_physicians_210721.csv'), stringsAsFactors=FALSE))
dprep2 <- as.data.table(read.csv(paste0(pkg.dir,'/inst/data/MOH_CNES_ventilators_210628.csv'), stringsAsFactors=FALSE))

# Latest Luciana file
dh <- as.data.table(read_excel('~/Downloads/imper_sixth.xlsx'), stringsAsFactors=FALSE)

# rename columns of interest in the same way as first dataset:
# setnames(dh, c('capital','competen','physicians','physiciansb','ICUBR', 'ventilators', 'nurses', 'tech_nurses', 'physiothera', 'isola_beds',  'intermed_beds'),
#              c('capital','competen','physician','physb','ICUBr','vent','nurse','tec_nurse','physicalterapist','isola','IU'))
setnames(dh, c('capital','competen','physicians','physiciansb', 'ventilators_inuse', 'nurses','tech_nurse', 'physiothera'),
         c('capital','competen','physician','physicians.specialist','vent','nurse', 'tec_nurse','physiotherapist'),
         skip_absent = T)


#save as csv in /data:
write.csv(dh, paste0(pkg.dir, 'inst/data/IPEA_ICUbeds_physicians_210824.csv'))

