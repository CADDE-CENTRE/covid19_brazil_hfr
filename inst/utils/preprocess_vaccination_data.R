####### LOG
# Now apply to all states
# Not a log change but a suggestion: use wget to download dados, so don t have to rename!

# https://opendatasus.saude.gov.br/dataset/covid-19-vacinacao/resource/ef3bd0b8-b605-474b-9ae5-c97390c197a8

# https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/PNI/vacina/uf/2021-10-07/uf%3DSE/part-00000-b9798f4e-eb5b-4362-9aa9-e7f42d77f38b.c000.csv
download_original <- FALSE
pkg.dir <- '~/git/covid19_brazil_hfr/'

if(download_original)
{
  options(timeout = Inf)
  tmp <- c("AP", "AM", 'BA', "GO", 'MA', "MG", 'PB','PR', 'RJ', 'RN', 'RS', 'RO', 'SC', 'SP')
  for(state in tmp)
  {
    # May need to modify the last '/part....csv'
    url = paste0('https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/PNI/vacina/uf/',
                 Sys.Date()-1,
                 '/uf%3D', state, '/part-00000-b9798f4e-eb5b-4362-9aa9-e7f42d77f38b.c000.csv')
    dest = paste0(pkg.dir, 'inst/data/dados_', state)
    download.file(url, dest)
  }
}



library(data.table)

# Helper function which takes vectors of strings and returns vectors of 1, 2, NAs
# This is quite slow, would be better if could gsub
# However probs with portuguese chars, for example can t write grepl(Unica, x)
extract_dose <- function(vec){
  i <- 1
  for (x in vec){
    if(grepl('1', x) | grepl('nica', x)){vec[i] <- 1} else if (grepl('2', x)){ vec[i] <- 2 }
    i <- i +1
  }
  return(as.integer(vec))
}

clean_vacc_data <- function(id, city, data.dir, 
                            pkg.dir = '~/git/covid19_brazil_hfr/') {
  
  city <- toupper(city)
  
  message("Reading data")
  dv <- fread(file.path(data.dir, sprintf("dados_%s.csv", tolower(id))),
              sep = ";",
              select = c("paciente_id", "paciente_datanascimento", "paciente_enumsexobiologico",
                         "paciente_endereco_nmmunicipio", 
                         "estabelecimento_municipio_nome", "vacina_dataaplicacao",
                         "vacina_descricao_dose", "vacina_nome"))
  setnames(dv, 
           c('paciente_id','paciente_datanascimento','paciente_enumsexobiologico','paciente_endereco_nmmunicipio','estabelecimento_municipio_nome','vacina_dataaplicacao','vacina_descricao_dose', 'vacina_nome'),
           c('patient','date_of_birth','sex','loc_res_label','loc_label','Date','dose','vaccine')
  )
  dv <- subset(dv, loc_label== city | loc_res_label == city)
  dv[,vaccine := plyr::revalue(vaccine, replace = c("Covid-19-AstraZeneca" = "Covishield",
                                                    "Covid-19-Coronavac-Sinovac/Butantan" = "Sinovac",
                                                    "Vacina covid-19 - Ad26.COV2.S - Janssen-Cilag" = "Janssen",
                                                    "Vacina covid-19 - BNT162b2 - BioNTech/Fosun Pharma/Pfizer" = "Pfizer",
                                                    "Vacina Covid-19 - Covishield" = "Covishield")
  )]
  dv[, state:= toupper(id)]
  dv[, Date2:= as.Date(gsub('([0-9]{4}-[0-9]{2}-[0-9]{2}).*','\\1',Date), format='%Y-%m-%d')]
  stopifnot(dv[, all(!is.na(Date2))])
  dv[, Date := NULL]
  tmp2 <- unique(subset(dv, select=patient))
  tmp2[, patient2 := paste0(id, seq_len(nrow(tmp2)))]
  dv <- merge(dv, tmp2, by='patient')
  dv[,patient:=NULL]
  dv[, dose2 := extract_dose(dose)] 
  set(dv, dv[, which(dose2=='')], 'dose2', NA_character_)
  set(dv, NULL, 'dose2', dv[, as.integer(dose2)]) # shouldnt be necessary 
  dv[, dose:=NULL]
  dv <- subset(dv, !is.na(dose2))
  dv <- unique(dv)
  setkey(dv, patient2, Date2,dose2)
  
  message("Check: >2 doses")
  # CHECK: for people with more than 2 vax, decide a single first dose and a single second dose
  # (there is still an issue with ppl changing year, but solved later)
  tmp2 <- dv[, list(nobs = length(Date2)), by='patient2']
  tmp2 <- subset(tmp2, nobs>2)$patient2
  tmp2 <- dv[patient2 %in% tmp2]
  setkey(tmp2, patient2, Date2,dose2)
  tmp2 <- tmp2[,list(sorted = !is.unsorted(dose2)), by = 'patient2']
  tmp2 <- tmp2[ sorted == FALSE, patient2,]
  tmp2 <- dv[patient2 %in% tmp2]
  # 'Dirty trick': set as dose1 the first first dose, and dose 2 the last second dose
  tmp2 <- tmp2[dose2 == 1,  `:=` (vaccine = vaccine[which.min(Date2)], Date2 = min(Date2)), by = "patient2"]
  tmp2 <- tmp2[dose2 == 2,  `:=` (vaccine = vaccine[which.max(Date2)], Date2 = max(Date2)), by = "patient2"]
  tmp2 <- unique(tmp2)
  setkey(tmp2, patient2, Date2,dose2)
  dv <- dv[!patient2 %in% unique(tmp2[,patient2]),,]
  dv <- rbind(dv, tmp2)
  setkey(dv, patient2, Date2,dose2)
  
  message("Check: Mixed doses")
  # some people supposedly had different vaccines: if they got second dose make first dose coherent.
  # If they got 2 second doses, consider first one in chronological order .
  tmp2 <- dv[ , list( nvax = length(unique(vaccine))) ,by=c('patient2')]
  tmp2 <- tmp2[nvax > 1, patient2]
  tmp2 <- dv[patient2 %in% tmp2 & dose2 == 2,
             {
               z <- which.min(Date2)
               list(vax2 = vaccine[z])
             }, by = 'patient2']
  dv <- merge(dv, tmp2, by = 'patient2', all.x = TRUE)
  dv[!is.na(vax2), vaccine := vax2]
  set(dv, NULL, 'vax2', NULL)
  dv <- unique(dv)
  setkey(dv, patient2, Date2,dose2)
  
  message("Check: Multiple same doses")
  # solve observations where same person has multiple administrations of the same dose
  dv <- dv[, 
           {
             z <- which.min(Date2)[1]
             list(Date2=Date2[z], loc_res_label=loc_res_label[z], loc_label=loc_label[z], vaccine=vaccine[z])
           }, 
           by=c('state','patient2','date_of_birth','sex','dose2')]
  tmp2 <- dv[, list(nobs = length(Date2)), by='patient2']
  tmp2 <- subset(tmp2, nobs > 2)
  if (nrow(tmp2) > 0) {
    tmp2 <- subset(dv, patient2 %in% tmp2$patient2)
    #	manual check --> all have 1 different date_of_birth or different sex
    tmp3 <- tmp2[, list(nage = length(unique(date_of_birth)), nsex = length(unique(sex))), by = 'patient2']
    tmp3[, check:= (nage + nsex) > 2]
    stopifnot(tmp3[,all(check)]) # substituted the manual check, even if there may be problems with nobs = 4 etc.
  
    dob_select <- function(x) {
      dt <- as.Date(names(x))
      dt <- dt[which.max(x)]
      dt
    }
    tmp3 <- tmp2[,
                 {
                   z <- table(date_of_birth)
                   z2 <- table(sex)
                   list(date_of_birth = dob_select(z),
                        sex= unlist(dimnames(z2))[which.max(z2)] )
                 },by='patient2']
    tmp2 <- merge(tmp2, tmp3, by=c('patient2','date_of_birth','sex'))
    rm(tmp3)
    dv <- rbind(subset(dv, !patient2 %in% tmp2$patient2), tmp2)
    tmp2 <- dv[, list(nobs = length(Date2)), by='patient2']
    stopifnot(nrow(subset(tmp2, nobs>2))==0)
  }

  message("Check: Wrong dose order")
  # Now we have everyone with at most 2 doses. If doses aren't in correct order, swap them
  setkey(dv, patient2, Date2)
  tmp2 <- dv[ ,
              { z1 <- which(dose2 == 1); z2 <- which(dose2 == 2)
              list( correct_order = fifelse(z1 < z2, 1  , 0 ))
              }
              , by = 'patient2']
  tmp2 <- tmp2[correct_order == 0, patient2]
  tmp2 <- dv[patient2 %in% tmp2,,]
  if(nrow(tmp2) >= 2){
    tmp2 <- tmp2[,dose2 := 1:2, by = 'patient2']
    dv <- rbind(subset(dv, !patient2%in%tmp2$patient2), tmp2)
  }
  
  message("Check: Two doses one day")
  # remain issue of two different doses the same day.
  tmp2 <- dv[, list(ndoses = length(unique(dose2))) ,by = c('patient2', 'Date2')]
  tmp2 <- tmp2[ndoses > 1, patient2,]
  tmp2 <- dv[patient2 %in% tmp2]
  tmp2 <- tmp2[dose2 == 1,,]
  dv <- rbind(subset(dv, !patient2%in%tmp2$patient2), tmp2)
  # if only second dose, make it first:
  tmp2 <- dv[, list(ndoses = length(unique(dose2))) ,by = c('patient2')]
  tmp2 <- tmp2[ndoses == 1, patient2,]
  tmp2 <- dv[patient2 %in% tmp2]
  tmp2[dose2 == 2, dose2:=1]
  dv <- rbind(subset(dv, !patient2%in%tmp2$patient2), tmp2)
  max_dv_p2 <- max(dv$patient2)

  message("Check: One dose is first dose")
  tmp2 <- dv[, .N ,by = 'patient2'][N == 1, patient2]
  stopifnot(all(dv[patient2 %in% tmp2, dose2] == 1))
  rm(tmp2)
  
  dv[, age := round(lubridate::time_length(difftime(as.Date("2021-01-01"),
                                                  date_of_birth), 
                                         "years"))]
  dv[, date_of_birth := NULL]
  
  message("Saving data")
  ## save and remove - load when needed later
  saveRDS(dv, 
          file = file.path(pkg.dir,'inst/data', sprintf('vacc_%s.rds', toupper(id))))
}

# ids_cities <- list(am = "manaus", ap = "macapa", ba = "salvador", go = "goiania", 
#                    ma = "sao luis", mg = "belo horizonte", pb = "joao pessoa", 
#                    pr = "curitiba", rj = "rio de janeiro", rn = "natal", 
#                    ro = "porto velho", rs = "porto alegre", sc = "florianopolis",
#                    sp = "sao paulo")

# 1: 
# 1: 
# ids_cities <- list(ac = "rio branco",am = "manaus", ap = "macapa", al = "maceio")
# ids_cities <- list(ba = "salvador",ce = "fortaleza", df = "brasilia", es = "vitoria", go = "goiania") # maybe keep vitoria last
# ids_cities <- list(ma = "sao luis", mt = "cuiaba", ms = "campo grande", mg = "belo horizonte", pa = "belem")
# ids_cities <- list(pb = "joao pessoa", pr = "curitiba", pe = "recife", pi = "teresina", mt = "cuiaba")
# ids_cities <- list(rj = "rio de janeiro", rn = "natal", rs = "porto alegre", ro = "porto velho", rr = "boa vista")
# ids_cities <- list(sc = "florianopolis", sp = "sao paulo", se  = "aracaju", to = "palmas")
  
Map(
  clean_vacc_data,
  names(ids_cities),
  ids_cities,
  MoreArgs = list(
    data.dir = "/home/andrea/Downloads",
    pkg.dir = "~/git/covid19_brazil_hfr"
  )
)

# 14 states only: 
tmp <- c("AP", "AM", 'BA', "GO", 'MA', "MG", 'PB',
         'PR', 'RJ', 'RN', 'RS', 'RO', 'SC', 'SP')
tmp <- paste0("vacc_", tmp, '.rds')
# tmp <- grep("vacc_", list.files("~/git/covid19_brazil_hfr/inst/data"), value = TRUE)
dv <- rbindlist(lapply(tmp, FUN = function(x) readRDS(file.path("~/git/covid19_brazil_hfr/inst/data", x))))

# CHANGE THIS WITH DATE = DATE of DOWNLOAD FROM SITE
outfile <- 'saudegovbr_vaccinations_210805.rds'

# make pretty and save
setnames(dv, colnames(dv), gsub('2','',colnames(dv)))
set(dv, NULL, 'loc_label', dv[, tolower(loc_label)])
set(dv, NULL, 'loc_res_label', dv[, tolower(loc_res_label)])
set(dv, dv[, which(loc_res_label=='')], 'loc_res_label', NA_character_)
saveRDS(dv, file=file.path('~/git/covid19_brazil_hfr/inst/data',outfile))
