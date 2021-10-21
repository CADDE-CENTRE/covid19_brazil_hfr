library(openxlsx)
library(lubridate)
library(ggplot2)
library(tidyr)
library(repr)
library(ggpubr)
library(readr)
library(tidyverse)
options(repr.plot.width=12, repr.plot.height=12)

####################################################################################
### Download new SIVEP dataset
### Check https://opendatasus.saude.gov.br/dataset/bd-srag-2021 to see the date of the latest release

# Previous was 30
RELEASE_DATE_i = "20-09-2021" 
print(RELEASE_DATE_i)

if(1) 
{
  sivepdir <- "~/Documents/P1Brazil/SIVEP/"
  pkg.dir <- "~/git/covid19_brazil_hfr/inst/"
}



options(timeout = 60*60)

# 2020 releases
url = paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2020/INFLUD-", RELEASE_DATE_i)
url = paste0(url, ".csv")
dest = paste0(sivepdir, "INFLUD-", RELEASE_DATE_i)
dest = paste0(dest, ".csv")
download.file(url, dest)

# 2021 releases
print(RELEASE_DATE_i)
url = paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2021/INFLUD21-", RELEASE_DATE_i)
url = paste0(url, ".csv", sep='')
dest = paste0(sivepdir, "INFLUD21-", RELEASE_DATE_i)
dest = paste0(dest, ".csv")
download.file(url, dest)


####################################################################################
if(0)
{
#### WATCH OUT!!!! If you want to add cities, remember to remove special Portuguese characters
CITIES <- toupper(c('Macapa', 'Manaus', 'Salvador', 'Goiania', 'Sao Luis', 'Campo Grande', 
            'Belo Horizonte', 'Joao Pessoa', 'Curitiba', 'Rio de Janeiro', 'Natal',
            'Porto Alegre', 'Porto Velho', 'Florianopolis', 'Sao Paulo', 'Palmas'))
print(CITIES)
}
# load the raw SIVEP data
RELEASE_DATE_i = RELEASE_DATE_i  # or "26-04-2021"

SIVEP_ORIGINAL20 <- paste0(sivepdir, "INFLUD-",RELEASE_DATE_i,".csv")
SIVEP_ORIGINAL21 <- paste0(sivepdir, "/INFLUD21-",RELEASE_DATE_i,".csv")

df_SIVEP_ORIGINAL20 <- read_csv2(SIVEP_ORIGINAL20)
df_SIVEP_ORIGINAL21 <- read_csv2(SIVEP_ORIGINAL21)

if(0)
{
# this is to make sure none of the CITIES were omitted
cities_sivep <- unique(df_SIVEP_ORIGINAL21$ID_MUNICIP)
for(i in 1:length(CITIES))
{
  city = CITIES[i]
  # print(city)
  print(city %in% cities_sivep)
}

# filter out CITIES
df_SIVEP_ORIGINAL20 = df_SIVEP_ORIGINAL20[which(df_SIVEP_ORIGINAL20$ID_MUNICIP %in% CITIES),] 
df_SIVEP_ORIGINAL21 = df_SIVEP_ORIGINAL21[which(df_SIVEP_ORIGINAL21$ID_MUNICIP %in% CITIES),]
}

# filter out only columns of interest
cols_of_interest20 = c("ID_MUNICIP", "ID_MN_RESI", "EVOLUCAO", "DT_SIN_PRI", "DT_EVOLUCA", "HOSPITAL",
                     "DT_INTERNA", "NU_IDADE_N", "CS_SEXO", "DT_NASC", 
                     "CS_GESTANT", "DT_ENTUTI", "DT_SAIDUTI", "SUPORT_VEN",
                     "CLASSI_FIN",
                     "OBESIDADE", "RENAL", "PNEUMOPATI", "DIABETES", "HEPATICA", 
                     "CARDIOPATI", "HEMATOLOGI", "FATOR_RISC", "CO_UNI_NOT")

cols_of_interest21 = c("ID_MUNICIP", "ID_MN_RESI", "EVOLUCAO", "DT_SIN_PRI", "DT_EVOLUCA", "HOSPITAL",
                      "DT_INTERNA", "NU_IDADE_N", "CS_SEXO", "DT_NASC", 
                      "CS_GESTANT", "DT_ENTUTI", "DT_SAIDUTI", "SUPORT_VEN",
                      "CLASSI_FIN",
                      "OBESIDADE", "RENAL", "PNEUMOPATI", "DIABETES", "HEPATICA", 
                      "CARDIOPATI", "HEMATOLOGI", "FATOR_RISC", "CO_UNI_NOT",
                      "VACINA_COV", "DOSE_1_COV", "DOSE_2_COV", "LAB_PR_COV")

df_SIVEP_ORIGINAL20 = df_SIVEP_ORIGINAL20[,cols_of_interest20]
df_SIVEP_ORIGINAL21 = df_SIVEP_ORIGINAL21[,cols_of_interest21]

# merge 2020 and 2021 data
df_SIVEP_ORIGINAL20$VACINA_COV <- NA
df_SIVEP_ORIGINAL20$DOSE_1_COV <- NA
df_SIVEP_ORIGINAL20$DOSE_2_COV <- NA
df_SIVEP_ORIGINAL20$LAB_PR_COV <- NA
df_SIVEP_original <- rbind(df_SIVEP_ORIGINAL20,df_SIVEP_ORIGINAL21)

rm(df_SIVEP_ORIGINAL20)
rm(df_SIVEP_ORIGINAL21)


# # get only suspected and confirmed
# df_SIVEP_original =  df_SIVEP_original[which(df_SIVEP_original$CLASSI_FIN %in% c(4,5)),] #4 is suspected, 5 is confirmed

# unpack hospital classification from Carlos updated for SP 
cnes <- readRDS(paste0(pkg.dir,"/data/hospital_classification_updated_210709.rds"))
cnes <- subset(cnes, select = c("CO_UNI_NOT", "CO_CONVENIO"))

# # filter out only columns of interest
# df_SIVEP_original = df_SIVEP_original[,c("ID_MUNICIP", "ID_MN_RESI", "EVOLUCAO", "DT_SIN_PRI", "DT_EVOLUCA",
#                                          "DT_INTERNA", "NU_IDADE_N", "CS_SEXO", "DT_NASC", 
#                                          "CS_GESTANT", "DT_ENTUTI", "DT_SAIDUTI", "SUPORT_VEN",
#                                          "CLASSI_FIN",
#                                          "OBESIDADE", "RENAL", "PNEUMOPATI", "DIABETES", "HEPATICA", 
#                                          "CARDIOPATI", "HEMATOLOGI", "FATOR_RISC", "CO_UNI_NOT")]

# merge with cnes database
df_SIVEP_original <- df_SIVEP_original %>% left_join(cnes, by = "CO_UNI_NOT") 

df_sivep_AM = df_SIVEP_original
if(0)
{
# filter out to have selected cities only
df_sivep_AM = df_SIVEP_original[which(df_SIVEP_original$ID_MUNICIP %in% CITIES),]
# print out the cities
print(unique(df_sivep_AM$ID_MUNICIP))
}

# rm(df_SIVEP_original)

# transform to dates
df_sivep_AM$DT_EVOLUCA = dmy(df_sivep_AM$DT_EVOLUCA)  # date of outcome: death or cure
df_sivep_AM$DT_SIN_PRI = dmy(df_sivep_AM$DT_SIN_PRI)  # date of symptoms
df_sivep_AM$DT_INTERNA = dmy(df_sivep_AM$DT_INTERNA )  # date of admission
df_sivep_AM$DT_NASC = dmy(df_sivep_AM$DT_NASC )  # date of birth
df_sivep_AM$DT_ENTUTI = dmy(df_sivep_AM$DT_ENTUTI )  # date of entering ICU
df_sivep_AM$DT_SAIDUTI = dmy(df_sivep_AM$DT_SAIDUTI )  # date of leaving ICU
df_sivep_AM$DOSE_1_COV = dmy(df_sivep_AM$DOSE_1_COV )  # date of 1st dose of COVID vaccine
df_sivep_AM$DOSE_2_COV = dmy(df_sivep_AM$DOSE_2_COV )  # date of 2nd dose of COVID vaccine

df_sivep_AM_filtered <- df_sivep_AM

# rmeove all dates that are too early
# df_sivep_AM_filtered = subset(df_sivep_AM, (DT_INTERNA >= dmy("26-02-2020") | (is.na(DT_INTERNA))))
# df_sivep_AM_filtered = subset(df_sivep_AM, (DT_ENTUTI >= dmy("26-02-2020") | (is.na(DT_ENTUTI))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_INTERNA >= dmy("01-12-2019") | (is.na(DT_INTERNA))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_EVOLUCA >= dmy("01-12-2019") | (is.na(DT_EVOLUCA))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_SIN_PRI >= dmy("01-12-2019") | (is.na(DT_SIN_PRI))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_INTERNA >= dmy("01-12-2019") | (is.na(DT_INTERNA))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_ENTUTI >= dmy("01-12-2019") | (is.na(DT_ENTUTI))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_SAIDUTI >= dmy("01-12-2019") | (is.na(DT_SAIDUTI))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DOSE_1_COV >= dmy("01-12-2019") | (is.na(DOSE_1_COV))))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DOSE_2_COV >= dmy("01-12-2019") | (is.na(DOSE_2_COV))))
# let's do it later
# df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_INTERNA <= dmy(RELEASE_DATE_i) | (is.na(DT_INTERNA))))
# df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_EVOLUCA <= dmy(RELEASE_DATE_i) | (is.na(DT_EVOLUCA))))
# df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_SIN_PRI <= dmy(RELEASE_DATE_i) | (is.na(DT_SIN_PRI))))
# df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_NASC <= dmy(RELEASE_DATE_i) | (is.na(DT_NASC))))
# df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_ENTUTI <= dmy(RELEASE_DATE_i) | (is.na(DT_ENTUTI))))
# df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DT_SAIDUTI <= dmy(RELEASE_DATE_i) | (is.na(DT_SAIDUTI))))


# add age where DOB given
df_sivep_AM_filtered$AgeCalculated = floor(as.numeric(df_sivep_AM_filtered$DT_SIN_PRI - df_sivep_AM_filtered$DT_NASC)/365)
df_sivep_AM_filtered$AgeCalculated <- ifelse(is.na(df_sivep_AM_filtered$AgeCalculated),
                                             df_sivep_AM_filtered$NU_IDADE_N, df_sivep_AM_filtered$AgeCalculated)
df_sivep_AM_filtered$NU_IDADE_N = df_sivep_AM_filtered$AgeCalculated

# get pregnant yes/no
df_sivep_AM_filtered$Pregnant <- ifelse(df_sivep_AM_filtered$CS_GESTANT %in% c(1,2,3,4),
                                        TRUE, FALSE)

# get ventilated yes/no
# df_sivep_AM_filtered$Ventilator <- ifelse(df_sivep_AM_filtered$SUPORT_VEN %in% c(1,2),
#                                           TRUE, FALSE)
df_sivep_AM_filtered$Ventilator <- ifelse(df_sivep_AM_filtered$SUPORT_VEN %in% c(1,2),
                                          TRUE, ifelse(is.na(df_sivep_AM_filtered$SUPORT_VEN), NA, FALSE))

# get SRAG covid/other
df_sivep_AM_filtered$SRAG <- ifelse(df_sivep_AM_filtered$CLASSI_FIN %in% c(4,5),
                                    'COVID', 'Other')

# get confirmed yes/no
df_sivep_AM_filtered$Confirmed <- ifelse(df_sivep_AM_filtered$CLASSI_FIN %in% c(5),
                                         TRUE, FALSE)

# comorbidities
# 1: True, 2: False, 9: NA
df_sivep_AM_filtered$Comorb_Obesity <- ifelse(df_sivep_AM_filtered$OBESIDADE %in% c(1),
                                              TRUE, ifelse(df_sivep_AM_filtered$OBESIDADE %in% c(9, NA), NA, FALSE))

df_sivep_AM_filtered$Comorb_Renal <- ifelse(df_sivep_AM_filtered$RENAL %in% c(1),
                                            TRUE, ifelse(df_sivep_AM_filtered$RENAL %in% c(9, NA), NA, FALSE))

df_sivep_AM_filtered$Comorb_Pneumon <- ifelse(df_sivep_AM_filtered$PNEUMOPATI %in% c(1),
                                              TRUE, ifelse(df_sivep_AM_filtered$PNEUMOPATI %in% c(9, NA), NA, FALSE))

df_sivep_AM_filtered$Comorb_Diabetes <- ifelse(df_sivep_AM_filtered$DIABETES %in% c(1),
                                               TRUE, ifelse(df_sivep_AM_filtered$DIABETES %in% c(9, NA), NA, FALSE))

df_sivep_AM_filtered$Comorb_Hepat <- ifelse(df_sivep_AM_filtered$HEPATICA %in% c(1),
                                            TRUE, ifelse(df_sivep_AM_filtered$HEPATICA %in% c(9, NA), NA, FALSE))

df_sivep_AM_filtered$Comorb_Cardio <- ifelse(df_sivep_AM_filtered$CARDIOPATI %in% c(1),
                                             TRUE, ifelse(df_sivep_AM_filtered$CARDIOPATI %in% c(9, NA), NA, FALSE))

df_sivep_AM_filtered$Comorb_Hemat <- ifelse(df_sivep_AM_filtered$HEMATOLOGI %in% c(1),
                                            TRUE, ifelse(df_sivep_AM_filtered$HEMATOLOGI %in% c(9, NA), NA, FALSE))

df_sivep_AM_filtered$Comorb_RiskFactor <- ifelse(df_sivep_AM_filtered$FATOR_RISC %in% c('S'),
                                                 TRUE, ifelse(df_sivep_AM_filtered$FATOR_RISC %in% c(9, NA), NA, FALSE))

# outcome
# change 9 to NA
df_sivep_AM_filtered$EVOLUCAO <- ifelse(df_sivep_AM_filtered$EVOLUCAO %in% c(9),
                                        NA, df_sivep_AM_filtered$EVOLUCAO)

# remove the DT_EVOLUCA date for those who were did not die (DT_EVOLUCA is date of outcome, EVOLUCAO=2 is death)
# df_sivep_AM_filtered = within(df_sivep_AM_filtered, DT_EVOLUCA[EVOLUCAO != 2] <- NA)
df_sivep_AM_filtered = within(df_sivep_AM_filtered, DT_EVOLUCA[EVOLUCAO %in% c(1,NA)] <- NA)

# change the outcome classes
df_sivep_AM_filtered = within(df_sivep_AM_filtered, EVOLUCAO[EVOLUCAO %in% c(1)] <- 'Cure')
df_sivep_AM_filtered = within(df_sivep_AM_filtered, EVOLUCAO[EVOLUCAO %in% c(2)] <- 'Death_Covid')
df_sivep_AM_filtered = within(df_sivep_AM_filtered, EVOLUCAO[EVOLUCAO %in% c(3)] <- 'Death_Other')

# change the 'Death_Other' for classes 4 and 5 to 'Death_Covid
df_sivep_AM_filtered$EVOLUCAO[(df_sivep_AM_filtered$EVOLUCAO == 'Death_Other') & (df_sivep_AM_filtered$CLASSI_FIN %in% c(4,5))] <- 'Death_Covid'

table(df_sivep_AM_filtered$EVOLUCAO, useNA="always")

# vaccine
# 1: True, 2: False, 9: NA
df_sivep_AM_filtered$VACINA_COV <- ifelse(df_sivep_AM_filtered$VACINA_COV %in% c(1),
                                                 TRUE, ifelse(df_sivep_AM_filtered$VACINA_COV %in% c(9, NA), NA, FALSE))


# drop  the unnecessary columns
df_sivep_AM_filtered =  df_sivep_AM_filtered[ , -which(names(df_sivep_AM_filtered) %in% c("AgeCalculated",
                                                                                          "DT_NASC",
                                                                                          "CS_GESTANT", "SUPORT_VEN",
                                                                                          #"CLASSI_FIN", 
                                                                                          "OBESIDADE", "RENAL", "PNEUMOPATI", 
                                                                                          "DIABETES", "HEPATICA", 
                                                                                          "CARDIOPATI", "HEMATOLOGI",
                                                                                          "FATOR_RISC", "CO_UNI_NOT"))]
# change the fields in the HOSPITAL column
df_sivep_AM_filtered$HOSPITAL <- ifelse(df_sivep_AM_filtered$HOSPITAL %in% c(1),
                                            TRUE, ifelse(df_sivep_AM_filtered$HOSPITAL %in% c(9, NA), NA, FALSE))

# rename the columns
df_sivep_AM_filtered = df_sivep_AM_filtered %>% 
  rename(
    Municip = ID_MUNICIP,
    Municip_resid = ID_MN_RESI,
    Admission = HOSPITAL,
    Outcome = EVOLUCAO,
    DateSymptoms = DT_SIN_PRI,
    DateAdmission = DT_INTERNA,
    DateDeath = DT_EVOLUCA,
    DateICUstart = DT_ENTUTI,
    DateICUend = DT_SAIDUTI,
    Age = NU_IDADE_N,
    Sex = CS_SEXO,
    Hospital_type = CO_CONVENIO,
    Vaccinated = VACINA_COV,
    Date1stDose = DOSE_1_COV,
    Date2ndDose = DOSE_2_COV,
    Vaccine_type = LAB_PR_COV
  )

# add admission to death column
df_sivep_AM_filtered$AdminToDeath = df_sivep_AM_filtered$DateDeath - df_sivep_AM_filtered$DateAdmission
df_sivep_AM_filtered$OnsetToDeath = df_sivep_AM_filtered$DateDeath - df_sivep_AM_filtered$DateSymptoms
df_sivep_AM_filtered$OnsetToAdmin = df_sivep_AM_filtered$DateAdmission - df_sivep_AM_filtered$DateSymptoms
df_sivep_AM_filtered$OnsetToICU = df_sivep_AM_filtered$DateICUstart - df_sivep_AM_filtered$DateSymptoms

##### remove couple more obvious errors
# df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (OnsetToAdmin >= 0) | (is.na(OnsetToAdmin)))

# remove were admission missing but ICU admission given -- we don't know what happened
df_sivep_AM_filtered <- df_sivep_AM_filtered[!((is.na(df_sivep_AM_filtered$DateAdmission)) & !(is.na(df_sivep_AM_filtered$DateICUstart))),]

# reported dates must be after Dec 2019  or NAN
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateICUstart >= dmy("01-12-2019")) | (is.na(DateICUstart)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateICUend >= dmy("01-12-2019")) | (is.na(DateICUend)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateAdmission >= dmy("01-12-2019")) | (is.na(DateAdmission)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateDeath >= dmy("01-12-2019")) | (is.na(DateDeath)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateSymptoms >= dmy("01-12-2019")) | (is.na(DateSymptoms)))



############
# fix the incorrect Admission time (day mixed with month) where negative admission to death time
swap_day_and_month <- function(date_entry)
{
  year <- year(date_entry)
  month <- month(date_entry)
  day <- day(date_entry)
  
  # if (day < 13)
  # {
  #   correct_date <- make_date(year=year, month=day, day = month)
  #   print('swapping')
  # }
  # else
  # {
  #   correct_date <- date_entry
  # }

  correct_date <- make_date(year=year, month=day, day = month)
  print('swapping')

  
  return(correct_date)
}

cache <- df_sivep_AM_filtered  # saving just in case
df_sivep_AM_filtered <- cache

# get rows where h2d is negative
df_neg <- filter(df_sivep_AM_filtered, AdminToDeath < 0)
# remove them from the main dataframe
df_sivep_AM_filtered <- filter(df_sivep_AM_filtered, is.na(df_sivep_AM_filtered$AdminToDeath) | (AdminToDeath >= 0))
# apply the function to fix the errors in df_neg
df_neg['DateAdmission'] <- lapply(df_neg['DateAdmission'], swap_day_and_month)
df_neg$AdminToDeath = df_neg$DateDeath - df_neg$DateAdmission
df_neg$OnsetToDeath = df_neg$DateDeath - df_neg$DateSymptoms
df_neg$OnsetToAdmin = df_neg$DateAdmission - df_neg$DateSymptoms
df_neg$OnsetToICU = df_neg$DateICUstart - df_neg$DateSymptoms
# bind the corrected df to the main one
df_sivep_AM_filtered <- rbind(df_sivep_AM_filtered, df_neg)

#############
# reported dates must be less than year 2022 date or NAN
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateICUstart < dmy("11-12-2021")) | (is.na(DateICUstart)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateICUend < dmy("11-12-2021")) | (is.na(DateICUend)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateAdmission < dmy("11-12-2021")) | (is.na(DateAdmission)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateDeath < dmy("11-12-2021")) | (is.na(DateDeath)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateSymptoms < dmy("11-12-2021")) | (is.na(DateSymptoms)))

######
# for all dates of admission above RELEASE date, try to swap month and day
df_SIVEP_swapped <- df_sivep_AM_filtered
df_admission_wrong <- filter(df_SIVEP_swapped, DateAdmission > dmy(RELEASE_DATE_i))
df_SIVEP_swapped <- filter(df_SIVEP_swapped, is.na(df_SIVEP_swapped$DateAdmission) | (DateAdmission <= dmy(RELEASE_DATE_i)))
df_admission_wrong['DateAdmission'] <- lapply(df_admission_wrong['DateAdmission'], swap_day_and_month)
df_admission_wrong$AdminToDeath = df_admission_wrong$DateDeath - df_admission_wrong$DateAdmission
df_admission_wrong$OnsetToAdmin = df_admission_wrong$DateAdmission - df_admission_wrong$DateSymptoms
df_SIVEP_swapped <- rbind(df_SIVEP_swapped, df_admission_wrong)
df_sivep_AM_filtered <- df_SIVEP_swapped


# reported dates must be less than Release date or NAN
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateICUstart <= dmy(RELEASE_DATE_i)) | (is.na(DateICUstart)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateICUend <= dmy(RELEASE_DATE_i)) | (is.na(DateICUend)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateAdmission <= dmy(RELEASE_DATE_i)) | (is.na(DateAdmission)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateDeath <= dmy(RELEASE_DATE_i)) | (is.na(DateDeath)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (DateSymptoms <= dmy(RELEASE_DATE_i)) | (is.na(DateSymptoms)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (Date1stDose <= dmy(RELEASE_DATE_i)) | (is.na(Date1stDose)))
df_sivep_AM_filtered = subset(df_sivep_AM_filtered, (Date2ndDose <= dmy(RELEASE_DATE_i)) | (is.na(Date2ndDose)))

change_vaccine_type <- function(vaccines){
  
  # AstraZeneca/fiocruz/covishield
  vaccines <- gsub(".*ASTRA.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*AZ.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  
  vaccines <- gsub(".*ASTRAZENE.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*ENECA.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*ZENICA.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*ASTE.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*AZASGTRENICA.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*ZANECA.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)

  vaccines <- gsub(".*OXFOR.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*OXO.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*OXF.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*OFX.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  
  vaccines <- gsub(".*FORD.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*CRUZ*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*FIOCR*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*FIO*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T) 
  
  vaccines <- gsub(".*SHIEL.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*SHILD.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T) 
  
  # coronavac/sinovac/Butant?
  vaccines <- gsub(".*CORONAVAC.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*CORONAV.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*CORONA VAC.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*CORONOV.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*COROV.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*CORNA.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*CORONACAV.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*CORONAC.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  
  vaccines <- gsub(".*CORONV.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*AVAC.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*OVAC.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*NIVAC.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*COVONAVA.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  
  vaccines <- gsub(".*SINO.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*SIVO.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  
  vaccines <- gsub(".*BUTAN.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*BUTAT.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*BUNTA.*", "SINOVAC", vaccines, fixed = FALSE, ignore.case=T) 
  
  # Janssen 
  vaccines <- gsub(".*JANSS.*", "JANSSEN", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*JASS.*", "JANSSEN", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*ONSON.*", "JANSSEN", vaccines, fixed = FALSE, ignore.case=T)
  vaccines <- gsub(".*JANS.*", "JANSSEN", vaccines, fixed = FALSE, ignore.case=T)
  
  # Sputnik/gamaleya
  vaccines <- gsub(".*SPUT.*", "SPUTNIK", vaccines, fixed = FALSE, ignore.case=T) 
  vaccines <- gsub(".*GAMALEYA.*", "SPUTNIK", vaccines, fixed = FALSE, ignore.case=T) 
  
  # pfizer
  vaccines <- gsub(".*PFI.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*IZER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*YZER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*ZIER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*AZER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*IFER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*AZER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*PZIR.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*EZER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  vaccines <- gsub(".*ZER.*", "PFIZER", vaccines, fixed = FALSE, ignore.case=T) # fixed = True is when no regex used
  

  # did not take a vaccine
  vaccines <- gsub(".*NAO TOMOU.*", "NO VACCINE", vaccines, fixed = FALSE, ignore.case=T)
  
  return(vaccines)
}

change_astra_type <- function(vaccines){
  vaccines <- gsub(".*ASTRA.*", "ASTRAZENECA", vaccines, fixed = FALSE, ignore.case=T)
  return(vaccines)
}

# df_test <- df_sivep_AM_filtered
df_sivep_AM_filtered$Vaccine_type <- change_vaccine_type(df_sivep_AM_filtered$Vaccine_type)
df_sivep_AM_filtered$Vaccine_type <- change_astra_type(df_sivep_AM_filtered$Vaccine_type)

df_sivep_AM_filtered$Vaccine_type <- ifelse(df_sivep_AM_filtered$Vaccine_type %in% c('PFIZER', 'ASTRAZENECA', 'SINOVAC', 'JANSSEN', 'SPUTNIK'),
                                            df_sivep_AM_filtered$Vaccine_type, NA)

table(df_sivep_AM_filtered$Vaccine_type, useNA="always")

# df_test <- df_sivep_AM_filtered
# df_test$Vaccine_type <- change_vaccine_type(df_test$Vaccine_type) 
# table(df_test$Vaccine_type, useNA="always")
# 
# df_test$Vaccine_type <- change_astra_type(df_test$Vaccine_type) 
# 
# df_test$Vaccine_type <- ifelse(df_test$Vaccine_type %in% c('PFIZER', 'ASTRAZENECA', 'SINOVAC', 'JANSSEN', 'SPUTNIK'),
#                                df_test$Vaccine_type, NA)
# 
# table(df_test$Vaccine_type, useNA="always")


######################
# drop  the unnecessary columns because file needs to be less than 100 MB
df_sivep_AM_filtered =  df_sivep_AM_filtered[ , -which(names(df_sivep_AM_filtered) %in% c("AdminToDeath",
                                                                                          "OnsetToDeath",
                                                                                          "OnsetToAdmin", "OnsetToICU"))]

######################

# save to csv and rds
if(0){
  write_csv(path = paste0(pkg.dir, "/data/SIVEP_hospital_",
                          tail(RELEASE_DATE_i,1), ".csv"), df_sivep_AM_filtered)
}


saveRDS(df_sivep_AM_filtered, file = paste0(pkg.dir,"/data/SIVEP_hospital_",
                                            tail(RELEASE_DATE_i,1), "-all.RDS"))

