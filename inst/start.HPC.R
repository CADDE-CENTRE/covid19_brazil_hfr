#pkg.dir <- system.file(package = "covid19.P1.IFR" )

#	function to make PBS header
make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal\nsource activate covid19AgeModel", hpc.q=NA, hpc.array=1 )
{	
	pbshead <- "#!/bin/sh"
	tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
	pbshead <- paste(pbshead, tmp, sep = "\n")
	tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")	
	pbshead <- paste(pbshead, tmp, sep = "\n")
	pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
	if(hpc.array>1)
	{
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
	}				
	if(!is.na(hpc.q))
	{
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	}		
	pbshead	<- paste(pbshead, hpc.load, sep = "\n")
	pbshead
}


## olli s args for HFR model
if(0)
{			
	loc_labels <- c('belo-horizonte','curitiba', "florianopolis", "goiania", "joao-pessoa", "macapa", 'manaus', "natal", 'porto-alegre', "porto-velho", 'rio-de-janeiro', 'salvador', 'sao-paulo', "sao-luis")
	
	args <- list()	
	args$stanmodel = 'age_hfr_210719d.stan'
	args$source_dir= '/rds/general/user/or105/home/libs/covid19.P1.IFR'
	args$preproc_file= 'inst/scripts/preprocessing.210719d2.R'
	args$cmdstan_file= 'inst/scripts/HFR.fit.210719d.cmdstan.R'
	args$postproc_oneloc_file = 'inst/scripts/HFR.postprocessing.oneloc.210719.R'
	args$postproc_allloc_file = 'inst/scripts/HFR.postprocessing.allloc.210719.R'
	args$postproc_mainfigs = 'inst/scripts/HFR.postprocessing.main.figures.R'
	args$postproc_tables = 'inst/scripts/HFR.postprocessing.tables.R'
	args$postproc_knit_rmd = 'inst/scripts/HFR.postprocessing.knit.report.R'
	args$postproc_rmd = 'inst/scripts/HFR.postprocessing.make.report.Rmd'		
	args$threshold_check = 14
	args$fr.exclude.P1 = 0
	args$use.only.hosp.deaths = 1
	args$genomic.data.type <- 'GISAID_darlan'
	args$fr.predictors.transformation <- 'value_cs'	
	args$keep.hospital.type <- 'Private_Hospitals,Public_Hospitals,Unknown'
	
	# all except excess deaths
	#args$fr.predictors <- "r.deaths.ooh.04,r.nicu.adm.02,r.nhosp.adm.04,r.nhosp.adm.75plus.04,icu.per.ventilator.02,icu.per.beds.02,icu.per.physicians.all.04,icu.per.physicians.specialist.04,hosps.per.beds.02,hosps.per.physicians.all.02,p.res.hosps.40"
	#args$fr.predictors.names <- "4-wk_Out_of_hospital_deaths,2-wk_ICU_admissions,4-wk_SARI_admissions,4-wk_SARI_admissions_among_75+,2-wk_ICU_admissions_per_ventilator,2-wk_ICU_admissions_per_ICU_bed,4-wk_ICU_admissions_per_physician,4-wk_ICU_admissions_per_critical_care_specialist,2-wk_SARI_admissions_per_ICU_bed,2-wk_SARI_admissions_per_physician,4-wk_prop_residents_in_SARI_admissions"	
	#args$out_dir= '/rds/general/project/ratmann_covid19/live/BrazilP1/hfr-fit-210619b3-noexcess'
	
	# only demand per resource
	#args$fr.predictors <- "icu.per.ventilator.02,icu.per.beds.02,icu.per.physicians.all.04,icu.per.physicians.specialist.04,hosps.per.beds.02,hosps.per.physicians.all.02,p.res.hosps.40"
	#args$fr.predictors.names <- "2-wk_ICU_admissions_per_ventilator,2-wk_ICU_admissions_per_ICU_bed,4-wk_ICU_admissions_per_physician,4-wk_ICU_admissions_per_critical_care_specialist,2-wk_SARI_admissions_per_ICU_bed,2-wk_SARI_admissions_per_physician,4-wk_prop_residents_in_SARI_admissions"
	#args$out_dir= '/rds/general/project/ratmann_covid19/live/BrazilP1/hfr-fit-210619d2-resources'
	
	# only demand per resource version3
	args$fr.predictors <- "icu.per.ventilator.02,icu.per.beds.02,icu.per.nurse.02,icu.per.tech.nurse.02,icu.per.physiotherapist.02,icu.per.physicians.all.02,icu.per.physicians.specialist.02,hosps.per.beds.02,hosps.per.IU.beds.02,hosps.per.physicians.all.02,hosps.per.ventilator.02,hosps.per.nurse.02"
	args$fr.predictors.names <- "2-wk_ICU_admissions_per_ventilator,2-wk_ICU_admissions_per_ICU_bed,2-wk_ICU_admissions_per_nurse,2-wk_ICU_admissions_per_technical_nurse,2-wk_ICU_admissions_per_physiotherapist,2-wk_ICU_admissions_per_physician,2-wk_ICU_admissions_per_critical_care_specialist,2-wk_SARI_admissions_per_bed,2-wk_SARI_admissions_per_IU_bed,2-wk_SARI_admissions_per_physician,2-wk_SARI_admissions_per_ventilator,2-wk_SARI_admissions_per_nurse"
	args$out_dir= '/rds/general/project/ratmann_covid19/live/BrazilP1/hfr-fit-210719d2-resources3'
		
	args$out_dir <- paste0(args$out_dir,'-ex',args$fr.exclude.P1,'-tr',gsub('_','',args$fr.predictors.transformation),'-hd',args$use.only.hosp.deaths)
}

## Andrea s args for HFR model
if(1)
{	
  loc_labels <- c('belo-horizonte','curitiba', "florianopolis", "goiania", "joao-pessoa", "macapa", 'manaus', "natal", 'porto-alegre', "porto-velho", 'rio-de-janeiro', 'salvador', 'sao-paulo', "sao-luis")
  
  args <- list()	
  args$stanmodel = 'age_hfr_210719d.stan'
  args$source_dir= '/rds/general/user/ab1820/home/git/covid19.P1.IFR'
  args$preproc_file= 'inst/scripts/HFR.preprocessing.R'
  args$cmdstan_file= 'inst/scripts/HFR.fit.cmdstan.R'
  args$postproc_oneloc_file = 'inst/scripts/HFR.postprocessing.oneloc.R'
  args$postproc_allloc_file = 'inst/scripts/HFR.postprocessing.allloc.R'
  args$postproc_mainfigs = 'inst/scripts/HFR.postprocessing.main.figures.R'
  args$postproc_tables = 'inst/scripts/HFR.postprocessing.tables.R'
  args$postproc_knit_rmd = 'inst/scripts/HFR.postprocessing.knit.report.R'
  args$postproc_rmd = 'inst/scripts/HFR.postprocessing.make.report.Rmd'		
  args$threshold_check = length(loc_labels)
  args$fr.exclude.P1 = 0
  args$use.only.hosp.deaths = 1
  args$genomic.data.type <- 'GISAID_darlan'
  args$fr.predictors.transformation <- 'value_cs'	
  args$keep.hospital.type <-  'Private_Hospitals,Public_Hospitals,Unknown'
  
  
  # only demand per resource version 8
  args$fr.predictors <- "icu.per.ventilator.02,icu.per.beds.02,icu.per.nurse.02,icu.per.tech.nurse.02,icu.per.physiotherapist.02,icu.per.physicians.all.02,icu.per.intensivists.02,hosps.per.saribeds.02,hosps.per.saribedsvent.02,hosps.per.physicians.all.02,hosps.per.ventilator.02,hosps.per.nurse.02"
  args$fr.predictors.names <- gsub(' ','_',"2-wk ICU admissions per ventilator,2-wk ICU admissions per ICU bed,2-wk ICU admissions per nurse,2-wk ICU admissions per nurse assistant,2-wk ICU admissions per physiotherapist,2-wk ICU admissions per physician,2-wk ICU admissions per intensivist,2-wk SARI admissions per critical bed, 2-wk SARI admissions per critical bed with ventilator,2-wk SARI admissions per physician,2-wk SARI admissions per ventilator,2-wk SARI admissions per nurse")
  # tmp <- gsub('^(.*?)preprocessing.(.*?).R$','\\2',args$preproc_file)
  args$out_dir= paste0('/rds/general/project/ratmann_covid19/live/BrazilP1/hfr-fit-', '210719h4','-resources8')  
  
  # first one is CLASS 4,5 + UNKNOWN, second only uses CLASS 5 (uncomment 1 exclusively)
  args$out_dir <- paste0(args$out_dir,'-ex',args$fr.exclude.P1,'-tr',gsub('_','',args$fr.predictors.transformation),'-hd',args$use.only.hosp.deaths)
  # args$out_dir <- paste0(args$out_dir,'-ex',args$fr.exclude.P1,'-tr',gsub('_','',args$fr.predictors.transformation),'-hd',args$use.only.hosp.deaths, '-onlyconfirmed')

  # If do not want to estimate outcomes among ICU and unknown outcomes:
  # args$out_dir <- paste0(args$out_dir, '-nocens')
  
  # If want to set 'unknown' to discharged:
  args$out_dir <- paste0(args$out_dir, '-unkdis')

}

## Charlie's args for HFR model
if(0)
{			
  
  args <- list()	
  
  args$genomic.data.type <- 'GISAID_darlan'
  #args$genomic.data.type <- 'genomica_fiocruz'
  #args$genomic.data.type <- 'random_sample|GISAID_alfonso_lutz'
  
  if (args$genomic.data.type == 'GISAID_darlan') {
    loc_labels <- c('belo-horizonte','curitiba', "florianopolis", "goiania", "joao-pessoa", "macapa", 
                    'manaus', "natal", 'porto-alegre', "porto-velho", 'rio-de-janeiro', 'salvador', 'sao-paulo', "sao-luis")
    args$threshold_check = 14
    gen_name <- "GISAID"
  } else if (args$genomic.data.type == 'random_sample|GISAID_alfonso_lutz') {
    loc_labels <- c('belo-horizonte', 'manaus', 'sao-paulo') # 'porto-alegre', 
    args$threshold_check = 3
    gen_name <- "random"
  } else if (args$genomic.data.type == 'genomica_fiocruz') {
    loc_labels <- c('belo-horizonte','curitiba', "florianopolis", "goiania", "joao-pessoa", "macapa", 
                    'manaus', "natal", 'porto-alegre', "porto-velho", 'rio-de-janeiro', 'salvador', 'sao-paulo', "sao-luis")
    args$threshold_check = 14
    gen_name <- "fiocruz"
  } else {
    stop("incorrect genomic.data.type specified")
  }
  
  args$stanmodel = 'age_hfr_210719d.stan' 
  args$source_dir= '/rds/general/user/cw1716/home/covid19.P1.IFR'  
  args$preproc_file= 'inst/scripts/preprocessing.210719h2.R'
  args$cmdstan_file= 'inst/scripts/HFR.fit.210619d.cmdstan.R'
  args$postproc_oneloc_file = 'inst/scripts/HFR.postprocessing.oneloc.R'
  args$postproc_allloc_file = 'inst/scripts/HFR.postprocessing.allloc.R'
  args$postproc_mainfigs = 'inst/scripts/HFR.postprocessing.main.figures.R'
  args$postproc_tables = 'inst/scripts/HFR.postprocessing.tables.R'
  args$postproc_knit_rmd = 'inst/scripts/HFR.postprocessing.knit.report.R'
  args$postproc_rmd = 'inst/scripts/HFR.postprocessing.make.report.Rmd'  
  args$fr.exclude.P1 = 0
  args$use.only.hosp.deaths = 1
  args$fr.predictors.transformation <- 'value_cs'
  args$keep.hospital.type <- 'Private_Hospitals,Public_Hospitals,Unknown'

  args$threshold_check = length(loc_labels) # Not really needed, but why not?

  # only demand per resource version 8
  args$fr.predictors <- "icu.per.ventilator.02,icu.per.beds.02,icu.per.nurse.02,icu.per.tech.nurse.02,icu.per.physiotherapist.02,icu.per.physicians.all.02,icu.per.intensivists.02,hosps.per.saribeds.02,hosps.per.saribedsvent.02,hosps.per.physicians.all.02,hosps.per.ventilator.02,hosps.per.nurse.02"
  args$fr.predictors.names <- gsub(' ','_',"2-wk ICU admissions per ventilator,2-wk ICU admissions per ICU bed,2-wk ICU admissions per nurse,2-wk ICU admissions per nurse assistant,2-wk ICU admissions per physiotherapist,2-wk ICU admissions per physician,2-wk ICU admissions per intensivist,2-wk SARI admissions per critical bed, 2-wk SARI admissions per critical bed with ventilator,2-wk SARI admissions per physician,2-wk SARI admissions per ventilator,2-wk SARI admissions per nurse")
  
  # Output directory name
  tmp <- gsub('^(.*?)preprocessing.(.*?).R$','\\2',args$preproc_file)
  args$out_dir= paste0('/rds/general/project/ratmann_covid19/live/BrazilP1/hfr-fit-', tmp,'-resources8')  
  args$out_dir <- paste0(args$out_dir,'-ex',args$fr.exclude.P1,'-tr',gsub('_','',args$fr.predictors.transformation),'-hd',args$use.only.hosp.deaths, "-", gen_name)
}


stopifnot(!dir.exists(args$out_dir))
dir.create(args$out_dir)	

########################################
#	prepare commands for cmdstan run and post processing for that location 
########################################
cmds <- vector('list',length(loc_labels))
for(i in seq_along(loc_labels))
{
	cmd2 <- ''		
	cmd2 <- paste0(cmd2,"CWD=$(pwd)\n")
	cmd2 <- paste0(cmd2,"echo $CWD\n")	
	tmpdir.prefix <- paste0('cvd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),'_',i)
	tmpdir <- paste0("$CWD/",tmpdir.prefix)
	cmd2 <- paste0(cmd2,"mkdir -p ",tmpdir,'\n')
	cmd2 <- paste0(cmd2,"cp ",file.path(args$out_dir,paste0(basename(args$out_dir),'_allloc_workspace.rda'))," ",tmpdir,'\n')
	tmp <- paste0('Rscript ', 
			file.path(args$source_dir, args$cmdstan_file), 
			' -outdir ', tmpdir,'',		
			' -loc_label ', loc_labels[i], '',
			' -file_prefix ',basename(args$out_dir), '',
			' -pkgdir ', args$source_dir, '',
			' -fr.exclude.P1 ', args$fr.exclude.P1, '',
			' -use.only.hosp.deaths ', args$use.only.hosp.deaths
		)
	cmd2 <- paste0(cmd2, tmp, '\n')
	cmd2 <- paste0(cmd2, tmp, '\n')	# note: not a mistake. we re-run by default to catch unexpected stan failure
	cmd2 <- paste0(cmd2, tmp, '\n')	# note: not a mistake. we re-run by default to catch unexpected stan failure
	cmd2 <- paste0(cmd2, tmp, '\n') # note: not a mistake. we re-run by default to catch unexpected stan failure
	tmp <- paste0('Rscript ', 
	        file.path(args$source_dir, args$postproc_oneloc_file), 
	        ' -outdir ', tmpdir,'',		
	        ' -loc_label ', loc_labels[i], '',
	        ' -file_prefix ',basename(args$out_dir), '',
	        ' -pkgdir ', args$source_dir
		)
	cmd2 <- paste0(cmd2, tmp, '\n')
	
	# PostProc. allloc
	cmd2 <- paste0(cmd2, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_',gsub('-|_| ','',loc_labels[i])),'* ', args$out_dir,'\n')
	cmd2 <- paste0(cmd2,"cp ",file.path(args$out_dir,'*workspace.rda')," ",tmpdir,'\n')
	cmd2 <- paste0(cmd2,"cp ",file.path(args$out_dir,'*fit.rds')," ",tmpdir,'\n')
	tmp <- paste0('Rscript ', 
			file.path(args$source_dir, args$postproc_allloc_file), 
			' -outdir ', tmpdir,'',		
			' -file_prefix ',basename(args$out_dir), '',
			' -threshold_check ', args$threshold_check, '',
			' -use.only.hosp.deaths ', args$use.only.hosp.deaths
		)
	cmd2 <- paste0(cmd2, tmp, '\n')	
	cmd2 <- paste0(cmd2, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_allloc* ', args$out_dir,'\n'))
	cmd2 <- paste0(cmd2, 'cp -R -n --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_*rds ', args$out_dir,'\n'))
	
	# PostProc. main figures
	tmp <- paste0('Rscript ', 
	              file.path(args$source_dir, args$postproc_mainfigs), 
	              ' -outdir ', tmpdir,'',		
	              ' -file_prefix ',basename(args$out_dir), '',
	              ' -pkgdir ', args$source_dir, '',
	              ' -threshold_check ', args$threshold_check, ''
	)
	cmd2 <- paste0(cmd2, tmp, '\n')	
	
	cmd2 <- paste0(cmd2, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_allloc_figure* ', args$out_dir,'\n'))
	cmd2 <- paste0(cmd2, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/',paste0('*Brizzi* ', args$out_dir,'\n'))
	cmd2 <- paste0(cmd2, 'cp -R -n --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_*rds ', args$out_dir,'\n'))

	# make tables
	tmp <- paste0('Rscript ', 
	              file.path(args$source_dir, args$postproc_tables), 
	              ' -outdir ', tmpdir,'',		
	              ' -file_prefix ',basename(args$out_dir), '',
	              ' -pkgdir ', args$source_dir, '',
	              ' -threshold_check ', args$threshold_check, ''
	)
	cmd2 <- paste0(cmd2, tmp, '\n')	
	cmd2 <- paste0(cmd2, 'cp -R -n --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_*rds ', args$out_dir,'\n'))
	
	
	# KnitRMD
	cmd2 <- paste0(cmd2, 'chmod -R g+rw ', args$out_dir,'\n')
	tmp <- paste0('Rscript ', 
			file.path(args$source_dir, args$postproc_knit_rmd), 
			' -outdir ', args$out_dir,'',		
			' -file_prefix ',basename(args$out_dir), '',
			' -script_rmd ',args$postproc_rmd, '',
			' -pkgdir ', args$source_dir, '',
			' -threshold_check ', args$threshold_check
		)
	cmd2 <- paste0(cmd2, tmp, '\n')	
	cmd2 <- paste0(cmd2, 'chmod -R g+rw ', args$out_dir,'\n')
	cmd2 <- paste0(cmd2,"cd $CWD\n")
	cmds[[i]] <- paste0(i,')\n',cmd2,';;\n')
}
cmd2 <- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')			


########################################
#	prepare PBS header for cmdstan run and post processing for that location
########################################
pbshead <- make.PBS.header(	hpc.walltime=71, 
		hpc.select=1, 
		hpc.nproc=4, 
		hpc.mem= "96gb", 
		hpc.load= "module load anaconda3/personal\nsource activate covid19.P1.IFR\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE", 
		hpc.q=NA,
		hpc.array= length(cmds) )

########################################
#	put everything together and write to file
########################################

cmd2 <- paste(pbshead,cmd2,sep='\n')
jobfile.stan <- gsub(':','',paste("bstan",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
jobfile.stan <- file.path(args$out_dir, jobfile.stan)
cat(cmd2, file=jobfile.stan)

########################################
#	prepare command for preprocessingjobfile.stan <- file.path(args$out_dir, jobfile.stan)
########################################

cmd1 <- ''				
cmd1 <- paste0(cmd1,"CWD=$(pwd)\n")
cmd1 <- paste0(cmd1,"echo $CWD\n")	
tmpdir.prefix <- paste0('cvd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),'_pre')
tmpdir <- paste0("$CWD/",tmpdir.prefix)
cmd1 <- paste0(cmd1,"mkdir -p ",tmpdir,'\n')		
tmp <- paste0('Rscript ', 
		file.path(args$source_dir, args$preproc_file), 
		' -outdir ', tmpdir,'',		
		' -file_prefix ',basename(args$out_dir), '',
		' -pkgdir ', args$source_dir,'',
		' -stanmodel ', args$stanmodel,'',
		' -genomic_data_type ',args$genomic.data.type,'',
		' -fr_predictors ',args$fr.predictors,'', 
		' -fr_predictors_names ',args$fr.predictors.names,'', 
		' -fr_predictors_transformation ',args$fr.predictors.transformation,'',
		' -keep.hospital.type ', "'",args$keep.hospital.type, "' "
	)
cmd1 <- paste0(cmd1, tmp, '\n')
cmd1 <- paste0(cmd1,"mkdir -p ",args$out_dir,'\n')
cmd1 <- paste0(cmd1, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ', args$out_dir,'\n')
cmd1 <- paste0(cmd1, 'chmod -R g+rw ', args$out_dir,'\n')
cmd1 <- paste0(cmd1, 'cd ', dirname(jobfile.stan),'\n')
cmd1 <- paste0(cmd1, 'qsub ', basename(jobfile.stan),'\n')

########################################
#	prepare PBS header for cmdstan run and post processing for that location
########################################
pbshead <- make.PBS.header(hpc.walltime=71, 
		hpc.select=1, 
		hpc.nproc=4, 
		hpc.mem= "96gb", 
		hpc.load= "module load anaconda3/personal\nsource activate covid19.P1.IFR\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE", 
		hpc.q=NA
	)

########################################
#	put everything together and submit job
########################################

cmd1 <- paste(pbshead,cmd1,sep='\n')
jobfile.preproc <- gsub(':','',paste("bpre",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
jobfile.preproc <- file.path(args$out_dir, jobfile.preproc)
cat(cmd1, file=jobfile.preproc)
cmd <- paste("qsub", jobfile.preproc)
cat(cmd)
cat(system(cmd, intern= TRUE))
