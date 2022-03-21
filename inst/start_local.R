# Specify output directory and run type
nocens <- 0;  unkdis <- 0; fiocruz <- 0; random <- 0

args <- list()	
args$source_dir <- '~/git/covid19_brazil_hfr/inst'
args$pkgdir <- dirname(args$source_dir)
# CAREFUL: DO NOT use underscores in the args$out_dir as this will cause a bug in HFR.postporcessing.allloc, in the search for available outputs
args$out_dir <- '~/Documents/P1Brazil/hfr-analyses'

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if ("-outdir" %in% args_line)
{
  tmp <- which(args_line == "-outdir")
  args$out_dir <- args_line[[tmp + 1]]
  if(grepl('^-',args$out_dir)){stop('Output directory starts with "-" ')}
  if(grepl('_', basename(args$out_dir))){stop('Cannot include underscore in output directory name')}
  
  args_line[[tmp + 1]] <- NULL; args_line[[tmp]] <- NULL 
}

if(length(args_line)>3)
{
  stop("At most 3 options different that -outdir can be chosen")
}
if(length(args_line) > 0) 
{
  nocens <- if("-nocens" %in% args_line){1}else{0}
  unkdis <- if("-unkdis" %in% args_line){1}else{0}
  fiocruz <- if("-fiocruz" %in% args_line){1}else{0}
  random <- if("-random" %in% args_line){1}else{0}
}
if(random & fiocruz)
{
  stop('Can only chose one of the two following options:\n-fiocruz\n-random')
}

# Add options nickname to output directory to inform HFR.preprocessing.R of type of wanted analysis
if(nocens){ args$out_dir <- paste0(args$out_dir, '-nocens')}
if(unkdis){args$out_dir <- paste0(args$out_dir, '-unkdis')}
if(fiocruz){args$out_dir <- paste0(args$out_dir, '-fiocruz')}
if(random){args$out_dir <- paste0(args$out_dir, '-random')}
 

# Specify locations and pandemic load indices variable names.

# For starters, just include manaus
# loc_labels <- 'manaus'
loc_labels <- c('belo-horizonte','curitiba', "florianopolis", "goiania", "joao-pessoa", "macapa", 'manaus', "natal", 'porto-alegre', "porto-velho", 'rio-de-janeiro', 'salvador', 'sao-paulo', "sao-luis")

args$fr.predictors <- "icu.per.ventilator.02,icu.per.beds.02,icu.per.nurse.02,icu.per.tech.nurse.02,icu.per.physiotherapist.02,icu.per.physicians.all.02,icu.per.intensivists.02,hosps.per.saribeds.02,hosps.per.saribedsvent.02,hosps.per.physicians.all.02,hosps.per.ventilator.02,hosps.per.nurse.02"
args$fr.predictors.names <- gsub(' ','_',"2-wk ICU admissions per ventilator,2-wk ICU admissions per ICU bed,2-wk ICU admissions per nurse,2-wk ICU admissions per nurse assistant,2-wk ICU admissions per physiotherapist,2-wk ICU admissions per physician,2-wk ICU admissions per intensivist,2-wk SARI admissions per critical bed, 2-wk SARI admissions per critical bed with ventilator,2-wk SARI admissions per physician,2-wk SARI admissions per ventilator,2-wk SARI admissions per nurse")

# Specify script locations
args$stanmodel = 'age_hfr_210719d.stan'
args$cmdstan_file= 'scripts/HFR.fit.cmdstan.R'
args$preproc_file= 'scripts/HFR.preprocessing.R'
args$postproc_oneloc_file = 'scripts/HFR.postprocessing.oneloc.R'
args$postproc_allloc_file = 'scripts/HFR.postprocessing.allloc.R'
args$postproc_knit_rmd = 'scripts/HFR.postprocessing.knit.report.R'
args$postproc_rmd = 'scripts/HFR.postprocessing.make.report.Rmd'		
# Add tables and main figures?
args$tables_file = 'scripts/HFR.postprocessing.tables.R'
args$mainfigures_file = 'scripts/HFR.postprocessing.main.figures.R'
args$fr.exclude.P1 = 0
args$use.only.hosp.deaths = 1
args$genomic.data.type <- 'GISAID_darlan'
args$fr.predictors.transformation <- 'value_cs'	
args$keep.hospital.type <-  'Private_Hospitals,Public_Hospitals,Unknown'
args$threshold_check = length(loc_labels)
# args$out_dir <- paste0(args$out_dir,'-ex',args$fr.exclude.P1,'-tr',gsub('_','',args$fr.predictors.transformation),'-hd',args$use.only.hosp.deaths, '-onlyconfirmed')

# If want to perform sensitivity runs with different data:
if(grepl('-fiocruz', args$out_dir))
{
  args$genomic.data.type <- 'genomica_fiocruz'
}
if(grepl("-random", args$out_dir))
{
  args$genomic.data.type <- "'random_sample|GISAID_alfonso_lutz'"
  loc_labels <- c('belo-horizonte', 'manaus', 'sao-paulo') 
  args$threshold_check = 3
}

stopifnot(!dir.exists(args$out_dir))
dir.create(args$out_dir)	

########################################
#	prepare command for preprocessingj
########################################

cmd1 <- ''				
# cmd1 <- paste0(cmd1,"CWD=$(pwd)\n")
# cmd1 <- paste0(cmd1,"echo $CWD\n")	
# tmpdir.prefix <- paste0('cvd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),'_pre')
# tmpdir <- paste0("$CWD/",tmpdir.prefix)
# cmd1 <- paste0(cmd1,"mkdir -p ",tmpdir,'\n')		
tmp <- paste0('Rscript ', 
              file.path(args$preproc_file), 
              ' -outdir ', args$out_dir,'',		
              ' -file_prefix ',basename(args$out_dir), '',
              ' -pkgdir ', args$pkgdir,'',
              ' -stanmodel ', args$stanmodel,'',
              ' -genomic_data_type ',args$genomic.data.type,'',
              ' -fr_predictors ',args$fr.predictors,'', 
              ' -fr_predictors_names ',args$fr.predictors.names,'', 
              ' -fr_predictors_transformation ',args$fr.predictors.transformation,'',
              ' -keep.hospital.type ', "'",args$keep.hospital.type, "' "
)
cmd1 <- paste0(cmd1, tmp, '\n')
# cmd1 <- paste0(cmd1,"mkdir -p ",args$out_dir,'\n')
# cmd1 <- paste0(cmd1, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ', args$out_dir,'\n')
# cmd1 <- paste0(cmd1, 'chmod -R g+rw ', args$out_dir,'\n')
# cmd1 <- paste0(cmd1, 'cd ', dirname(jobfile.stan),'\n')
# cmd1 <- paste0(cmd1, 'qsub ', basename(jobfile.stan),'\n')

########################################
#	prepare commands for cmdstan run and post processing for that location 
########################################
cmds <- vector('list',length(loc_labels))
for(i in seq_along(loc_labels))
{
  cmd2 <- ''		
  # cmd2 <- paste0(cmd2,"CWD=$(pwd)\n")
  # cmd2 <- paste0(cmd2,"echo $CWD\n")	
  # tmpdir.prefix <- paste0('cvd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),'_',i)
  # tmpdir <- paste0("$CWD/",tmpdir.prefix)
  # cmd2 <- paste0(cmd2,"mkdir -p ",tmpdir,'\n')
  # cmd2 <- paste0(cmd2,"cp ",file.path(args$out_dir,paste0(basename(args$out_dir),'_allloc_workspace.rda'))," ",tmpdir,'\n')
  tmp <- paste0('Rscript ', 
                file.path(args$source_dir, args$cmdstan_file), 
                ' -outdir ', args$out_dir,'',		
                ' -loc_label ', loc_labels[i], '',
                ' -file_prefix ',basename(args$out_dir), '',
                ' -pkgdir ', args$pkgdir, '',
                ' -fr.exclude.P1 ', args$fr.exclude.P1, '',
                ' -use.only.hosp.deaths ', args$use.only.hosp.deaths
  )
  cmd2 <- paste0(cmd2, tmp, '\n')
  cmd2 <- paste0(cmd2, tmp, '\n')	# note: not a mistake. we re-run by default to catch unexpected stan failure
  cmd2 <- paste0(cmd2, tmp, '\n')	# note: not a mistake. we re-run by default to catch unexpected stan failure
  cmd2 <- paste0(cmd2, tmp, '\n') # note: not a mistake. we re-run by default to catch unexpected stan failure
  
  tmp <- paste0('Rscript ', 
                file.path(args$source_dir, args$postproc_oneloc_file), 
                ' -outdir ', args$out_dir,'',		
                ' -loc_label ', loc_labels[i], '',
                ' -file_prefix ',basename(args$out_dir), '',
                ' -pkgdir ', args$pkgdir
  )
  cmd2 <- paste0(cmd2, tmp, '\n')
  # cmd2 <- paste0(cmd2, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_',gsub('-|_| ','',loc_labels[i])),'* ', args$out_dir,'\n')
  # cmd2 <- paste0(cmd2,"cp ",file.path(args$out_dir,'*workspace.rda')," ",tmpdir,'\n')
  # cmd2 <- paste0(cmd2,"cp ",file.path(args$out_dir,'*fit.rds')," ",tmpdir,'\n')
  cmds[[i]] <- paste0(cmd2,'\n')
}
cmd2 <- (paste0(cmds, collapse=''))	

if(length(loc_labels) > 1)
{
  
  tmp <- paste0('Rscript ', 
                file.path(args$source_dir, args$postproc_allloc_file), 
                ' -outdir ', args$out_dir,'',		
                ' -file_prefix ',basename(args$out_dir), '',
                ' -threshold_check ', args$threshold_check, '',
                ' -use.only.hosp.deaths ', args$use.only.hosp.deaths
  )
  cmd2 <- paste0(cmd2, tmp, '\n')	
  # cmd2 <- paste0(cmd2, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/',paste0(basename(args$out_dir),'_allloc* ', args$out_dir,'\n'))
  # cmd2 <- paste0(cmd2, 'chmod -R g+rw ', args$out_dir,'\n')
  tmp <- paste0('Rscript ', 
                file.path(args$source_dir, args$postproc_knit_rmd), 
                ' -outdir ', args$out_dir,'',		
                ' -file_prefix ',basename(args$out_dir), '',
                ' -script_rmd ',args$postproc_rmd, '',
                ' -pkgdir ', args$pkgdir, '',
                ' -threshold_check ', args$threshold_check
  )
  cmd2 <- paste0(cmd2, tmp, '\n')	
  # cmd2 <- paste0(cmd2, 'chmod -R g+rw ', args$out_dir,'\n')
  # cmd2 <- paste0(cmd2,"cd $CWD\n")
  tmp <- paste0('Rscript ', 
                file.path(args$source_dir, args$mainfigures_file), 
                ' -outdir ', args$out_dir,'',		
                ' -file_prefix ',basename(args$out_dir), '',
                ' -pkgdir ', args$pkgdir, '',
                ' -threshold_check ', args$threshold_check
  )
  cmd2 <- paste0(cmd2, tmp, '\n')	
  

  
  tmp <- paste0('Rscript ', 
                file.path(args$source_dir, args$tables_file), 
                ' -outdir ', args$out_dir,'',		
                ' -file_prefix ',basename(args$out_dir), '',
                ' -pkgdir ', args$pkgdir, ''
                 # ' -threshold_check ', args$threshold_check
  )
  cmd2 <- paste0(cmd2, tmp, '\n')	
}


########################################
#	put everything together and write to file
########################################

cmd2 <- paste(cmd1,cmd2,sep='\n')
jobfile.stan <- gsub(':','',paste("bstan",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
jobfile.stan <- file.path(args$source_dir, jobfile.stan)
cat(cmd2, file=jobfile.stan)

########################################
# Submit analysis
#######################################
tmp <- paste0('mv ', jobfile.stan, ' ', args$out_dir)
cmd2 <- paste0(cmd2, tmp, '\n')

tmp <- file.path(args$out_dir, basename(jobfile.stan))
tmp <- paste0(tmp,'.log.txt' )

system(paste0('sh ', jobfile.stan, ' | tee ',tmp))

