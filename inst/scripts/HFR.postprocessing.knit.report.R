# postprocessing.knit.report.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n HFR.postprocessing.knit.report.R \n \n -------------------------------- \n")

suppressMessages(library(rmarkdown, quietly = TRUE))

if(1)
{  
	pkg.dir2 <- '~/git/covid19_brazil_hfr/'
	out.dir2 <- 'test'
	file.prefix2 <- basename(out.dir2)  
	script.rmd <- '/scripts/HFR.postprocessing.make.report.Rmd'
	threshold_check <- 1
}


# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-outdir')	
	stopifnot(args_line[[3]]=='-file_prefix')
	stopifnot(args_line[[5]]=='-script_rmd')
	stopifnot(args_line[[7]]=='-pkgdir')
	stopifnot(args_line[[9]]=='-threshold_check')
	out.dir2 <- args_line[[2]]
	file.prefix2 <- args_line[[4]]
	script.rmd <- args_line[[6]]
	pkg.dir2 <- args_line[[8]]
	threshold_check <- as.integer(args_line[[10]])
} 

####################################################################
# Check that at least threshold number of chains are done running
####################################################################

tmp <- list.files(out.dir2, pattern='_fit.rds')
cat('\ntry to bugfix: start\n')
print(out.dir2)
print(tmp)
cat('\ntry to bugfix: end\n')
if(length(tmp) < threshold_check)
{
	stop("Not enough runs have completed yet!")
}

####################################################################
# source functions
####################################################################

tmp <- list.files(file.path(pkg.dir2,'inst/R'), full.names=TRUE)
cat('sourcing code from files\n', tmp)
invisible(lapply(tmp, source))

####################################################################
# source workspace
####################################################################

tmp <- paste0(file.path(out.dir2,file.prefix2), '_allloc_workspace.rda')
cat('\nAttempting to load ', tmp)
load(file = tmp)

####################################################################
# execute RMD script
####################################################################

##	
cat(paste0("\n ----------- create report ----------- \n"))
cat(paste0("\nRMD script ",file.path(pkg.dir2, '/inst/', script.rmd)))

if(!grepl('^\\/', out.dir2)){
  out.dir2 <- file.path(getwd(), out.dir2)
}

cat(paste0("\noutput ",paste0(file.path(out.dir2,file.prefix2), '_allloc_report.html')))

Sys.setenv(RSTUDIO_PANDOC= "/usr/lib/rstudio/bin/pandoc")

rmarkdown::render( file.path(pkg.dir2,'inst', script.rmd), 
		output_file= paste0(file.path(out.dir2,file.prefix2), '_allloc_report.html'), 
		params = list(
				stanModelFile = args$file.stan.model,
				job_dir= out.dir2,
				job_tag= file.prefix2,				
				selected.locs = args$selected.locs
		))

cat(paste("\n -------------------------------- \n \n post-processing-knit-report.R \n \n -------------------------------- \n"))
