make_health_care_predictors <- function(dhcadm, dicu, dooh, dhcnres, dhcb, dw, dr, dpop, ddates, dhr, args)
{
	dhca <- copy(dhcadm)
	
	# icu admissions (all cause) in selected hospital type
	dhca <- merge(dhca, dicu, by=c('week','loc_label'), all.x=TRUE)
	
	# out of hospital deaths
	dhca <- merge(dhca, subset(dooh, select=c(loc_label,week,deaths.ooh)), by=c('week', 'loc_label'), all.x=TRUE)
	
	# proportion of SARI admissions among non residents
	dhca <- merge(dhca, subset(dhcnres,select=c(loc_label,week,p.res.hosps)), by=c('week', 'loc_label'))
	
	# equipment and personnel
	tmp <- subset(dhcb, select=c(loc_label, week, 
		ventilator,
		physicians.all, 
		physicians.specialist,
		intensivists,
		nurse, 
		tech.nurse,
		physiotherapist,
		icubeds.all, 
		icu.per.physicians.specialist,
		icu.per.intensivists,
		icu.per.physicians.all,
		icu.per.nurse, 
		icu.per.tech.nurse,
		icu.per.physiotherapist,
		hosps.per.physicians.all, 
		hosps.per.saribeds,
		hosps.per.saribedsvent,
		hosps.per.allbeds,
		hosps.per.ventilator, 
		hosps.per.nurse, 
		hosps.per.tech.nurse, 
		hosps.per.physiotherapist))
	dhca <- merge(dhca, tmp, by=c('week','loc_label'), all.x=TRUE)
	
	# exp decay with time
	dhca <- merge(dhca, unique(subset(dw, select=c(week, week.start))), by=c('week'))
	tmp <- unique(subset(dhca, select=c(week,week.start)))
	setkey(tmp, week.start)
	tmp[, time.expdecay := c( rep(1, 19-min(tmp$week)), exp( log(0.7)/(45-19) * c(0:(45-19))), rep(.7, max(tmp$week)-45) )]
	dhca <- merge(dhca, tmp, by=c('week','week.start'), all.x=TRUE)
	
	# excess deaths from registry
	cat('\nmake hospital collapse predictor excess deaths ...')
	tmp <- dr[, list(reg.deaths.excess=sum(reg.deaths.all.excess)), by=c('loc_label','week')]
	dhca <- merge(dhca, tmp, by=c('week','loc_label'), all.x=TRUE)
	
	# set missing values to zero
	dhca <- suppressWarnings(melt(dhca,id.vars=c('loc_label','week','week.start')))
	set(dhca, dhca[,which(is.na(value))], 'value', 0.)
	dhca <- dcast.data.table(dhca, loc_label+week+week.start~variable,value.var='value')
	
	# standardise by pop size
	tmp <- subset(dpop, age.label%in%c('75-79','80-84','85-89','90+'))
	tmp <- tmp[, list(pop75 = sum(pop), tpop=tpop[1]), by='loc_label']
	dhca <- merge(dhca, tmp, by=c('loc_label'))
	dhca[, r.nhosp.adm := nhosp.adm/tpop*1e5]
	dhca[, r.nicu.adm := nicu.adm/tpop*1e5]
	dhca[, r.deaths.ooh := deaths.ooh/tpop*1e5]
	dhca[, r.reg.deaths.excess := reg.deaths.excess/tpop*1e5]
	dhca[, r.nhosp.adm.75plus := nhosp.adm.75plus/pop75*1e5]    
	set(dhca, NULL, c('pop75','tpop'), NULL)
	set(dhca, NULL, c('nhosp.adm.75plus','deaths.ooh','reg.deaths.excess'), NULL)
	
	# rollover future weeks, padding data for last two / four weeks
	dhca <- melt(dhca, id.vars=c('loc_label','week','week.start'))
	tmp <- dhca[, list( week = min(week) - 1:4, value = value[week==min(week)] ), by=c('loc_label','variable')]
	tmp <- rbind(dhca, tmp, fill=TRUE)	
	tmp2 <- dhca[, list( week = max(week) + 1:4, value = value[week==max(week)] ), by=c('loc_label','variable')]
	tmp <- rbind(tmp, tmp2, fill=TRUE)	
	setkey(tmp, loc_label, variable, week)
	tmp2 <- tmp[variable%in%c('nhosp.adm','nicu.adm','r.nhosp.adm','r.nicu.adm','r.deaths.ooh','r.reg.deaths.excess','r.nhosp.adm.75plus'),
			{
				value.20 <- data.table::frollsum(value, c(3), algo="exact", align="right")
				value.40 <- data.table::frollsum(value, c(5), algo="exact", align="right")
				value.22 <- data.table::frollsum(value, c(3), algo="exact", align="center")
				value.44 <- data.table::frollsum(value, c(5), algo="exact", align="center")		
				value.02 <- data.table::frollsum(value, c(3), algo="exact", align="left")
				value.04 <- data.table::frollsum(value, c(5), algo="exact", align="left")
				list( week=week, value.20=value.20, value.40=value.40, value.22=value.22, value.44=value.44, value.02=value.02, value.04=value.04 )		
			}, 
			by=c('loc_label','variable')]	
	tmp3 <- tmp[variable%in%c('physicians.all','physicians.specialist', 'intensivists','nurse','tech.nurse','physiotherapist','icubeds.all','ventilator','p.res.hosps','time.expdecay'),
			{
				value.20 <- data.table::frollmean(value, c(3), algo="exact", align="right")
				value.40 <- data.table::frollmean(value, c(5), algo="exact", align="right")
				value.22 <- data.table::frollmean(value, c(3), algo="exact", align="center")
				value.44 <- data.table::frollmean(value, c(5), algo="exact", align="center")		
				value.02 <- data.table::frollmean(value, c(3), algo="exact", align="left")
				value.04 <- data.table::frollmean(value, c(5), algo="exact", align="left")
				list( week=week, value.20=value.20, value.40=value.40, value.22=value.22, value.44=value.44, value.02=value.02, value.04=value.04 )		
			}, 
			by=c('loc_label','variable')]	
	tmp <- rbind(tmp2,tmp3)
	dhca <- merge(dhca, tmp, by=c('loc_label','variable','week'))
	setnames(dhca, 'value', 'value.00')
	
	# add predictors that are ratios depending on ICU admission and rolled variables
	tmp <- subset(dhca, variable=='nicu.adm')	
	tmp <- merge(tmp, subset(dhcb, select=c(loc_label,week,icubeds.all,ventilator,physicians.specialist,intensivists,physicians.all,nurse,tech.nurse,physiotherapist)), by=c('loc_label','week'))
	tmp <- melt(tmp, id.vars=c('loc_label','variable','week','week.start','ventilator','icubeds.all','physicians.specialist','intensivists','physicians.all','nurse','tech.nurse','physiotherapist'), variable.name='frollsum.type')
	tmp[, icu.per.beds := value / icubeds.all]
	tmp[, icu.per.physicians.specialist := value / physicians.specialist]
	tmp[, icu.per.intensivists := value / intensivists]
	tmp[, icu.per.physicians.all := value / physicians.all]
	tmp[, icu.per.ventilator := value / ventilator]	
	tmp[, icu.per.nurse := value / nurse]
	tmp[, icu.per.tech.nurse := value / tech.nurse]
	tmp[, icu.per.physiotherapist := value / physiotherapist]	
	tmp <- melt(tmp, 
		id.vars=c('loc_label','week','week.start','frollsum.type'),
		measure.vars=c('icu.per.beds','icu.per.physicians.specialist','icu.per.intensivists','icu.per.physicians.all','icu.per.ventilator','icu.per.nurse','icu.per.tech.nurse','icu.per.physiotherapist'))
	tmp <- dcast.data.table(tmp, loc_label+week+week.start+variable~frollsum.type, value.var='value')
	dhca <- rbind(dhca, tmp)
	
	# add predictors that are ratios depending on hosp admission and rolled variables
	tmp <- subset(dhca, variable=='nhosp.adm')
	tmp <- merge(tmp, subset(dhcb, select=c(loc_label,week,icubeds.all,saribeds.all, saribedsvent.all,beds.all,ventilator,physicians.specialist,intensivists,physicians.all,nurse,tech.nurse,physiotherapist)), by=c('loc_label','week'))
	tmp <- melt(tmp, id.vars=c('loc_label','variable','week','week.start','ventilator','icubeds.all','saribeds.all','saribedsvent.all','beds.all','physicians.specialist','intensivists','physicians.all','nurse','tech.nurse','physiotherapist'), variable.name='frollsum.type')
	tmp[, hosps.per.physicians.all := value / physicians.all]
	tmp[, hosps.per.saribeds := value / saribeds.all]
	tmp[, hosps.per.saribedsvent := value / saribedsvent.all]
	tmp[, hosps.per.allbeds := value / beds.all]
	tmp[, hosps.per.nurse := value / nurse]
	tmp[, hosps.per.ventilator := value / ventilator]
	tmp <- melt(tmp, id.vars=c('loc_label','week','week.start','frollsum.type'),measure.vars=c('hosps.per.physicians.all','hosps.per.saribeds', 'hosps.per.saribedsvent','hosps.per.allbeds','hosps.per.nurse','hosps.per.ventilator'))
	tmp <- dcast.data.table(tmp, loc_label+week+week.start+variable~frollsum.type, value.var='value')	
	dhca <- rbind(dhca, tmp)
	
	dhca <- subset(dhca, !variable%in%c('nicu.adm','nhosp.adm','physicians.all','physicians.specialist', 'intensivists','nurse','tech.nurse','physiotherapist','icubeds.all','interm.beds.all','interm.beds.all','ventilator'))
	
	# We actually want correlations with age-std hfr. Do here or in text? 
	# determine correlation with all-age hfr
	dhca <- merge(dhca, dhfrt, by=c('loc_label','week.start'), all.x=TRUE)
	dhcc <- dhca[!is.na(hfr), 
			list(
					cor.00=cor(value.00,hfr), 
					cor.20=cor(value.20,hfr), 
					cor.40=cor(value.40,hfr),
					cor.22=cor(value.22,hfr),
					cor.44=cor(value.44,hfr),
					cor.02=cor(value.02,hfr),
					cor.04=cor(value.04,hfr)
			), 
			by=c('loc_label','variable')]
	
	# count which rollsum wins for each variable
	dhcc <- melt(dhcc, id.vars=c('loc_label','variable'), variable.name='frollsum.type') 
	tmp <- dhcc[, list(cor.max= frollsum.type[which.max(value)] ), by=c('variable','loc_label')]
	tmp <- tmp[, list(value=length(loc_label)), by=c('variable','cor.max')]
	setkey(tmp, variable, cor.max)
	dhccw <- tmp[, list(frollsum.type=cor.max[which.max(value)]), by='variable']
	set(dhccw, NULL, 'frollsum.type', dhccw[,gsub('cor','value',frollsum.type)])
	dhca <- melt(dhca, id.vars=c('loc_label','variable', 'week', 'week.start', 'biweek', 'hfr'), variable.name='frollsum.type')
	set(dhca, NULL, 'variable', dhca[, paste0(variable,gsub('value','',frollsum.type))])
	set(dhca, NULL, 'frollsum.type', NULL)
	set(dhccw, NULL, 'variable', dhccw[, paste0(variable,gsub('value','',frollsum.type))])
	set(dhccw, NULL, 'frollsum.type', NULL)
	set(dhcc, NULL, 'variable', dhcc[, paste0(variable,gsub('cor','',frollsum.type))])
	set(dhcc, NULL, 'frollsum.type', NULL)
	
	cat('\ncorrelation predictors <-> HFR\n')
	tmp <- dhcc[, list(avg_abs_cor_with_hfr=mean(abs(value)), avg_cor_with_hfr=mean(value)), by='variable']
	tmp <- tmp[order(-avg_abs_cor_with_hfr)]
	print(tmp, n=200)
	
	cat('\nhighest correlating predictors are\n', dhccw$variable)
	cat('\nselected  predictors are\n', args$fr.predictors)
		
	# plot predictor against hfr
	cat('\nplotting predictors ...')
	z <- sort(unique(c(dhccw$variable)))
	z <- z[!z%in%args$fr.predictors]
	z <- c(z, args$fr.predictors)
	
	names(z) <- gsub('2', '3',names(z))
		
	if(args$make.preprocessing.plots)
	{			
		for(i in seq_along(z))
		{
			tmp <- subset(dhca, variable==z[i])
			tmp2 <- subset(dhcc, variable==z[i])
			tmp3 <- tmp[, range(value)]
			tmp4 <- tmp[, range(hfr, na.rm=TRUE)*c(1,1.2)]
			p <- ggplot(tmp, aes(x=week.start)) +
					geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
					geom_step(data=subset(tmp, !is.na(hfr)), aes(y=( hfr - tmp4[1] ) / diff(tmp4) * diff(tmp3) + tmp3[1]), color = 'black') +
					geom_point(data=subset(tmp, !is.na(hfr)), aes(y=( hfr - tmp4[1] ) / diff(tmp4) * diff(tmp3) + tmp3[1]), color = 'black') +
					geom_step(aes(y=value, colour=change_city_label(loc_label))) +
					geom_point(aes(y=value, colour=change_city_label(loc_label))) +
					geom_text(data=tmp2, aes(x= min(tmp$week.start)+70, y= max(tmp$value)*.96, label= paste0('r=',round(value, d=2))), hjust=0.5, vjust=0.5, colour='black') +
					scale_x_date(breaks='2 months',expand=c(0,0)) +
					scale_y_continuous(expand=c(0,0)) +
					scale_y_continuous(sec.axis= sec_axis(~ (.-tmp3[1])/diff(tmp3)*diff(tmp4)+tmp4[1], name = "in-hospital fatality rate", labels = scales::percent), expand=c(0,0)) +
					scale_colour_manual(values = args$city.palette(length(unique(tmp$loc_label)))) +
					facet_wrap(~change_city_label(loc_label), ncol=5) +
					labs(x='', y=ifelse(is.na(names(z)[i]) | names(z)[i]=='', z[i], names(z)[i]), colour='') +
					guides(colour='none') +
					theme_bw() +
					theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
							strip.background = element_blank(),
							legend.position='bottom',
							axis.line.y.right = element_line(color = 'black'), 
							axis.ticks.y.right = element_line(color = 'black'),
							axis.text.y.right = element_text(color = 'black'))			
			ggsave(file=paste0(args$out.base,'_hfr_predictors_',z[i],'_bycity.pdf'),p,w=12,h=8)
			ggsave(file=paste0(args$out.base,'_hfr_predictors_',z[i],'_bycity.png'),p,w=12,h=8)
		}		
	}
	
	#	add reference week with lowest HFR for all locations to dhca
	dhca <- merge(dhca, subset(ddates, select = c('loc_label', 'w1.start', 'w2.start')), by = 'loc_label')
	dhca <- merge(dhca, dhr, by='loc_label')
	
	#	determine weeks in initial phase, non.P1 phase, and P.1 phase
	dhca[, week.cat := 'with.P1']
	set(dhca, dhca[, which(week.start<w2.start)], 'week.cat', 'before.P1')
	set(dhca, dhca[, which(week.start<w1.start)], 'week.cat', 'initial')
	set(dhca, NULL, c('w1.start', 'w2.start'), NULL)
	# there should not be any more initial since hosps.start is on / after w1.start
	# initial no longer excluded in below
	stopifnot( !nrow(subset(dhca, week.cat=='initial')) ) 
	
	#	after reference week is determined, set excess deaths to non-negative
	tmp <- dhca[, which(grepl('excess',variable))]
	set(dhca, tmp, 'value', dhca[tmp, pmax(0, value)])
	stopifnot(!nrow(subset(dhca, value<0)))
	
	return(list(dhca=dhca, dhcc=dhcc, dhccw=dhccw))
}






make_health_care_predictors.b2 <- function(dhcadm, dicu, dhcv, dooh, dhcnres, dhcb, dw, dr, dpop, ddates, dhr, args)
{
	dhca <- copy(dhcadm)
	
	# icu admissions (all cause) in selected hospital type
	dhca <- merge(dhca, dicu, by=c('week','loc_label'), all.x=TRUE)
	
	# ventilator predictors
	dhca <- merge(dhca, subset(dhcv, select=c(loc_label, week, ventilator_total, r.ventilator)), by=c('week','loc_label'), all.x=TRUE)
	
	# out of hospital deaths
	dhca <- merge(dhca, subset(dooh, select=c(loc_label,week,deaths.ooh)), by=c('week', 'loc_label'), all.x=TRUE)
	
	# proportion of SARI admissions among non residents
	dhca <- merge(dhca, dhcnres, by=c('week', 'loc_label'))
	
	# beds and physicians
	dhca <- merge(dhca, subset(dhcb, select=c(loc_label, week, physicians.all, physicians.specialist,  icubeds.all, icu.per.physicians.specialist, icu.per.physicians.all, hosps.per.physicians.all, hosps.per.saribeds, hosps.per.saribedsvent)), by=c('week','loc_label'), all.x=TRUE)
	
	# exp decay with time
	dhca <- merge(dhca, unique(subset(dw, select=c(week, week.start))), by=c('week'))
	tmp <- unique(subset(dhca, select=c(week,week.start)))
	setkey(tmp, week.start)
	tmp[, time.expdecay := c( rep(1, 19-min(tmp$week)), exp( log(0.7)/(45-19) * c(0:(45-19))), rep(.7, max(tmp$week)-45) )]
	dhca <- merge(dhca, tmp, by=c('week','week.start'), all.x=TRUE)
	
	# excess deaths from registry
	cat('\nmake hospital collapse predictor excess deaths ...')
	tmp <- dr[, list(reg.deaths.excess=sum(reg.deaths.all.excess)), by=c('loc_label','week')]
	dhca <- merge(dhca, tmp, by=c('week','loc_label'), all.x=TRUE)
	
	# set missing values to zero
	for(x in c('nhosp.adm','nhosp.adm.75plus','nicu.adm','reg.deaths.ooh','reg.deaths.excess'))
	{
		set(dhca, which(is.na(dhca[[x]])), x, 0.)
	}
	
	# standardise by pop size
	tmp <- subset(dpop, age.label%in%c('75-79','80-84','85-89','90+'))
	tmp <- tmp[, list(pop75 = sum(pop), tpop=tpop[1]), by='loc_label']
	dhca <- merge(dhca, tmp, by=c('loc_label'))
	dhca[, r.nhosp.adm := nhosp.adm/tpop*1e5]
	dhca[, r.nicu.adm := nicu.adm/tpop*1e5]
	dhca[, r.deaths.ooh := deaths.ooh/tpop*1e5]
	dhca[, r.reg.deaths.excess := reg.deaths.excess/tpop*1e5]
	dhca[, r.nhosp.adm.75plus := nhosp.adm.75plus/pop75*1e5]    
	set(dhca, NULL, c('pop75','tpop'), NULL)
	set(dhca, NULL, c('nhosp.adm.75plus','deaths.ooh','reg.deaths.excess'), NULL)
	
	# rollover future weeks, padding data for last two / four weeks
	dhca <- melt(dhca, id.vars=c('loc_label','week','week.start'))
	tmp <- dhca[, list( week = min(week) - 1:4, value = value[week==min(week)] ), by=c('loc_label','variable')]
	tmp <- rbind(dhca, tmp, fill=TRUE)	
	tmp2 <- dhca[, list( week = max(week) + 1:4, value = value[week==max(week)] ), by=c('loc_label','variable')]
	tmp <- rbind(tmp, tmp2, fill=TRUE)	
	setkey(tmp, loc_label, variable, week)
	tmp2 <- tmp[!variable%in%c('physicians.all','physicians.specialist','icubeds.all','ventilator_total','p.res.hosps','p.nonres.hosps','time.expdecay'),
			{
				value.20 <- data.table::frollsum(value, c(3), algo="exact", align="right")
				value.40 <- data.table::frollsum(value, c(5), algo="exact", align="right")
				value.22 <- data.table::frollsum(value, c(3), algo="exact", align="center")
				value.44 <- data.table::frollsum(value, c(5), algo="exact", align="center")		
				value.02 <- data.table::frollsum(value, c(3), algo="exact", align="left")
				value.04 <- data.table::frollsum(value, c(5), algo="exact", align="left")
				list( week=week, value.20=value.20, value.40=value.40, value.22=value.22, value.44=value.44, value.02=value.02, value.04=value.04 )		
			}, 
			by=c('loc_label','variable')]	
	tmp3 <- tmp[variable%in%c('physicians.all','physicians.specialist','icubeds.all','ventilator_total','p.res.hosps','p.nonres.hosps','time.expdecay'),
			{
				value.20 <- data.table::frollmean(value, c(3), algo="exact", align="right")
				value.40 <- data.table::frollmean(value, c(5), algo="exact", align="right")
				value.22 <- data.table::frollmean(value, c(3), algo="exact", align="center")
				value.44 <- data.table::frollmean(value, c(5), algo="exact", align="center")		
				value.02 <- data.table::frollmean(value, c(3), algo="exact", align="left")
				value.04 <- data.table::frollmean(value, c(5), algo="exact", align="left")
				list( week=week, value.20=value.20, value.40=value.40, value.22=value.22, value.44=value.44, value.02=value.02, value.04=value.04 )		
			}, 
			by=c('loc_label','variable')]	
	tmp <- rbind(tmp2,tmp3)
	dhca <- merge(dhca, tmp, by=c('loc_label','variable','week'))
	setnames(dhca, 'value', 'value.00')
	
	# add predictors that are ratios depending on ICU admission and rolled variables
	tmp <- subset(dhca, variable=='nicu.adm')
	tmp <- merge(tmp, subset(dhcv, select=c(loc_label,week,ventilator_total)), by=c('loc_label','week'))
	tmp <- merge(tmp, subset(dhcb, select=c(loc_label,week,icubeds.all,physicians.specialist,physicians.all)), by=c('loc_label','week'))
	tmp <- melt(tmp, id.vars=c('loc_label','variable','week','week.start','ventilator_total','icubeds.all','physicians.specialist','physicians.all'), variable.name='frollsum.type')
	tmp[, icu.per.beds := value / icubeds.all]
	tmp[, icu.per.physicians.specialist := value / physicians.specialist]
	tmp[, icu.per.physicians.all := value / physicians.all]
	tmp[, icu.per.ventilator := value / ventilator_total]
	tmp <- melt(tmp, id.vars=c('loc_label','week','week.start','frollsum.type'),measure.vars=c('icu.per.beds','icu.per.physicians.specialist','icu.per.physicians.all','icu.per.ventilator'))
	tmp <- dcast.data.table(tmp, loc_label+week+week.start+variable~frollsum.type, value.var='value')
	dhca <- subset(dhca, !variable%in%c('icu.per.beds','icu.per.physicians.specialist','icu.per.physicians.all','icu.per.ventilator'))
	dhca <- rbind(dhca, tmp)
	
	# add predictors that are ratios depending on hosp admission and rolled variables
	tmp <- subset(dhca, variable=='nhosp.adm')
	tmp <- merge(tmp, subset(dhcb, select=c(loc_label,week,icubeds.all,physicians.all)), by=c('loc_label','week'))
	tmp <- melt(tmp, id.vars=c('loc_label','variable','week','week.start','icubeds.all','physicians.all'), variable.name='frollsum.type')
	tmp[, hosps.per.physicians.all := value / physicians.all]
	tmp[, hosps.per.saribeds := value / saribeds.all] 
	tmp[, hosps.per.saribedsvent := value / saribedsvent.all]
	tmp <- melt(tmp, id.vars=c('loc_label','week','week.start','frollsum.type'),measure.vars=c('hosps.per.physicians.all','hosps.per.saribeds','hosps.per.saribedsvent'))
	tmp <- dcast.data.table(tmp, loc_label+week+week.start+variable~frollsum.type, value.var='value')
	dhca <- subset(dhca, !variable%in%c('hosps.per.physicians.all','hosps.per.saribeds', 'hosps.per.saribedsvent'))
	dhca <- rbind(dhca, tmp)
	
	dhca <- subset(dhca, !variable%in%c('nicu.adm','ventilator_total','nhosp.adm','physicians.all','physicians.specialist','icubeds.all'))
	
	# determine correlation with all-age hfr
	dhca <- merge(dhca, dhfrt, by=c('loc_label','week.start'), all.x=TRUE)
	dhcc <- dhca[!is.na(hfr), 
			list(
					cor.00=cor(value.00,hfr), 
					cor.20=cor(value.20,hfr), 
					cor.40=cor(value.40,hfr),
					cor.22=cor(value.22,hfr),
					cor.44=cor(value.44,hfr),
					cor.02=cor(value.02,hfr),
					cor.04=cor(value.04,hfr)
			), 
			by=c('loc_label','variable')]
	
	# count which rollsum wins for each variable
	dhcc <- melt(dhcc, id.vars=c('loc_label','variable'), variable.name='frollsum.type') 
	tmp <- dhcc[, list(cor.max= frollsum.type[which.max(value)] ), by=c('variable','loc_label')]
	tmp <- tmp[, list(value=length(loc_label)), by=c('variable','cor.max')]
	setkey(tmp, variable, cor.max)
	dhccw <- tmp[, list(frollsum.type=cor.max[which.max(value)]), by='variable']
	set(dhccw, NULL, 'frollsum.type', dhccw[,gsub('cor','value',frollsum.type)])
	dhca <- melt(dhca, id.vars=c('loc_label','variable', 'week', 'week.start', 'biweek', 'hfr'), variable.name='frollsum.type')
	set(dhca, NULL, 'variable', dhca[, paste0(variable,gsub('value','',frollsum.type))])
	set(dhca, NULL, 'frollsum.type', NULL)
	set(dhccw, NULL, 'variable', dhccw[, paste0(variable,gsub('value','',frollsum.type))])
	set(dhccw, NULL, 'frollsum.type', NULL)
	set(dhcc, NULL, 'variable', dhcc[, paste0(variable,gsub('cor','',frollsum.type))])
	set(dhcc, NULL, 'frollsum.type', NULL)
	
	cat('\ncorrelation predictors <-> HFR\n')
	tmp <- dhcc[, list(avg_abs_cor_with_hfr=mean(abs(value)), avg_cor_with_hfr=mean(value)), by='variable']
	tmp <- tmp[order(-avg_abs_cor_with_hfr)]
	print(tmp, n=200)
	
	cat('\nhighest correlating predictors are\n', dhccw$variable)
	cat('\nselected  predictors are\n', args$fr.predictors)
	
	# plot predictor against hfr
	cat('\nplotting predictors ...')
	z <- sort(unique(c(dhccw$variable)))
	z <- z[!z%in%args$fr.predictors]
	z <- c(z, args$fr.predictors)
	
	if(args$make.preprocessing.plots)
	{			
		for(i in seq_along(z))
		{
			tmp <- subset(dhca, variable==z[i])
			tmp2 <- subset(dhcc, variable==z[i])
			tmp3 <- tmp[, range(value)]
			tmp4 <- tmp[, range(hfr, na.rm=TRUE)*c(1,1.2)]
			p <- ggplot(tmp, aes(x=week.start)) +
					geom_vline(data=ddates, aes(xintercept=w2.start), colour='grey50') +
					geom_step(data=subset(tmp, !is.na(hfr)), aes(y=( hfr - tmp4[1] ) / diff(tmp4) * diff(tmp3) + tmp3[1]), color = 'black') +
					geom_point(data=subset(tmp, !is.na(hfr)), aes(y=( hfr - tmp4[1] ) / diff(tmp4) * diff(tmp3) + tmp3[1]), color = 'black') +
					geom_step(aes(y=value, colour=change_city_label(loc_label))) +
					geom_point(aes(y=value, colour=change_city_label(loc_label))) +
					geom_text(data=tmp2, aes(x= min(tmp$week.start)+70, y= max(tmp$value)*.96, label= paste0('r=',round(value, d=2))), hjust=0.5, vjust=0.5, colour='black') +
					scale_x_date(breaks='2 months',expand=c(0,0)) +
					scale_y_continuous(expand=c(0,0)) +
					scale_y_continuous(sec.axis= sec_axis(~ (.-tmp3[1])/diff(tmp3)*diff(tmp4)+tmp4[1], name = "in-hospital fatality rate", labels = scales::percent), expand=c(0,0)) +
					scale_colour_manual(values = args$city.palette(length(unique(dhcv$loc_label)))) +
					facet_wrap(~change_city_label(loc_label), ncol=5) +
					labs(x='', y=ifelse(is.na(names(z)[i]) | names(z)[i]=='', z[i], names(z)[i]), colour='') +
					guides(colour='none') +
					theme_bw() +
					theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
							strip.background = element_blank(),
							legend.position='bottom',
							axis.line.y.right = element_line(color = 'black'), 
							axis.ticks.y.right = element_line(color = 'black'),
							axis.text.y.right = element_text(color = 'black'))			
			ggsave(file=paste0(args$out.base,'_hfr_predictors_',z[i],'_bycity.pdf'),p,w=12,h=8)
			ggsave(file=paste0(args$out.base,'_hfr_predictors_',z[i],'_bycity.png'),p,w=12,h=8)
		}		
	}
	
	#	add reference week with lowest HFR for all locations to dhca
	dhca <- merge(dhca, subset(ddates, select = c('loc_label', 'w1.start', 'w2.start')), by = 'loc_label')
	dhca <- merge(dhca, dhr, by='loc_label')
	
	#	determine weeks in initial phase, non.P1 phase, and P.1 phase
	dhca[, week.cat := 'with.P1']
	set(dhca, dhca[, which(week.start<w2.start)], 'week.cat', 'before.P1')
	set(dhca, dhca[, which(week.start<w1.start)], 'week.cat', 'initial')
	set(dhca, NULL, c('w1.start', 'w2.start'), NULL)
	# there should not be any more initial since hosps.start is on / after w1.start
	# initial no longer excluded in below
	stopifnot( !nrow(subset(dhca, week.cat=='initial')) ) 
	
	#	after reference week is determined, set excess deaths to non-negative
	tmp <- dhca[, which(grepl('excess',variable))]
	set(dhca, tmp, 'value', dhca[tmp, pmax(0, value)])
	stopifnot(!nrow(subset(dhca, value<0)))
	
	return(list(dhca=dhca, dhcc=dhcc, dhccw=dhccw))
}