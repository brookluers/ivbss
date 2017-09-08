library(dplyr)
library(ggplot2)
library(reshape2)
d <- read.csv('directions_sep1.txt', header=T)
maxdir <- 10
nplotdir <- 5
dspeed <- filter(d, varname %in% paste("Speed[", -(30:0), "]", sep=''))
dspeed <- left_join(dspeed, data.frame(varname = paste("Speed[", -(30:0), "]", sep=''), t = 30:0), by=c('varname'='varname'))
dspeedm <- melt(dspeed, id=c('varname','t'))
pcnames <- paste('pc',1:maxdir,sep='')
cdnames <-paste('cd',1:maxdir,sep='')
dspeedm$pc <- dspeedm$variable %in% pcnames
dspeedm$cd <- dspeedm$variable %in% cdnames
dspeedm$dircat <- with(dspeedm, ifelse(pc, "PC", ifelse(cd, "Cov. Dir.", "Mean Dir.")))
dspeedm <- left_join(dspeedm, data.frame(name=c('meandir',pcnames,cdnames), dir_ix=c(1,1:maxdir,1:maxdir)), by=c('variable'='name'))

drange <- filter(d, varname %in% paste("FcwRange[", -(30:0), "]", sep=''))
drange <- left_join(drange, data.frame(varname = paste("FcwRange[", -(30:0), "]", sep=''), t = 30:0), by=c('varname'='varname'))
drangem <- melt(drange, id=c('varname','t'))
pcnames <- paste('pc',1:maxdir,sep='')
cdnames <-paste('cd',1:maxdir,sep='')
drangem$pc <- drangem$variable %in% pcnames
drangem$cd <- drangem$variable %in% cdnames
drangem$dircat <- with(drangem, ifelse(pc, "PC", ifelse(cd, "Cov. Dir.", "Mean Dir.")))
drangem <- left_join(drangem, data.frame(name=c('meandir',pcnames,cdnames), dir_ix=c(1,1:maxdir,1:maxdir)), by=c('variable'='name'))


pspeed <- ggplot(filter(dspeedm,dir_ix<nplotdir | is.na(dir_ix)), aes(x=-t, y= value)) + geom_line(aes(group=variable,color=as.factor(dir_ix))) + facet_wrap(~dircat, scales='free',nrow=2) + scale_color_brewer(palette="YlOrRd",name='Direction index')+ggtitle("Speed")+theme(legend.position=c(1,0),legend.justification=c(1,0))
prange <- ggplot(filter(drangem,dir_ix<nplotdir | is.na(dir_ix)), aes(x=-t, y= value)) + geom_line(aes(group=variable,color=as.factor(dir_ix))) + facet_wrap(~dircat, scales='free',nrow=2) + scale_color_brewer(palette="YlOrRd",name='Direction index')+ggtitle("Range")+theme(legend.position=c(1,0),legend.justification=c(1,0))

ggsave('plotdir_speed.pdf',pspeed,width=9.5,height=5,units='in')
ggsave('plotdir_range.pdf',prange,width=9.5,height=5,units='in')

