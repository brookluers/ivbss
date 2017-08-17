s <- read.csv('/nfs/turbo/ivbss/LvFot/summary.txt',header=T)
s <- subset(s, Trip==1)[,c('Driver','TODTripStart')]
s$trip1starttime <- s$TODTripStart
s$TODTripStart <- NULL
write.csv(s, file='/scratch/stats_flux/luers/summary_trip1_starttime.txt',
	     quote=FALSE,row.names=FALSE)
