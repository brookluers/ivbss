s <- read.csv('summary.txt',header=T)
s <- subset(s, Trip==1)[,c('Driver','TODTripStart')]
s$trip1starttime <- s$TODTripStart
s$TODTripStart <- NULL
write.csv(s, file='summary_trip1_starttime.txt',quote=FALSE,row.names=FALSE)
