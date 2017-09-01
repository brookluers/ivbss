
library(ggplot2)
ndir = 1
na <- c("meandir0",
  paste("cd", 0:(ndir-1), sep=''),
  paste("pc", 0:(ndir-1), sep=''))

nneighbor <- 100
w <- nneighbor %/% 2
maxId <- 2
ycounts <- matrix(nrow=maxId,ncol=3)
ycounts[,1] <- 1:maxId
colnames(ycounts)<-c('Driver','ny1','ny0')
for (cid in 1:maxId){
    s <- read.csv(sprintf("smproj_multi_%03d.txt", cid), header=T)
    final <- data.frame(Driver=rep(cid, nrow(s)))
    for (cn in na){
    	temp <- s[order(s[,cn]), c(cn, "Brake")]
    	tnrow <- nrow(temp)
    	temp <- cbind(temp, phat = do.call('c',
	     lapply(sapply(1:tnrow, function(x) return(max(1,x - w):min(tnrow, x + w))),
    	     function(ix) return(mean(temp[ix,"Brake"])))))
    	final <- setNames(cbind(final, temp[,c(cn,'phat')]),
	      	 	 c(names(final), cn, paste('phat_',cn,sep='')))
    }
    ycounts[cid,'ny1'] <- sum(s$Brake)
    ycounts[cid,'ny0'] <- nrow(s) - ycounts[cid,'ny1']
    write.table(final, file=paste('phat_driver',cid,'.txt',sep=''),
    		quote=F,sep=',',row.names=FALSE)
}
write.table(ycounts, file='ycounts.txt', quote=F,sep=',',row.names=FALSE)
