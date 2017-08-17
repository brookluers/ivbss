library(ggplot2)
maxID = 2

lowc = "#fff7bc"
highc = "#d95f0e"
fscale = scale_fill_gradient(low = lowc, high=highc,name="Pr(brake | B^t x)")
h <- read.csv("hmap_multi_cd1_cd2_001.txt", header=T)
p_cd <- ggplot(h, aes(x=cd1,y=cd2)) + geom_tile(aes(fill=ypred)) +
  xlab("cov dir 1") + ylab("cov dir2") +
  fscale
p_m <- ggplot(h, aes(x=meandir,y=cd1)) + geom_tile(aes(fill=ypred)) + 
    xlab("mean direction") + ylab("cov dir 1") +
    fscale

for (curID in 1:maxID){
    h <- read.csv(sprintf("hmap_multi_cd1_cd2_%03d.txt", curID), header=T)
    ggsave(sprintf("hmap_cd1_cd2_%03d.pdf", curID), 
    		p_cd %+% h +
		ggtitle(paste("Driver", curID)),
		width = 6, height=5.5, units='in')
    h <- read.csv(sprintf("hmap_multi_meandir_cd1_%03d.txt", curID), header=T)
    ggsave(sprintf("hmap_meandir_cd1_%03d.pdf", curID),
    	    p_m %+% h +
	    ggtitle(paste("Driver", curID)),
	    width = 6, height=5.5, units='in')
}
