library(ggplot2)
maxID = 24


cat("\n")
cat("-------R: Plotting heatmaps-------\n")
lowc = "#fff7bc"
highc = "#d95f0e"
fscale = scale_fill_gradient(low = lowc, high=highc,name="Pr(brake | B^t x)")

p_cd <- ggplot()  + geom_tile(aes(x=cd1,y=cd2, fill=ypred)) +
  xlab("cov dir 1") + ylab("cov dir2") #+
  #fscale
p_m <- ggplot() + geom_tile(aes(x=meandir,y=cd1,fill=ypred)) +
    xlab("mean direction") + ylab("cov dir 1") #+
    #fscale

for (curID in 1:maxID){
    h <- read.csv(sprintf("/scratch/stats_flux/luers/hmap_pca_cd1_cd2_%03d.txt", curID), header=T)
    ggsave(sprintf("hmap_cd1_cd2_pca_%03d.pdf", curID), 
    		p_cd %+% h +
		ggtitle(paste("Driver", curID)),
		width = 6, height=5.5, units='in')
    h <- read.csv(sprintf("/scratch/stats_flux/luers/hmap_pca_meandir_cd1_%03d.txt", curID), header=T)
    ggsave(sprintf("hmap_meandir_cd1_pca_%03d.pdf", curID),
    	    p_m %+% h +
	    ggtitle(paste("Driver", curID)),
	    width = 6, height=5.5, units='in')
}
