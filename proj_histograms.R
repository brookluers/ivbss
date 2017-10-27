library(dplyr)
library(ggplot2)

d <- read.csv('proj_histograms.txt', header=T)

maxID <- max(d$Driver)

mt <- theme(panel.background=element_rect(fill='white',color='white'),
      panel.grid.major.y=element_line(color='lightgrey',linetype='dotted'),
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      strip.background=element_rect(fill='white',color='white'))

# big overlapping histogram

ggsave("fig/proj_hist_global.png",
	ggplot(d) + 
	  geom_ribbon(aes(x=bin_lwr, ymin=0, ymax=hdens, group=Driver),alpha=I(1/10)) +
	  facet_wrap(~Dir, scales='free') + mt + 
	  ggtitle("All drivers") + xlab("projected value") + ylab("density (within driver)"),
	width=9.5, height=4, units='in', dpi=300)

# Histograms for each driver
for (i in 1:maxID) {
    ggsave(sprintf('fig/proj_hist_%03d.png', i),
    	    ggplot(filter(d, Driver==i)) + geom_ribbon(aes(x=bin_lwr, ymin=0, ymax=hdens)) + 
	    	ggtitle(sprintf("Driver %d", i)) +  facet_wrap(~Dir,scales='free')+ mt + 
		xlab("projected value") + ylab("density"),
    	   width=9, height=3.5, units='in', dpi=300)
}