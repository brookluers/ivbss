library(ggplot2)
h <- read.csv("hmap_001.txt", header=T)
p <- ggplot(h, aes(x=meandir,y=cd1)) + geom_tile(aes(fill=ypred)) +
  xlab("mean direction") + ylab("cov. dir1")
ggsave("hmap_001.pdf", p, width=6,height=5.5,units='in')