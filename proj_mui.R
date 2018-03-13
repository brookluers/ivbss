library(Matrix)
library(tidyverse)

dirdf <- read_csv("directions_mui.txt")
ndir_proj <- 5
dirmat <- Matrix(select(dirdf, paste('doc_mui_', 1:ndir_proj, sep='')) %>% as.matrix)
vnames <- dirdf$varname
i <- 1
idstr <- sprintf("%03d", i)
lagdat <- read_csv(paste("data/lagdat_small_", idstr, ".txt",sep=''))
Xmat <- Matrix(select(lagdat, vnames) %>% as.matrix)
write_csv(bind_cols(select(lagdat, Driver, Trip, Time, Brake_1sec),
          as_tibble(as.matrix(Xmat %*% dirmat))),
          paste('projdat_doc_mui_', idstr, '.txt', sep=''))
