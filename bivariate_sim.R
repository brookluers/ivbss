library(ggplot2)
set.seed(1)
m0 <- 0 # x1 coordinate mean for y==1
m1 <- 0 # x1 coordinate mean for y==0
rho0 <- 0.6 # y==0 correlation of 0.6
rho1 <- 0 # y==1 correlation of 0
theta0 <- acos(rho0)
theta1 <- acos(rho1)
n <- 10000
xneg <- rnorm(n, mean=m0, sd=1) # x1 coordinate for y==0
xtilde_neg <- rnorm(n)
xpos <- rnorm(n, mean=m1, sd=1) # x1 coordinate for y==1
xtilde_pos <- rnorm(n)
xtilde_pos <- xtilde_pos / sd(xtilde_pos)
Xneg <- cbind(xneg, xtilde_neg)
Xpos <- cbind(xpos, xtilde_pos)
Xneg_c <- scale(Xneg, center=TRUE, scale=FALSE)
Xpos_c <- scale(Xpos, center=TRUE, scale=FALSE)
Id <- diag(n) # identity matrix
Pneg <- tcrossprod(qr.Q(qr(Xneg_c[,1,drop=FALSE])))
X2neg_c <- cbind(Xneg_c[,1], (Id - Pneg) %*% Xneg_c[, 2])
X2neg_c <- X2neg_c %*% diag(1/sqrt(colSums(X2neg_c^2)))
xneg_cor <- X2neg_c[,2] + (1 / tan(theta0)) * X2neg_c[,1]
Ppos <- tcrossprod(qr.Q(qr(Xpos_c[,1,drop=FALSE])))
X2pos_c <- cbind(Xpos_c[,1], (Id - Ppos) %*% Xpos_c[,2])
X2pos_c <- X2pos_c %*% diag(1/sqrt(colSums(X2pos_c^2)))
xpos_cor <- X2pos_c[,2] + (1 / tan(theta1)) * X2pos_c[,1]
xneg <- xneg / sd(xneg)
xneg_cor <- xneg_cor / sd(xneg_cor)
xpos_cor <- xpos_cor / sd(xpos_cor)
xpos <- xpos / sd(xpos)


d <- data.frame(rbind(cbind(xpos, xpos_cor, 1),
	cbind(xneg, xneg_cor, 0)))
names(d) <- c('x1','x2','y')
d[d$y==1,]$x1 <- d[d$y==1,]$x1 - mean(d[d$y==1,]$x1)
d[d$y==1,]$x2 <- d[d$y==1,]$x2 - mean(d[d$y==1,]$x2)
d[d$y==0,]$x1 <- d[d$y==0,]$x1 - mean(d[d$y==0,]$x1)
d[d$y==0,]$x2 <- d[d$y==0,]$x2 - mean(d[d$y==0,]$x2)

cat("Mean vector for Y==1:\n")
cat(apply(subset(d, y==1)[,1:2],2,mean), sep=', ')
cat('\n\n')
cat("Mean vector for Y==0:\n")
cat(apply(subset(d, y==0)[,1:2],2,mean), sep=', ')
cat('\n\n')
cat("Sample covariance for Y==1:\n")
cat('[')
cat(cov(subset(d,y==1)[,1:2])[1:2],sep='\t')
cat('\n')
cat(cov(subset(d,y==1)[,1:2])[3:4],sep='\t')
cat("]")
cat('\n\n')
cat("Sample covariance for Y==0:\n")
cat("[")
cat(cov(subset(d,y==0)[,1:2])[1:2],sep='\t')
cat('\n')
cat(cov(subset(d,y==0)[,1:2])[3:4],sep='\t')
cat("]")
cat("\n\n")
cat("Marginal sample covariance:\n")
cat('[')
cat(cov(d[,1:2])[1:2], sep='\t')
cat('\n')
cat(cov(d[,1:2])[3:4], sep='\t')
cat("]")
cat('\n')
write.table(d, "bivariate_doc.txt", sep=',', quote=F,row.names=F)
ggsave('plotpoints.pdf',ggplot(d) + geom_point(aes(x=x1,y=x2,color=as.factor(y)),alpha=I(1/3)) + guides(color='none'),width=5,height=5,units='in')
