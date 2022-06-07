cmd_args=commandArgs()
INPUT<-cmd_args[3]

filecluster<-read.table(INPUT,header=T,sep=" ",dec=".")
beta<-as.matrix(filecluster[,2:3])
gamma<-as.matrix(filecluster[,4])

D <- diag(length(gamma))
S <- t(beta)%*%solve(D, beta) 
alpha <- solve(S, t(beta) %*% solve(D, gamma))
alpha<-as.vector(alpha)

round(alpha[1],digits=5)
