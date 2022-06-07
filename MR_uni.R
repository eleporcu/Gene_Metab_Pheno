cmd_args=commandArgs()
INPUT<-cmd_args[3]

mr_TwoSampleMR <- function(zx, zy, se_zy){
  res <- summary(lm(zy ~ -1 + zx, weights = 1/se_zy^2))
  beta <- res$coef["zx","Estimate"]
  se <- res$coef["zx","Std. Error"]/min(1,res$sigma) #sigma is the residual standard error
  return (c(beta, se))
}


out<-c("gene","alpha_ORIGINAL","SE_ORIGINAL","P_ORIGINAL","N_ORIGINAL","P_het_ORIGINAL","alpha","SE","P","N","P_het","N_outlier")

filecluster<-read.table(INPUT,header=T,sep=" ",dec=".")
beta<-as.vector(filecluster[,2])
gamma<-as.vector(filecluster[,(length(filecluster[1,])-1)])
N_gwas<-filecluster[,length(filecluster[1,])]
N_gwas<-as.vector(N_gwas)
N_eQTLs<-as.vector(filecluster[,3])
zrev = (abs(beta) - abs(gamma)) /sqrt((1/filecluster[,3]) + (1/N_gwas) )
exclude<-which(zrev < -2)

if(length(exclude)>0) {
beta<- beta[-exclude]
gamma<- gamma[-exclude]
N_gwas<-N_gwas[-exclude]
N_eQTLs<-N_eQTLs[-exclude]
}

if (length(beta)==1) {print("ONLY 1 SNP in the model. Exit.");q();}

res.mr_two = mr_TwoSampleMR(beta, gamma, 1/sqrt(N_gwas))
alpha = res.mr_two[1]
se = res.mr_two[2]
p = 2*pnorm(abs(alpha/se), lower.tail = F)

alphaORIGINAL<-alpha
seORIGINAL<-se
pvalORIGINAL<-p

Nstart=length(gamma)
d<-gamma - alpha[1]*beta
var_d_vec<- (1/N_gwas)+se^2*beta^2+(1/N_eQTLs)*as.vector(alpha^2)+(1/N_eQTLs)*as.vector(se^2)
var_d<-diag(as.vector(var_d_vec))


z<-t(as.matrix(d))%*%solve(as.matrix(var_d))%*%as.matrix(d)
N<-length(d)
p_hetORIGINAL<- 1-pchisq(z,N-1)

d<-abs(d)
outlier<-which(d==max(d))[1]
p_het<-p_hetORIGINAL

while (p_het<0.05 & N>3) {

gamma<-as.vector(gamma[-outlier])
beta<-as.vector(beta[-outlier])
N_gwas<-as.vector(N_gwas[-outlier])
N_eQTLs<-as.vector(N_eQTLs[-outlier])
N=length(beta)

res.mr_two = mr_TwoSampleMR(beta, gamma, 1/sqrt(N_gwas))
alpha = res.mr_two[1]
se = res.mr_two[2]
p = 2*pnorm(abs(alpha/se), lower.tail = F)

d<-gamma - alpha[1]*beta
var_d_vec<- (1/N_gwas)+se^2*beta^2+(1/N_eQTLs)*as.vector(alpha^2)+(1/N_eQTLs)*as.vector(se^2)
var_d<-diag(as.vector(var_d_vec))


z<-t(as.matrix(d))%*%solve(as.matrix(var_d))%*%as.matrix(d)
N<-length(d)
p_het<- 1-pchisq(z,N-1)
d<-abs(d)
outlier<-which(d==max(d))[1]

}  #end while

line<-c(colnames(filecluster)[2],signif(alphaORIGINAL,5),signif(seORIGINAL,5),signif(pvalORIGINAL,5),Nstart,signif(p_hetORIGINAL,5),signif(alpha,5),signif(se,5),signif(p,5),N,signif(p_het,5),Nstart-N)
out<-rbind(out,line)

write.table(out,file=paste(colnames(filecluster)[2],".alpha",sep=""),quote=F,col.names=F,row.names=F)


warnings()
