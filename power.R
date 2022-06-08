NmQTL<-8000
Ngwas<-300000
NeQTL<-32000
total_effect<- 0.035
NSNPs_eQTL<-6
h2_GE<-0.064
NSNPs_mQTL<-5
h2_ME<-0.043

mr_TwoSampleMR <- function(zx, zy, se_zy){
  res <- summary(lm(zy ~ -1 + zx, weights = 1/se_zy^2))
  beta <- res$coef["zx","Estimate"]
  se <- res$coef["zx","Std. Error"]/min(1,res$sigma) #sigma is the residual standard error
  return (c(beta, se))
}

direct<-seq(0.01,0.99,0.02)
proportion<-seq(log10(10),log10(0.1),-0.04) 

for (j in 1:length(direct)) {
direct_eff<-direct[j]
power<-c()
powerTWMR<-c()
for (k in 1:length(proportion)) {
sigma<-10^proportion[k]
countTWMR<-0
countGENETICMAP<-0
countMetabWMR<-0

for (l in 1:500) {
beta_eQTLs<-rnorm(NSNPs_eQTL,0,sqrt(h2_GE/NSNPs_eQTL))+rnorm(NSNPs_eQTL,0,sqrt(1/NeQTL))
scaled_beta_eQTLs<- sqrt(h2_GE)*beta_eQTLs/sqrt(sum(beta_eQTLs^2))
beta_eQTLs<-scaled_beta_eQTLs
beta_GWAS<- beta_eQTLs*total_effect+rnorm(NSNPs_eQTL,0,sqrt(1/Ngwas))

res.mr_two = mr_TwoSampleMR(beta_eQTLs, beta_GWAS, rep(1/sqrt(Ngwas),NSNPs_eQTL))
alpha = res.mr_two[1]
se = res.mr_two[2]
pval = 2*pnorm(abs(alpha/se), lower.tail = F)

if(pval<0.05/(500)) {countTWMR=countTWMR+1}
for (m in 1:80) {


beta_mQTLs<-beta_eQTLs*sigma*sqrt((1-direct_eff)*total_effect/sigma)+rnorm(NSNPs_eQTL,0,sqrt(1/NmQTL))
res.mr_two = mr_TwoSampleMR(beta_eQTLs, beta_mQTLs, rep(1/sqrt(NmQTL),NSNPs_eQTL))
alpha = res.mr_two[1]
se = res.mr_two[2]
pval = 2*pnorm(abs(alpha/se), lower.tail = F)

if(pval<0.05/(500*80)) {countGENETICMAP=countGENETICMAP+1}
}
}
for (n in 1:80) {

beta_mQTLs<-rnorm(NSNPs_mQTL,0,sqrt(h2_ME/NSNPs_mQTL))+rnorm(NSNPs_mQTL,0,sqrt(1/NmQTL))
scaled_beta_mQTLs<- sqrt(h2_ME)*beta_mQTLs/sqrt(sum(beta_mQTLs^2))
beta_mQTLs<-scaled_beta_mQTLs

beta_GWAS<- beta_mQTLs*sqrt((1-direct_eff)*total_effect/sigma)+rnorm(NSNPs_mQTL,0,sqrt(1/Ngwas))
res.mr_two = mr_TwoSampleMR(beta_mQTLs, beta_GWAS, rep(1/sqrt(Ngwas),NSNPs_mQTL))
alpha = res.mr_two[1]
se = res.mr_two[2]
pval = 2*pnorm(abs(alpha/se), lower.tail = F)

if(pval<0.05/80) {countMetabWMR=countMetabWMR+1}
}
power<-c(power,countTWMR/500-(countGENETICMAP/(500*80))*(countMetabWMR/80))

}


if (j==1) {powermatrix<-power}
else{powermatrix<-data.frame(powermatrix,power)}
print(j)
}

 
write.table(powermatrix,file="powermatrix.txt",sep=" ",dec=".",col.names=F,row.names=F,quote=F)


