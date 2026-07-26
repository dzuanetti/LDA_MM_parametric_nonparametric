#
##################
# Avalia o número de grupos
##################
#
#caminho<-"/Users/Daiane/Downloads/rstudio-export_Edvaldo/" - menos simulações
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/BBC_sport_DP/"
#
Sj<-scan(file=paste(caminho,"Sj_BBC_sport.txt",sep=""))
alpha<-scan(file=paste(caminho,"alpha_BBC_sport.txt",sep=""))
#
iter<-5000
#
library(coda)
plot(alpha,type='l')
log_v<-mcmc(alpha)
effectiveSize(log_v)
geweke.diag(log_v)
mean(alpha)
#
m<-737
S<-matrix(Sj,ncol=m,nrow=iter,byrow=TRUE)
#
K<-NULL
for (i in 1:nrow(S)) K[i]<-length(table(S[i,]))
#
k5<-which(K==5)
#
Sj.j<-S #matriz de agrupamento a posteriori
prob.eq<-matrix(0,nrow=ncol(Sj.j),ncol=ncol(Sj.j))
library(compiler)
enableJIT(3)
for (i in 1:ncol(Sj.j)){
	for (j in 1:ncol(Sj.j)){
		prob.eq[i,j]<-round(sum(Sj.j[,i]==Sj.j[,j])/length(Sj.j[,i]),5)*100}}
#
thresh<-0.50*100 # definindo os grupos finais
clust_f<-c(1,rep(0,(ncol(Sj.j)-1)))
for (i in 2:ncol(Sj.j)){
#for (i in 310:514){
	if (max(prob.eq[i,1:(i-1)])>thresh) clust_f[i]<-clust_f[which(prob.eq[i,1:(i-1)]==max(prob.eq[i,1:(i-1)]))[1]] else clust_f[i]<-max(clust_f[1:(i-1)]+1)}
#
table(S[iter,],clust_f)
#
thesing<-0.3*100 # juntando outliers que aparecem pelo menos 30% das vezes juntos
singl<-which(clust_f %in% which(table(clust_f)==1))
prob.eq.sin<-matrix(prob.eq[singl,],ncol=ncol(prob.eq))
for (i in 1:nrow(prob.eq.sin)){
	prob.eq.sin[i,singl[i]]<-0
	if (max(prob.eq.sin[i,])>thesing) clust_f[singl[i]]<-clust_f[which(prob.eq.sin[i,]==max(prob.eq.sin[i,]))[1]]}
while (length(table(clust_f))<max(clust_f)){ # exclude empty clusters
	categr<-as.numeric(as.character(data.frame(table(clust_f))[,1]))
	categd<-seq(1:length(table(clust_f)))
	dif<-which(categr!=categd)
	clust_f[which(clust_f>dif[1])]<-clust_f[which(clust_f>dif[1])]-1}
#
caminho<-"/Users/Daiane/Downloads/bbcsport/"
#
S_verd<-scan(file=paste(caminho,"bbcsport_label.txt",sep=""))
table(S_verd,clust_f)
#
cat('',clust_f,file=paste(caminho,"Grupos_MM_DP.txt",sep=""))
#
############ Calculando a verossimilhança marginal
#
caminho<-"/Users/Daiane/Downloads/bbcsport/"
dados<-scan(paste(caminho,"bbcsport.txt",sep=""))
y_ap<-matrix(dados,nrow=737,ncol=4613,byrow=TRUE)
summary(apply(y_ap,1,sum))
#
library(MGLM)
beta0<-rep(1,ncol(y_ap))
somaWs<-NULL
for (k in 1:5) somaWs<-rbind(somaWs,apply(y_ap[clust_f==k,],2,sum))
logmargvero<-0
for (i in 1:nrow(y_ap)) logmargvero<-logmargvero+ddirmn(y_ap[i,], beta0+somaWs[clust_f[i],])

