library(ggplot2)
library(tibble)
library(tidyr)
library(readr)
library(dplyr)
library(tidytext)
library(topicmodels)
library(compiler)

################# sports data set #################
#caminho<-"/Users/Daiane/Downloads/bbcsport/"
caminho<-"/home/mariella/Edvaldo/"
dados<-scan(paste(caminho,"bbcsport.txt",sep=""))
#
dados_ap<-matrix(dados,nrow=737,ncol=4613,byrow=TRUE)
#
summary(apply(dados_ap,1,sum))
dados_ap = data.frame(dados_ap)
#
### correlation matriz among words
bow_ap <- matrix(unlist(dados_ap), nrow = nrow(dados_ap))
coorre_ap<-matrix(0,ncol(bow_ap),ncol(bow_ap))
#
enableJIT(3)
for (i in 1:nrow(coorre_ap)){
	for (j in i:ncol(coorre_ap)){
		coorre_ap[i,j]<-sum(bow_ap[,i]>0 & bow_ap[,j]>0)+0.001
		if (j>i) coorre_ap[j,i]<-coorre_ap[i,j]}}
probs_ocor_ap<-coorre_ap/nrow(bow_ap)
#
################# LDA with a fixed number of topics #################
#
### Fitting the model and several criteria
start_time <- Sys.time()
set.seed(100)
enableJIT(3)
for (K_2 in c(2:20,25,30,35,40,45,50,75,100,125,150,175,200,250,300,400,500)){
	cat('\n',K_2)
	ap_lda2_ap <- LDA(dados_ap, k = K_2, method="Gibbs")
	#	
	### Selection criteria 
	#
	#Likelihood
	vero2_lda <- ap_lda2_ap@loglikelihood  
	#
	#Perplexity
	perp2_lda <- exp(-vero2_lda/sum(dados_ap))
	#
	#BIC
	BIC <- vero2_lda-((nrow(dados_ap)+ncol(dados_ap)-1)*K_2-nrow(dados_ap))/2*log(sum(dados_ap))
	#
	#auxiliar adjust to calculate the other criteria
	### TOP 10 TERMS in each topic
	ap_top_terms2 <- terms(ap_lda2_ap, 10)
	ap_top_terms_2 = t(ap_top_terms2)
	ap_top_terms_2_v2 = matrix(unlist(ap_top_terms_2), nrow = nrow(ap_top_terms_2))
	ap_top_terms_2_v3<-matrix(0,nrow(ap_top_terms_2_v2),ncol(ap_top_terms_2_v2))
	for (i in 1:K_2){
		for (j in 1:10){
			ap_top_terms_2_v3[i,j] = which(colnames(dados_ap) == ap_top_terms_2_v2[i,j])}}
	pa_top_terms2 = ap_top_terms_2_v3
	top_termos<-pa_top_terms2
	soma1<-0
	soma2<-0
	soma3<-0
	for (topicos in 1:K_2){
		for (i in 1:(length(top_termos[topicos,])-1)){
			for (j in (i+1):length(top_termos[topicos,])){
        		soma1<-soma1+log((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]])/((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,i]])*(probs_ocor_ap[top_termos[topicos,j],top_termos[topicos,j]])))                                                                                                                                                                          
				soma2<-soma2+(log((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]])/((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,i]])*(probs_ocor_ap[top_termos[topicos,j],top_termos[topicos,j]])))/-log(probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]]))
        		soma3<-soma3+log((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]])/((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,i]])))}}}
	pmi_lda<-soma1/(K_2*45)
	npmi_lda<-soma2/(K_2*45)
	lcp_lda<-soma3/(K_2*45)
#
	cat('',vero2_lda,file=paste(caminho,"vero_LDA_sport.txt",sep=""),append=T)
	cat('',perp2_lda,file=paste(caminho,"perp_LDA_sport.txt",sep=""),append=T)
	cat('',BIC,file=paste(caminho,"BIC_LDA_sport.txt",sep=""),append=T)	
	cat('',pmi_lda,file=paste(caminho,"pmi_LDA_sport.txt",sep=""),append=T)	
	cat('',npmi_lda,file=paste(caminho,"npmi_LDA_sport.txt",sep=""),append=T)	
	cat('',lcp_lda,file=paste(caminho,"lcp_LDA_sport.txt",sep=""),append=T)}
#
##### calculating sBIC
#
learning_coefficient <- function(n_docs, n_words, super_model, sub_model) {
  # Compute learning coefficient and its multiplicity
  # using formulas from Hayashi (2021)
  #
  # N = n_docs: number of documents
  # M = n_words: vocabulary size
  # H = super_model: number of topics in a model of a higher order
  # r = sub_model: number of topics in a model of a lower order
  #
  # Returns:
  # learn_coeff = λ_Hr
  # m = multiplicity of λ_Hr

  N <- n_docs
  M <- n_words
  H <- super_model
  r <- sub_model

  if ((N + r) <= (M + H) && (M + r) <= (N + H) && (H + r) <= (M + H)) {
    if ((M + N + H + r - 1) %% 2 != 0) {
      learn_coeff <- (2 * (H + r) * (M + N) - (M - N)^2 - (H + r)^2) / 8 - N / 2
      m <- 1} else {
      learn_coeff <- (2 * (H + r) * (M + N) - (M - N)^2 - (H + r)^2 + 1) / 8 - N / 2
      m <- 2}} else if ((M + H) < (N + r)) {
    learn_coeff <- (M * H + N * r - H * r - N) / 2
    m <- 1} else if ((N + H) < (M + r)) {
    learn_coeff <- (N * H + M * r - H * r - N) / 2
    m <- 1} else {
	learn_coeff <- (M * N - N) / 2
    m <- 1}
  return(c(learn_coeff,m))
}
#
#### calculating log(LKk)
#
log_LKk<-function(loglikKuper, paramKuper, Ntot) {
	logLKk<-loglikKuper-paramKuper[1,]*log(Ntot)+(paramKuper[2,]-1)*log(log(Ntot))
	return(logLKk)}
#
#######
#
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/BBC_sport_DP/LDA_param_sport/"
topic_sequence <- c(2:20,25,30,35,40,45,50,75,100,125,150,175,200,250,300,400,500)
N <- nrow(dados_ap) # Number of documents
V <- ncol(dados_ap) # Vocabulary size
log_liks<-scan(file=paste(caminho,"vero_LDA_sport.txt",sep=""))
#
# https://github.com/VikaNa/sBIC
#
###
#
param_sBIC<-list()
enableJIT(3)
for (k in 1:length(topic_sequence)){
	lambda<-NULL
	m<-NULL
	for (j in 1:k){
		param<-learning_coefficient(N, V, topic_sequence[k], topic_sequence[j])
		lambda[j]<-param[1]
		m[j]<-param[2]}
	param_sBIC[[k]]<-rbind(lambda,m)}
#
Ntot<-sum(dados_ap)
#
logLKk<-list()
for (i in 1:length(topic_sequence)) logLKk[[i]]<-log_LKk(log_liks[i], param_sBIC[[i]], Ntot)
#
mn<-max(sapply(logLKk,max))
#
library(gmp)
library(Rmpfr)
#
LKk<-list()
for (i in 1:length(topic_sequence)){
	Li<-mpfr(rep(0,length(logLKk[[i]])), precBits = 128)
	for (j in 1:length(logLKk[[i]])){
    	val<-mpfr(logLKk[[i]][j] - mn, precBits = 128)
    	Li[j] <- exp(val)}
    LKk[[i]]<-Li}
#
L<-list(LKk[[1]][1])
sBIC<-c(log(LKk[[1]][1]) + mn)
#
enableJIT(3)
for (i in 2:length(topic_sequence)){
	b<-Reduce(`+`,L[1:(i-1)])-LKk[[i]][i]
	c<-mpfr(0, precBits = 128)
	for (j in 1:(i-1)) c<-c-LKk[[i]][j]*L[[j]]
    prec <- 50
    mpfr_prec <- prec * 3.321928  # bits from decimal digits approx
    temp <- mpfr(0, precBits = mpfr_prec) # Increase precision if needed
    repeat{
    	b_mpfr <- as(mpfr(b,precBits = mpfr_prec), "mpfr")
    	c_mpfr <- as(mpfr(c,precBits = mpfr_prec), "mpfr")
    	discriminant <- b_mpfr^2 - 4 * c_mpfr
		if (discriminant < 0) stop("Negative discriminant encountered.")
		temp<-(-b_mpfr+sqrt(discriminant))/2
		val_ln <- log(max(temp, mpfr(0, precBits = mpfr_prec))) + mn
		if (!is.infinite(val_ln)) break
		prec <- prec + 50
		mpfr_prec <- prec * 3.321928}
    L[[i]] <- temp
    sBIC[i] <- log(max(L[[i]], mpfr(0, precBits = mpfr_prec))) + mn}
cat('',as.numeric(sBIC),file=paste(caminho,"sBIC_LDA_sport.txt",sep=""))
cat('',mn,file=paste(caminho,"mn_LDA_sport.txt",sep=""))
#
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)
#
#### analize the results
#
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/BBC_sport_DP/LDA_param_sport/"
log_liks<-scan(file=paste(caminho,"vero_LDA_sport.txt",sep=""))
perp<-scan(file=paste(caminho,"perp_LDA_sport.txt",sep=""))
pmi<-scan(file=paste(caminho,"pmi_LDA_sport.txt",sep=""))
npmi<-scan(file=paste(caminho,"npmi_LDA_sport.txt",sep=""))
lcp<-scan(file=paste(caminho,"lcp_LDA_sport.txt",sep=""))
#BIC<-scan(file=paste(caminho,"BIC_LDA_sport.txt",sep="")), esse está errado, vou calcular pela log
BIC<-log_liks-((nrow(dados_ap)+ncol(dados_ap)-1)*topic_sequence-nrow(dados_ap))/2*log(sum(dados_ap))
sBIC2<-scan(file=paste(caminho,"sBIC_LDA_sport.txt",sep=""))
#
plot(topic_sequence,log_liks,type='l')
K<-topic_sequence[which(log_liks==max(log_liks))]
K
plot(topic_sequence,perp,type='l')
K<-topic_sequence[which(perp==min(perp))]
K
plot(topic_sequence,pmi,type='l')
K<-topic_sequence[which(pmi==max(pmi))]
K
plot(topic_sequence,npmi,type='l')
K<-topic_sequence[which(npmi==max(npmi))]
K
plot(topic_sequence,lcp,type='l')
K<-topic_sequence[which(lcp==max(lcp))]
K
plot(topic_sequence,BIC,type='l')
K<-topic_sequence[which(BIC==max(BIC))]
K
plot(topic_sequence,sBIC2,type='l')
K<-topic_sequence[which(sBIC2==max(sBIC2))]
K
