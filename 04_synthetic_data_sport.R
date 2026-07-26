#
##### Simula dados na estrutura do bbcsport
#
caminho<-"/Users/Daiane/Downloads/bbcsport/"
#caminho<-"/home/mariella2/Edvaldo/"
dados<-scan(paste(caminho,"bbcsport.txt",sep=""))
y_verd<-matrix(dados,nrow=737,ncol=4613,byrow=TRUE)
summary(apply(y_verd,1,sum))
S_verd<-scan(file=paste(caminho,"bbcsport_label.txt",sep=""))
#
#library(lsa)
retirar<-which(apply(y_verd,2,sum)<10)
y_verd<-y_verd[,-retirar] # tenho que tirar se não da palavra com frequência zero e não roda o NMF e talvez outros
phi<-NULL
for (k in 1:5) phi<-rbind(phi,apply(y_verd[S_verd==k-1,],2,sum)/sum(y_verd[S_verd==k-1]))
apply(phi,1,sum)
#cos<-NULL
#for (k in 1:4){
#	for (j in (k+1):5){
#		cos<-rbind(cos, c(k,j,cosine(phi[k,], phi[j,])))}} # todos são razoavelmente diferentes
#
dim(phi)
#
##################################### NMF
#
library(NMF)
library(compiler)
topic_sequence<-c(2:20,25,30)
m<-nrow(y_verd)
#
replicas<-20
set.seed(100)
for (rep in 1:replicas){
dados<-matrix(0,nrow=m,ncol=ncol(y_verd))
for (i in 1:m) dados[i,]<-rmultinom(1,size=sum(y_verd[i,]),phi[S_verd[i]+1,])
#
y_ap=as.matrix(dados)
filtra<-apply(y_ap,2,sum)
if (sum(filtra==0)>0) y_ap<-y_ap[,-which(filtra==0)]
#
### correlation matriz among words
#
bow_ap <- matrix(unlist(y_ap), nrow = nrow(y_ap))
coorre_ap<-matrix(0,ncol(bow_ap),ncol(bow_ap))
#
enableJIT(3)
for (i in 1:nrow(coorre_ap)){
	for (j in i:ncol(coorre_ap)){
		coorre_ap[i,j]<-sum(bow_ap[,i]>0 & bow_ap[,j]>0)+0.001
		if (j>i) coorre_ap[j,i]<-coorre_ap[i,j]}}
probs_ocor_ap<-coorre_ap/nrow(bow_ap)
#
topicof<-NULL
enableJIT(3)
for (K_2 in topic_sequence){
	cat('\n', K_2, rep)
#	set.seed(100)
	nmf_model <- nmf(y_ap, rank = K_2, method = "brunet", nrun = 3)
	#
	### select the most frequent words for each topic to selection criteria
	#
	basis_matrix <- basis(nmf_model)
	topico<-NULL
	for (i in 1:nrow(basis_matrix)) topico[i]<-which(basis_matrix[i,]==max(basis_matrix[i,]))
	topicof<-rbind(topicof,topico)
	#
	posicoes_2 <- matrix(0,nrow = K_2, ncol = 10)
	for (k in 1:K_2){
		posic<-apply(matrix(y_ap[topico==k,],ncol=ncol(y_ap)),2,sum)
		posicoes_2[k,]<-order(posic, decreasing = TRUE)[1:10]}
	soma1_2<-0
	soma2_2<-0
	soma3_2<-0
	for (topicos in 1:K_2){
		for (i in 1:(length(posicoes_2[topicos,])-1)){
			for (j in (i+1):length(posicoes_2[topicos,])){
				soma1_2<-soma1_2+log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])*(probs_ocor_ap[posicoes_2[topicos,j],posicoes_2[topicos,j]])))                                                                                                                                                                          
				soma2_2<-soma2_2+(log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])*(probs_ocor_ap[posicoes_2[topicos,j],posicoes_2[topicos,j]])))/-log(probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]]))
				soma3_2 <-soma3_2+log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])))}}}
	pmi_nmf_2<-soma1_2/(K_2*45)
	npmi_nmf_2<-soma2_2/(K_2*45)
	lcp_nmf_2<-soma3_2/(K_2*45)
	cat('',pmi_nmf_2,file=paste(caminho,"pmi_NMF_dados1_",rep,".txt",sep=""),append=T)	
	cat('',npmi_nmf_2,file=paste(caminho,"npmi_NMF_dados1_",rep,".txt",sep=""),append=T)	
	cat('',lcp_nmf_2,file=paste(caminho,"lcp_NMF_dados1_",rep,".txt",sep=""),append=T)}
#
pmi<-scan(file=paste(caminho,"pmi_NMF_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(pmi==max(pmi))]
cat('',K,file=paste(caminho,"K_pmi_NMF_dados1.txt",sep=""),append=T)
cat('',topicof[topic_sequence==K,],file=paste(caminho,"Sj_pmi_NMF_dados1_",rep,".txt",sep=""),append=T)
#
npmi<-scan(file=paste(caminho,"npmi_NMF_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(npmi==max(npmi))]
cat('',K,file=paste(caminho,"K_npmi_NMF_dados1.txt",sep=""),append=T)
cat('',topicof[topic_sequence==K,],file=paste(caminho,"Sj_npmi_NMF_dados1_",rep,".txt",sep=""),append=T)
#
lcp<-scan(file=paste(caminho,"lcp_NMF_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(lcp==max(lcp))]
cat('',K,file=paste(caminho,"K_lcp_NMF_dados1.txt",sep=""),append=T)
cat('',topicof[topic_sequence==K,],file=paste(caminho,"Sj_lcp_NMF_dados1_",rep,".txt",sep=""),append=T)}
#	
################################## run parametric LDA
#
library(dplyr)
library(topicmodels)
library(compiler)
library(gmp)
library(Rmpfr)
#
# calculating sBIC
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
# calculating log(LKk)
#
log_LKk<-function(loglikKuper, paramKuper, Ntot) {
	logLKk<-loglikKuper-paramKuper[1,]*log(Ntot)+(paramKuper[2,]-1)*log(log(Ntot))
	return(logLKk)}
#
topic_sequence<-c(2:20,25,30,35,40,45,50)
#
replicas<-20
set.seed(100)
for (rep in 1:replicas){
dados<-matrix(0,nrow=nrow(y_verd),ncol=ncol(y_verd))
for (i in 1:nrow(y_verd)) dados[i,]<-rmultinom(1,size=sum(y_verd[i,]),phi[S_verd[i]+1,])
#
y_ap=as.matrix(dados)
filtra<-apply(y_ap,2,sum)
if (sum(filtra==0)>0) y_ap<-y_ap[,-which(filtra==0)]
#
### correlation matriz among words
#
bow_ap <- matrix(unlist(y_ap), nrow = nrow(y_ap))
coorre_ap<-matrix(0,ncol(bow_ap),ncol(bow_ap))
#
enableJIT(3)
for (i in 1:nrow(coorre_ap)){
	for (j in i:ncol(coorre_ap)){
		coorre_ap[i,j]<-sum(bow_ap[,i]>0 & bow_ap[,j]>0)+0.001
		if (j>i) coorre_ap[j,i]<-coorre_ap[i,j]}}
probs_ocor_ap<-coorre_ap/nrow(bow_ap)

dados_ap = data.frame(dados)
#set.seed(100)
Sf<-NULL
enableJIT(3)
for (K_2 in topic_sequence){
	cat('\n',rep, K_2)
	ap_lda2_ap <- LDA(dados_ap, k = K_2, method="Gibbs")
	Sj<-NULL
	for (i in 1:nrow(dados_ap)) Sj[i]<-which(posterior(ap_lda2_ap)$topics[i,]==max(posterior(ap_lda2_ap)$topics[i,]))[1]
	Sf<-rbind(Sf,Sj)
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
	cat('',vero2_lda,file=paste(caminho,"vero_LDA_dados1_",rep,".txt",sep=""),append=T)
	cat('',perp2_lda,file=paste(caminho,"perp_LDA_dados1_",rep,".txt",sep=""),append=T)
	cat('',BIC,file=paste(caminho,"BIC_LDA_dados1_",rep,".txt",sep=""),append=T)	
	cat('',pmi_lda,file=paste(caminho,"pmi_LDA_dados1_",rep,".txt",sep=""),append=T)	
	cat('',npmi_lda,file=paste(caminho,"npmi_LDA_dados1_",rep,".txt",sep=""),append=T)	
	cat('',lcp_lda,file=paste(caminho,"lcp_LDA_dados1_",rep,".txt",sep=""),append=T)}
#
N <- nrow(dados_ap) # Number of documents
V <- ncol(dados_ap) # Vocabulary size
log_liks<-scan(file=paste(caminho,"vero_LDA_dados1_",rep,".txt",sep=""))
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
Ntot<-sum(dados_ap)
#
logLKk<-list()
for (i in 1:length(topic_sequence)) logLKk[[i]]<-log_LKk(log_liks[i], param_sBIC[[i]], Ntot)
#
mn<-max(sapply(logLKk,max))
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
cat('',as.numeric(sBIC),file=paste(caminho,"sBIC_LDA_dados1_",rep,".txt",sep=""))
#
pmi<-scan(file=paste(caminho,"pmi_LDA_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(pmi==max(pmi))]
cat('',K,file=paste(caminho,"K_pmi_LDA_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_pmi_LDA_dados1_",rep,".txt",sep=""),append=T)
#
npmi<-scan(file=paste(caminho,"npmi_LDA_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(npmi==max(npmi))]
cat('',K,file=paste(caminho,"K_npmi_LDA_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_npmi_LDA_dados1_",rep,".txt",sep=""),append=T)
#
lcp<-scan(file=paste(caminho,"lcp_LDA_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(lcp==max(lcp))]
cat('',K,file=paste(caminho,"K_lcp_LDA_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_lcp_LDA_dados1_",rep,".txt",sep=""),append=T)
#
vero<-scan(file=paste(caminho,"vero_LDA_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(vero==max(vero))]
cat('',K,file=paste(caminho,"K_vero_LDA_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_vero_LDA_dados1_",rep,".txt",sep=""),append=T)
#
perp<-scan(file=paste(caminho,"perp_LDA_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(perp==min(perp))]
cat('',K,file=paste(caminho,"K_perp_LDA_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_perp_LDA_dados1_",rep,".txt",sep=""),append=T)
#
bic<-scan(file=paste(caminho,"BIC_LDA_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(bic==max(bic))]
cat('',K,file=paste(caminho,"K_bic_LDA_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_bic_LDA_dados1_",rep,".txt",sep=""),append=T)
#
sbic<-scan(file=paste(caminho,"sBIC_LDA_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(sbic==max(sbic))]
cat('',K,file=paste(caminho,"K_sbic_LDA_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_sbic_LDA_dados1_",rep,".txt",sep=""),append=T)	
}		
#
#################################### run MM parametrico
#
library(compiler)
#
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  val}
#
rDiric<-function(gama){
  X<-rgamma(length(gama),gama,1)
  Y<-X/sum(X)
  return(Y)}
#
a0<-1
b0<-1
iter<-250
topic_sequence<-c(2:20)
m<-nrow(y_verd)
#
replicas<-20
set.seed(100)
for (rep in 1:replicas){
dados<-matrix(0,nrow=m,ncol=ncol(y_verd))
for (i in 1:m) dados[i,]<-rmultinom(1,size=sum(y_verd[i,]),phi[S_verd[i]+1,])
#
y_ap=as.matrix(dados)
filtra<-apply(y_ap,2,sum)
if (sum(filtra==0)>0) y_ap<-y_ap[,-which(filtra==0)]
#
### correlation matriz among words
#
bow_ap <- matrix(unlist(y_ap), nrow = nrow(y_ap))
coorre_ap<-matrix(0,ncol(bow_ap),ncol(bow_ap))
#
enableJIT(3)
for (i in 1:nrow(coorre_ap)){
	for (j in i:ncol(coorre_ap)){
		coorre_ap[i,j]<-sum(bow_ap[,i]>0 & bow_ap[,j]>0)+0.001
		if (j>i) coorre_ap[j,i]<-coorre_ap[i,j]}}
probs_ocor_ap<-coorre_ap/nrow(bow_ap)
#
#set.seed(100)
Sf<-NULL
for (K_2 in topic_sequence){
	cat('\n', K_2, rep)
	#starting point
	S_2_ap<-sample(1:K_2,nrow(y_ap),replace=T)
	S_tot_2_ap<-S_2_ap
	n_k_2_ap<-NULL
	for (k in 1:K_2) n_k_2_ap[k]<-sum(S_2_ap==k)
	posteriori_theta_2_ap <-matrix(rDiric(n_k_2_ap+a0),nrow=K_2,ncol=1) 
	posteriori_theta_tot_2_ap<-posteriori_theta_2_ap
	m_2_ap =  matrix(data=NA,nrow=K_2,ncol=ncol(y_ap))
	for (k in 1:K_2) m_2_ap[k,]<-apply(as.matrix(y_ap[S_2_ap==k,],ncol=ncol(y_ap)),2,sum)
	posteriori_phi_k_2_ap<-NULL
	for (k in 1:K_2) posteriori_phi_k_2_ap<-rbind(posteriori_phi_k_2_ap,rDiric(m_2_ap[k,] + b0))
	posteriori_phi_tot_2_ap = posteriori_phi_k_2_ap
	### GIBBS
	enableJIT(3)
	for (it in 1:iter){ 
		#Atualiza os Ss
		for (i in 1:nrow(y_ap)){
			log_probs1_ap<-NULL
			for (top in 1:K_2) log_probs1_ap[top]<-dmultinom(y_ap[i,], size = sum(y_ap[i,]),prob=posteriori_phi_k_2_ap[top,], log = TRUE)
			log_probs_ap<-log(c(posteriori_theta_2_ap))+log_probs1_ap
			log_probs_ap<-log_probs_ap-max(log_probs_ap)
			probs_ap<-exp(log_probs_ap)
			probs_ap<-probs_ap/sum(probs_ap)
			S_2_ap[i]<-rDiscreta(probs_ap)}
		S_tot_2_ap<-rbind(S_tot_2_ap,S_2_ap)
  		#Atualiza theta
		for (k in 1:K_2) n_k_2_ap[k]<-sum(S_2_ap==k)
		posteriori_theta_2_ap<-rDiric(n_k_2_ap+a0)
		posteriori_theta_tot_2_ap<-cbind(posteriori_theta_tot_2_ap,posteriori_theta_2_ap)
		#Atualiza phi
		for (k in 1:K_2){ 
			m_2_ap[k,]<-apply(as.matrix(y_ap[S_2_ap==k,],ncol=ncol(y_ap)),2,sum)
			posteriori_phi_k_2_ap[k,] = rDiric(m_2_ap[k,] + b0)}}
	#Theta from the last iteration
	Sf<-rbind(Sf,S_2_ap) 
	theta_2 = posteriori_theta_tot_2_ap[,iter+1]
	#
	### Selection criteria
	#
	soma1_2<-0
	soma2_2<-0
	soma3_2<-0
	soma4_2<-0
	soma5_2<-0
	posicoes_2 <- matrix(0,nrow = K_2, ncol = 10)
	for (topicos in 1:K_2){
		soma4_2<-soma4_2+n_k_2_ap[topicos]*log(theta_2[topicos])
		for (i in 1:ncol(y_ap)) soma5_2<-soma5_2+m_2_ap[topicos,i]*log(posteriori_phi_k_2_ap[topicos,i])
		for (i in 1:(length(posicoes_2[topicos,])-1)){
			for (j in (i+1):length(posicoes_2[topicos,])){
				posicoes_2[topicos,]<-rev(tail(order(posteriori_phi_k_2_ap[topicos,]), 10))      
				soma1_2<-soma1_2+log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])*(probs_ocor_ap[posicoes_2[topicos,j],posicoes_2[topicos,j]])))                                                                                                                                                                          
				soma2_2<-soma2_2+(log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])*(probs_ocor_ap[posicoes_2[topicos,j],posicoes_2[topicos,j]])))/-log(probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]]))
				soma3_2 <-soma3_2+log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])))}}}
	pmi_mm_2<-soma1_2/(K_2*45)
	npmi_mm_2<-soma2_2/(K_2*45)
	lcp_mm_2<-soma3_2/(K_2*45)
	vero_mm_2<-soma4_2+soma5_2
	perp_mm_2<-exp(-vero_mm_2/sum(y_ap))
	BIC<-vero_mm_2-(ncol(y_ap)*K_2-1)/2*log(sum(y_ap))
	cat('',vero_mm_2,file=paste(caminho,"vero_MM_dados1_",rep,".txt",sep=""),append=T)
	cat('',perp_mm_2,file=paste(caminho,"perp_MM_dados1_",rep,".txt",sep=""),append=T)
	cat('',BIC,file=paste(caminho,"BIC_MM_dados1_",rep,".txt",sep=""),append=T)	
	cat('',pmi_mm_2,file=paste(caminho,"pmi_MM_dados1_",rep,".txt",sep=""),append=T)	
	cat('',npmi_mm_2,file=paste(caminho,"npmi_MM_dados1_",rep,".txt",sep=""),append=T)	
	cat('',lcp_mm_2,file=paste(caminho,"lcp_MM_dados1_",rep,".txt",sep=""),append=T)}
#
pmi<-scan(file=paste(caminho,"pmi_MM_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(pmi==max(pmi))]
cat('',K,file=paste(caminho,"K_pmi_MM_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_pmi_MM_dados1_",rep,".txt",sep=""),append=T)
#
npmi<-scan(file=paste(caminho,"npmi_MM_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(npmi==max(npmi))]
cat('',K,file=paste(caminho,"K_npmi_MM_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_npmi_MM_dados1_",rep,".txt",sep=""),append=T)
#
lcp<-scan(file=paste(caminho,"lcp_MM_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(lcp==max(lcp))]
cat('',K,file=paste(caminho,"K_lcp_MM_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_lcp_MM_dados1_",rep,".txt",sep=""),append=T)
#
vero<-scan(file=paste(caminho,"vero_MM_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(vero==max(vero))]
cat('',K,file=paste(caminho,"K_vero_MM_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_vero_MM_dados1_",rep,".txt",sep=""),append=T)
#
perp<-scan(file=paste(caminho,"perp_MM_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(perp==min(perp))]
cat('',K,file=paste(caminho,"K_perp_MM_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_perp_MM_dados1_",rep,".txt",sep=""),append=T)
#
bic<-scan(file=paste(caminho,"BIC_MM_dados1_",rep,".txt",sep=""))
K<-topic_sequence[which(bic==max(bic))]
cat('',K,file=paste(caminho,"K_bic_MM_dados1.txt",sep=""),append=T)
cat('',Sf[topic_sequence==K,],file=paste(caminho,"Sj_bic_MM_dados1_",rep,".txt",sep=""),append=T)	
}	
#
#################################### RUN MM DP
#
library(compiler)
library(MGLM)
#
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  val}
###
################
### gera valores da posteriori de alpha
#
gera_eta_alpha<-function(alpha,a_alph_prior,b_alph_prior,K,m){
	eta_aux<-rbeta(1,alpha+1,m)
	aux_prob<-(a_alph_prior+K-1)/(m*(b_alph_prior-log(eta_aux)))
	prob_alpha<-aux_prob/(1+aux_prob)
	unif<-runif(1)
	if (unif<=prob_alpha) alpha<-rgamma(1,a_alph_prior+K,b_alph_prior-log(eta_aux)) else alpha<-rgamma(1,a_alph_prior+K-1,b_alph_prior-log(eta_aux))
	return(alpha)}
#
m<-nrow(y_verd)
amostrasfin<-100
burnin<-20
saltos<-3
AmostrasTotal<-burnin+amostrasfin*saltos
alpha<-1
a_alph_prior<-b_alph_prior<-0.01 # parâmetros da dist gama para o parâmetro de concentração
#
replicas<-20
set.seed(100)
for (rep in 1:replicas){
dados<-matrix(0,nrow=m,ncol=ncol(y_verd))
for (i in 1:m) dados[i,]<-rmultinom(1,size=sum(y_verd[i,]),phi[S_verd[i]+1,])
#
y_ap=as.matrix(dados)
filtra<-apply(y_ap,2,sum)
if (sum(filtra==0)>0) y_ap<-y_ap[,-which(filtra==0)]
#
Nt<-nrow(y_ap)
#
# Initialize Sj and posterior parameters
#
#set.seed(100)
beta0<-rep(1,ncol(y_ap)) # G0
Sj<-sample(1:20,Nt,replace=TRUE)
nk<-table(Sj)
K<-length(nk)
betapost<-matrix(0,nrow=ncol(y_ap),ncol=K)
#
for (k in 1:K) betapost[,k]<-apply(y_ap[Sj==k,],2,sum)+beta0
#
# Gibbs sampling of cluster indicator
#
enableJIT(3)
for (int in 1:AmostrasTotal){
	cat('\n',rep,int,K,alpha)
	for (i in 1:Nt){
		nk<-numeric(K)
		for (k in 1:K) 	nk[k]<-sum(Sj[-i]==k)
		betapostnew<-betapost
		if (nk[Sj[i]]>0) betapostnew[,Sj[i]]<-betapostnew[,Sj[i]]-y_ap[i,]
		prob<-numeric(K)
		for (k in 1:K) prob[k]<-log(nk[k])+ddirmn(y_ap[i,], betapostnew[,k])
		prob<-c(prob,log(alpha)+ddirmn(y_ap[i,], beta0))
		prob<-prob-max(prob)
		prob<-exp(prob)		
		Snew<-rDiscreta(prob/sum(prob))
#
		if (Snew != Sj[i] & Snew <= K){
			Sj[i]<-Snew
			nk[Sj[i]]<-nk[Sj[i]]+1
			betapost<-betapostnew
			betapost[,Sj[i]]<-betapost[,Sj[i]]+y_ap[i,]}
#
		if (Snew != Sj[i] & Snew > K){
			Sj[i]<-Snew
			nk<-c(nk,1)
			betapost<-betapostnew
			betapost<-cbind(betapost,y_ap[i,]+beta0)}
#
		while (length(table(Sj))<max(Sj)){ # exclude empty clusters
			categr<-as.numeric(as.character(data.frame(table(Sj))[,1]))
			categd<-seq(1:length(table(Sj)))
			dif<-which(categr!=categd)
			Sj[which(Sj>dif[1])]<-Sj[which(Sj>dif[1])]-1
			betapost<-matrix(c(betapost[,-dif[1]]),nrow=ncol(y_ap))
			K<-ncol(betapost)}
#
		if (length(table(Sj))<K){
			betapost<-matrix(c(betapost[,-K]),nrow=ncol(y_ap))
			K<-ncol(betapost)}
#
		K<-ncol(betapost)}
	#
	alpha<-gera_eta_alpha(alpha,a_alph_prior,b_alph_prior,K,Nt)
	#
	if (int>burnin & int%%saltos==0){
		cat('',Sj,file=paste(caminho,"Sj_BBC_dados1_",rep,".txt",sep=""),append=T)
		cat('',alpha,file=paste(caminho,"Alpha_BBC_dados1_",rep,".txt",sep=""),append=T)}
}
Sj<-scan(file=paste(caminho,"Sj_BBC_dados1_",rep,".txt",sep=""))
S<-matrix(Sj,ncol=m,nrow=length(Sj)/m,byrow=TRUE)
Sj.j<-S #matriz de agrupamento a posteriori
prob.eq<-matrix(0,nrow=ncol(Sj.j),ncol=ncol(Sj.j))
enableJIT(3)
for (i in 1:ncol(Sj.j)){
	for (j in 1:ncol(Sj.j)){
		prob.eq[i,j]<-round(sum(Sj.j[,i]==Sj.j[,j])/length(Sj.j[,i]),5)*100}}
#
thresh<-0.50*100 # definindo os grupos finais
clust_f<-c(1,rep(0,(ncol(Sj.j)-1)))
for (i in 2:ncol(Sj.j)){
	if (max(prob.eq[i,1:(i-1)])>thresh) clust_f[i]<-clust_f[which(prob.eq[i,1:(i-1)]==max(prob.eq[i,1:(i-1)]))[1]] else clust_f[i]<-max(clust_f[1:(i-1)]+1)}
#
if (sum(table(clust_f)==1)>0){
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
		clust_f[which(clust_f>dif[1])]<-clust_f[which(clust_f>dif[1])]-1}}
K<-length(table(clust_f))
somaWs<-NULL
for (k in 1:K) somaWs<-rbind(somaWs,apply(y_ap[clust_f==k,],2,sum))
logmargvero<-0
for (i in 1:nrow(y_ap)) logmargvero<-logmargvero+ddirmn(y_ap[i,], beta0+somaWs[clust_f[i],])
cat('',K,file=paste(caminho,"K_MM_DP_dados1.txt",sep=""),append=T)
cat('',logmargvero,file=paste(caminho,"logmargvero_MM_DP_dados1.txt",sep=""),append=T)
cat('',clust_f,file=paste(caminho,"Sj_MM_DP_dados1_",rep,".txt",sep=""),append=T)	
}
#
############################# RUN LDA HDP
#
library(hdp)
library(compiler)
hdp_setdata2<-function (hdp, dpindex, data) 
{
    if (class(hdp) != "hdpState") 
        stop("hdp must have class hdpState")
    if (!validObject(hdp)) 
        stop("input hdp is not a valid hdpState object")
    if (any(dpindex < 1) | any(dpindex > hdp@numdp) | any(dpindex%%1 != 
        0) | any(duplicated(dpindex))) {
        stop("dpindex must be positive integers no greater than\n         numdp(hdp) with no duplicates")
    }
#    if (!class(data) %in% c("matrix", "data.frame")) {
#        stop("data must be data.frame or matrix")
#    }
    if (nrow(data) != length(dpindex)) 
        stop("nrow(data) must equal length(dpindex)")
    if (ncol(data) != numcateg(hdp)) 
        stop("ncol(data) must equal numcateg(hdp)")
    if (any(data%%1 != 0) | any(data < 0)) {
        stop("data must contain non-negative integer values")
    }
    datass <- apply(data, 1, function(x) rep(1:ncol(data), x))
    HELDOUT <- 0L
    for (jj in 1:length(dpindex)) {
        if (hdp@dpstate[dpindex[jj]] != HELDOUT) {
            stop("Cannot set data for DPs that are not held out")
        }
        hdp@dp[[dpindex[jj]]]@numdata <- length(datass[[jj]])
        hdp@dp[[dpindex[jj]]]@datass <- datass[[jj]]
    }
    if (!validObject(hdp)) 
        warning("Not a valid hdpState object.")
    return(hdp)
}
#
m<-nrow(y_verd)
#
replicas<-20
set.seed(100)
for (rep in 1:replicas){
dados<-matrix(0,nrow=m,ncol=ncol(y_verd))
for (i in 1:m) dados[i,]<-rmultinom(1,size=sum(y_verd[i,]),phi[S_verd[i]+1,])
#
y_ap=as.matrix(dados)
filtra<-apply(y_ap,2,sum)
if (sum(filtra==0)>0) y_ap<-y_ap[,-which(filtra==0)]
#
hdp<-hdp_init(ppindex=0, cpindex=1, hh=rep(1, ncol(y_ap)), alphaa=rep(0.01, 100), alphab=rep(0.01, 100))
hdp<-hdp_adddp(hdp, nrow(y_ap), 1, 1)
hdp<-hdp_setdata2(hdp, 2:(nrow(y_ap)+1), as.matrix(y_ap))
hdp<-dp_activate(hdp, 1:(nrow(y_ap)+1), initcc=20) #initcc é o número de grupos iniciais
hdp<-hdp_posterior(hdp, burnin=1500, n=500, space=5, cpiter=2)
#
a<-hdp_extract_components(hdp)
b<-a@clust_categ_counts
grupo_pal<-NULL
enableJIT(3)
for (i in 1:ncol(y_ap)){
	teste<-which(b[[500]][i,]==max(b[[500]][i,]))
	grupo_pal[i]<-teste[1]}
grupo_texto<-NULL
enableJIT(3)
for (i in 1:nrow(y_ap)){
	pal<-NULL
	for (j in 1:ncol(y_ap)) pal<-c(pal,rep(grupo_pal[j],y_ap[i,j]))
	grupo_texto[i]<-names(table(pal))[which.max(table(pal))]}
K<-length(table(grupo_texto))
cat('\n',rep, K)
cat('',mean(hdp@lik[seq(3000,4000,5)]),file=paste(caminho,"loglikmedia_LDA_HDP_dados1.txt",sep=""),append=T)
cat('',quantile(hdp@lik[seq(3000,4000,5)],0.5),file=paste(caminho,"loglikmediana_LDA_HDP_dados1.txt",sep=""),append=T)
cat('',K,file=paste(caminho,"K_LDA_HDP_dados1.txt",sep=""),append=T)
cat('',grupo_texto,file=paste(caminho,"Sj_LDA_HDP_dados1_",rep,".txt",sep=""),append=T)
}


