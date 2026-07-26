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
################# sports data set #################
#caminho<-"/Users/Daiane/Downloads/bbcsport/"
caminho<-"/home/mariella/Edvaldo/"
dados<-scan(paste(caminho,"bbcsport.txt",sep=""))
#
dados_ap<-matrix(dados,nrow=737,ncol=4613,byrow=TRUE)
#
summary(apply(dados_ap,1,sum))
y_ap = as.matrix(dados_ap)
#dados_ap = data.frame(dados)
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
################# MM with a fixed number of topics #################
#
### Fitting the model and several criteria
a0<-1
b0<-1
iter<-500
start_time <- Sys.time()
set.seed(100)
for (K_2 in c(2:20,25,30,35,40,45,50,75,100,125,150,175,200,250,300,400,500)){
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
		cat('\n', K_2, it)
		#Atualiza os Ss
		for (i in 1:nrow(dados_ap)){
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
	BIC<-vero_mm_2-(ncol(dados_ap)*K_2-1)/2*log(sum(dados_ap))
	cat('',vero_mm_2,file=paste(caminho,"vero_MM_sport.txt",sep=""),append=T)
	cat('',perp_mm_2,file=paste(caminho,"perp_MM_sport.txt",sep=""),append=T)
	cat('',BIC,file=paste(caminho,"BIC_MM_sport.txt",sep=""),append=T)	
	cat('',pmi_mm_2,file=paste(caminho,"pmi_MM_sport.txt",sep=""),append=T)	
	cat('',npmi_mm_2,file=paste(caminho,"npmi_MM_sport.txt",sep=""),append=T)	
	cat('',lcp_mm_2,file=paste(caminho,"lcp_MM_sport.txt",sep=""),append=T)}
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) 		
#
##### calculating sBIC
#
learning_coefficient_1 <- function(n_docs, n_words, super_model, sub_model) {
  # Compute learning coefficient and its multiplicity
  # using formulas from Drton and Plummer (2017)
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

  learn_coeff <- 1/2*(H*(M-1)+r-1)
  m<-1
  return(c(learn_coeff,m))
}
learning_coefficient_2 <- function(n_docs, n_words, super_model, sub_model) {
  # Compute learning coefficient and its multiplicity
  # using formulas from Drton and Plummer (2017)
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

  learn_coeff <- 1/2*(H*r*(M-1)-1)
  m<-1
  return(c(learn_coeff,m))
}
#
#### calculating log(LKk)
#
log_LKk<-function(loglikKuper, paramKuper, N) {
	logLKk<-loglikKuper-paramKuper[1,]*log(N)+(paramKuper[2,]-1)*log(log(N))
	return(logLKk)}
#
#######
#
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/BBC_sport_DP/MM_param_sport/"
topic_sequence <- c(2:20,25,30,35,40,45,50,75,100,125,150,175,200,250,300,400,500)
N <- nrow(dados_ap) # Number of documents
V <- ncol(dados_ap) # Vocabulary size
log_liks<-scan(file=paste(caminho,"vero_MM_sport.txt",sep=""))
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
		param<-learning_coefficient_1(N, V, topic_sequence[k], topic_sequence[j])
		lambda[j]<-param[1]
		m[j]<-param[2]}
	param_sBIC[[k]]<-rbind(lambda,m)}
#
ntot<-nrow(dados_ap)
#
logLKk<-list()
for (i in 1:length(topic_sequence)) logLKk[[i]]<-log_LKk(log_liks[i], param_sBIC[[i]], ntot)
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
	cat('\n', i)
	b<-Reduce(`+`,L[1:(i-1)])-LKk[[i]][i]
	c<-mpfr(0, precBits = 128)
	for (j in 1:(i-1)) c<-c-LKk[[i]][j]*L[[j]]
    prec <- 5000
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
cat('',as.numeric(sBIC),file=paste(caminho,"sBIC_MM_sport.txt",sep=""))
cat('',mn,file=paste(caminho,"mn_MM_sport.txt",sep=""))
#
#### fiz a segunda maneira e chamei de sBIC_MM_sporte_2
#
#
#### analize the results
#
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/BBC_sport_DP/MM_param_sport/"
log_liks<-scan(file=paste(caminho,"vero_MM_sport.txt",sep=""))
perp<-scan(file=paste(caminho,"perp_MM_sport.txt",sep=""))
pmi<-scan(file=paste(caminho,"pmi_MM_sport.txt",sep=""))
npmi<-scan(file=paste(caminho,"npmi_MM_sport.txt",sep=""))
lcp<-scan(file=paste(caminho,"lcp_MM_sport.txt",sep=""))
BIC<-scan(file=paste(caminho,"BIC_MM_sport.txt",sep=""))
sBIC2<-scan(file=paste(caminho,"sBIC_MM_sport.txt",sep=""))
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
