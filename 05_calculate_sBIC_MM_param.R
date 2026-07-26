library(compiler)
library(gmp)
library(Rmpfr)
#
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
caminho<-"/Users/Daiane/Downloads/bbcsport/"
#caminho<-"/home/mariella2/Edvaldo/"
dados<-scan(paste(caminho,"bbcsport.txt",sep=""))
y_verd<-matrix(dados,nrow=737,ncol=4613,byrow=TRUE)
summary(apply(y_verd,1,sum))
S_verd<-scan(file=paste(caminho,"bbcsport_label.txt",sep=""))
#
#library(lsa)
retirar<-which(apply(y_verd,2,sum)<10)
y_verd<-y_verd[,-retirar] # tenho que tirar se não da palavra com frequência zero e não roda o NMF e talvez 
topic_sequence <- c(2:20)
N <- nrow(y_verd) # Number of documents
V <- ncol(y_verd) # Vocabulary size
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/Simulacao/Dados_2/MM_param_dados2/"
#
for (rep in 1:20){
log_liks<-scan(file=paste(caminho,"vero_MM_dados2_",rep,".txt",sep=""))
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
ntot<-nrow(y_verd)
#
logLKk<-list()
for (i in 1:length(topic_sequence)) logLKk[[i]]<-log_LKk(log_liks[i], param_sBIC[[i]], ntot)
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
	cat('\n', rep, i)
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
K<-topic_sequence[which(sBIC==max(sBIC))]
cat('',K,file=paste(caminho,"K_sBIC_MM_dados2.txt",sep=""),append=T)
cat('',as.numeric(sBIC),file=paste(caminho,"sBIC_MM_dados2_",rep,".txt",sep=""))
#
param_sBIC<-list()
enableJIT(3)
for (k in 1:length(topic_sequence)){
	lambda<-NULL
	m<-NULL
	for (j in 1:k){
		param<-learning_coefficient_2(N, V, topic_sequence[k], topic_sequence[j])
		lambda[j]<-param[1]
		m[j]<-param[2]}
	param_sBIC[[k]]<-rbind(lambda,m)}
#
ntot<-nrow(y_verd)
#
logLKk<-list()
for (i in 1:length(topic_sequence)) logLKk[[i]]<-log_LKk(log_liks[i], param_sBIC[[i]], ntot)
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
	cat('\n', rep, i)
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
K<-topic_sequence[which(sBIC==max(sBIC))]
cat('',K,file=paste(caminho,"K_sBIC2_MM_dados2.txt",sep=""),append=T)
cat('',as.numeric(sBIC),file=paste(caminho,"sBIC2_MM_dados2_",rep,".txt",sep=""))
}
