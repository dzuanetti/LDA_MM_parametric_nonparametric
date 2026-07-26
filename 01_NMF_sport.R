# Load libraries
library(NMF)
#caminho<-"/Users/Daiane/Downloads/bbcsport/"
caminho<-"/home/mariella2/Edvaldo/"
dados<-scan(paste(caminho,"bbcsport.txt",sep=""))
#
dados_ap<-matrix(dados,nrow=737,ncol=4613,byrow=TRUE)
#
summary(apply(dados_ap,1,sum))
y_ap = as.matrix(dados_ap)
##
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
#### run NMF
#
start_time <- Sys.time()
enableJIT(3)
for (K_2 in c(2:20,25,30,35,40,45,50,75,100,125,150,175,200,250,300,400,500)){
	cat('\n', K_2)
	set.seed(100)
	nmf_model <- nmf(y_ap, rank = K_2, method = "brunet", nrun = 3)
	#
	### select the most frequent words for each topic to selection criteria
	#
	basis_matrix <- basis(nmf_model)
	topico<-NULL
	for (i in 1:nrow(basis_matrix)) topico[i]<-which(basis_matrix[i,]==max(basis_matrix[i,]))
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
	cat('',pmi_nmf_2,file=paste(caminho,"pmi_NMF_sport.txt",sep=""),append=T)	
	cat('',npmi_nmf_2,file=paste(caminho,"npmi_NMF_sport.txt",sep=""),append=T)	
	cat('',lcp_nmf_2,file=paste(caminho,"lcp_NMF_sport.txt",sep=""),append=T)}
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) 		
#
#
#### analize the results
#
topic_sequence<-c(2:20,25,30,35,40,45,50,75,100,125,150,175,200,250,300,400,500)
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/BBC_sport_DP/NMF_sport/"
pmi<-scan(file=paste(caminho,"pmi_NMF_sport.txt",sep=""))
npmi<-scan(file=paste(caminho,"npmi_NMF_sport.txt",sep=""))
lcp<-scan(file=paste(caminho,"lcp_NMF_sport.txt",sep=""))
#
plot(topic_sequence,pmi,type='l')
K<-topic_sequence[which(pmi==max(pmi))]
K
plot(topic_sequence,npmi,type='l')
K<-topic_sequence[which(npmi==max(npmi))]
K
plot(topic_sequence,lcp,type='l')
K<-topic_sequence[which(lcp==max(lcp))]
K
