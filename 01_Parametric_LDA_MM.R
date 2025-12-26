library(ggplot2)
library(tibble)
library(tidyr)
library(readr)
library(dplyr)
library(tidytext)
library(topicmodels)
library(LDA)

################# ASSOCIATED PRESS data set #################

data("AssociatedPress", package = "topicmodels")

ap_td <- tidy(AssociatedPress) 

ap_dtm <- ap_td %>%
  anti_join(stop_words, by = c(term = "word")) %>%
  cast_dtm(document, term, count)

dados = as.matrix(ap_dtm)

dados_ap = data.frame(dados)

################# LDA exemple - K = 2 #################

### Fitting the model
ap_lda2_ap <- LDA(dados_ap, k = 2, method="Gibbs")

### TOP 10 TERMS in each topic
ap_top_terms2 <- terms(ap_lda2_ap, 10)  

### correlation matriz among words
bow_ap <- matrix(unlist(dados_ap), nrow = nrow(dados_ap))
coorre_ap<-matrix(0,ncol(bow_ap),ncol(bow_ap))

for (i in 1:nrow(coorre_ap)){
  for (j in i:ncol(coorre_ap)){
    coorre_ap[i,j]<-sum(bow_ap[,i]>0 & bow_ap[,j]>0)+0.001
    if (j>i) coorre_ap[j,i]<-coorre_ap[i,j]}}
#Na diagonal principal temos o número de documentos em que cada palavra 
ocorre e fora da diagonal principal temos o número de documentos em que o 
par de palavras (i,j) acontece

probs_ocor_ap<-coorre_ap/nrow(bow_ap) 
#Aqui de fato é a probabilidade de cada palavra e dos pares de palavras

### Selection criteria 

#Likelihood
vero2_lda <- ap_lda2_ap@loglikelihood  

#Perplexity
perp2_lda <- exp(-vero2_ap/sum(dados_ap))

#auxiliar adjust to calculate the other criteria 
ap_top_terms_2 = t(ap_top_terms2)
ap_top_terms_2_v2 = matrix(unlist(ap_top_terms_2), nrow = nrow(ap_top_terms_2))
ap_top_terms_2_v3<-matrix(0,nrow(ap_top_terms_2_v2),ncol(ap_top_terms_2_v2))

for (i in 1:K_2){
  for (j in 1:10){
    ap_top_terms_2_v3[i,j] = which(colnames(dados_ap) == ap_top_terms_2_v2[i,j])
    
  }}
pa_top_terms2 = ap_top_terms_2_v3


pmi_lda<-NULL
npmi_lda<-NULL
lcp_lda<-NULL
modelos_ap<-c(2) #vetor de ks
for (mod in 1:length(modelos_ap)){
  K<-modelos_ap[mod]
  top_termos<-get(paste("pa_top_terms",eval(modelos_ap[mod]),sep=""))
  
  soma1<-0
  soma2<-0
  soma3<-0
  for (topicos in 1:K){
    for (i in 1:(length(top_termos[topicos,])-1)){
      for (j in (i+1):length(top_termos[topicos,])){
        
        soma1<-soma1+
          log((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]])/
                ((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,i]])*
                   (probs_ocor_ap[top_termos[topicos,j],top_termos[topicos,j]])))                                                                                                                                                                          
        
        
        soma2<-soma2+
          (log((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]])/
                 ((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,i]])*
                    (probs_ocor_ap[top_termos[topicos,j],top_termos[topicos,j]])))
             /-log(probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]]))
        
        soma3 <-soma3+
          log((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,j]])/
                ((probs_ocor_ap[top_termos[topicos,i],top_termos[topicos,i]])))
        
        
      }}}
  pmi_lda[mod]<-soma1/(K*45)
  npmi_lda[mod]<-soma2/(K*45)
  lcp_lda[mod]<-soma3/(K*45)
}
pmi_lda
npmi_lda
lcp_lda

################# MM exemplo - K = 2 #################

#Simula valores de uma distribuição discreta para simular as variáveis S

rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  val}

#Simula vetores da distribuição Dirichlet
rDiric<-function(gama){
  X<-rgamma(length(gama),gama,1)
  Y<-X/sum(X)
  return(Y)}

y_ap<-as.matrix(dados_ap) #documentos

#Valores dos hiperparâmetros das distribuições a priori (alfa e beta)
a0<-1
b0<-1


K_2<-2 #fixar o número de tópicos

#Chute inicial (aleatório) para o vetor S
set.seed(nrow(dados_ap))
S_2_ap<-sample(1:K_2,nrow(y_ap),replace=T)
S_tot_2_ap<-S_2_ap
table(S_2_ap)
n_k_2_ap<-NULL
for (k in 1:K_2) n_k_2_ap[k]<-sum(S_2_ap==k)

posteriori_theta_2_ap <-matrix(rDiric(n_k_2_ap+a0),nrow=K_2,ncol=1) 
#posteriori Dirichlet correta da probabilidade de cada tópico
posteriori_theta_tot_2_ap<-posteriori_theta_2_ap

#Primeiro passo das probabilidades de cada palavra dentro de cada tópico da 
posteriori Dirichlet dela. 
m_2_ap =  matrix(data=NA,nrow=K_2,ncol=ncol(dados_ap))
for (k in 1:K_2) m_2_ap[k,]<-apply(as.matrix(y_ap[S_2_ap==k,],
ncol=ncol(y_ap)),2,sum)

posteriori_phi_k_2_ap<-NULL
for (k in 1:K_2) posteriori_phi_k_2_ap<-rbind(posteriori_phi_k_2_ap,
rDiric(m_2_ap[k,] + b0))
posteriori_phi_tot_2_ap = posteriori_phi_k_2_ap

### GIBBS

library(compiler)

iter<-500
for (it in 1:iter){ 
  cat('\n', K_2, it)
  #Atualiza os Ss
  enableJIT(3)
  for (i in 1:nrow(dados_ap)){
    log_probs1_ap<-NULL
    for (top in 1:K_2) {
      log_probs1_ap[top]<-dmultinom(y_ap[i,], size = sum(y_ap[i,]),
      prob=posteriori_phi_k_2_ap[top,], log = TRUE)}
    log_probs_ap<-log(c(posteriori_theta_2_ap))+log_probs1_ap
    log_probs_ap<-log_probs_ap-max(log_probs_ap)
    probs_ap<-exp(log_probs_ap)
    probs_ap<-probs_ap/sum(probs_ap)
    S_2_ap[i]<-rDiscreta(probs_ap)}
  S_tot_2_ap<-rbind(S_tot_2_ap,S_2_ap)
  
  #Atualiza theta
  for (k in 1:K_2){
    n_k_2_ap[k]<-sum(S_2_ap==k)}
  posteriori_theta_2_ap<-rDiric(n_k_2_ap+a0)
  posteriori_theta_tot_2_ap<-cbind(posteriori_theta_tot_2_ap, 
  posteriori_theta_2_ap)
  
  #Atualiza phi
  for (k in 1:K_2){ 
    m_2_ap[k,]<-apply(as.matrix(y_ap[S_2_ap==k,],ncol=ncol(y_ap)),2,sum)
    posteriori_phi_k_2_ap[k,] = rDiric(m_2_ap[k,] + b0)}
  
}

#Theta da última iteração 
theta_2 = posteriori_theta_tot_2_ap[,501]

### Selection criteria

soma1_2<-0
soma2_2<-0
soma3_2<-0
posicoes_2 <- matrix (0,nrow = K_2, ncol = 10)

for (topicos in 1:K_2){
  for (i in 1:(length(posicoes_2[topicos,])-1)){
    for (j in (i+1):length(posicoes_2[topicos,])){
      
      posicoes_2[topicos,]<-rev(tail(order(posteriori_phi_k_2_ap[topicos,]), 10))
      
      soma1_2<-soma1_2+
        log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/
              ((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])*
                 (probs_ocor_ap[posicoes_2[topicos,j],posicoes_2[topicos,j]])))                                                                                                                                                                          
      
      
      soma2_2<-soma2_2+
        (log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/
               ((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])*
                  (probs_ocor_ap[posicoes_2[topicos,j],posicoes_2[topicos,j]])))/
           -log(probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]]))
      
      soma3_2 <-soma3_2+
        log((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,j]])/
              ((probs_ocor_ap[posicoes_2[topicos,i],posicoes_2[topicos,i]])))
      
      
    }}}
pmi_mm_2<-soma1_2/(K_2*45)
npmi_mm_2<-soma2_2/(K_2*45)
lcp_mm_2<-soma3_2/(K_2*45)


soma4_2<-0
soma5_2<-0

for (topicos in 1:K_2){
  
  soma4_2<-soma4_2+
    n_k_2_ap[topicos]*log(theta_2[topicos])
  
  for (i in 1:ncol(y_ap)){
    
    soma5_2<-soma5_2+
      m_2_ap[topicos,i]*log(posteriori_phi_k_2_ap[topicos,i])
    
  }}

vero_mm_2 = soma4_2+soma5_2
perp_mm_2 = exp(-vero_2/sum(y_ap))
