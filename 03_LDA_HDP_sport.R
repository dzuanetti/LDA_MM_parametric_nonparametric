##
###### Lendo dados
##
caminho<-"/Users/Daiane/Downloads/bbcsport/"
#caminho<-"/home/daiane/Edvaldo/"
dados<-scan(paste(caminho,"bbcsport.txt",sep=""))
#
y_ap<-matrix(dados,nrow=737,ncol=4613,byrow=TRUE)
#
summary(apply(y_ap,1,sum))
#
#
#http://127.0.0.1:13369/doc/html/Search?objects=1&port=13369
#https://github.com/nicolaroberts/hdp

library(hdp)
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
start_time <- Sys.time()
hdp<-hdp_init(ppindex=0, cpindex=1, hh=rep(1, ncol(y_ap)), alphaa=rep(0.01, 100), alphab=rep(0.01, 100))
hdp<-hdp_adddp(hdp, nrow(y_ap), 1, 1)
hdp<-hdp_setdata2(hdp, 2:(nrow(y_ap)+1), as.matrix(y_ap))
hdp<-dp_activate(hdp, 1:(nrow(y_ap)+1), initcc=1, seed=100)
hdp<-hdp_posterior(hdp, burnin=1500, n=500, space=5, cpiter=2, seed=100)
#
caminho<-"/Users/Daiane/Dropbox/Artigo_Edvaldo/BBC_sport_DP/"		
cat('',hdp@lik,file=paste(caminho,"lik_LDA_HDP_sport.txt",sep=""))
mean(hdp@lik[seq(3000,4000,5)])
quantile(hdp@lik[seq(3000,4000,5)],0.5)
#
plot_lik(hdp)
numcluster(hdp)
plot_numcluster(hdp)
#
#hdp_multi_chain(hdp)
#
a<-hdp_extract_components(hdp)
slotNames(a)
a@numcluster
b<-a@clust_categ_counts #quantas vezes cada palavra caiu em cada grupo por iteração.
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)
#
grupo_pal<-NULL
for (i in 1:ncol(y_ap)){
	teste<-which(b[[500]][i,]==max(b[[500]][i,]))
	grupo_pal[i]<-teste[1]}
#
grupo_texto<-NULL
for (i in 1:nrow(y_ap)){
	pal<-NULL
	for (j in 1:ncol(y_ap)) pal<-c(pal,rep(grupo_pal[j],y_ap[i,j]))
	nk<-NULL
	for (k in 1:6) nk[k]<-sum(pal==k)
	grupo_texto<-rbind(grupo_texto,round(nk/sum(nk)*100,1))}
#
grupo_fim<-NULL
set.seed(100)
for (i in 1:nrow(grupo_texto)){
	grupo_fim[i]<-rDiscreta(grupo_texto[i,]/100)
}

grupo_texto<-NULL
for (i in 1:nrow(y_ap)){
	pal<-NULL
	for (j in 1:ncol(y_ap)) pal<-c(pal,rep(grupo_pal[j],y_ap[i,j]))
	grupo_texto[i]<-names(table(pal))[which.max(table(pal))]}
#
caminho<-"/Users/Daiane/Downloads/bbcsport/"
#
S_verd<-scan(file=paste(caminho,"bbcsport_label.txt",sep=""))
table(S_verd,grupo_fim)
table(S_verd,grupo_texto) # esse aqui ficou melhor, vou manter



"seed"               "settings"           "hdp"                "lik"                "numcluster"        
"cp_values"          "clust_categ_counts" "clust_dp_counts"    "numcomp"            "prop.ex"           
"comp_cos_merge"     "comp_categ_counts"  "comp_dp_counts"     "comp_categ_distn"   "comp_dp_distn"   


