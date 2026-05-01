## 
# STA3 - TP5- Clustering (suite)
library(ggplot2)
library(FactoMineR)
library(factoextra)

rm(list=objects());graphics.off()
# setwd("à personnaliser")

#---chargement des données

load("seeds.RData")  # fichier binaire R
variety = seeds$variety
seeds = seeds[,-8]       # enlever la variable quali
head(seeds)
n = nrow(seeds); p = ncol(seeds)

###################
# CAH 
###################

########
###---a- cah: méthode "complete" par défaut
########
res1 = hclust(dist(seeds))
res1
# Cluster method   : complete 
# Distance         : euclidean 
# Number of objects: 210

########
###---b- dendrogramme
########
par(mfrow=c(1,2))
plot(res1,cex=0.5,hang=-1) # dendrogramme
abline(h=6)        # 4 classes si on coupe à la hauteur 6

#accéder à la classification
cutree(res1,h=6) # classification à une hauteur donnée
cutree(res1,k=4) # idem

#éboulis des hauteurs

names(res1)
head(res1$height) 
head(rev(res1$height))

barplot(rev(res1$height)[1:15],main="diagramme des hauteurs")
abline(h=6)       # 4 classes

# indice de silhouette
library(cluster)
sil = silhouette(cutree(res1,k=4), dist(seeds))
mean(sil[,3])         # indice moyen 0.37
fviz_silhouette(sil)  # visualisation


########
###---c- représentation graphique
########
library(FactoMineR)

# on fait la visu "à la main"
res.pca2 = PCA(cbind(seeds,grp=as.factor(variety)), 
               quali.sup = 8, graph=FALSE)
plot(res.pca2,habillage=8)
plot(res.pca2,axes=2:3,habillage=8)

# visu des classifications de 2 à 5 groupes
 
# on sauvegarde l'environnement avant le changement de marge
old.par <- par(no.readonly = TRUE)
# pour profiter au maximum de l'espace disponible
par (mfrow=c(2,2), mai=c(0.7,0.7,0.4,0.2))

for (i in 2:5) {
  grp = cutree(res1,h=rev(res1$height)[i])
  plot(res.pca2$ind$coord[,1],res.pca2$ind$coord[,2],
       col=grp, pch=grp,
       main=paste("nb de groupes=",i),
       xlim=c(-4,5),ylim=c(-4,4),cex=0.5,xlab="Dim 1",ylab="Dim 2")
}

# puis on utilise fviz_Cluster

# elle demande d'utiliser hcut qui appelle hclust + cutree
hc.cut = hcut(seeds, k = 4, hc_method = "complete")
sum(hc.cut$cluster!=cutree(res1,k=4)) #petite vérification

# les 4 graphes. On remarque que l'axe 1 est en miroir par
# rapport à notre programmation: 
# suivant le type de méthode utilisé pour calculer les directions
# propres u ou -u sont choisis
# cela ne change pas l'interprétation

plt =lapply(2:5, function(k) fviz_cluster(hcut(seeds, k = k, hc_method = "complete"),
                                          main="complete"))
cowplot::plot_grid(plt[[1]],plt[[2]],plt[[3]],plt[[4]],nrow=2,ncol=2)

# dendrogram avec groupes:

fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE)

########
###---d- autre distance
########
res2 = hclust(dist(seeds),method="ward.D2")

pltw =lapply(2:4, function(k) fviz_cluster(hcut(seeds, k = k, hc_method = "ward.D2"),
                                           main="ward.D2"))
cowplot::plot_grid(plt[[1]],plt[[2]],plt[[3]],
                   pltw[[1]],pltw[[2]],pltw[[3]],nrow=2,ncol=3)

# plot alluvial
library(alluvial)
alluvial(data.frame(complete=cutree(res1,k=4),
                    wardD2=cutree(res2,k=4)), 
         freq = rep(1,n), #col=2, 
         col = c("red", "blue", "green","orange")[cutree(res1,k=4)])


#comparaison numérique des partition avec l'index de Rand
mclust::adjustedRandIndex(cutree(res2,k=4), cutree(res1,k=4))


########
###---  centrer réduire?
########
seedscr = data.frame(scale(seeds))
rescr = hclust(dist(seedscr),method="ward.D2")
par(mfrow=c(1,2))
barplot(rev(rescr$height)[1:15],main="diagramme des hauteurs")
abline(h=15) 

# l'éboulis  définit trois groupes
plot(rescr,cex=0.5) # dendrogramme
abline(h=15)         # 3 classes

mclust::adjustedRandIndex(cutree(res2,k=3), cutree(rescr,k=3)) #0.68
par(mfrow=c(1,1))
alluvial(data.frame(brut=cutree(res2,k=3),
                    reduit=cutree(rescr,k=3)), 
         freq = rep(1,n), #col=2, 
         col = c("red", "blue", "green","orange")[cutree(res2,k=3)])

# La standardisation évite que les variables à grande échelle domine pour la contruction des clusters. 
mclust::adjustedRandIndex(cutree(rescr,k=3),variety)#0.79
mclust::adjustedRandIndex(cutree(res2,k=3),variety) #0.71

########
###---e- Clustering de variables
########
# on choisit une dissemblance adaptée aux variables
dd = as.dist((1 - cor(seeds))/2)
par(mfrow=c(1,1))
plot(hclust(dd,method="ward.D2")) 

########
### Rmixmod
########

########
###---a- modèle de mélange
########
library(Rmixmod)
set.seed(123)
res.mixmod = mixmodCluster(seeds, 3,
                           models=mixmodGaussianModel (list=c("Gaussian_pk_Lk_Ck")))

?mixmodGaussianModel 
# "Gaussian_pk_Lk_Ck": Volume, shape et orientation sont libres
# d'où le nombre de paramètres
7 * 3  + 3* 7*8/2  + 2     #107 = 7 * 3 (mu) + 3* 7*8/2 (var) + 2 (pi)

# res.mixmod est une nouvelle forme d'objet (S4)
# qui s'interroge de façon différente
bm = res.mixmod@bestResult
bm@model  #"Gaussian_pk_Lk_Ck"
bm@parameters@nbFreeParam  # 107 
bm@parameters@mean
bm@parameters@variance
bm@parameters@proportions

########
###---b- histogramme des probabilités conditionnelles
########

# on vérifie que la somme sur les lignes fait 1
apply(bm@proba,1,sum)

# on prend les plus grandes pour chaque observations
hist(apply(bm@proba,1,max)) # données en général bien séparées

########
###---c- les fonctions plot
########

plot(res.mixmod, main=bm@model )        # projection sur les axes
plotCluster(bm, seeds, main=bm@model)   # la visu sur deux axes canoniques
fviz_cluster(list(data=seeds,cluster=bm@partition))

########
###---d- comparaison
########

sil2 = silhouette(bm@partition, dist(seeds))
mean(sil2[,3])         # indice moyen 0.44
fviz_silhouette(sil2)  # visualisation

alluvial(data.frame(ward.d2=c(3,2,1)[cutree(res2,k=3)],
                    gmm=bm@partition), 
         freq = rep(1,n), #col=2, 
         col = c("red", "blue", "green","orange")[cutree(res2,k=3)])

mclust::adjustedRandIndex(bm@partition,variety) # 0.72

########
###---e- model par défaut
########

res.mixmod0 = mixmodCluster(seeds, 3)
bm0 = res.mixmod0@bestResult
bm0@model  #"Gaussian_pk_Lk_C"
7 * 3  +  (7*8/2 +2)  + 2
bm0@parameters@nbFreeParam #53

# on vérifie la proportionalité de variances
bm0@parameters@variance[[1]]/bm0@parameters@variance[[2]]
bm0@parameters@variance[[3]]/bm0@parameters@variance[[2]]

sil0 = silhouette(bm0@partition, dist(seeds))
mean(sil0[,3])         # indice moyen 0.38
fviz_silhouette(sil0)
mclust::adjustedRandIndex(bm0@partition,variety) # 0.55

### choix parmi tous les modèles à trois classes
res.mixmod_all = mixmodCluster(seeds, 3, 
                               models=mixmodGaussianModel(family = "all"))
bm_all@criterion # BIC est le critère de choix de modèle
bm_all = res.mixmod_all@bestResult
bm_all@model #"Gaussian_p_L_Dk_A_Dk"
bm_all@parameters@nbFreeParam  #91

sil_all = silhouette(bm_all@partition, dist(seeds))
mean(sil_all[,3])         # indice moyen 0.45
fviz_silhouette(sil_all)
mclust::adjustedRandIndex(bm_all@partition,variety) # 0.62

# ici, le centrage réduction ne change rien, car il est pris en compte dans le modèle
res.mixmod_cr = mixmodCluster(seedscr, 3, models=mixmodGaussianModel(family = "all"))
bm_cr = res.mixmod_cr@bestResult
sil_cr = silhouette(bm_cr@partition, dist(seeds))
mean(sil_cr[,3])         # indice moyen 0.45
fviz_silhouette(sil_cr)
mclust::adjustedRandIndex(bm_cr@partition,variety) # 0.62