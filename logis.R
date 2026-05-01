#STA1
setwd("mon_repertoire")
library(ggplot2)
rm(list=objects());graphics.off()

########################## Détecteur de spam  ###################
#### TREVOR et HASTIE http://www-stat.stanford.edu/~tibs/ElemStatLearn/

#commandes du TP précedent

spam = read.table("spam.data")
dim(spam) # 4601 x 58
p = ncol(spam)-1          # nombre de variables explicatives
n = nrow(spam)
names(spam)[p+1] = "Y"

###################################
### Q1 - classification 
###################################

# seuil de la règle de Bayes = 0.5
res = glm(Y~.,data=spam,family=binomial)         # estimation
predproba = predict(res,type="response")         # prédiction en proba

glm.pred = ifelse(predproba>.5,1,0)

# Rem
sum(abs(glm.pred-as.numeric(predproba>.5)))

# table de classification
table(Y=spam$Y , glm.pred )
1-mean(glm.pred==spam$Y) #0.06868072= taux de mauvais classement
mean(glm.pred!=spam$Y)

# on peut tracer les faux positifs et faux négatifs en fonction du seuil
seuil = seq(0,1,0.01)

# nombre de faux positifs et faux négatifs en fonction du seuil
# sapply pour transformer la liste en vecteur
nbspam = sum(spam$Y)   # nombre de spams
FP = sapply(seuil,function(s)  sum( predproba>=s & spam$Y==0))
FN = sapply(seuil,function(s)  sum( predproba< s & spam$Y==1))
VP = nbspam - FN
VN = (n-nbspam) - FP

## alternative:
# on peut calculer à partir de la table de classification pour différents seuils
# mais c'est un peu plus compliqué car table ne renvoie que deux valeurs pour s=0 et s=1
table(Y=spam$Y,pred=predproba>0)
sap = sapply(seq(0.01,0.99,0.01),
             function(x){c(table(Y=spam$Y,pred=predproba>x))})

# ajout des cas spéciaux s=0 et s=1 
# (pb avec table qui ne renvoie alors que deux valeurs)

tc = data.frame(seuil,rbind(c(0,0,n-nbspam, nbspam), # s = 0
                            t(sap),  # VN, FN, FP, VP pour les seuils de 0.01 à 0.99
                            c( n-nbspam,nbspam,0,0))) # s = 1
names(tc)[2:5] = c("VN","FN","FP","VP")
FP = tc$FP; FN = tc$FN; VP = tc$VP; VN = tc$VN

## Fin de l'alternative

# le plot
matplot(seuil,cbind(FN,FP,FN+FP),
        type="l",ylab="erreur",
        lty=c(2,3,1),col=c("red","brown","black"))

legend("top",c("FN","FP","totale"),lty=c(2,3,1),col=c("red","brown","black"))


# avec ggplot2 pour gérer des séries chronologiques
library(zoo)
plt = autoplot(zoo(data.frame(FN,FP,tot=FN+FP),
                   order.by=seuil),
         facet = NULL) +  geom_line()

#le seuil minimum (estimé)
wm = which.min(FN+FP)                # index dans seuil
me = (FN+FP)[wm]                     # erreur min 
c(seuil[wm], (FN+FP)[wm]/n)          # s=0.43 mean err=0.065

plt + geom_point(aes(x=seuil[wm],y=(FN+FP)[wm]),  pch=19,col=2)
# le seuil théorique est à 0.5. 
# noter la forme de la courbe autour entre 0.4 et 0.6

VP[wm]/nbspam      # sensibilité = 0.916 on arrête 92% de spams
FP[wm]/(n-nbspam)  # 1-specificité = 5% des courriers ordinaires passent en spam

###################################
### Q2- ROC 
###################################

## à la main

#seuil s
seuil = seq(0,1,.01)

# on recalcule (pour bien comprendre 1-specificité et sensibilité)
cv_ROC= sapply(seuil,function(s)  {c(
        #1-specificité = FPR = FP/NEG=taux faux positifs
  x_ROC = sum( predproba>=s & spam$Y==0)/sum(spam$Y==0), 
        #sensibilité = TPR= VP/POS = Taux de vrais positifs 
  y_ROC = sum( predproba>=s & spam$Y==1)/sum(spam$Y==1))})

plot(cv_ROC[1,],cv_ROC[2,],type="l")
lines(c(0,1),c(0,1),lty=2)          # règle aléatoire

# on peut aussi utiliser les calculs de la question précédente
x_ROC = FP / (n-nbspam) #1-specificité  
y_ROC = VP / nbspam     #sensibilité

sum(x_ROC!=cv_ROC[1,])

### en utilisant library(ROCR)
library(ROCR)
pred = prediction(predproba,spam$Y)  

# sort un objet de type prediction (qui calcul les TP, FN, etc en fonction du seuil)
#objet S4, syntaxe particulière pour accéder aux composantes
class(pred) 
mode(pred) # "S4"
slotNames(pred)
# [1] "predictions" "labels"      "cutoffs"     "fp"          "tp"          "tn"          "fn"         
# [8] "n.pos"       "n.neg"       "n.pos.pred"  "n.neg.pred" 

# accéder aux FP
slot(pred,"fp")       # les informations à disposition
class(pred@fp)        # c'est une liste qui contient un vecteur des prédictions
head(unlist(pred@fp)) # le vecteur des faux positifs

### courbe ROC
ROC = performance(pred,"sens","fpr")  # prépare les infos pour la courbe ROC
plot(ROC, xlab="", main="courbes ROC")               
lines(c(0,1),c(0,1),lty=2)          # règle aléatoire 

#  lines(cv_ROC[1,],cv_ROC[2,],type="l",col=2) 
# toutes petites différences, dues à la discrétisation du seuil

### AUC
perf = performance(pred, "auc")
slotNames(perf)
(AUC=round(unlist(perf@y.values),4) ) #AUC = 0.9774


###################################
###  blocage de 95% des spams

a = min(which(unlist(ROC@y.values) >= 0.95)) # numéro de l'obs de sensibilité>95%: 1867
s = unlist(ROC@alpha.values)[a]              # seuil= s =0.2793256
ROC@x.values[[1]][a]                      # tx de faux positif :  1 - beta = 0.108

pred95 = ifelse(predproba>s,1,0)
table(Y=spam$Y, pred95 = pred95)
#    pred95
# Y      0    1
#   0 2486  302
#   1   91 1722

# on vérifie :
1722/(1722+91)  # [1] 0.9498069 95% des spams bloqués
302/(302+2486)  # 0.1083214: faux positifs
mean(pred95==spam$Y) #0.9145838 et une erreur totale plus forte

###################################
###  Q3 modèle plus simple avec les variables les plus significatives
###################################

res2 = glm(Y~V5+V6+V7+V8+V16+V17+V21+V23+V25+V27+V45+V46+V52+V53+V56+V57,
           family=binomial,data=spam)

glm.pred2 = ifelse(predict(res2,type="response")>.5,1,0)

table(Y=spam$Y, glm.pred2 )
mean(glm.pred2!=spam$Y) # 0.08367746 # légèrement moins bon

# modèle M2: erreur supérieure, mais attention, il y a moins de variables
# on ne peut donc prendre une décision définitive (il faudrait avoir un échantillon test)

pred2 = prediction(predict(res2,type="response"),spam$Y)
ROC2 = performance(pred2,"sens","fpr")
plot(ROC2,xlab="", main="courbes ROC", col = 2,
     add =TRUE)     # pour superposer au graphe précédent           

(AUC2 =  round (performance(pred2, 'auc')@y.values[[1]] , 4) ) # 0.9689


####################
### avec LDA
####################
library(MASS)
res.lda = lda(Y ~ ., spam)

#proba
res.lda$prior       # on retrouve l'estimation des proportions

# aggregate calcule des statistiques par niveau [première colonne: groupe]
mb_k = sapply(spam[,-(p+1)], function(x) tapply(x,spam$Y, mean) )
mb_k

res.lda$means       # on retrouve l'estimation des moyennes

# variance totale
var(spam-mb_k[spam$Y,])

# variance par groupe
V_k = by(as.matrix(spam[,-(p+1)]),list(as.factor(spam$Y)),cov) 

# performance (sur l'apprentissage)
1-mean(spam$Y==predict(res.lda, spam)$class)  # 0.111
pred.lda = predict(res.lda, spam)$posterior[,2]  # pour p(Y=1|x)

(AUC.lda =  round (performance(prediction(pred.lda,spam$Y), 'auc')@y.values[[1]] , 4) ) # 0.9543
ROC.lda = performance(prediction(pred.lda,spam$Y),"sens","fpr")
plot(ROC.lda,xlab="", main="courbes ROC", col = 3,  add =TRUE) 

####################
### avec QDA
####################
res.qda = qda(Y ~ ., spam)

1-mean(spam$Y==predict(res.qda, spam)$class)  # 0.1671376

pred.qda = predict(res.qda, spam)$posterior[,2]
(AUC.qda =  round (performance(prediction(pred.qda,spam$Y), 'auc')@y.values[[1]] , 4) ) #  0.9509

ROC.qda = performance(prediction(pred.qda,spam$Y),"sens","fpr")
plot(ROC.qda,xlab="", main="courbes ROC", col = 4,  add =TRUE) 

legend("bottomright",legend= c(paste("AUC M: ",AUC),
                               paste("AUC M2: ",AUC2),
                               paste("AUC lda: ",AUC.lda),
                               paste("AUC qda: ",AUC.qda)), 
       col=1:4, lty=1)

####################
### avec ggplot
####################

### Q1 autres méthodes
ggplot(data=data.frame(FN,FP), aes(x=seuil))+
  geom_line(aes(y=FN),col="red",lty=2)+
  geom_line(aes(y=FP),col="orange",lty=3)+
  geom_line(aes(y=FP+FN),col="black",lty=1)+
  geom_point(aes(x=seuil[wm], y=me),col="red")

# si on veut une légende, c'est plus compliqué... il faut faire un df spécifique
# on peut aussi utiliser la fonction melt
tcforggplot = data.frame(seuil=rep(seuil,3), 
                         err = c(FN,FP,FN+FP),
                         type= rep(c("FN","FP","FN+FP"),rep(length(seuil),3)) )
ggplot(data=tcforggplot, aes(x=seuil, y=err, color=type, linetype=type))+ 
  geom_line()+
  scale_color_manual(name = "erreur",
                     values = c("red","black","orange"),
                     labels = c("faux neg", "total", "faux pos"))+
  scale_linetype_manual(name = "erreur",
                        values = c(2,1,3),
                        labels = c("faux neg", "total", "faux pos"))+
  geom_point(aes(x=seuil[wm], y=me),col="red")+
  theme(legend.position="top")

### Q2 alternatives
# se créer le df qui va bien
df = data.frame( FPR = c(unlist(ROC@x.values),unlist(ROC2@x.values),
                         unlist(ROC.lda@x.values),unlist(ROC.qda@x.values),
                         unlist(ROC@x.values)),
                 TPR = c(unlist(ROC@y.values),unlist(ROC2@y.values),
                         unlist(ROC.lda@y.values),unlist(ROC.qda@y.values),
                         unlist(ROC@x.values)),
                type = factor(rep(1:5,c(length(unlist(ROC@y.values)),
                                 length(unlist(ROC2@x.values)), 
                                 length(unlist(ROC.lda@y.values)),
                                 length(unlist(ROC.qda@y.values)),
                                 length(unlist(ROC@y.values)))))
)
ggplot(df,aes(FPR,TPR,color=type,linetype=type))+geom_line()+
  labs(title= "ROC curve", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)" )+
  scale_color_manual(name = "",values=1:5,labels=c("M", "M2", "lda","qda","aléa"))+
  scale_linetype_manual(name = "",values=1:5,labels=c("M", "M2", "lda","qda","aléa"))