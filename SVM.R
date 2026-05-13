# STA3- C. Keribin
# TP12- SVM 
###################
library(ggplot2)
library(e1071)

rm(list=objects()); graphics.off()
# setwd("Mon_Repertoire")

##########################
### 1- Support Vector classifier: un cas séparable
#########################

df = readRDS("SVM-ex1.RDS")
plot(df$x1,df$x2,col=as.numeric(df$y),pch=as.numeric(df$y))
# GGally::ggpairs(df,aes(color=y))

###
### 1-2
###
res.svc = svm(y~., data=df , kernel ="linear", scale =FALSE)# 
summary(res.svc)
c = coef(res.svc)
# (Intercept)          x1          x2 
# 1.7114266  -0.3296334  -0.8325440

res.svc$index            # index des vecteurs supports
# 2 14 15

res.svc$SV               # coordonnées des vecteurs supports
df[res.svc$index,]       # on vérifie

res.svc$coefs            #  coefficients times the training labels: alpha_i y_i
sum(res.svc$coefs)       # 0 , OK !
res.svc$rho              # rho est The negative intercept = -1.711427

# on peut retrouver les coefficients de la droite de séparation
c           
# (Intercept)          x1          x2 
#   1.7114266  -0.3296334  -0.8325440

c(b0 =-res.svc$rho,                                # 1.7114266         
  b1 = sum(res.svc$coefs*df$x1[res.svc$index]),   # -0.3296334 OK !
  b2 = sum(res.svc$coefs*df$x2[res.svc$index]) )   # -0.8325440 OK !


###
### 1-3 plot
###

plsvm = function(res,df,marge=FALSE){
  n = nrow(df)
  pch = c("o","s")[1:n %in% res$index+1]
  c = coef(res)

  gg = ggplot(df,aes(x1, x2, color=y)) +
    geom_point(aes(pch=pch),size=3) + 
    geom_abline(intercept=-c[1]/c[3], slope=-c[2]/c[3]) 
  
  if (marge) gg = gg +
    geom_abline(intercept=-(c[1] + 1)/c[3], slope=-c[2]/c[3], linetype=2)+
    geom_abline(intercept=-(c[1] - 1)/c[3], slope=-c[2]/c[3], linetype=2)
  gg
  
  # plot(df[,1:2],col=as.integer(df$y) , pch=pch)
  # abline(-c[1]/c[3], -c[2]/c[3], col = "red")
  # abline(-(c[1] + 1)/c[3], -c[2]/c[3], col = "blue",lty=2)
  # abline(-(c[1] - 1)/c[3], -c[2]/c[3], col = "blue",lty=2)
  # 
}

plsvm(res.svc,df,marge=TRUE)
plot(res.svc, df)  # attention, échange x1 et x2 !

##########################
### 2- SVC cas non séparé
#########################

df = readRDS("SVM-ex2.RDS")

res.svc = svm(y~., data=df , kernel ="linear", cost=10, scale =FALSE) 
summary(res.svc)
#   SVM-Type:  C-classification 
# SVM-Kernel:  linear 
# cost:  10 
# 
# Number of Support Vectors:  7
# ( 4 3 )

# plot(res.svc,df)
plsvm(res.svc,df,marge=TRUE)

###
### 2-2 classification
###
cbind(res.svc$decision.values,                     # fonction de décision yf(x)
      as.matrix(cbind(1,df[1:2]))%*% coef(res.svc))# OK !

data.frame(fit=res.svc$fitted, dv =  res.svc$decision.values[,1])

# mais pourquoi si res.svc$decision.values <0 alors +1 ?

coeff(res.svc)
# (Intercept)          x1          x2 
#    1.506707   -1.555547   -1.611534

# coef(res.svc) donne un vecteur orthogonal
# au plan séparateur orienté vers les valeurs négatives
# donc valeur de la fonction de décision est cohérente

# ce choix dépend de la première valeur lue dans le fichier: ici -1 
# donc la frontière de décicion va se référer à ce choix, cf l'entête  $decision.values

head(res.svc$decision.values)
#           -1/1     #  la fonction de décision est calculée pour -1 (par rapport à 1)
# 1  1.000221752
# 2 -0.039398515
# 3  2.686402236
# 4  2.231079764

# ce choix est explicité dans la sortie labels 
res.svc$labels   # 1 2 => -1 / 1 car as.numeric (factor(c(-1,1))) # 1 2

# donc, si on veut retrouver les prédictions
data.frame(fitted = res.svc$fitted,  dv = ifelse(res.svc$decision.values<  0, 1, -1) )

# matrice de confusion sur l'échantillon d'apprentissage
table(pred=predict(res.svc), vrai=df$y)   # 3 mal classées
#       vrai
# pred -1 1
#   -1  8 1
#    1  2 9

##
## Rem: 
## si on met les +1 en premier dans le fichier
##

dfbis = df[c(11:20,1:10),]
res.svc.bis = svm(y~., data=dfbis , kernel ="linear", cost=10, scale =FALSE) 
summary(res.svc.bis)
# ( 3 4 )     # mais ceci est défini par rapport à labels (et non par rapport à levels)
# Levels: 
#   -1 1
res.svc.bis$labels  # 2 1 (soit 1 -1)

coef(res.svc.bis)   
# (Intercept)          x1          x2 
# -1.507109    1.555244    1.612603 

# le vecteur orthogonal a bien changé de direction.
# valeurs opposées de la fonction de décision car la référence a changé
res.svc.bis$decision.values
#            1/-1    # donc si decision.values>0 alors 1
# 1  -0.999451801
# 2   0.039776602
# 3  -2.686471751
# 4  -2.234092914

## pour les curieu.se.s
# voir le commentaire dans le code source svm.cpp de cette méthode 
# // Labels are ordered by their first occurrence in the training set. 

###
### 2-3 choix du coût
###

### quand le coût/poids d'une erreur diminue, le nombre de vecteurs support augmente

res.svc2 = svm(y~., data=df , kernel ="linear", cost = 0.1, scale =FALSE)
summary(res.svc2) # Number of Support Vectors:  16 ( 8 8 )
table(pred=predict(res.svc2), vrai=df$y)  # ici l'erreur d'apprentissage est meilleure
plsvm(res.svc2,df,marge=TRUE)

res.svc2 = svm(y~., data=df , kernel ="linear", cost = 0.01, scale =FALSE)
summary(res.svc2) # Number of Support Vectors:  20 ( 10 10 )  tout le jeu de données !
table(pred=predict(res.svc2), vrai=df$y)  # ici l'erreur d'apprentissage est moins bonne
plsvm(res.svc2,df)

# mais attention à un éventuel sur-apprentissage

### choix du paramètre par VC
set.seed (1)
tune.out = tune(svm ,y~.,data=df ,kernel ="linear",
                ranges =list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100) ))
summary(tune.out)

# le modèle sélectionné cost = 0.1
bestmod = tune.out$best.model
summary(bestmod)    # ( 8  8 )
bestmod$cost        # 0.1

# performance en test
dftest = readRDS("SVM-ex2-test.RDS")

table( vrai= dftest$y, pred = predict (bestmod ,dftest)) # 3 mal classifiés

# et sur un autre paramètre, c'est moins bien
res.svc0.01 = svm(y~., data=df , kernel ="linear", cost = 0.01, scale =FALSE)
ypred0.01 = predict (res.svc0.01 ,dftest)  
table( vrai= dftest$y, pred = ypred0.01)  # 6 mal classifiés


###
### courbe ROC
###
library(ROCR)
best.svc = svm(y~., data=df , kernel ="linear", cost = bestmod$cost, scale =FALSE, 
              probability=TRUE)
plsvm(best.svc,df,TRUE)

pred.svc = predict(best.svc,newdata=dftest, probability=TRUE)

roc.svc = attributes(pred.svc)$probabilities[,2]

plot( performance(prediction(roc.svc,dftest$y),"tpr","fpr") )


#######################
### SVM non linéaire
#######################
df = readRDS("SVM-ex3.RDS")
plot(df$x1,df$x2,col=as.numeric(df$y),pch=19)

# on sépare en apprentissage et test
set.seed(2022)
train =  sample(200,100)
res.svm = svm(y~., data=df [train ,], kernel ="radial", gamma =1, cost =1)

plot(res.svm , df[train ,])


summary(res.svm)
table( vrai= df$y[train], pred = res.svm$fitted) # peu d'erreur d'apprentissage

# si le coût est plus grand, 
# la frontière est plus irrégulière, et l'erreur d'apprentissage est plus faible
res.svm = svm(y~., data=df [train ,], kernel ="radial", gamma =1, cost =1e5)

plot(res.svm , df[train ,])

summary(res.svm)
table( vrai= df$y[train], pred = res.svm$fitted)  # 0

#en CV pour le meilleur choix de gamma et cout
set.seed (1)
tune.out = tune(svm , y~., data=df[train ,], kernel ="radial",
                ranges =list(cost=c(0.1 ,1 ,10 ,100 ,1000),
                             gamma=c(0.5,1,2,3,4) ))
summary(tune.out)

# performance du meilleur modèle sur le test
bestmod = tune.out$best.model
best.svm = svm(y~., data=df[train,], kernel ="radial", 
               cost = bestmod$cost,gamma=bestmod$gamma,
               probability=TRUE)
pred.svm = predict(best.svm,newdata=df[-train,], probability=TRUE)
 
table( vrai= df$y[-train], pred = pred.svm)   # 11

roc.svm = attributes(pred.svm)$probabilities[,1]
plot( performance(prediction(roc.svm,df$y[-train]),"tpr","fpr") )
