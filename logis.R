## STA3 - TP7- GHQ
#  setwd("mon répertoire")
rm(list=objects()) ; graphics.off()
library(ggplot2)

# lecture des données et étude descriptive sommaire
df = read.table("GHQ.data",header=TRUE)
summary(df)
str(df)

# on souhaite modéliser la probabilité d'être atteint de la maladie
# en fonction de deux variables (sex quantitative, ghq quatitative)
# Données groupées: K=17 groupes de taille nk=cases+noncases différentes
# d'où la création d'une variable taille
df$tot = df$noncases+df$cases

# une variable proportion observée
df$prop = df$cases/df$tot


###################
### 1- régression linéaire
###################
### P_k= a+b ghq_k + eps_k indep eps_k iid gaussien

res.lm = lm(prop~ghq,data=df);summary(res.lm)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.18938    0.08857   2.138   0.0494 *  
#  ghq        0.11233    0.01667   6.739 6.65e-06 ***
#Residual standard error: 0.2094 on 15 degrees of freedom
#Multiple R-squared: 0.7517,  Adjusted R-squared: 0.7352 
#F-statistic: 45.41 on 1 and 15 DF,  p-value: 6.648e-06 

# cette F-statistique est la stat de test de significativité globale
# H0: modèle iid (coeff associé à ghq=0, ie b=0) contre H1: modèle d'étude (b<>0)
# (SCR(iid) - SCR(modèle))/( SCR(modèle)/différence de ddl)

# comme c'est une régression simple, F est égale au carré de la stat de Student
# de nullité du coefficient b
(6.739)^2 # > qf(.95,1,15) donc le coefficient est significatif au risque (de première espèce 5%)

# ggplot (attention, si on met col et pch dans ggplot,  deux régressions sont faites)
plt = ggplot(data=df)+ aes(x=ghq, y=prop)+
   geom_point(aes(x=ghq, y=prop, col=sex, pch=sex), size=3)+
   stat_smooth(method="lm",formula=y~x, se=FALSE,col=1)+
  labs(y="proportion de cas")
plt
# la probabilité d'être atteint alors qu'on a un score > 6 est >1

# on voit un (léger) pattern dans les résidus
library(MASS)
plot(df$ghq,studres(res.lm)) # motif dans les résidus

par(mfrow=c(2,2))            # d'autres représentations
plot(res.lm)
par(mfrow=c(1,1))

# glm fait aussi du lm...
summary(glm(prop~ghq,data=df))  
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    0.18938    0.08857   2.138   0.0494 *  
#   ghq          0.11233    0.01667   6.739 6.65e-06 ***
# (Dispersion parameter for Gaussian family taken to be 0.04386396) => RSE^2
# 
# Null deviance: 2.64999  on 16  degrees of freedom
# Residual deviance: 0.65796  on 15  degrees of freedom

# on fait le lien entre les sorties de lm 
# et celles de glm (par défaut, fait de la régression linéaire)

# le param de dispersion (glm) = RSE^2 où RSE est une sortie de lm
c(0.65796/15)                           # 0.043864
# et sa racine carré (RSE), estimation de sigma, l'écart-type du bruit
sqrt(c(0.65796/15, 0.0438639) )         # 0.2094372  
# remarquons que le calcul (Null deviance - Residual deviance)/(Residual deviance/15) 
# permet de retrouver la stat de significativité globale de la régression
# 
(2.64999 - 0.65796)/( 0.65796  / 15)    # 45.41382  => F stat


###################
### 2- régression logistique en données groupées
###################
### cases_k~B(noncases, pi_k) indep avec logit(pi_k)= a+b ghq_k

res.glm = glm(cbind(cases,noncases)~ghq,
            family=binomial,data=df)
summary(res.glm)
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.4536     0.5633  -6.131 8.73e-10 ***
# ghq           1.4402     0.2979   4.834 1.34e-06 ***
#   ---
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 83.176  on 16  degrees of freedom
# Residual deviance:  5.744  on 15  degrees of freedom
# AIC: 22.22
# 
# Number of Fisher Scoring iterations: 6

###################
### 3- Représentation graphique
###################
  

# tracé des prédictions
# construire une fonction
f = function(x,a,b)    1/(1+exp(-a-b*x))  
par(mfrow=c(1,1))
plt = plt + stat_function(fun=function(x)f(x,a=res.glm$coef[1],b= res.glm$coef[2]),
                          col="red",lwd=1, lty=2)
plt

# la loi binomiale est évidemment adaptée pour modéliser un phénomène binaire.
# le modèle de régression logistique permet de modéliser une probabilité
# comprise entre 0 et 1, dépend de covariables explicatives.
# visuellement, l'ajustement est meilleur.

###################
### 4- ICs
###################
# loi asymptotique de l'EMV
# (vcov)^{-1/2} (thetachap-theta) -> loi N(0,Id)
# l'approximpation gaussienne est utilisée à distance finie
alpha = 0.05
q = qnorm(1-alpha/2)  

# IC de l'intercept
a = res.glm$coef[1]
s.est = sqrt(vcov(res.glm)[1,1])
ICint = c(est=a,min=a-q*s.est,max=a+q*s.est)  #l'IC
# est.(Intercept) min.(Intercept) max.(Intercept) 
#       -3.453574       -4.557608       -2.349539

# IC de la "pente", coefficient associé à ghq
b = res.glm$coef[2]
s.est = sqrt(vcov(res.glm)[2,2])
ICghq = c(est=b,min=b-q*s.est,max=b+q*s.est)  #l'IC
#   est.ghq   min.ghq   max.ghq 
# 1.4402323 0.8562711 2.0241935 


confint(res.glm)
# avec confint lesintervalles de confiance des composantes du paramètre
# sont calculés à partir de la vraisemblance profilée
# d'où la légère différence de résultat 
#  (Intercept) -4.7369172 -2.485148
#  ghq          0.9406541  2.124828


###################
### 5- IC de l'espérance de la réponse
###################
## cas ghq=0: on prend simplement la fonction de lien inverse des bornes de l'IC calculé à Q4
# delta méthode inutile ici car cas d'une fonction strictement monotone

1/(1+exp(-ICint))

# Rem:  fonction inverse est déjà codée! 
binomial()$linkinv(eta=ICint) 

## pour ghq = 1:5, 
# étape 1 : on calcule l'IC du régresseur linéaire
alpha = 0.05; q = qnorm(1-alpha/2)

# à la main
A = cbind(rep(1,5),1:5)
c = res.glm$coef
est = A%*%c
vest = diag(A%*%vcov(res.glm)%*%t(A))
s.est = sqrt(vest)
# 0.3780032 0.3821204 0.5715661 0.8275880 1.1048341

ICxtheta = data.frame (est=est, min=est-q*s.est, max=est+q*s.est)
#      est        min        max
#1 -2.0133413 -2.7542140 -1.2724687
#2 -0.5731090 -1.3220512  0.1758331
#3  0.8671233 -0.2531257  1.9873723
#4  2.3073556  0.6853128  3.9293983
#5  3.7475879  1.5821528  5.9130229

# étape 2: et on prend la fonction de lien inverse pour ramener dans l'échelle de la probabilité
pred = 1/(1+exp(- ICxtheta)) 
#       est       min       max
#1 0.1178093 0.0598491 0.2188350
#2 0.3605197 0.2104772 0.5438454
#3 0.7041468 0.4370543 0.8794649
#4 0.9094844 0.6649234 0.9807234
#5 0.9769684 0.8295092 0.9973033


# la représentation graphique avec ggplot
# la commande suivante ne fonctionne, pas, à cause de data=df dans ggplot
plt + geom_segment(aes(x=1:5, y=pred$min, xend=1:5, yend=pred$max),
                   col="purple",lty=2)

# alternative 1 = ne pas mettre data dans ggplot
ggplot()+ 
  geom_point(aes(x=df$ghq, y=df$prop, col=df$sex, pch=df$sex), size=3)+
  stat_smooth(method="lm",se=TRUE)+
  stat_function(fun=function(x)f(x,a=res.glm$coef[1],b= res.glm$coef[2]),
                col="red", lwd=1, lty=2)+
  geom_segment(aes(x=1:5, y=pred$min, xend=1:5, yend=pred$max),col="purple")

# alternative 2 voir plus loin

###################
### 6- avec la commande predict
###################
pred.link = predict(res.glm,newdata=data.frame(ghq=1:5),type="link", se.fit=TRUE)
cbind(s.est, pred.link$se.fit)   # idem !                                 # std err du régresseur linéaire
data.frame(min = 1/(1+exp(-pred.link$fit+q*pred.link$se.fit)), 
           max = 1/(1+exp(-pred.link$fit-q*pred.link$se.fit)) )  # idem

# predict permet aussi d'exprimer dans l'échelle de la proba
pred.prob = predict(res.glm,newdata=data.frame(ghq=1:5),type="response", se.fit=TRUE) # std err de la delta méthode

cbind(pred.prob$fit, pred$est)  # OK pour l'estimation ponctuelle
cbind(pred.prob$se.fit, s.est *pred.prob$fit*(1-pred.prob$fit)) # OK pour std err par delta methode

# et alternative 2 pour la visu: on fait un nouveau df
# pred contient déjà l'IC des prédictions faites en utilisant la fonction de lien inverse sur les bornes
pred$ghq = 1:5

# on peut comparer avec la méthode delta 
pred$minDelta= pred.prob$fit-q* pred.prob$se.fit
pred$maxDelta= pred.prob$fit+q* pred.prob$se.fit
plt + 
   geom_segment(data=pred,aes(x=ghq, y=max, xend=ghq, yend=min),col="red",lwd=1) +
   geom_segment(data=pred,aes(x=ghq, y=maxDelta, xend=ghq, yend=minDelta),col="black",lty=2,lwd=1.1)+
   stat_smooth(data=df,method="glm",se=TRUE, method.args = list(family=binomial())) 

# la visu des IC n'est pas bonne...
# et d'ailleurs, il y a un pb avec stat_smooth
# Warning message:
#    In eval(family$initialize) : non-integer #successes in a binomial glm!

# c'est parce que  stat_smooth ne fait pas le glm qu'on croit...
# stat_smooth fait la régression suivante
# (il considère 17 observations individuelles au lieu de sum_k n_k=120 en tout)
# donc un rapport de std error de l'ordre de sqrt(120/17) = 2.65

res.glm2 = glm(prop~ghq, family=binomial,data=df)
summary(res.glm2)
#             Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -3.4757     2.0117  -1.728   0.0840 . # std err bcp plus grande donc IC + larges
# ghq           1.4668     0.7742   1.895   0.0581 .

summary(res.glm)
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.4536     0.5633  -6.131 8.73e-10 ***
# ghq           1.4402     0.2979   4.834 1.34e-06 ***
   
# on peut arranger cela en mettant un poids sur les observations
res.glm3 = glm(prop~ghq, family=binomial,weights=tot,data=df)
summary(res.glm3)   # idem res.glm

# les bons ICs
plt + 
   geom_segment(data=pred,aes(x=ghq, y=max, xend=ghq, yend=min),color="red",lwd=1) +
   geom_segment(data=pred,aes(x=ghq, y=maxDelta, xend=ghq, yend=minDelta),color="black",lty=2,lwd=1.1)+
   stat_smooth(data=df,aes(x=ghq,y=prop, weight=tot),method="glm",se=TRUE, 
               method.args = list(family=binomial()),fill="blue" )


########################### pour les curieu.ses
# alternative 3 
# on crée un fonction qui calcule le min
fICm = function(x){
   est = c(1,x)%*% res.glm$coef 
   s.est = sqrt(diag(c(1,x)%*%vcov(res.glm)%*%c(1,x)))
   min=est-q*s.est
   1/(1+exp(-min))
} 
# donc il faut vectoriser l'argument pour l'appel dans segments
FICm= Vectorize(fICm)
FICm(1:5)

# idem pour borne sup
fICM = function(x){
   1/(1+exp(- c(1,x)%*% res.glm$coef -q*sqrt(diag(c(1,x)%*%vcov(res.glm)%*%c(1,x)))  ))
} 
FICM= Vectorize(fICM)
FICM(1:5)
df$min=FICm(df$ghq)
df$max=FICM(df$ghq)

plt = plt +   geom_segment(aes(x=ghq, y=FICm(ghq),xend=ghq, yend=FICM(ghq)),col="violet")

# et tant qu'à faire, les courbes d'IC
plt = plt +
   stat_function(fun=FICm, col="violet",  lty=2)+
   stat_function(fun=FICM, col="violet",  lty=2)

# et un remplissage
plt + geom_ribbon(aes(ymin=FICm(ghq), ymax=FICM(ghq)),fill="lightblue", alpha=0.7) 


###### avec plot
# lin
col = as.numeric(factor(df$sex))

plot(df$ghq,df$prop,main="données GHQ",xlab="score ghq", ylab="proportion de cas",
     col=col,pch=col)
abline(res.lm)   # ajustement linéaire

# logistique
plot(df$ghq,df$prop,main="données GHQ",xlab="score ghq", ylab="proportion de cas")
abline(res.lm) 
curve(f(x,res.glm$coef[1], res.glm$coef[2]),from=0, to=10,col=2,add=TRUE)
legend("bottomright", c("observations","lm","glm"),
       col=c(1,1,2),lty=0:2,pch=c(1,-1,-1))

#IC
plot(df$ghq,df$prop,main="données GHQ")
curve(predict(res.glm,newdata=data.frame(ghq=x),type="response"),lty=2,col=2,add=TRUE,
      from=0, to=10)
segments(1:5,pred$min,1:5,pred$max,col=2) 
segments(1:5, pred.prob$fit-q* pred.prob$se.fit, 
         1:5, pred.prob$fit+q* pred.prob$se.fit ,col="black", lty=2,lwd=2)
legend("bottomright", c("observation","prediction", "IC", "IC par delta"),
       col=c("black","red","red","black"),lty=c(0,2,1,2),pch=c(1,-1,-1,-1),cex=0.7)