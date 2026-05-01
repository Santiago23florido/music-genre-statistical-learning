##################
# STA203 - Project
##################
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(cluster)
library(MASS)

rm(list=objects())
graphics.off(); seeds=read.table("data/Music_2026.txt", header=TRUE, sep=";")
dir.create("output/part1", recursive=TRUE, showWarnings=FALSE)
dir.create("output/part1/hclust_pca", recursive=TRUE, showWarnings=FALSE)
#################
# PART I
#################

####
#Question 1
####
summary(seeds)        

#Genre Proportion Graphic Representation          
summary(as.factor(seeds$GENRE))
gp = ggplot(seeds) +
  aes(x=GENRE, fill=as.factor(GENRE)) +
  geom_bar(aes(y=after_stat(count/sum(count)))) +
  xlab("Music genre") +
  ylab("Proportion of tracks") +
  labs(fill="Music genre")
ggsave("output/part1/genre_proportions.png", gp)
prop.table(table(seeds$GENRE))

n = nrow(seeds)
p = ncol(seeds)

#Log-transformation diagnostics on positive variables
numeric_seeds = seeds[, sapply(seeds, is.numeric), drop=FALSE]
positive_seeds = numeric_seeds[, vapply(numeric_seeds, function(x) all(x > 0), logical(1)), drop=FALSE]

#Skewness calculation
skewness_pop = function(x){
  sx = sqrt(mean((x - mean(x))^2))
  if (sx == 0) return(0)
  mean((x - mean(x))^3)/(sx^3)
}

# Heteroscedasticity score based on the association between level and spread With LLM suggested implementation
# Source for the association measure: NIST/SEMATECH rank correlation documentation
# Source for the spread-versus-level diagnostic idea: NIST/SEMATECH spread-location plot
heteroscedasticity_score = function(x, n_bins=10){
  breaks = unique(quantile(x, probs=seq(0, 1, length.out=n_bins + 1), names=FALSE, type=8))
  if (length(breaks) < 4) return(0)
  groups = cut(x, breaks=breaks, include.lowest=TRUE, labels=FALSE)
  mean_by_group = tapply(x, groups, mean)
  sd_by_group = tapply(x, groups, sd)
  valid = is.finite(mean_by_group) & is.finite(sd_by_group)
  if (sum(valid) < 3) return(0)
  score = suppressWarnings(cor(mean_by_group[valid], sd_by_group[valid], method="spearman"))
  if (!is.finite(score)) return(0)
  score
}

#Normality test on log transformed function
shapiro_log_test = function(x, max_n=5000){
  y = log(x)
  if (length(y) > max_n){
    y = y[seq_len(max_n)]
  }
  test = shapiro.test(y)
  c(statistic=unname(test$statistic), p_value=test$p.value)
}

#Characteristic extraction for asymmetric identification
log_candidates = data.frame(
  variable = names(positive_seeds),
  mean_value = vapply(positive_seeds, mean, numeric(1)),
  median_value = vapply(positive_seeds, median, numeric(1)),
  skewness = vapply(positive_seeds, skewness_pop, numeric(1)),
  heteroscedasticity = vapply(positive_seeds, heteroscedasticity_score, numeric(1))
)

log_candidates$mean_median_ratio = log_candidates$mean_value / log_candidates$median_value

top_n = min(30, nrow(log_candidates))

top_ratio_variables = log_candidates$variable[order(log_candidates$mean_median_ratio, decreasing=TRUE)[1:top_n]]
top_skewness_variables = log_candidates$variable[order(log_candidates$skewness, decreasing=TRUE)[1:top_n]]
top_heteroscedasticity_variables = log_candidates$variable[order(log_candidates$heteroscedasticity, decreasing=TRUE)[1:top_n]]

#Data preparation for ploting for asymmetric identification
mean_median_plot_data = log_candidates[match(top_ratio_variables, log_candidates$variable), ]
mean_median_plot_data = mean_median_plot_data[order(mean_median_plot_data$mean_median_ratio, decreasing=TRUE), ]
skewness_plot_data = log_candidates[match(top_skewness_variables, log_candidates$variable), ]
skewness_plot_data = skewness_plot_data[order(skewness_plot_data$skewness, decreasing=TRUE), ]
heteroscedasticity_plot_data = log_candidates[match(top_heteroscedasticity_variables, log_candidates$variable), ]
heteroscedasticity_plot_data = heteroscedasticity_plot_data[order(heteroscedasticity_plot_data$heteroscedasticity, decreasing=TRUE), ]


#Ploting for asymmetric identification
mean_median_plot = ggplot(mean_median_plot_data) +
  aes(x=reorder(variable, mean_median_ratio), y=mean_median_ratio) +
  geom_col(fill="steelblue") +
  coord_flip() +
  scale_y_log10() +
  xlab("Positive variable") +
  ylab("Mean / median ratio (log scale)")
ggsave("output/part1/log_candidates_mean_median_ratio.png", mean_median_plot)


skewness_plot = ggplot(skewness_plot_data) +
  aes(x=reorder(variable, skewness), y=skewness) +
  geom_col(fill="firebrick") +
  coord_flip() +
  scale_y_log10() +
  xlab("Positive variable") +
  ylab("Skewness coefficient (log scale)")
ggsave("output/part1/log_candidates_skewness.png", skewness_plot)

heteroscedasticity_plot = ggplot(heteroscedasticity_plot_data) +
  aes(x=reorder(variable, heteroscedasticity), y=heteroscedasticity) +
  geom_col(fill="purple4") +
  coord_flip() +
  xlab("Positive variable") +
  ylab("Heteroscedasticity score")
ggsave("output/part1/log_candidates_heteroscedasticity.png", heteroscedasticity_plot)


#Candidates common to all the interest variables for asymmetric identification
common_candidates = Reduce(intersect, list(top_ratio_variables, top_skewness_variables, top_heteroscedasticity_variables))
common_candidates

#Normality test on log transformed
common_normality = t(vapply(positive_seeds[common_candidates], shapiro_log_test, numeric(2)))
common_normality = data.frame(
  variable = rownames(common_normality),
  shapiro_w = common_normality[, "statistic"],
  p_value = common_normality[, "p_value"],
  row.names = NULL
)

common_normality = common_normality[order(common_normality$shapiro_w, decreasing=TRUE), ]
common_normality

#Normality test on log transformed Ploting
common_normality_plot = ggplot(common_normality) +
  aes(x=reorder(variable, shapiro_w), y=shapiro_w) +
  geom_col(fill="darkgreen") +
  coord_flip() +
  xlab("Common candidate variable") +
  ylab("Shapiro-Wilk statistic on log(x)")
ggsave("output/part1/common_candidates_log_normality.png", common_normality_plot)

# Shared preprocessing for the retained variables
seeds_filtered = seeds
seeds_filtered[, c("PAR_SC_V", "PAR_ASC_V")] = log(seeds_filtered[, c("PAR_SC_V", "PAR_ASC_V")])
seeds_filtered = seeds_filtered[, -unique(which(abs(cor(seeds_filtered[,-p])) > 0.99 & upper.tri(cor(seeds_filtered[,-p])), arr.ind=TRUE)[,"col"])]
seeds_filtered = seeds_filtered[, !names(seeds_filtered) %in% c("PAR_ASE_M", "PAR_ASE_MV", "PAR_SFM_M", "PAR_SFM_MV")]
p_filtered = ncol(seeds_filtered)

# PCA-specific scaling
seeds_pca = seeds_filtered
seeds_pca[,-p_filtered] = scale(seeds_pca[,-p_filtered], center=TRUE, scale=TRUE)
p_pca = ncol(seeds_pca)


####
#Question 2
####

res = PCA(seeds_pca[,-p_pca],ncp=40, graph=FALSE)
#Proper values study

round(res$eig,4) # variance of each axe
sum(res$eig[,1])

#Inertie percetage of each component
inertie_percentage = ggplot() + aes(x=1:length(res$eig[,2]), y=res$eig[,2]) + geom_col() + 
  geom_hline(yintercept=100/p_pca, lty=2) +
  ggtitle("Percentage of Inertia Explained by Each Principal Component") +
  xlab("Principal component") +
  ylab("Explained inertia (%)")
ggsave("output/part1/inertie_percentage.png", inertie_percentage)


#Correlation of variables study
V = res$var
#Quality of representation
V$cos2

#First principal plane
#LLM suggestion for use of fviz_pca_ind to adequatly represent the genre in the  PCA graphs
plt1 = factoextra::fviz_pca_ind(res, axes = c(1,2), habillage = seeds_pca[,p_pca], label = "none")
plt2 = plot(res, axes = c(1,2), choix = "var")
pca_first_plane = cowplot::plot_grid(plt1, plt2, ncol = 2, nrow = 1)
ggsave("output/part1/pca_first_plane.png", pca_first_plane, width = 14, height = 6)

#Second principal plane
plt3 = factoextra::fviz_pca_ind(res, axes = c(3,4), habillage = seeds_pca[,p_pca], label = "none")
plt4 = plot(res, axes = c(3,4), choix = "var")
pca_second_plane = cowplot::plot_grid(plt3, plt4, ncol = 2, nrow = 1)
ggsave("output/part1/pca_second_plane.png", pca_second_plane, width = 14, height = 6)

#Inertias and quality for both planes
planes_inertia = sum(res$eig[1:4,2])
planes_cos2 = mean(rowSums(res$ind$cos2[,1:4]))

c( planes_inertia = planes_inertia, planes_cos2 = planes_cos2)

#which(res$eig[,2] > 100 / p)

#Candidates to select a  4 6 10 34

#Inertias and quality for Candidate variables
for (k in c(4, 6, 10, 34)) {
  cat(
    "k =", k,
    " | inertia =", sum(res$eig[1:k,2]),
    " | mean cos2 =", mean(rowSums(res$ind$cos2[, 1:k, drop=FALSE])),
    "\n"
  )
}

####
#Question 3
####

genre = as.integer(as.factor(seeds_pca[,p_pca]))
k_grp = nlevels(as.factor(seeds_pca[,p_pca]))

#Hclust estimation for all ACP resultant candidates 4 6 10 34
dir.create("output/part1/hclust_pca/not_normalized", recursive=TRUE, showWarnings=FALSE)
dir.create("output/part1/hclust_pca/normalized", recursive=TRUE, showWarnings=FALSE)

#Hclust estimation for all ACP resultant candidates 4 6 10 34
#ggplot +ggsave  rendering to expensive LLM suggestion to implement png + plot suggested
for (mode in c("not_normalized", "normalized")) {
  for (d in c(4, 6, 10, 34)) {
    X_acp = res$ind$coord[, 1:d, drop=FALSE]
    if (mode == "normalized") {
      X_acp = scale(X_acp, center=TRUE, scale=TRUE)
    }
    dd = dist(X_acp)
    
    hc = hclust(dd, method="ward.D2")
    
    #Plot the first heights 
    height_count = min(15, length(hc$height))
    height_profile = rev(hc$height)[1:height_count]
    
    #Select k by maximizing the average silhouette width indice de silhouette based on (Rousseeuw et al. 1987) LLM suggested implementation of the function
    k_candidates = 2:min(10, nrow(X_acp) - 1)
    mean_silhouette_by_k = sapply(k_candidates, function(k){
      grp_k = cutree(hc, k=k)
      mean(silhouette(grp_k, dd)[,3])
    })
    k_selected = k_candidates[which.max(mean_silhouette_by_k)]
    h_selected = rev(hc$height)[k_selected]
    
    grp = cutree(hc, k=k_selected)
    
    sil_ward = silhouette(grp, dd)
    sil_genre = silhouette(genre, dd)
    
    ward_mean = mean(sil_ward[,3])
    genre_mean = mean(sil_genre[,3])
    
    #Plot the first heights to support the choice of k
    png(sprintf("output/part1/hclust_pca/%s/heights_pca_%02d_components.png", mode, d), width=1200, height=700, res=150)
    barplot(height_profile, main=sprintf("Height diagram for Ward clustering on first %d PCs (%s)", d, mode))
    abline(h=h_selected, lty=2, col="red")
    dev.off()
    
    #Plot the mean silhouette criterion used to choose k
    png(sprintf("output/part1/hclust_pca/%s/mean_silhouette_pca_%02d_components.png", mode, d), width=1200, height=700, res=150)
    bar_mid = barplot(mean_silhouette_by_k, names.arg=k_candidates,
                      main=sprintf("Average silhouette width by k on first %d PCs (%s)", d, mode),
                      xlab="Number of clusters k", ylab="Average silhouette width")
    abline(v=bar_mid[which.max(mean_silhouette_by_k)], lty=2, col="red")
    dev.off()
    
    png(sprintf("output/part1/hclust_pca/%s/dendrogram_pca_%02d_components.png", mode, d), width=1200, height=700, res=150)
    plot(hc, cex=0.5, main=sprintf("Ward dendrogram on first %d PCs (%s)", d, mode))
    abline(h=h_selected, lty=2, col="red")
    rect.hclust(hc, k=k_selected, border=2:(k_selected+1))
    dev.off()
    
    #Build silhouettes for Ward clustering and true genre labels
    plt_sil_ward = fviz_silhouette(sil_ward) +
      ggtitle(sprintf("Ward silhouette on first %d PCs", d),
              subtitle=sprintf("%s | selected k = %d | mean = %.3f", mode, k_selected, ward_mean))
    plt_sil_genre = fviz_silhouette(sil_genre) +
      ggtitle(sprintf("GENRE silhouette on first %d PCs", d),
              subtitle=sprintf("%s | k = %d | mean = %.3f", mode, k_grp, genre_mean))
    plt_sil = cowplot::plot_grid(plt_sil_ward, plt_sil_genre, ncol=2, nrow=1)
    ggsave(sprintf("output/part1/hclust_pca/%s/silhouette_pca_%02d_components.png", mode, d), plt_sil, width=14, height=6)
    
    cat(
      mode,
      " | PCA dimensions =", d,
      " | selected k =", k_selected,
      " | mean Ward silhouette =", ward_mean,
      " | mean GENRE silhouette =", genre_mean,
      "\n"
    )
  }
}

####
#Question 4
####

set.seed(103)
train = sample(c(TRUE, FALSE), n, replace=TRUE, prob=c(2/3, 1/3))


#################
# PART II
#################

#Restrict to Classical and Jazz for binary classification
binary_genres = c("Classical", "Jazz")
binary_train = seeds_filtered[train & seeds_filtered[,p_filtered] %in% binary_genres, ]
binary_test = seeds_filtered[!train & seeds_filtered[,p_filtered] %in% binary_genres, ]

#Check the requested sample sizes
c(nrow(binary_train), nrow(binary_test))

####
#Question 1
####
#Directory for output storage
dir.create("output/part2/models", recursive=TRUE, showWarnings=FALSE)

binary_train$Y = factor(binary_train$GENRE, levels=c("Classical", "Jazz"))
binary_test$Y = factor(binary_test$GENRE, levels=c("Classical", "Jazz"))

# Full logistic model with all retained variables
train_ModT = binary_train[, !names(binary_train) %in% "GENRE"]
test_ModT = binary_test[, !names(binary_test) %in% "GENRE"]

ModT = glm(Y ~ ., family=binomial, data=train_ModT)
summary(ModT)

coef_ModT = summary(ModT)$coefficients

# Significant variables in ModT
vars_Mod1 = rownames(coef_ModT)[coef_ModT[,4] < 0.05]
vars_Mod1 = vars_Mod1[vars_Mod1 != "(Intercept)"]

vars_Mod2 = rownames(coef_ModT)[coef_ModT[,4] < 0.20]
vars_Mod2 = vars_Mod2[vars_Mod2 != "(Intercept)"]

# Reduced datasets for Mod1 and Mod2
train_Mod1 = binary_train[, c("Y", vars_Mod1), drop=FALSE]
test_Mod1 = binary_test[, c("Y", vars_Mod1), drop=FALSE]

train_Mod2 = binary_train[, c("Y", vars_Mod2), drop=FALSE]
test_Mod2 = binary_test[, c("Y", vars_Mod2), drop=FALSE]

# Reduced logistic models
Mod1 = glm(Y ~ ., family=binomial, data=train_Mod1)
Mod2 = glm(Y ~ ., family=binomial, data=train_Mod2)

summary(Mod1)
summary(Mod2)

# Stepwise AIC logistic model
ModAIC = stepAIC(ModT)

summary(ModAIC)

ModAIC$anova

# Save fitted models
saveRDS(ModT, "output/part2/models/ModT.rds")
saveRDS(Mod1, "output/part2/models/Mod1.rds")
saveRDS(Mod2, "output/part2/models/Mod2.rds")
saveRDS(ModAIC, "output/part2/models/ModAIC.rds")

# Save retained variable names
writeLines(
  c(
    paste("ModT :", paste(names(train_ModT)[names(train_ModT) != "Y"], collapse=" + ")),
    paste("Mod1 :", paste(vars_Mod1, collapse=" + ")),
    paste("Mod2 :", paste(vars_Mod2, collapse=" + ")),
    paste("ModAIC :", deparse(formula(ModAIC)))
  ),
  "output/part2/models/model_formulas.txt"
)

# Save summaries
capture.output(summary(ModT), file="output/part2/models/ModT_summary.txt")
capture.output(summary(Mod1), file="output/part2/models/Mod1_summary.txt")
capture.output(summary(Mod2), file="output/part2/models/Mod2_summary.txt")
capture.output(summary(ModAIC), file="output/part2/models/ModAIC_summary.txt")
capture.output(ModAIC$anova, file="output/part2/models/ModAIC_anova.txt")
