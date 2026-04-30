##################
# STA203 - Project
##################
library(ggplot2)
library(FactoMineR)
library(factoextra)

rm(list=objects())
graphics.off(); seeds=read.table("data/Music_2026.txt", header=TRUE, sep=";")
dir.create("output/part1", recursive=TRUE, showWarnings=FALSE)
#################
# PART I
#################

#Question 1
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

ncol(seeds)

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

# PCA
#Apply log and center and scale data
seeds[, c("PAR_SC_V", "PAR_ASC_V")] = log(seeds[, c("PAR_SC_V", "PAR_ASC_V")])
seeds[, sapply(seeds, is.numeric)] = scale(seeds[, sapply(seeds, is.numeric)], center=TRUE, scale=TRUE)

n = nrow(seeds)
p = ncol(seeds)

res = PCA(seeds[,-p], graph=FALSE)

round(res$eig,4) # variance of each axe
sum(res$eig[,1])

#Inertie percetage of each component
inertie_percentage = ggplot() + aes(x=1:length(res$eig[,2]), y=res$eig[,2]) + geom_col() + 
  geom_hline(yintercept=100/p, lty=2) +
  ggtitle("Percentage of Inertia Explained by Each Principal Component") +
  xlab("Principal component") +
  ylab("Explained inertia (%)")
ggsave("output/part1/inertie_percentage.png", inertie_percentage)

