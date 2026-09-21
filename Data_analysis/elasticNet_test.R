# https://glmnet.stanford.edu/articles/glmnet.html
library(glmnet)

PP_mass_cluster_freq_1 # cluster freq Mass 2001
PP_mass_cluster_freq_3 # cluster freq Mass 2007
ggCPP_intermed_gene_presence_absence_consensus_2001 # gene-presence absence Mass 2001
ggCPP_intermed_gene_presence_absence_consensus # maybe better because other has na for some

# calculate gene freqs, based on freq1, freq3?
PP_mass_cluster_freq_1_genes_mtx <- sapply(ggCPP_intermed_gene_presence_absence_consensus[-1,-1], as.integer) * PP_mass_cluster_freq_1
PP_mass_cluster_freq_3_genes_mtx <- sapply(ggCPP_intermed_gene_presence_absence_consensus[-1,-1], as.integer) * PP_mass_cluster_freq_3

PP_mass_cluster_freq_1_genes <- rowSums(PP_mass_cluster_freq_1_genes_mtx)
PP_mass_cluster_freq_3_genes <- rowSums(PP_mass_cluster_freq_3_genes_mtx)


fit <- glmnet(PP_mass_cluster_freq_1_genes_mtx, PP_mass_cluster_freq_3_genes)

plot(fit)
print(fit)
coef(fit, s = 0.1)

fit2 <- glmnet(t(PP_mass_cluster_freq_1_genes_mtx), PP_mass_cluster_freq_3)
plot(fit2)
print(fit2)
coef(fit2, s = 0.1)

fit3 <- glmnet(t(ggCPP_intermed_gene_presence_absence_consensus[-1,-1]), PP_mass_cluster_freq_3)
plot(fit3)
print(fit3)
coef(fit3, s = 0.1)

cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus[-1,-1]), PP_mass_cluster_freq_3)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")

cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus[-1,-1]), PP_mass_cluster_freq_3-PP_mass_cluster_freq_1)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")

cvfit <- cv.glmnet(t(PP_mass_cluster_freq_1_genes_mtx), PP_mass_cluster_freq_3)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
predict(cvfit, newx = x[1:5,], s = "lambda.min")

# add VT info to feature vector?
ggCPP_intermed_gene_presence_absence_consensus_withVT <- rbind(ggCPP_intermed_gene_presence_absence_consensus, c("VT",PP_mass_VT))

# fit elastic net to explain pre-post vacc changes in GPSC frequency
cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus_withVT[-1,-1]), PP_mass_cluster_freq_3 - PP_mass_cluster_freq_1)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))
# 4
(ggCPP_intermed_gene_presence_absence_consensus_withVT)[which(as.vector(coef(cvfit, s = "lambda.min")) != 0),1]
# "0" "group_4121" "group_6291" "VT" 
# "0" "group_4121" "group_1826" "group_6291" "group_7700" "VT"
# "0" "group_4121"  "group_8047"  "group_6841"  "group_8293"  "group_7584" "group_3305"  "group_12369" "group_11908" "group_9711"  "group_3143"  "group_4827" "group_6009"  "group_782"   "group_1826"  "group_6291"  "group_7700"  "VT"  
# "0" "group_4121" "group_8047" "group_7584" "group_3305" "group_6009" "group_1826" "group_6291" "VT"
# "0" "group_4121" "group_8047" "group_7584" "group_4827" "group_6009" "group_782" "group_1826" "group_6291" "VT"

# this does not look to bad!
# check delta
ggC_delta_ranking2[which(as.vector(coef(cvfit, s = "lambda.min")) != 0)[-c(1,length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0)))]-1]
#1371 1490 2277 2591 3147 3350 3588 3850 
#955 1985 1963 1984 1961 1920  553 1373
# delta ranking would probably only consider gene under NFDS which <= max(ggC_delta_ranking2) * 0.35 = 712
# but maybe I don't necessarily expect this to match anyway?

cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus_withVT[-1,-1]), PP_mass_cluster_freq_3)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))

# check whether pre-vacc population has predictive power for post (yes! but the vacc loose pred power)
cvfit <- cv.glmnet(cbind(t(ggCPP_intermed_gene_presence_absence_consensus_withVT[-1,-1]),PP_mass_cluster_freq_1), PP_mass_cluster_freq_3)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))


cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus_withVT[-1,-1]), PP_mass_cluster_freq_3)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))
# 11 / 18 (there seems to be quite some variability)

ggCPP_intermed_gene_presence_absence_consensus_withVT_genesNeg <- ggCPP_intermed_gene_presence_absence_consensus_withVT
ggCPP_intermed_gene_presence_absence_consensus_withVT_genesNeg[-1, 1 + which(PP_mass_VT==1)] <- apply(ggCPP_intermed_gene_presence_absence_consensus_withVT_genesNeg[-1, 1 + which(PP_mass_VT==1)], c(1,2),as.numeric) * (-1)

cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus_withVT_genesNeg[-1,-1]), PP_mass_cluster_freq_3)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))
# 11 / 18 (there seems to be quite some variability)

ggCPP_intermed_gene_presence_absence_consensus_clusterFreq <- t(apply(ggCPP_intermed_gene_presence_absence_consensus[-1,-1], c(1,2), as.numeric)) * as.numeric(PP_mass_cluster_freq_1)
ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq <- cbind(ggCPP_intermed_gene_presence_absence_consensus_clusterFreq, PP_mass_VT)
rownames(ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq) <- ggCPP_intermed_gene_presence_absence_consensus[1,-1] 
colnames(ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq) <- c(ggCPP_intermed_gene_presence_absence_consensus[-1,1], "VT")

cvfit <- cv.glmnet(ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq, PP_mass_cluster_freq_3)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))
# 25
colnames(ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq)[which(as.vector(coef(cvfit, s = "lambda.min")) != 0)-1]

#run1 <- colnames(ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq)[which(as.vector(coef(cvfit, s = "lambda.min")) != 0)-1]
#run2 <- colnames(ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq)[which(as.vector(coef(cvfit, s = "lambda.min")) != 0)-1]
#run3 <- colnames(ggCPP_intermed_gene_presence_absence_consensus_withVT_clusterFreq)[which(as.vector(coef(cvfit, s = "lambda.min")) != 0)-1]

intersect(intersect(run1, run2), run3)
# "group_6841" "group_8486" "group_3992" "group_3099" "group_7491" "group_6539" "group_1442"
# "group_3731" "group_5081" "group_1498" "group_8593" "group_8396" "group_7871" "group_4871"
# "group_6291" "VT"


# repeat for UK
UK_cluster_freqs[[1]] # cluster freq UK 1st year
UK_cluster_freqs[[7]] # cluster freq UK last year
UK_ggCaller_intermed_consensus # intermed consensus pangenome matrix UK
UK_VT # PCV7
UK_VT2 # PCV13

UK_ggCaller_intermed_consensus_clusterFreq <- t(apply(UK_ggCaller_intermed_consensus[-1,-1], c(1,2), as.numeric)) * as.numeric(UK_cluster_freqs[[1]])
UK_ggCaller_intermed_consensus_withVT_clusterFreq <- cbind(UK_ggCaller_intermed_consensus_clusterFreq, UK_VT,UK_VT2)
rownames(UK_ggCaller_intermed_consensus_withVT_clusterFreq) <- UK_ggCaller_intermed_consensus[1,-1] 
colnames(UK_ggCaller_intermed_consensus_withVT_clusterFreq) <- c(UK_ggCaller_intermed_consensus[-1,1], "VT_PCV7","VT_PCV13")

cvfit <- cv.glmnet(UK_ggCaller_intermed_consensus_withVT_clusterFreq, UK_cluster_freqs[[7]])
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))
# 16
colnames(UK_ggCaller_intermed_consensus_withVT_clusterFreq)[which(as.vector(coef(cvfit, s = "lambda.min")) != 0)-1]

UK_ggCaller_intermed_consensus_withVT <- rbind(UK_ggCaller_intermed_consensus, c("VT_PCV7",UK_VT), c("VT_PCV13",UK_VT2))

cvfit <- cv.glmnet(cbind(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]),UK_cluster_freqs[[1]]), UK_cluster_freqs[[7]])
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))
# 16
cbind(t(UK_ggCaller_intermed_consensus_withVT),c("preVT",UK_cluster_freqs[[1]]))[1,which(as.vector(coef(cvfit, s = "lambda.min")) != 0)]



cvfit <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs[[7]] - UK_cluster_freqs[[1]])
#cvfit <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs[[7]])
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(fit3, s = 0.1)) != 0))

fit3 <- glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs[[7]] - UK_cluster_freqs[[1]])
plot(fit3)
print(fit3)
coef(fit3, s = 0.1)

fit3 <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs[[7]] - UK_cluster_freqs[[1]])
plot(fit3) # I think this plot is saying me that 1 variable is the best?
print(fit3)
coef(fit3, s = "lambda.min")
length(which(as.vector(coef(fit3, s = "lambda.min")) != 0))
# 1
intersect(which(as.vector(coef(cvfit, s = "lambda.min")) != 0), which(as.vector(coef(fit3, s = "lambda.min")) != 0))
# and an overlap of 1

fit3 <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs_winter[[7]])
plot(fit3)
print(fit3)
coef(fit3, s = 0.1)
length(which(as.vector(coef(fit3, s = "lambda.min")) != 0))
# 6 
intersect(which(as.vector(coef(cvfit, s = "lambda.min")) != 0), which(as.vector(coef(fit3, s = "lambda.min")) != 0))
# only one overlap?


fit3 <- cv.glmnet(t(cbind(UK_ggCaller_intermed_consensus_withVT[-1,-1], UK_ggCaller_intermed_consensus_withVT[-1,-1])), c(UK_cluster_freqs_winter[[7]],UK_cluster_freqs_winter[[1]]))
plot(fit3)
print(fit3)
coef(fit3, s = "lambda.min")
length(which(as.vector(coef(fit3, s = "lambda.min")) != 0))
# 1

# vt1 vs pre
fit3 <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs_winter[[4]] - UK_cluster_freqs_winter[[1]])
plot(fit3) # I think this plot is saying me that 0 variables is the best?
print(fit3)
coef(fit3, s = "lambda.min")
length(which(as.vector(coef(fit3, s = "lambda.min")) != 0))
# but then it outputs 75 coeffs
intersect(which(as.vector(coef(cvfit, s = "lambda.min")) != 0), which(as.vector(coef(fit3, s = 0.1)) != 0))
# and an overlap of 10

# vt1 vs vt2
fit3 <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs_winter[[7]] - UK_cluster_freqs_winter[[4]])
plot(fit3) # I think this plot is saying me that 0 variables is the best?
print(fit3)
fit3$lambda.min
coef(fit3, s = "lambda.min")
length(which(as.vector(coef(fit3, s = "lambda.min")) != 0))
intersect(which(as.vector(coef(cvfit, s = 0.1)) != 0), which(as.vector(coef(fit3, s = 0.1)) != 0))



UK_ggCaller_intermed_consensus_withVTandpreVT <- rbind(UK_ggCaller_intermed_consensus, c("VT_PCV7",UK_VT), c("VT_PCV13",UK_VT2), c("pre_VT",UK_cluster_freqs_winter[[1]]))
# alpha = 0.5 = elastic net
fit3 <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVTandpreVT[-1,-1]), UK_cluster_freqs_winter[[7]] - UK_cluster_freqs_winter[[1]])
plot(fit3) # I think this plot is saying me that 0 variables is the best?
print(fit3)
fit3$lambda.min
coef(fit3, s = "lambda.min")
length(which(as.vector(coef(fit3, s = "lambda.min")) != 0))
# 42
#intersect(which(as.vector(coef(cvfit, s = 0.1)) != 0), which(as.vector(coef(fit3, s = 0.1)) != 0))
UK_ggCaller_intermed_consensus_withVTandpreVT[which(as.vector(coef(fit3, s = "lambda.min")) != 0)+1,1]


ggCPP_intermed_gene_presence_absence_consensus_withVTandpreVT <- rbind(ggCPP_intermed_gene_presence_absence_consensus, c("VT",PP_mass_VT),c("preVT",PP_mass_cluster_freq_1))
# fit elastic net to explain pre-post vacc changes in GPSC frequency
cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus_withVTandpreVT[-1,-1]), PP_mass_cluster_freq_3 - PP_mass_cluster_freq_1)
cvfit <- cv.glmnet(t(ggCPP_intermed_gene_presence_absence_consensus_withVT[-1,-1]), PP_mass_cluster_freq_3 - PP_mass_cluster_freq_1)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
length(which(as.vector(coef(cvfit, s = "lambda.min")) != 0))
# 6
which(as.vector(coef(cvfit, s = "lambda.min")) != 0)

genes1 <- which(as.vector(coef(cvfit, s = "lambda.min")) != 0)
ggCPP_intermed_gene_presence_absence_consensus_withVTandpreVT[which(as.vector(coef(cvfit, s = "lambda.min")) != 0)+1,1]

# plot tanglegram / alluvial plot
library(ggplot2)
library(ggalluvial)
library(scales)
intermed_genes_UKUS <- intersect(names(ggC_delta_data2),names(UK_delta_data))
intermed_genes_UKUS_df <- data.frame(cbind(1:length(intermed_genes_UKUS),names(sort(ggC_delta_data2[intermed_genes_UKUS])), names(sort(UK_delta_data[intermed_genes_UKUS])), rep(1,length(intermed_genes_UKUS))))
colnames(intermed_genes_UKUS_df) <- c("Position","US","UK","test_col")
ggplot(as.data.frame(intermed_genes_UKUS_df), aes(y = Position, axis1 = US, axis2 = UK)) +
  geom_alluvium(aes(fill = test_col), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("UC Berkeley admissions and rejections, by sex and department")

ggplot(intermed_genes_UKUS_df,
       aes(axis1 = US, axis2 = UK)) +
  geom_alluvium(aes(fill = US), width = 0.2) +
  geom_stratum(width = 0.2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Vector 1", "Vector 2")) +
  theme_minimal()

df <- data.frame(
  value = names(sort(ggC_delta_data2[intermed_genes_UKUS])),
  pos1 = match(names(sort(ggC_delta_data2[intermed_genes_UKUS])), names(sort(ggC_delta_data2[intermed_genes_UKUS]))),
  pos2 = match(names(sort(ggC_delta_data2[intermed_genes_UKUS])), names(sort(UK_delta_data[intermed_genes_UKUS]))),
  NFDS_US = c(rep(1,500), rep(0,length(names(sort(ggC_delta_data2[intermed_genes_UKUS])))-500))
)

ggplot(df) +
  geom_segment(aes(x = 1, xend = 2, y = pos1, yend = pos2, color = NFDS_US), size = 1, alpha = 0.4) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  labs(x = "Vector 1 → Vector 2", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
