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

# repeat for UK
UK_cluster_freqs[[1]] # cluster freq UK 1st year
UK_cluster_freqs[[7]] # cluster freq UK last year
UK_ggCaller_intermed_consensus # intermed consensus pangenome matrix UK
UK_VT # PCV7
UK_VT2 # PCV13

UK_ggCaller_intermed_consensus_withVT <- rbind(UK_ggCaller_intermed_consensus, c("VT_PCV7",UK_VT), c("VT_PCV13",UK_VT2))

cvfit <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs[[7]] - UK_cluster_freqs[[1]])
cvfit <- cv.glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs[[7]])
plot(cvfit)
print(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")

fit3 <- glmnet(t(UK_ggCaller_intermed_consensus_withVT[-1,-1]), UK_cluster_freqs[[7]] - UK_cluster_freqs[[1]])
plot(fit3)
print(fit3)
coef(fit3, s = 0.1)




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
