#------------------------------------------------------------------------
# Contact Zoe Ji (chji1@bwh.harvard.edu) or Lauren Breithaupt (lbreithaupt@mgh.harvard.edu) for any bug report.
# Last update: 6.16.2023
#------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
#------------------------------------------------------------------------
load("~/npx_baseline.Rdata") ### our baseline dataset that contains protein expression (npx values) and clinical data. 
load("~/protein.RData") ### the list of proteins that we are going to analyze (N=73)
#------------------------------------------------------------------------
# Main Model 1: account for age.
#------------------------------------------------------------------------

# (1.1) AN/atyp-AN vs. HC
npx_baseline$group <- factor(npx_baseline$group, levels = c("HC", "ED")) ### AN/atyp-AN grup was combined and labeled as ED group.
group2_main1 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ group + ageyrs, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_main1) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_main1 <- group2_main1 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (1.2) AN vs. atyp-AN vs. HC
### we created a low_wt variables to categorize AN (lw_wt) and atyp-AN (non_low_wt)
npx_baseline$low_wt <- ifelse(npx_baseline$group=="HC", "hc", npx_baseline$low_wt)
npx_baseline$low_wt <- factor(npx_baseline$low_wt, levels = c("hc", "low_wt", "non_low_wt"))
npx_baseline$low_wt1 <- relevel(npx_baseline$low_wt, ref = "low_wt")

group3_main1 <- lapply(protein, function(i){
  model1 <- lm(npx_baseline[,i] ~ low_wt + ageyrs, data=npx_baseline)
  ### change reference group to get the third comparison
  model2 <- lm(npx_baseline[,i] ~ low_wt1 + ageyrs, data=npx_baseline)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (1.3) Association with BMI
bmi_main1 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ bmiz + ageyrs, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_main1) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_main1 <- bmi_main1 %>% 
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Main 2: account for age, smoking status (Y/N), body temperature, 
# antihistamine use (Y/N), any psychiatric medication use(Y/N), psychiatric comorbidity (Y/N)
#------------------------------------------------------------------------
npx_baseline$smokerm <- factor(npx_baseline$smokerm, levels=c(0,1), labels=c("No", "Yes"))
npx_baseline$temperature <- as.numeric(npx_baseline$temperature)
npx_baseline$antihis <- factor(npx_baseline$antihis, levels=c(0, 1), labels=c("No", "Yes"))
npx_baseline$psych_med_new <- factor(npx_baseline$psych_med_new, levels=c(FALSE, TRUE), labels=c("No", "Yes"))
npx_baseline$ksads_any <- factor(npx_baseline$ksads_any, levels=c(FALSE, TRUE), labels=c("No", "Yes"))

# (2.1) AN/atyp-AN vs. HC
group2_main2 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ group + ageyrs + smokerm + temperature + 
                antihis + psych_med_new + ksads_any, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_main2) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_main2 <- group2_main2 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (2.2) AN vs. atyp-AN vs. HC
group3_main2 <- lapply(protein, function(i){
  model1 <- lm(npx_baseline[,i] ~ low_wt + ageyrs + smokerm + temperature + 
                 antihis + psych_med_new + ksads_any, data=npx_baseline)
  ### change reference group to get the third comparison
  model2 <- lm(npx_baseline[,i] ~ low_wt1 + ageyrs + smokerm + temperature + 
                 antihis + psych_med_new + ksads_any, data=npx_baseline)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (2.3) Association with BMI
bmi_main2 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ bmiz + ageyrs + smokerm + temperature + 
                antihis + psych_med_new + ksads_any, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_main2) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_main2 <- bmi_main2 %>% 
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Secondary Model 1: removed subjects who were smoking, account for age
#------------------------------------------------------------------------
npx1 <- npx_baseline %>% dplyr::filter(smokerm!="Yes")

# (3.1) AN/atyp-AN vs. HC
group2_sup1 <- lapply(protein, function(i){
  model <- lm(npx1[,i] ~ group + ageyrs, data=npx1)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_sup1) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_sup1 <- group2_sup1 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (3.2) AN vs. atyp-AN vs. HC
group3_sup1 <- lapply(protein, function(i){
  model1 <- lm(npx1[,i] ~ low_wt + ageyrs, data=npx1)
  model2 <- lm(npx1[,i] ~ low_wt1 + ageyrs, data=npx1)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (3.3) Association with BMI
bmi_sup1 <- lapply(protein, function(i){
  model <- lm(npx1[,i] ~ bmiz + ageyrs, data=npx1)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_sup1) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_sup1 <- bmi_sup1 %>% 
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Secondary Model 2: removed subjects whose body temperature>99, account for age
#------------------------------------------------------------------------
npx2 <- npx_baseline %>% dplyr::filter(temperature <= 99)

# (4.1) AN/atyp-AN vs. HC
group2_sup2 <- lapply(protein, function(i){
  model <- lm(npx2[,i] ~ group + ageyrs, data=npx2)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_sup2) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_sup2 <- group2_sup2 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (4.2) AN vs. atyp-AN vs. HC
group3_sup2 <- lapply(protein, function(i){
  model1 <- lm(npx2[,i] ~ low_wt + ageyrs, data=npx2)
  model2 <- lm(npx2[,i] ~ low_wt1 + ageyrs, data=npx2)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (4.3) Association with BMI
bmi_sup2 <- lapply(protein, function(i){
  model <- lm(npx2[,i] ~ bmiz + ageyrs, data=npx2)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_sup2) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_sup2 <- bmi_sup2 %>% 
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Secondary Model 3: account for age and number of psychiatric diagnoses
#------------------------------------------------------------------------

# (5.1) AN/atyp-AN vs. HC
group2_sup3 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ group + ageyrs + ksads_count, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_sup3) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_sup3 <- group2_sup3 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (5.2) AN vs. atyp-AN vs. HC
group3_sup3 <- lapply(protein, function(i){
  model1 <- lm(npx_baseline[,i] ~ low_wt + ageyrs + ksads_count, data=npx_baseline)
  model2 <- lm(npx_baseline[,i] ~ low_wt1 + ageyrs + ksads_count, data=npx_baseline)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (5.3) Association with BMI
bmi_sup3 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ bmiz + ageyrs + ksads_count, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_sup3) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_sup3 <- bmi_sup3 %>% 
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Secondary Model 4: account for age and anxiety t score
#------------------------------------------------------------------------
npx_baseline$anxiety_t_score <- as.numeric(npx_baseline$anxiety_t_score)

# (6.1) AN/atyp-AN vs. HC
group2_sup4 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ group + ageyrs + anxiety_t_score, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_sup4) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_sup4 <- group2_sup4 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (6.2) AN vs. atyp-AN vs. HC
group3_sup4 <- lapply(protein, function(i){
  model1 <- lm(npx_baseline[,i] ~ low_wt + ageyrs + anxiety_t_score, data=npx_baseline)
  model2 <- lm(npx_baseline[,i] ~ low_wt1 + ageyrs + anxiety_t_score, data=npx_baseline)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (6.3) Association with BMI
bmi_sup4 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ bmiz + ageyrs + anxiety_t_score, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_sup4) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_sup4 <- bmi_sup4 %>% 
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Secondary Model 5: account for age and antidepressants use (Y/N)
#------------------------------------------------------------------------
npx_baseline$psych_med_antidepress <- factor(npx_baseline$psych_med_antidepress,
                                             levels = c(0,1),
                                             labels = c("No", "Yes"))

# (7.1) AN/atyp-AN vs. HC
group2_sup5 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ group + ageyrs + psych_med_antidepress, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_sup5) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_sup5 <- group2_sup5 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (7.2) AN vs. atyp-AN vs. HC
group3_sup5 <- lapply(protein, function(i){
  model1 <- lm(npx_baseline[,i] ~ low_wt + ageyrs + psych_med_antidepress, data=npx_baseline)
  model2 <- lm(npx_baseline[,i] ~ low_wt1 + ageyrs + psych_med_antidepress, data=npx_baseline)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (7.3) Association with BMI
bmi_sup5 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ bmiz + ageyrs + psych_med_antidepress, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_sup5) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_sup5 <- bmi_sup5 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Secondary Model 6: account for age and antihistamine+NSAIDs use (Y/N)
#------------------------------------------------------------------------
npx_baseline$antihis_NSAID <- factor(npx_baseline$antihis_NSAID, levels=c(0, 1), labels=c("No", "Yes"))

# (8.1) AN/atyp-AN vs. HC
group2_sup6 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ group + ageyrs + antihis_NSAID, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(group2_sup6) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
group2_sup6 <- group2_sup6 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (8.2) AN vs. atyp-AN vs. HC
group3_sup6 <- lapply(protein, function(i){
  model1 <- lm(npx_baseline[,i] ~ low_wt + ageyrs + antihis_NSAID, data=npx_baseline)
  model2 <- lm(npx_baseline[,i] ~ low_wt1 + ageyrs + antihis_NSAID, data=npx_baseline)
  test <- cbind(i, rbind(summary(model1)$coef[2:3,], summary(model2)$coef[3,]), 
                c("AN-HC", "AtypAN-HC", "AtypAN-AN")) %>% as.data.frame()
  colnames(test) <- c("Protein", "Coefficient", "Std", "t.value", "p.value", "Comp")
  rownames(test) <- NULL
  test
}) %>% do.call(rbind,.) %>% as.data.frame %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

# (8.3) Association with BMI
bmi_sup6 <- lapply(protein, function(i){
  model <- lm(npx_baseline[,i] ~ bmiz + ageyrs + antihis_NSAID, data=npx_baseline)
  c(i, summary(model)$coef[2,])
}) %>% do.call(rbind,.) %>% as.data.frame

colnames(bmi_sup6) <- c("Protein", "Coefficient", "Std", "t.value", "p.value")
bmi_sup6 <- bmi_sup6 %>%
  dplyr::mutate(Coefficient = as.numeric(Coefficient),
                p.value = as.numeric(p.value),
                FDR = p.adjust(p.value, method="BH")) %>%
  dplyr::arrange(p.value)

#------------------------------------------------------------------------
# Visualization
#------------------------------------------------------------------------

#-----------------
# Volcano plots
#-----------------

# FDR vs. log2FC for AN/atyp-AN vs. HC
edhc_plot <- group2_main1 %>% dplyr::select(Protein, Coefficient, p.value, FDR)
colnames(edhc_plot) <- c("Protein", "log2FC", "PValue", "FDR") 

edhc_plot <- edhc_plot %>% 
  ### change from ED vs. HC to HC vs. ED
  mutate(log2FC = -log2FC) %>%
  mutate(
    Significance = case_when(
      FDR <= 0.3 & FDR > 0.1 ~ "FDR 0.3", 
      FDR <= 0.1 & FDR >= 0.06 ~ "FDR 0.1",
      FDR < 0.06 ~ "FDR 0.05", 
      TRUE ~ "Unchanged")
  )

p1 <- ggplot(edhc_plot, aes(log2FC, -log(FDR,10))) +
  geom_point(aes(color = Significance), size = 2) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

top_protein1 <- edhc_plot %>% dplyr::filter(FDR <= 0.1)
arrowLab1 <- data.frame(lab = "Higher in AN/atyp-AN", x = -0.3, y = -0.25)
pp1 <- p1 + geom_segment(aes(x = 0, xend= -0.15, y= -0.25, yend= -0.25), 
                         arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(data = arrowLab1,aes(x=x,y=y,label = lab),size = 3) +
  geom_label_repel(data = top_protein1, mapping = aes(log2FC, -log(FDR,10), label = Protein), size = 2) +
  labs(title="A: HC vs. AN/atyp-AN")

# FDR vs. log2FC for AN vs. HC
anhc_plot <- group3_main1 %>% dplyr::filter(Comp == "AN-HC") %>% dplyr::select(Protein, Coefficient, p.value, FDR)
colnames(anhc_plot) <- c("Protein", "log2FC", "PValue", "FDR") 

anhc_plot <- anhc_plot %>% 
  mutate(log2FC = -log2FC) %>%
  mutate(
    Significance = case_when(
      FDR <= 0.3 & FDR > 0.1 ~ "FDR 0.3", 
      FDR <= 0.1 & FDR >= 0.06 ~ "FDR 0.1",
      FDR < 0.06 ~ "FDR 0.05", 
      TRUE ~ "Unchanged")
  )

p2 <- ggplot(anhc_plot, aes(log2FC, -log(FDR,10))) +
  geom_point(aes(color = Significance), size = 2) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

top_protein2 <- anhc_plot %>% dplyr::filter(FDR <= 0.1)
arrowLab2 <- data.frame(lab = "Higher in AN", x = -0.32, y = -0.25)
pp2 <- p2 + geom_segment(aes(x = 0, xend = -0.2, y= -0.25, yend= -0.25), 
                         arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(data = arrowLab2,aes(x=x,y=y,label = lab),size = 3) +
  geom_label_repel(data = top_protein2, 
                   mapping = aes(log2FC, -log(FDR,10), label = Protein), 
                   size = 2) +
  labs(title="B: HC vs. AN")

# FDR vs. log2FC for association with BMI
bmi_plot <- bmi_main1 %>% dplyr::select(Protein, Coefficient, p.value, FDR)
colnames(bmi_plot) <- c("Protein", "log2FC", "PValue", "FDR") 

bmi_plot <- bmi_plot %>% 
  mutate(
    Significance = case_when(
      FDR <= 0.3 & FDR > 0.1 ~ "FDR 0.3", 
      FDR <= 0.1 & FDR >= 0.06 ~ "FDR 0.1",
      FDR < 0.06 ~ "FDR 0.05", 
      TRUE ~ "Unchanged")
  )

p3 <- ggplot(bmi_plot, aes(log2FC, -log(FDR,10))) +
  geom_point(aes(color = Significance), size = 2) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5)))

top_protein3 <- bmi_plot %>% dplyr::filter(FDR <= 0.1)
pp3 <- p3 + geom_label_repel(data = top_protein3, 
                             mapping = aes(log2FC, -log(FDR,10), label = Protein), 
                             size = 2) +
  labs(title="C: Association with BMI")

ggarrange(
  pp1, pp2, pp3, nrow=3,
  common.legend = TRUE, legend = "right") 

#-----------------
# Heatmaps
#-----------------

# Heatmap for AN/atyp-AN vs. HC
edhc_plot <- edhc_plot %>% 
  mutate(Group = case_when(FDR > 0.3 ~ "No difference",
                           log2FC > 0 & FDR < 0.06 ~ "Positive (FDR < 0.05)",
                           log2FC > 0 & FDR <= 0.3 & FDR >= 0.06 ~ "Positive (FDR < 0.3)",
                           log2FC < 0 & FDR < 0.06 ~ "Negative (FDR < 0.05)",
                           log2FC < 0 & FDR <= 0.3 & FDR >= 0.06 ~ "Negative (FDR < 0.3)"))
edhc_plot$Group <- factor(edhc_plot$Group, 
                          levels = c("Positive (FDR < 0.05)", "Positive (FDR < 0.3)",
                                     "No difference", "Negative (FDR < 0.3)", 
                                     "Negative (FDR < 0.05)"))

h1 <- ggplot(edhc_plot, aes(y=Protein, x=1)) +
  geom_tile(aes(fill=Group)) +
  scale_fill_manual(values = alpha(c("#CA0300", "#E16519", "#ffffff", "#2986CC", "#0B5394"))) +
  xlab("HC vs. AN/atyp-AN") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Heatmap for AN vs. HC
anhc_plot <- anhc_plot %>% 
  mutate(Group = case_when(FDR > 0.3 ~ "No difference",
                           log2FC > 0 & FDR < 0.06 ~ "Positive (FDR < 0.05)",
                           log2FC > 0 & FDR <= 0.3 & FDR >= 0.06 ~ "Positive (FDR < 0.3)",
                           log2FC < 0 & FDR < 0.06 ~ "Negative (FDR < 0.05)",
                           log2FC < 0 & FDR <= 0.3 & FDR >= 0.06 ~ "Negative (FDR < 0.3)"))
anhc_plot$Group <- factor(anhc_plot$Group, 
                          levels = c("Positive (FDR < 0.05)", "Positive (FDR < 0.3)",
                                     "No difference", "Negative (FDR < 0.3)", 
                                     "Negative (FDR < 0.05)"))

h2 <- ggplot(anhc_plot, aes(y=Protein, x=1)) +
  geom_tile(aes(fill=Group)) +
  scale_fill_manual(values = alpha(c("#CA0300", "#E16519", "#ffffff", "#2986CC", "#0B5394"))) +
  xlab("HC vs. AN") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Heatmap for association with BMI
bmi_plot <- bmi_plot %>% 
  mutate(Group = case_when(FDR > 0.3 ~ "No difference",
                           log2FC > 0 & FDR < 0.06~ "Positive (FDR < 0.05)",
                           log2FC > 0 & FDR <= 0.3 & FDR >= 0.06 ~ "Positive (FDR < 0.3)",
                           log2FC < 0 & FDR < 0.06 ~ "Negative (FDR < 0.05)",
                           log2FC < 0 & FDR <= 0.3 & FDR >= 0.06 ~ "Negative (FDR < 0.3)"))
bmi_plot$Group <- factor(bmi_plot$Group, 
                         levels = c("Positive (FDR < 0.05)", "Positive (FDR < 0.3)",
                                    "No difference", "Negative (FDR < 0.3)", 
                                    "Negative (FDR < 0.05)"))

h3 <- ggplot(bmi_plot, aes(y=Protein, x=1)) +
  geom_tile(aes(fill=Group)) +
  scale_fill_manual(values = alpha(c("#CA0300", "#E16519", "#ffffff", "#2986CC", "#0B5394"))) +
  xlab("Association with BMI") +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggarrange(
  h1, h2, h3, labels = c("A", "B", "C"), nrow=1,
  common.legend = TRUE, legend = "right"
) 

#------------------------------------------------------------------------
# Since you made it to the end of our code - run the code below for a fun surprise!
#------------------------------------------------------------------------

plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))
rect(3, 3, 7, 7, col = "gray", border = NA)
rect(4, 7, 6, 9, col = "gray", border = NA)
points(4.5, 8.5, pch = 15, col = "white", cex = 2)
points(5.5, 8.5, pch = 15, col = "white", cex = 2)
polygon(c(5, 5.5, 4.5), c(8, 7.5, 7.5), col = "pink", border = NA)
polygon(c(3, 4, 4), c(8, 9, 7), col = "gray", border = NA)
polygon(c(7, 6, 7), c(8, 9, 7), col = "gray", border = NA)
segments(4, 7.5, 3, 7, lwd = 2)
segments(4, 7.5, 3, 8, lwd = 2)
segments(4, 7.5, 3, 9, lwd = 2)
segments(6, 7.5, 7, 7, lwd = 2)
segments(6, 7.5, 7, 8, lwd = 2)
segments(6, 7.5, 7, 9, lwd = 2)
text(5, 1, "Coding Cat! Meow!", col = "gray30", cex = 1.5, font = 2, adj = 0.5)
# Optional: Save the plot as an image file
# dev.copy(png, "cat.png")
# dev.off()
