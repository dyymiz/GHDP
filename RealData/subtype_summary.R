rm(list=ls())

library("readxl")
library("tidyverse")
library("dplyr")
library(ggplot2)
library(mgcv)
library(fields)
library(reshape2)
library(scales)
library(plyr)

load("C:/Users/ydyhf/Box Sync/DI/code/Real_data/Restart_2021/GHDP_un_sup/All stuff_LungCancerDatamRNAMethLUSCLUAD_infinite_scaled_selected_probs_restart_2021.all.RData")
load("C:/Users/ydyhf/OneDrive/Desktop/TCGA-Assembler-2021/LungCancerDataLUSCLUADClinical.restart.2021.all.RData")
#load("C:/Users/ydyhf/Box Sync/DI/code/Real_data/Restart_2021/GHDP_un_sup_non/All.Stuff.non.RData")

### TODO: the two lists are not complete
LUAD_subtype <- read_excel("C:/Users/ydyhf/Box Sync/DI/code/Real_data/LUAD tables.xlsx", sheet = 7, skip = 4)
LUSC_subtype <- read_excel("C:/Users/ydyhf/Box Sync/DI/code/Real_data/LUSC tables/LUSC table.xls", skip = 2)

# group singletons together
data_to_plot <- data
data_to_plot[[2]] <- exp(data_to_plot[[2]])/(1+exp(data_to_plot[[2]]))
rownames(data_to_plot[[1]]) = rownames(data_to_plot[[2]]) = colnames(data_to_plot[[1]]) = colnames(data_to_plot[[2]]) =  NULL
data_to_plot <- lapply(1:t, function(x) data_to_plot[[x]][order(All.Stuff$est_c.v[[1]]), 
                                                       order(All.Stuff$est_c.v.c[[1]][[x]])])

# lapply(1:t, function(x) image.plot(1:dim(data_to_plot[[x]])[2], 
#                                    1:dim(data_to_plot[[x]])[1], 
#                                    t(data_to_plot[[x]][order(All.Stuff$est_c.v[[1]]), 
#                                                        order(All.Stuff.non$est_c.v.c[[1]][[x]])] ),
#                                    col=grey(seq(1, 0, length = 256))))

### ggplot ###
fn.melt <- function(data){
  n=dim(data[[1]])[1]
  p=dim(data[[1]])[2]
  tmp <- cbind(xx=1, melt(data[[1]], id.vars=c(n,p)))
  for(xx in 2:t){
    n=dim(data[[xx]])[1]
    p=dim(data[[xx]])[2]
    C.melt <- cbind(xx, melt(data[[xx]], id.vars=c(n,p)))
  }
  C.melt <- rbind(tmp, C.melt)
  colnames(C.melt) <- (c('xx','n','p','value'))
  C.melt
}

order <- fn.melt(data_to_plot)

heat.data <- function(data){
  p <- ggplot(data,aes(p, n)) + 
    geom_tile(aes(fill = value)) +
    # facet_wrap(~xx)+
    # geom_hline(yintercept=c(28.5, 20), size=1) +
    geom_segment(aes(x = 0, y = 28.5, xend = max(p), yend = 28.5), size = 1) +
    geom_segment(aes(x = 0, y = 43.5, xend = max(p), yend = 43.5), size = 1) +
    geom_segment(aes(x = 0, y = 56.5, xend = max(p), yend = 56.5), size = 1) +
    scale_fill_gradient(low = "white",high = "black")+
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          legend.title=element_text(size=24), 
          axis.title=element_text(size=24),
          axis.text.y = element_text(angle = 0), 
          axis.title.y = element_text(angle = 0, vjust = 0.5), 
          axis.text = element_text(size = 22), 
          legend.text=element_text(size=22),
          legend.key.size = unit(1, "cm"))

  p
}

order_mRNA <- order %>% subset(xx == 1)
order_meth <- order %>% subset(xx == 2)

p <- heat.data(order_mRNA)
p <- p + xlab("p1")
ggsave("clrRNA.png", height = 8, width =12)
p

p <- heat.data(order_meth)
p <- p + xlab("p2")
ggsave("clrMeth.png", height = 8, width =12)
p
# lapply(1:t, function(x) image.plot(t(data[[x]][c(which(All.Stuff$est_c.v[[1]]%in% c(12)),which(All.Stuff$est_c.v[[1]]%in% c(5)),which(All.Stuff$est_c.v[[1]]%in% c(8)),which(All.Stuff$est_c.v[[1]]%in% c(13)), which(!All.Stuff$est_c.v[[1]] %in% c(12, 5, 8, 13)) ), order(All.Stuff.non$est_c.v.c[[1]][[x]]) ]),col=grey(seq(0, 1, length = 256))))

lapply(1:t, function(x) image.plot(t(data_to_plot[[x]][order(All.Stuff$est_c.v[[1]]), order(All.Stuff$est_c.v.c[[1]][[x]])] ),col=grey(seq(0, 1, length = 256))))


table(All.Stuff$est_c.v[[1]])
table(All.Stuff$est_c.v[[1]][1:45])
table(All.Stuff$est_c.v[[1]][46:133])


LUAD_subtype$`Tumor ID` <- paste0(gsub("-",".", LUAD_subtype$`Tumor ID`),".01")
LUSC_subtype$`Tumor ID` <- paste0(gsub("-",".", LUSC_subtype$`Tumor ID`),".01")
LUSC_subtype$`Tumor ID` <- gsub("LUSC","TCGA", LUSC_subtype$`Tumor ID`)

sum(clinical.slt$bcr_patient_barcode %in% LUAD_subtype$`Tumor ID`)
sum(clinical.slt$bcr_patient_barcode %in% LUSC_subtype$`Tumor ID`)

table(All.Stuff$est_c.v[[1]][names(All.Stuff$est_c.v[[1]]) %in% LUAD_subtype$`Tumor ID`])
table(All.Stuff$est_c.v[[1]][names(All.Stuff$est_c.v[[1]]) %in% LUSC_subtype$`Tumor ID`])

# data_cluster <- data.frame("Tumor ID" = names(All.Stuff$est_c.v[[1]]), cluster = All.Stuff$est_c.v[[1]], cancer = c(rep("LUSC", 45), rep("LUAD", 88))) 
data_cluster <- data.frame("Tumor ID" = names(All.Stuff$est_c.v[[1]]), cluster = All.Stuff$est_c.v[[1]], cancer = c(rep("LUSC", 68), rep("LUAD", 171))) 
colnames(data_cluster)[1] <- "Tumor ID"

LUAD_merge <- merge(data_cluster, LUAD_subtype, by = "Tumor ID")
LUSC_merge <- merge(data_cluster, LUSC_subtype, by = "Tumor ID")
dat <- table(data_cluster$cancer, data_cluster$cluster)
table(LUAD_merge$expression_subtype, LUAD_merge$cluster)
table(LUSC_merge$...101, LUSC_merge$cluster)

test <- fisher.test(dat)
test

dat <- matrix(c(13, 6, 20,66, 2, 4,3, 80, 78,0, 3, 73),nrow = 4, byrow = TRUE)
chisq.test(dat)

#### Survival curve ###
library(survival)
library("survminer")
lung <- merge(data_cluster, clinical.slt, by = "Tumor ID")
fit <- survfit(Surv(days, censor) ~ cluster, data = lung)

summary(fit,time=c(365*6))

p <- ggsurvplot(fit, data = lung, xlab = "Survival time (in days)", pval = T) 
p$plot <- p$plot + labs(
  title    = "Kaplan Meier Survival Curves"
)
p <- ggpar(
  p,
  font.legend = c(14),
  font.title    = c(22),
  font.x        = c(20),
  font.y        = c(20),
  font.xtickslab = c(18),
  font.ytickslab = c(18),
  legend = "top"
)
p
ggsave("surv.png", height = 6, width =10)

### Logrank test pairwise ###
survdiff(Surv(days, censor) ~ cluster, data = lung)


#### Expression difference in cluster 1,2 and 3 ####
### LUAD ###
LUAD_subtype_exp <-  LUAD_subtype %>% 
  mutate(STK11 = case_when(STK11...14 != "none" ~ "active",
                           STK11...14 == "none" ~ "nonactive"),
         KRAS = case_when(KRAS...15 != "none" ~ "active", 
                          KRAS...15 == "none" ~ "nonactive"), 
         EGFR = case_when(EGFR...17 != "none" ~ "active", 
                          EGFR...17 == "none" ~ "nonactive"), 
         NF1 = case_when(NF1...28 != "none" ~ "active", 
                         NF1...28 == "none" ~ "nonactive"),
         TP53 = case_when(TP53...13 != "none" ~ "active",
                          TP53...13 == "none" ~ "nonactive"), 
         CoMutation = ifelse(NF1 == "active" & TP53 == "active", "Yes", "No"))
a<-table(LUAD_subtype_exp$STK11, LUAD_subtype_exp$expression_subtype)
a[2, ]/colSums(a)
a<-table(LUAD_subtype_exp$KRAS, LUAD_subtype_exp$expression_subtype)
a[1, ]/colSums(a)
a<-table(LUAD_subtype_exp$NF1, LUAD_subtype_exp$expression_subtype)
a[1, ]/colSums(a)
a<-table(LUAD_subtype_exp$TP53, LUAD_subtype_exp$expression_subtype)
a[1, ]/colSums(a)
a<-table(LUAD_subtype_exp$CoMutation, LUAD_subtype_exp$expression_subtype)
a[2, ]/colSums(a)
a<-table(LUAD_subtype_exp$EGFR, LUAD_subtype_exp$expression_subtype)
a[1, ]/colSums(a)

LUAD_merge_exp <- LUAD_merge %>% 
  mutate(STK11 = case_when(STK11...14 != "none" ~ "active",
                           STK11...14 == "none" ~ "nonactive"),
         KRAS = case_when(KRAS...15 != "none" ~ "active", 
                          KRAS...15 == "none" ~ "nonactive"), 
         EGFR = case_when(EGFR...17 != "none" ~ "active", 
                          EGFR...17 == "none" ~ "nonactive"), 
         NF1 = case_when(NF1...28 != "none" ~ "active", 
                         NF1...28 == "none" ~ "nonactive"),
         TP53 = case_when(TP53...13 != "none" ~ "active",
                          TP53...13 == "none" ~ "nonactive"), 
         CoMutation = ifelse(NF1 == "active" & TP53 == "active", "Yes", "No"))

table(LUAD_merge$expression_subtype, LUAD_merge$cluster)
a <- table(LUAD_merge_exp$KRAS, LUAD_merge_exp$cluster)
a[1, ]/colSums(a)
a <- table(LUAD_merge_exp$STK11, LUAD_merge_exp$cluster)
a[2, ]/colSums(a)
a <- table(LUAD_merge_exp$NF1, LUAD_merge_exp$cluster)
a[1, ]/colSums(a)
a <- table(LUAD_merge_exp$TP53, LUAD_merge_exp$cluster)
a[1, ]/colSums(a)
a <- table(LUAD_merge_exp$CoMutation, LUAD_merge_exp$cluster)
a[2, ]/colSums(a)
a <- table(LUAD_merge_exp$EGFR, LUAD_merge_exp$cluster)
a[1, ]/colSums(a)
a <- table(LUAD_merge_exp$`methylation signature`, LUAD_merge_exp$cluster)
(a[1,] + a[2, ])/colSums(a)
a/matrix(rep(colSums(a),each=3), nrow = 3)

### LUSC ###
LUSC_subtype_exp <- LUSC_subtype %>% 
  mutate(KEAP1_exp = case_when(!is.na(KEAP1...15) ~ "active", 
                               is.na(KEAP1...15) ~ "nonactive"),
         NFE2L2_exp = case_when(!is.na(NFE2L2...18) ~ "active",
                                is.na(NFE2L2...18) ~ "nonactive"), 
         PTEN_exp = case_when(!is.na(PTEN...40) ~ "active",
                              is.na(PTEN...40) ~ "nonactive"), 
         TP53_exp = case_when(!is.na(TP53...14) ~ "active",
                              is.na(TP53...14) ~ "nonactive"),
         PIK3CA_exp = case_when(!is.na(PIK3CA...33) ~ "active", 
                                is.na(PIK3CA...33) ~ "nonactive"), 
         RB1_exp = case_when(!is.na(RB1...41) ~ "active", 
                             is.na(RB1...41) ~ "nonactive"), 
         PTEN_exp = case_when(!is.na(PTEN...40) ~ "active", 
                              is.na(PTEN...40) ~ "nonactive"),
         NF1_exp = case_when(!is.na(NF1...46) ~ "active", 
                             is.na(NF1...46) ~ "nonactive"))

table(LUSC_merge$...101, LUSC_merge$cluster)
## classical 
table(LUSC_subtype_exp$KEAP1_exp, LUSC_subtype_exp$...101)
table(LUSC_subtype_exp$KEAP1...59, LUSC_subtype_exp$...101)

# this one
table(LUSC_subtype_exp$NFE2L2_exp, LUSC_subtype_exp$...101)
table(LUSC_subtype_exp$NFE2L2...62, LUSC_subtype_exp$...101)

table(LUSC_subtype_exp$PTEN_exp, LUSC_subtype_exp$...101)
table(LUSC_subtype_exp$PTEN...83, LUSC_subtype_exp$...101)

# this one
as.vector(table(LUSC_subtype_exp$SOX2...72, LUSC_subtype_exp$...101))/table(LUSC_subtype_exp$...101)

table(LUSC_subtype_exp$TP53_exp, LUSC_subtype_exp$...101)

table(LUSC_subtype_exp$PIK3CA_exp, LUSC_subtype_exp$...101)
# this one
as.vector(table(LUSC_subtype_exp$PIK3CA...76, LUSC_subtype_exp$...101))/table(LUSC_subtype_exp$...101)

## primitive 
table(LUSC_subtype_exp$RB1_exp, LUSC_subtype_exp$...101)

table(LUSC_subtype_exp$PTEN_exp, LUSC_subtype_exp$...101)
as.vector(table(LUSC_subtype_exp$PTEN...83, LUSC_subtype_exp$...101))/table(LUSC_subtype_exp$...101)

## basel 
table(LUSC_subtype_exp$NF1_exp, LUSC_subtype_exp$...101)

############## use clusters ##################
LUSC_merge_exp <- LUSC_merge %>% 
  mutate(KEAP1_exp = case_when(!is.na(KEAP1...15) ~ "active", 
                               is.na(KEAP1...15) ~ "nonactive"),
         NFE2L2_exp = case_when(!is.na(NFE2L2...18) ~ "active",
                                is.na(NFE2L2...18) ~ "nonactive"), 
         PTEN_exp = case_when(!is.na(PTEN...40) ~ "active",
                              is.na(PTEN...40) ~ "nonactive"), 
         TP53_exp = case_when(!is.na(TP53...14) ~ "active",
                              is.na(TP53...14) ~ "nonactive"),
         PIK3CA_exp = case_when(!is.na(PIK3CA...33) ~ "active", 
                                is.na(PIK3CA...33) ~ "nonactive"), 
         RB1_exp = case_when(!is.na(RB1...41) ~ "active", 
                             is.na(RB1...41) ~ "nonactive"), 
         PTEN_exp = case_when(!is.na(PTEN...40) ~ "active", 
                              is.na(PTEN...40) ~ "nonactive"),
         NF1_exp = case_when(!is.na(NF1...46) ~ "active", 
                             is.na(NF1...46) ~ "nonactive"))
## classical 
table(LUSC_merge_exp$KEAP1_exp, LUSC_merge_exp$cluster)
table(LUSC_merge_exp$KEAP1...59, LUSC_merge_exp$cluster)

table(LUSC_merge_exp$NFE2L2_exp, LUSC_merge_exp$cluster)
table(LUSC_merge_exp$NFE2L2...62, LUSC_merge_exp$cluster)

table(LUSC_merge_exp$PTEN_exp, LUSC_merge_exp$cluster)
table(LUSC_merge_exp$PTEN...83, LUSC_merge_exp$cluster)

table(LUSC_merge_exp$SOX2...72, LUSC_merge_exp$cluster)

table(LUSC_merge_exp$TP53_exp, LUSC_merge_exp$cluster)

table(LUSC_merge_exp$PIK3CA_exp, LUSC_merge_exp$cluster)
table(LUSC_merge_exp$PIK3CA...76, LUSC_merge_exp$cluster)

## primitive 
table(LUSC_merge_exp$RB1_exp, LUSC_merge_exp$cluster)

table(LUSC_merge_exp$PTEN_exp, LUSC_merge_exp$cluster)
table(LUSC_merge_exp$PTEN...83, LUSC_merge_exp$cluster)

## basel 
table(LUSC_merge_exp$NF1_exp, LUSC_merge_exp$cluster)


