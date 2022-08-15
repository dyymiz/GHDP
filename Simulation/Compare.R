##### Compare cluster accuracy #####
library(tidyverse)
library(ggplot2)

#### load GHDP ####

t.l <- c(1, 2, 4, 5)
k.l <- c(3, 4, 5)
GHDP_mt <- matrix(NA, ncol = 5)
for(tau in t.l){
  for (k in k.l) {
    for(i in 1:100){
      path_file <- paste("YourPath",paste("GHDP_0_k", k, "tau", tau, "i", i, ".RData", sep = "_"), sep = "")
      print(path_file)
      load(path_file)
      GHDP_mt <- rbind(GHDP_mt, matrix(c(k, tau, 1 - unlist(GHDP$mean.taxicab.v.e)/2, unlist(GHDP$RI), i), ncol=5))
    }
  }
}

dim(GHDP_mt)
GHDP_df <- as.data.frame(GHDP_mt[-1 ,])
rownames(GHDP_df) <- NULL
colnames(GHDP_df) <- c("H", "tau", "taxi", "RI", "rep")
GHDP_df <- GHDP_df %>% mutate(method = "Proposed")
GHDP_df_1 <- c()
for(tau in t.l){
  for (k in k.l) {
    tmp <- GHDP_df[GHDP_df$H == k & GHDP_df$tau == tau, ] 
    tmp <- tmp[!tmp$RI %in% boxplot(tmp$RI)$out, ]
    GHDP_df_1 <- rbind(GHDP_df_1, tmp)
  }
}
