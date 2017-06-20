library(dplyr)
x <- read.delim("EGEUV1_SraRunTable.txt", header = TRUE, as.is = TRUE)
xp <- x %>% select(Experiment_s, Center_Name_s, population_s) %>%
  rename(Experiment = Experiment_s, population = population_s, 
         center = Center_Name_s)
rownames(xp) <- xp$Experiment
write.table(xp, file = "EGEUV1_pheno.txt", row.names = TRUE, col.names = TRUE,
            sep = "\t", quote = FALSE)
