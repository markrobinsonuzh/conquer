library(dplyr)
x <- read.delim("EMTAB2805_SraRunInfo.csv", header = TRUE, as.is = TRUE,
                sep = ",")
xp <- x[, "LibraryName", drop = FALSE] %>% 
  tidyr::separate(LibraryName, into = c("cell_cycle_stage", "cell"), sep = "_", 
                  remove = FALSE)
rownames(xp) <- xp$LibraryName
write.table(xp, file = "EMTAB2805_pheno.txt", row.names = TRUE, col.names = TRUE,
            sep = "\t", quote = FALSE)
