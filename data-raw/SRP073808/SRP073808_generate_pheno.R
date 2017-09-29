library(dplyr)
x <- read.delim("SRP073808_SraRunInfo.csv", header = TRUE, as.is = TRUE,
                sep = ",")
xp <- x[, c("Run", "LibraryName"), drop = FALSE]
xp$LibraryName <- gsub("_single-cell: In vitro cultured H7 human embryonic stem cells \\(WiCell\\) and H7-derived downstream early mesoderm progenitors", "", xp$LibraryName)
rownames(xp) <- xp$Run
write.table(xp, file = "SRP073808_pheno.txt", row.names = TRUE, col.names = TRUE,
            sep = "\t", quote = FALSE)
