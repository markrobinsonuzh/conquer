library(dplyr)

x <- read.delim("E-MTAB-3929.sdrf.txt", header = TRUE, as.is = TRUE)
xp <- x[, c("Source.Name", "Characteristics.individual.", "Characteristics.organism.part.", 
            "Characteristics.developmental.stage.", "Characteristics.treatment."), drop = FALSE]
rownames(xp) <- xp$Source.Name
write.table(xp, file = "EMTAB3929_pheno.txt", row.names = TRUE, col.names = TRUE,
            sep = "\t", quote = FALSE)

