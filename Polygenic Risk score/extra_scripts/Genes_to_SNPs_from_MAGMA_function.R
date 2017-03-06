GENES_to_snps <- scan(file = "0.5_CLOZUK_PGC_SNPs.genes.annot", what = "", sep = "\n")
y <- strsplit(GENES_to_snps, "[[:space:]]+")
names(y) <- sapply(y, '[[', 1)
y <- lapply(y, '[', -1)
y <- lapply(y, '[', -1)
a <- unlist(y)
b <- a[5:length(a)]

length(unique(b))

y[[1]] <- NULL
y[[1]] <- NULL