tax = read.csv("taxonomy.csv", header = FALSE)
counts = read.csv("counts.csv", header = TRUE)
merge(x = tax, y = counts, by.x = "V9", by.y ="spcies_id", all.y = TRUE)
apply(tax$V9, FUN=str_trim)

