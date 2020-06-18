tax = read.csv("data/taxonomy.csv", header = FALSE)
counts = read.csv("data/counts.csv", header = TRUE)

library(tidyverse)
tax = data.frame(apply(tax, 2, str_trim))

merged = merge(x = tax, y = counts, by.x = "V9", by.y ="spcies_id", all.y = TRUE)
merged = na.omit(merged)

by_group <- merged %>% group_by(V3)
res = by_group %>% summarise(
  disp = mean(prop),
  fwd = sum(nb_fwd),
  rev = sum(nb_rev)
)

comp_mat = as.matrix(cbind(res$disp, res$fwd))
comp_mat = t(comp_mat)
rownames(comp_mat) = res$V3
heatmap(comp_mat)
matcor = cor(comp_mat, )

comp_mat2 = as.matrix(res$disp)
rownames(comp_mat2) = res$V3
comp_mat2 = comp_mat2[-c(3)]
barplot(log(t(scale(comp_mat2))), las=2)

chi = t(as.matrix(res))
colnames(chi) = chi[1,]
chi = chi[-c(1,2),-c(1)]
chi = apply(chi, 2, as.numeric)
chisq.test(chi)
