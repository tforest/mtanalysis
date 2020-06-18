## Time-stamp: <2020-06-13 19:45:11 chl>

##
## Analysis of synteny.
## Input data = ~/tmp/synteny.txt
##

library(ape)
library(fpc)
library(ggplot2)
library(gggenes)

set.seed(101)

WD <- "~/tmp"

s <- scan(paste(WD, "synteny_oriented.txt", sep = "/"), what = "character")
meta <- read.csv("metadata.csv", header = FALSE)

meta$V9 <- trimws(meta$V9)

g <- list()
for (r in seq(along = s)) {
  values <- unlist(strsplit(s[r], ","))
  len <- length(values)
  tmp <- values[2:len]
  genes <- substr(tmp, 2, 10)
  pfx <- substr(tmp, 0, 1)
  ## HACK remove duplicates
  dup <- duplicated(genes)
  g[[r]] <- genes[!dup]
  attr(g[[r]], "strand") <- pfx[!dup]
  attr(g[[r]], "order") <- trimws(meta[meta$V9 %in% values[1], "V5"])
  names(g)[r] <- values[1]
}

genes <- unique(unlist(g))

ilen <- unlist(lapply(g, length))
table(ilen)

## remove low frequency item
ix <- names(which(ilen < 10 | ilen > 17))
g <- within(g, rm(list = ix))

## remove duplicates (FIX)
ix <- duplicated(names(g))
g <- g[!ix]

## fill an empty data frame with available genes as columns, where
## each cell is a number reflecting position number of each gene.
dg <- as.data.frame(setNames(replicate(length(genes), numeric(length(g)),
                                       simplify = FALSE), genes),
                    row.names = names(g))

for (k in seq(along = g)) {
  dg[k,] <- NA
  dg[k, g[[k]]] <- seq(1, length(g[[k]]))
}

# filter on order
o <- unlist(lapply(g, function(x) attributes(x)$order))
ix <- which(o %in% c("Sordariales", "Eurotiales", "Boletales"))
sel <- names(o[ix])
clr <- as.numeric(as.factor(o[ix]))

## replace identifers by species name + id
nams = character(length(sel))
for(i in seq(length(sel)))
  nams[i] <- trimws(paste((meta[which(sel[i] == meta[, 9]), 7]), sel[i]))

dg_sel = dg[sel,]

rownames(dg_sel) <- nams

## Manhattan distance, or Kendall's distance for ranks
dd <- dist(dg_sel, method = "manhattan")
hc <- hclust(dd)

## Cutting at 3 doesn't work well, but 4 is better, yet
## it would be nice to explain why NC 026920 -- NC 001329
## cluster together, and apart from other...
cl <- cutree(hc, 4)
table(cl, clr)

plot(as.phylo(hc), tip.color = clr, cex = .5)

## same as above, using model validation for selecting
## the optimal number of clusters.
## For an explanation of how it works, see, e.g.:
## https://aliquote.org/post/using-bootstrap-in-cluster-analysis/
hc.boot <- clusterboot(dd, B = 500, bootmethod = "boot",
                      clustermethod = disthclustCBI, k = 4,
                      method = "complete", cut = "number")

hc.boot
## The cluster with 4 taxons discussed above appears to be
## less reliable than other, so maybe this is just an artifact,
## although it will hardly be merged with any other cluster.

## dataviz based on {gggenes}
## Requires a data frame like this:
##    molecule  gene  start    end  strand direction
## 1   Genome5  genA 405113 407035 forward         1
## 2   Genome5  genB 407035 407916 forward         1

d <- data.frame(taxon = NA, gene = NA, start = NA, end = NA,
               strand = NA, direction = NA, order = NA)
start <- 1
for (j in seq(along = g)) {
  len <- length(g[[j]])
  stop <- start + len - 1
  st <- attr(g[[j]], "strand")
  ord <- attr(g[[j]], "order")
  ## HACK because of some mising information
  if (length(ord) == 0) ord <- ""
  d[start:stop, 1] <- names(g)[j]
  d[start:stop, 2] <- as.character(unlist(g[[j]]))
  d[start:stop, 3] <- seq(1, len)
  d[start:stop, 4] <- seq(1, len) + 1
  d[start:stop, 5] <- "forward"
  d[start:stop, 6] <- ifelse(st == "+", 1, 2)
  d[start:stop, 7] <- rep(ord, len)
  start <- start + len
}

d$taxon <- factor(d$taxon)
d$gene <- factor(d$gene)
d$order <- factor(d$order)

## Sanity check:
## sum(unlist(lapply(g, length))) == nrow(d)

ds <- subset(d, taxon %in% sel)

# names + IDs
nams = character(length(ds$taxon))
for(i in seq(length(ds$taxon)))
  nams[i] <- trimws(paste((meta[which(meta[, 9] == ds$taxon[i]), 7]), ds$taxon[i]))
ds$taxon <- nams

ds$taxon <- factor(ds$taxon, levels = unique(ds[order(ds$order), "taxon"]))
ds$gene <- droplevels(ds$gene)
ds$order <- droplevels(ds$order)
clr = nlevels(ds$gene)
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))



## all species, sorted by order
p <- ggplot(ds, aes(xmin = start, xmax = end,
                   y = taxon, fill = gene, forward = direction))
p <- p + geom_gene_arrow(arrowhead_height = unit(2, "mm"),
                        arrowhead_width = unit(1, "mm"),
                        arrow_body_height = unit(2, "mm")) +
  scale_fill_manual("", values = getPalette(clr)) +
  theme_genes() +
  theme(axis.text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(x = NULL, y = NULL)
p

## all species, facetted by order
p <- ggplot(ds, aes(xmin = start, xmax = end,
                   y = taxon, fill = gene, forward = direction))
p <- p + geom_gene_arrow(arrowhead_height = unit(2, "mm"),
                        arrowhead_width = unit(1, "mm"),
                        arrow_body_height = unit(2, "mm")) +
  scale_fill_manual("", values = getPalette(clr)) +
  facet_wrap(~ order, scales = "free", ncol = 1) +
  theme_genes() +
  theme(axis.text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(x = NULL, y = NULL)
p

## ggsave("~/Desktop/case2.pdf", width = 6, height = 6)

