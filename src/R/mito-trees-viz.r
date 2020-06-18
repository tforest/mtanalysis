## Time-stamp: <2020-06-09 12:41:58 chl>

##
## Analysis of cluster of orthologous genes.
## Depends on `mito-trees.r'
##

suppressPackageStartupMessages(library(heatmaply))
library(ade4)

set.seed(101)

PLOT <- FALSE

## helper functions

## take a matrix and return a data frame with proper
## colnames and rownames
fmt <- function(d, i) {
  out <- as.data.frame(d[[i]])
  rownames(out) <- colnames(out)
  return(out[ix, ix])
}

WD <- "~/tmp"

fs <- list.files(WD, "distances-*", full.names = TRUE)

f <- lapply(fs, read.table, header = TRUE)

ix <- intersect(intersect(colnames(f[[1]]), colnames(f[[2]])), colnames(f[[3]]))

ns <- tolower(unlist(strsplit(fs, "_"))[seq(2, 6, 2)])
for (i in seq_along(ns)) assign(ns[i], fmt(f, i))

## Heatmap using plotly backend. Add file = "outfile.png" to save as image.
if (PLOT) {
  p1 <- heatmaply(sordariales, dendrogram = "column", symm = TRUE)
  p2 <- heatmaply(boletales, dendrogram = "column", symm = TRUE)
  p3 <- heatmaply(eurotiales, dendrogram = "column", symm = TRUE)
  subplot(p1, p2, p3)
}

## Mantel correlation test
cat("\n------------ Sordariales vs. Boletales ------------\n")
print(mantel.rtest(as.dist(sordariales), as.dist(boletales), nrepet = 999))
cat("\n------------ Sordariales vs. Eurotiales -----------\n")
print(mantel.rtest(as.dist(sordariales), as.dist(eurotiales), nrepet = 999))
cat("\n------------ Eurotiales vs. Boletales -------------\n")
print(mantel.rtest(as.dist(eurotiales), as.dist(boletales), nrepet = 999))

# Clustering : partitionning aroung medoids
cluster::pam(as.dist(sordariales), 3, diss = TRUE)


# CMDscale
fit <- cmdscale(sordariales)

x <- fit[,1]
y <- fit[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     #main="Orthogroups pairwise distances MDS", 
     type="p", pch=15, col=c("green"))

fit <- cmdscale(eurotiales)

x <- fit[,1]
y <- fit[,2]

points(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="p", col=c("red"), pch=15)

fit <- cmdscale(boletales)

x <- fit[,1]
y <- fit[,2]

points(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="Metric MDS", type="p", pch=15, col=c("black"))

