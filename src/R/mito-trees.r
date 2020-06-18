## Time-stamp: <2020-06-08 19:47:42 chl>

##
## Analysis of cluster of orthologous genes.
## Input data are assumed to be located in the same directory
## (here, "~/tmp"), unzipped, with a common prefix ("ARBRES_").
## Output data (distance matrices) are written in the same
## directory.
##

library(phangorn)

DEBUG <- TRUE

## helper functions
tip_order <- function(x) unlist(lapply(strsplit(x$tip.label, "~"),
                               function(x) x[1]))
tip_name <- function(x) lapply(strsplit(x$tip.label, "~"), function(x) x[2])

species <- c("SORDARIALES", "BOLETALES", "EUROTIALES")

for (s in species) {
  ## NOTE: rename suffix _TREES to prefix ARBRES_
  wd <- paste("~/tmp/ARBRES", s, sep = "_")
  fs <- list.files(wd, pattern = "\\.nw$", full.names = TRUE)
  n <- length(fs)


  ts <- list()
  o <- vector("character", n)
  for (i in seq_along(fs)) {
    ts[[i]] <- read.tree(fs[i])
    o[i] <- unique(tip_order(ts[[i]]))
    ts[[i]]$tip.label <- tip_name(ts[[i]])
  }

  cc <- matrix(NA, nr = n, nc = n)
  for (i in 1:n) {
    for (j in 1:n) {
      matches <- unlist(intersect(ts[[i]]$tip.label, ts[[j]]$tip.label))
      len <- length(matches)
      if (len > 5) {
        values <- wRF.dist(keep.tip(ts[[i]], matches),
                          keep.tip(ts[[j]], matches),
                          normalize = TRUE)
        cc[i, j] <- values
      }
    }
  }

  if (DEBUG) cat(s, length(o), "-", paste(o, collapse = " "), "\n")
  colnames(cc) <- rownames(cc) <- o
  write.table(cc, file = paste("~/tmp/distances", s, sep = "_"),
              col.names = TRUE, row.names = FALSE)
}

## cc[upper.tri(cc, diag = TRUE)] <- NA

## This won't work since we need to work on common tips
## mts <- do.call(c, ts)
## wRF.dist(mts)
