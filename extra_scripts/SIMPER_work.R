install.packages('vegan')
install.packages('cluster')
library('vegan')
simper
simper2<-function (comm, group, permutations = 0, trace = FALSE, parallel = getOption("mc.cores"), 
                   ...) 
{browser()
  EPS <- sqrt(.Machine$double.eps)
  if (any(rowSums(comm, na.rm = TRUE) == 0)) 
    warning("you have empty rows: results may be meaningless")
  pfun <- function(x, comm, comp, i, contrp) {
    groupp <- group[perm[x, ]]
    ga <- comm[groupp == comp[i, 1], , drop = FALSE]
    gb <- comm[groupp == comp[i, 2], , drop = FALSE]
    n.a <- nrow(ga)
    n.b <- nrow(gb)
    for (j in seq_len(n.b)) {
      for (k in seq_len(n.a)) {
        mdp <- abs(ga[k, , drop = FALSE] - gb[j, , drop = FALSE])
        mep <- ga[k, , drop = FALSE] + gb[j, , drop = FALSE]
        contrp[(j - 1) * n.a + k, ] <- mdp/sum(mep)
      }
    }
    colMeans(contrp)
  }
  comm <- as.matrix(comm)
  comp <- t(combn(unique(as.character(group)), 2))
  outlist <- NULL
  P <- ncol(comm)
  nobs <- nrow(comm)
  perm <- getPermuteMatrix(permutations, nobs, ...)
  if (ncol(perm) != nobs) 
    stop(gettextf("'permutations' have %d columns, but data have %d rows", 
                  ncol(perm), nobs))
  nperm <- nrow(perm)
  if (nperm > 0) 
    perm.contr <- matrix(nrow = P, ncol = nperm)
  if (is.null(parallel)) 
    parallel <- 1
  hasClus <- inherits(parallel, "cluster")
  isParal <- hasClus || parallel > 1
  isMulticore <- .Platform$OS.type == "unix" && !hasClus
  if (isParal && !isMulticore && !hasClus) {
    parallel <- makeCluster(parallel)
  }
  for (i in seq_len(nrow(comp))) {
    group.a <- comm[group == comp[i, 1], , drop = FALSE]
    group.b <- comm[group == comp[i, 2], , drop = FALSE]
    n.a <- nrow(group.a)
    n.b <- nrow(group.b)
    contr <- matrix(ncol = P, nrow = n.a * n.b)
    for (j in seq_len(n.b)) {
      for (k in seq_len(n.a)) {
        md <- abs(group.a[k, , drop = FALSE] - group.b[j, 
                                                       , drop = FALSE])
        me <- group.a[k, , drop = FALSE] + group.b[j, 
                                                   , drop = FALSE]
        contr[(j - 1) * n.a + k, ] <- md/sum(me)
      }
    }
    average <- colMeans(contr)
    if (nperm > 0) {
      if (trace) 
        cat("Permuting", paste(comp[i, 1], comp[i, 2], 
                               sep = "_"), "\n")
      contrp <- matrix(ncol = P, nrow = n.a * n.b)
      if (isParal) {
        if (isMulticore) {
          perm.contr <- mclapply(seq_len(nperm), function(d) pfun(d, 
                                                                  comm, comp, i, contrp), mc.cores = parallel)
          perm.contr <- do.call(cbind, perm.contr)
        }
        else {
          perm.contr <- parSapply(parallel, seq_len(nperm), 
                                  function(d) pfun(d, comm, comp, i, contrp))
        }
      }
      else {
        perm.contr <- sapply(1:nperm, function(d) pfun(d, 
                                                       comm, comp, i, contrp))
      }
      p <- (rowSums(apply(perm.contr, 2, function(x) x >= 
                            average - EPS)) + 1)/(nperm + 1)
    }
    else {
      p <- NULL
    }
    overall <- sum(average)
    sdi <- apply(contr, 2, sd)
    ratio <- average/sdi
    ava <- colMeans(group.a)
    avb <- colMeans(group.b)
    ord <- order(average, decreasing = TRUE)
    cusum <- cumsum(average[ord]/overall)
    out <- list(species = colnames(comm), average = average, 
                overall = overall, sd = sdi, ratio = ratio, ava = ava, 
                avb = avb, ord = ord, cusum = cusum, p = p)
    outlist[[paste(comp[i, 1], "_", comp[i, 2], sep = "")]] <- out
  }
  if (isParal && !isMulticore && !hasClus) 
    stopCluster(parallel)
  attr(outlist, "permutations") <- nperm
  attr(outlist, "control") <- attr(perm, "control")
  class(outlist) <- "simper"
  outlist
}

data(agriculture)
## Example 1 in ref:
##  Dissimilarities using Euclidean metric and without standardization
d.agr <- daisy(agriculture, metric = "euclidean", stand = FALSE)
d.agr
as.matrix(d.agr)[,"DK"] # via as.matrix.dist(.)
## compare with
as.matrix(daisy(agriculture, metric = "gower"))