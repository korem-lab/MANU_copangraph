library(vegan)
library(dplyr)


get_x <- function(x_fl) {
  x <- read.csv(x_fl, row.names = 1, check.names = FALSE) %>%
    t()
  return(x)
}

get_outcome <- function(o_fl) {
  o <- read.csv(o_fl, row.names = 1)
  return(o)
}

set_X_order_to_outcome <- function(X, outcome) {
  return(X[rownames(outcome),])
}

get_adonis <- function(X, group) {
  props = c(0, 0.1)
  rkbm_ts = c(10, 50)
  methods = c('euclidean', 'bray')
  res_list = list()
  for (p in props) {
    for (t in rkbm_ts) {
      for (m in methods) {
        X_filt = X
        X_filt[X_filt < t] <- 0
        non_zero_prop = colSums(X_filt != 0)/nrow(X_filt)
        X_filt = X_filt[, non_zero_prop > p]
        print(dim(X))
        print(dim(X_filt))
        res = adonis2(X_filt ~ group, method=m, permutations=4999)
        res = cbind(res, rbkm_t = t, nonz_prop = p, method = m)
        res_list[[length(res_list)+1]] <- res
      }
    }
  }
  res_final = do.call(rbind, res_list)
  res_final = cbind(res_final, dof = sub("\\d+$", "", row.names(res_final)))
  rownames(res_final) = seq_len(nrow(res_final))
  print(res_final)
  return(res_final)
}


args <- commandArgs(trailingOnly =  TRUE)
print(args)

if (length(args) >= 2) {
    x <- get_x(args[1])
    outcome <- get_outcome(args[2])
    print(rownames(outcome))
    x <- set_X_order_to_outcome(x, outcome)
    print(rownames(x))
    adonis_results <- get_adonis(x, outcome$outcome)
    write.csv(adonis_results, paste0(sub("\\.csv$", "", args[1]), "_adonis.csv"), row.names = FALSE, quote = FALSE)
}