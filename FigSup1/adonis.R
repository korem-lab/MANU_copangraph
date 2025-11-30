options(warn = -1) # Suppress all warnings
library(vegan)

library(dplyr)


get_x <- function(x_fl, is_dist) {
  if (is_dist == "dist") {
    x <- read.csv(x_fl, row.names = 1, check.names = FALSE)
    return( as.dist(as.matrix(x)))
  } else {
    x <- read.csv(x_fl, row.names = 1, check.names = FALSE)
  }
  return(x)
}

get_outcome <- function(o_fl) {
  o <- read.csv(o_fl, row.names = 1, check.names = FALSE)
  return(o)
}

set_X_order_to_outcome <- function(X, outcome) {
  return(X[rownames(outcome),])
}

get_adonis <- function(X, group, is_factor, dist) {
  if (is_factor) {
    group = factor(group)
  }

  if (dist == "dist"){
    #print("Using dist")
    res = adonis2(X ~ group, method=m, permutations=4999)
    #print(res)
    return(res)
  }
  methods = c('jaccard', 'bray')
  res_list = list()
      for (m in methods) {
        X_filt = X
        non_zero_prop = colSums(X_filt != 0)/nrow(X_filt)
        res = adonis2(X_filt ~ group, method=m, permutations=9999)
        res = cbind(res, method = m)
        res_list[[length(res_list)+1]] <- res
      }
  res_final = do.call(rbind, res_list)
  res_final = cbind(res_final, dof = sub("\\d+$", "", row.names(res_final)))
  rownames(res_final) = seq_len(nrow(res_final))
  #print(res_final)
  return(res_final)
}


args <- commandArgs(trailingOnly =  TRUE)
#print(args)

if (length(args) >= 2) {
    x <- get_x(args[1], args[4])
    outcome <- get_outcome(args[2])
    if (args[3] == "factor") {
    #  print('outcome treated as factor')
      is_factor <- TRUE
    } else {
    #  print('outcome treated as continuous')
      is_factor <- FALSE
    }
    #print(rownames(outcome))
    #print(rownames(x))
    if (args[4] != "dist") {
      x <- set_X_order_to_outcome(x, outcome)
    }
    adonis_results <- get_adonis(x, outcome$outcome, is_factor, args[4])
    write.csv(adonis_results, paste0(sub("\\.csv$", "", args[1]), "_adonis.csv"), row.names = FALSE, quote = FALSE)
}
