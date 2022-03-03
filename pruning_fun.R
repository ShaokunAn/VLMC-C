pruneing_fun <- function(node_name, obs_seq, X, alpha=0.95, num_simu=100000){
  ##### pruning_fun compute a logical value to descide whether to prune a node away from the tree
  ##### it estimates covariance matrix with plug-in methods, computes its corresponding eigenvalues 
  ##### with asymptotic distribution, alpha*100% quantiles are computed based on sampling from the given distribution
  # input: 
  #   node_name: a (forward-order) string of given branch uw, u\in X, w \in bigcup X^m
  #   obs_seq: a string consisting a time series observational sequences
  #   X: the state space
  #   alpha: quantile of asymptotic distribution, default 0.95
  #   num_simu: number of samples from asymptotic distribtution to approximate cutoff value
  # output:
  #   IfPruned: a logical value to determine if prune the test branch, i.e. prune branch from uw to w
  
  
  # compute covariance matrix based on plug-in estimators
  set.seed(1)
  S <- obs_seq
  k <- nchar(node_name)
  
  seg_k <- substring(S, 1:(nchar(S)-k+1), k:nchar(S)) # cut observational sequence to strings with length k
  seg_km1 <- substring(S, 1:(nchar(S)-k+2), (k-1):nchar(S)) # cut observational sequence to strings with length k-1
  seg_kp1 <- substring(S, 1:(nchar(S)-k), (k+1):nchar(S)) # cut observational sequence to strings with length k+1
  
  w <- substr(node_name, 2, nchar(node_name))
  u <- substr(node_name, 1, 1)
  uw <- node_name
  
  if(length(which(seg_k==uw))==length(which(seg_km1==w))){
    warning("Elements in covariance matrix will all be 0 since Nuw=Nw! branch will be directly pruned away without testing. Enlarge THRESHOLD.GEN might help to remove this issue.")
    out <- list(ifpruned=TRUE, Delta_quan=1)
    return(out)
  }
  
  hat_Puw <- length(which(seg_k==uw))/nchar(S)
  hat_Pw <- length(which(seg_km1==w))/nchar(S)
  
  Nuwv <- c()
  Nwv <- c()
  l <- 1
  for(i in 1:length(X)){
    tt <- length(which(seg_kp1==paste(uw, X[i], sep='')))
    # only record the string whose number of occurence larger than 0
    # and consider their joint distributions
    if(tt>0){
      Nuwv[l] <- tt
      names(Nuwv)[l] <- X[i]
      Nwv[l] <- length(which(seg_k==paste(w, X[i], sep=''))) 
      names(Nwv)[l] <- X[i]
      l <- l+1
    }
  }
  hat_Puwv <- Nuwv/nchar(S)
  hat_Pwv <- Nwv/nchar(S)
  
  sigma <- matrix(0, length(Nuwv), length(Nuwv))
  for(i in 1:nrow(sigma)){
    for (j in 1:ncol(sigma)) {
      sigma[i,j] <- as.integer(i==j)/hat_Puwv[i]-1/hat_Puw-as.integer(i==j)/hat_Pwv[i]+1/hat_Pw
      sigma[i,j] <- sqrt(hat_Puwv[i]*hat_Puwv[j])*sigma[i,j]
    }
  }
  rownames(sigma) <- names(hat_Puwv)
  colnames(sigma) <- names(hat_Puwv)
  
  # compute its eigenvalues
  eigen_sigma <- eigen(sigma)
  eigen_sigma <- eigen_sigma$values
  eigen_sigma[eigen_sigma<0] <- 0
  
  # compute Delta value
  
  Delta <- 0
  # to ensure sum_x N(uwx)==N(uw) and sum_x N(wx)=N(w), the computation of N(uw) and N(w) need little change
  seg_Nuw <- substring(S, 1:(nchar(S)-k), k:(nchar(S)-1))
  Nuw <- length(which(seg_Nuw==uw))
  seg_Nw <- substring(S, 1:(nchar(S)-k+1), (k-1):(nchar(S)-1))
  Nw <- length(which(seg_Nw==w))
  # Nuwv same as defined above
  # Nwv same as defined above
  
  Nwv <- as.numeric(Nwv)
  Nuw <- as.numeric(Nuw)
  Nw <- as.numeric(Nw)
  
  Euwv <- Nwv*Nuw/Nw
  Delta <- sum((Nuwv-Euwv)^2/Euwv)
  
  # sampels from the given distribution to estimate its alpha*100% quantile
  gauss_mat <- matrix(0, ncol=num_simu, nrow=nrow(sigma)) # each row contains samples from a gauss distribution
  for(i in 1:nrow(sigma)){
    gauss_mat[i,] <- rnorm(num_simu)
  }
  gauss_mat <- gauss_mat^2
  sample_weight_chisq <- apply(eigen_sigma*gauss_mat, 2, sum)
  
  cutoff <- sort(sample_weight_chisq, decreasing = FALSE)[floor(alpha*num_simu)]
  
  # if Delta < cutoff, return TRUE and prune the node
  out_logi <- (Delta<cutoff)
  # compute the quantile of statistic Delta
  out_delta_quan <- c(Delta, sample_weight_chisq)
  out_delta_quan <- which(sort(out_delta_quan, decreasing = FALSE)==Delta)
  out_delta_quan <- out_delta_quan[1]/(length(sample_weight_chisq)+1)
  
  if(out_logi)
    print(paste("Testing branch", node_name, '...... and prune it'))
  else
    print(paste("Testing branch", node_name, '...... and keep it'))
  
  out <- list(ifpruned=out_logi, Delta_quan=out_delta_quan)
  return(out)
  
}