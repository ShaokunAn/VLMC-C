
rm(list=ls())
set.seed(100)

# Initialize 
N <- 100000 # the length of observation sequences, need to be large enough s.t. Nwu>0
R <- 1 # number of simulations
K <- 3 # K-order VLMC

# Set transfer probability ------------------------------------------------

X <- c("0", "1", "2", "3")
v1 <- c(1,2,3,4)
v2 <- c(2,5,1,3)
v3 <- c(1,2,1,4)
v4 <- c(3,2,5,1)
branches <- c("1", "2", "00", "10", "20", "30", "03", "13", "23", "033", "133", "233", "333")
trans_prob <- matrix(0, nrow = length(branches), ncol = length(X))
rownames(trans_prob) <- branches
colnames(trans_prob) <- X
trans_prob <- t(trans_prob)
trans_prob[, "00"] <- v1/sum(v1)
trans_prob[, c("03", "13", "23")] <- v2/sum(v2)
trans_prob[, c("10", "20", "30", "1", "2", "033", "133", "233")] <- v3/sum(v3)
trans_prob[, "333"] <- v4/sum(v4)
trans_prob <- t(trans_prob)


# generate observation sequence -------------------------------------------

# generate the initial states with length K
s_K <- seq(length(X)^K) # s_K is all possible states sequence of states sequence with length K
for (i in 0:(length(s_K)-1)){
  mod <- i
  tt1 <- i%%length(X)
  while(mod%/%length(X)>0){
    if(mod%/%length(X)<length(X)){
      tt1 <- c(mod%/%length(X), tt1)
      break()
    }else{
      mod <- mod%/%length(X)  
      tt1 <- c(mod%%length(X),tt1)
    }
  }
  tt1 <- as.character(tt1)
  tt1 <- paste(tt1, collapse = "")
  tt1 <- paste("000000", tt1, sep="")
  tt1 <- substr(tt1, start = nchar(tt1)-K+1, stop = nchar(tt1))
  s_K[i+1] <- tt1
}

s0 <- sample(s_K, R, replace = TRUE)
S_all <- s0
pos <- K # record the number of generated samples
s_k <- seq(R) # s_k stores the last k(maximum to K) states of generated samples of R sequences
for (n in 1:N){
  
  if(n %% 100 ==0){
    print(paste("generate ",n, "-th sample", sep=""))
  }
  
  flag <- rep(0, R)
  flag_old <- flag
  for(k in 1:K){
    s_k[flag==0] <- substr(S_all[flag==0], start = (pos-k+1), stop = pos)
    flag[flag==0] <- as.integer(s_k[flag==0]%in%branches) # update flag
    
    if(sum(flag-flag_old)==0)
      next()
    
    tt_trans <- trans_prob[s_k[flag-flag_old==1], ] # geneate new sample for sequences[flag-flag_old==1]
    for(j in 1:sum(flag-flag_old)){
      if(sum(flag-flag_old)==1){
        s_new <- sample(X, 1, replace = TRUE, prob = tt_trans)
      }else{
        s_new <- sample(X, 1, replace = TRUE, prob = tt_trans[j,])  
      }
      S_all[which(flag-flag_old==1)[j]] <- paste(S_all[which(flag-flag_old==1)[j]], s_new, sep='')
    }
    if(sum(flag)==R){
      pos <- pos+1
      break()
    }
    flag_old <- flag
    if(k==K){
      stop("ERROR: Not all branches are founded!")
    }
  }
}

saveRDS(S_all, "./simulate_data.RDS")

