#! /usr/bin/R
rm(list=ls())
set.seed(100)
# N <- 1000000 # the length of observation sequences, need to be large enough s.t. Nwu>0
# N <- 500000
N <- 100000
# R <- 1000 # number of simulations
R <- 10
# R <- 5
# K <- 3 # K-order VLMC
K <- 5

U <- "0" #default u in wu w.l.o.g.

# Set transfer probability ------------------------------------------------

X <- c("0", "1", "2", "3")
# v1 <- c(1,2,3,4)
# v2 <- c(2,4,1,3)
# v3 <- c(1,2,1,3)
# v4 <- c(3,2,5,1)

# v1 <- c(1,2,3,10)
# v2 <- c(2,10,1,3)
# v3 <- c(1,2,1,10)
# v4 <- c(3,2,10,1)

v1 <- c(1,2,3,4)
v2 <- c(2,5,1,3)
v3 <- c(1,2,1,4)
v4 <- c(3,2,5,1)
branches <- c("1", "2", "00", "10", "20", "30", "03", "13", "23", "133", "233", "333", "0033", "2033", "3033", "01033", "11033", "21033", "31033")
trans_prob <- matrix(0, nrow = length(branches), ncol = length(X))
rownames(trans_prob) <- branches
colnames(trans_prob) <- X
trans_prob <- t(trans_prob)
trans_prob[, "00"] <- v1/sum(v1)
trans_prob[, c("03", "13", "23", "0033", "21033")] <- v2/sum(v2)
trans_prob[, c("10", "20", "30", "1", "2", "133", "233", "2033", "31033")] <- v3/sum(v3)
trans_prob[, c("333", "3033", "01033", "11033")] <- v4/sum(v4)
trans_prob <- t(trans_prob)


# generate observation sequence -------------------------------------------

# generate the initial states with length K
s_K <- seq(length(X)^K) #s_K is all possible states sequence of states sequence with length K
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
s_k <- seq(R) #s_k stores the last k(maximum to K) states of generated samples of R sequences
for (n in 1:N){
  print(paste("generate ",n, "-th sample", sep=""))
  
  flag <- rep(0, R)
  flag_old <- flag
  for(k in 1:K){
    s_k[flag==0] <- substr(S_all[flag==0], start = (pos-k+1), stop = pos)
    flag[flag==0] <- as.integer(s_k[flag==0]%in%branches) #update flag
    
    if(sum(flag-flag_old)==0)
      next()
    
    tt_trans <- trans_prob[s_k[flag-flag_old==1], ] #geneate new sample for sequences[flag-flag_old==1]
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

saveRDS(S_all, "./new_S_all_simudata_K=5_N=1e5.RDS")


rm(list=ls())
library(VLMC)
library(data.tree)
# source("plot results.R")
source("pruning_fun.R")
source("predict_fun.R")
# S_all <- readRDS('./S_all_simudata.RDS')
S_all <- readRDS("./new_S_all_simudata_K=5_N=1e5.RDS")
X <- as.character(0:3)
THRESHOLD.GEN <- max(2, round(sqrt(nchar(S_all[1]))/length(X)))
M <- 10


# Generate results of vlmc-corrected --------------------------------------
# time0 <- Sys.time()
si_seq <- c(1:M)
tau_list <- list()
ALPHA_list <- list()

for(si in si_seq){
  # print(paste("S",si, sep=""))
  S <- S_all[si] # given observational sequences
  
  # perform substring(strings, 1:?,?:?) during each iteration. 
  # in each iteration, we can construt a deeper level for the current level, consisting children of 
  # nodes in the current level. WIDTH-FIRST!
  k <- 1
  # t_begin <- Sys.time()
  tau0 <- Node$new("x1")
  tt_record_subtree_old <- list()
  while (1){
    
    tt_seg <- substring(S, 1:(nchar(S)-k+1), k:nchar(S))
    tt_seg_unique <- sort(unique(tt_seg))
    tt_seg_count <- table(tt_seg)
    tt_seg_count <- tt_seg_count[tt_seg_count>=THRESHOLD.GEN]
    
    if(length(tt_seg_count)==0) # no more strings with N(w)>=2
      break()
    
    if(k>1){
      tt_tree_dataframe <- data.frame(from=substr(names(tt_seg_count), 2, k), # record sequence in w1w2...w_{t-1}w_t forward direction of time
                                      to=names(tt_seg_count),
                                      num=as.integer(tt_seg_count)
      ) # the new added state(the initial state of the state string) is put on top
      tt_tree_dataframe$pathString <- paste(tt_tree_dataframe$from, tt_tree_dataframe$to, sep="/")
      
      tt_unique_row <- unique(tt_tree_dataframe$from)
      l <- 1
      # add new nodes via following for loop
      tt_record_subtree_new <- list() # record new added nodes in current iteration
      for(i in 1:length(tt_unique_row)){
        
        tt_subtree_data <- tt_tree_dataframe[tt_tree_dataframe[,1]==tt_unique_row[i],] # compute 
        # tt_subtree <- Node$new(tt_unique_row[i])
        # construct each child of subtree 
        for(j in 1:nrow(tt_subtree_data)){
          tt_record_subtree_new[[l]] <- tt_record_subtree_old[[tt_subtree_data[j,1]]]$AddChild(tt_subtree_data[j,2], num=tt_subtree_data[j, "num"])
          names(tt_record_subtree_new)[l] <- tt_record_subtree_new[[l]]$name
          l <- l+1
        }
      }
      tt_record_subtree_old <- tt_record_subtree_new
      
    }else{
      for(i in 1:length(tt_seg_count)){
        tt_record_subtree_old[[i]] <- tau0$AddChild(names(tt_seg_count)[i], num=as.integer(tt_seg_count[i]))
        names(tt_record_subtree_old)[i] <- tt_record_subtree_old[[i]]$name
      }
    }
    k <- k+1
  }
  # t_initialtree <- Sys.time()
  # print(paste("time for constructing initial tree: ", t_initialtree-t_begin, sep=""))
  
  # prune from initial tree
  tau <- tau0
  tau$Set(tested_flag=FALSE) # all nodes are set as flag_node=FALSE, including internal nodes
  tau$Set(len=as.integer(tau$Get("level")-1)) # set depth for all nodes including internal nodes
  tau$Set(Delta_quan=-1) # set Delta_quan to record the quantile of Delta for kept branches
  # Delta_quan of kept nodes (leaves) are their quantiles while for untested nodes(not leaves) are noted by -1
  tau0_maxlen <- max(tau$Get('len'))
  tau0_leafcount <- tau$leafCount
  tau0_totalcount <- tau$totalCount
  print("tau0:")
  print(paste('maximum length of contexts:', tau0_maxlen))
  print(paste('number of leaf counts:', tau0_leafcount))
  print(paste('number of total counts:', tau0_totalcount))
  
  # t_prunestart <- Sys.time()
  iter <- 1
  MAX_HEIGHT <- tau$height - 1
  ALPHA_seq <- c()
  while (1) {
    # test branches with same length during each iteration
    # leafcount may not decrease when prune away one node, since if its structure is a-b, cut off b, a turns to new leaf
    # and the leafcount does not change
    print(paste("left leave branches:", tau$leafCount, "iteration:", iter))
    tt_len <- tau$Get("len", filterFun=function(x) !x$tested_flag & x$isLeaf)
    tt_max_len <- max(tt_len)
    tt_testing_nodes <- tt_len[tt_len==tt_max_len]
    
    
    ALPHA <- 1-0.05/(tau$leafCount*MAX_HEIGHT) # ensure the nodes with same length (from same level) use same quantiles
    ALPHA_seq[iter] <- ALPHA
    # if set alpha=1-0.05/tau$leafCount in function as before, the order of hypothesis test (pruning process)
    # will affect results since the nodes in same level use various alpha values (various cutoff values)
    
    for(i in 1:length(tt_testing_nodes)){
      tt_testing_node <- names(tt_testing_nodes)[i]
      tt_testing_result <- pruneing_fun(node_name = tt_testing_node,
                                        obs_seq = S, X=X, alpha=ALPHA)
      
      # if tt_testing_logi==TRUE, prune the node, which equals to prune the branch
      if(tt_testing_result$ifpruned){
        # mark the testing node to update its state, i.e. {prune} or {kept and update its tested_flag=TRUE}
        tt_prune_pos <- rep(FALSE, tau$totalCount)
        tt_leaves_name <- tau$Get("name")
        tt_prune_pos[which(tt_leaves_name==tt_testing_node)] <- TRUE
        tau$Set(prune_tmp=tt_prune_pos) # data tree Get and Set values in same order
        Prune(tau, function(x) !x$prune_tmp)
      }else{
        #if not prune the node, set the its tested_flag=TRUE, in case of duplicate test
        # recored the quantile of its corresponding statistic Delta 
        tt_flag <- tau$Get('tested_flag')
        tt_flag[tt_testing_node] <- TRUE
        tau$Set(tested_flag=tt_flag)
        tt_quan <- tau$Get("Delta_quan")
        tt_quan[tt_testing_node] <- tt_testing_result$Delta_quan
        tau$Set(Delta_quan=tt_quan)
      }
      
    }
    tt_flag <- tau$Get("tested_flag", filterFun = isLeaf)
    if(all(tt_flag) | tau$count==0){
      # print("Pruning process done!")
      break()
    }
    iter <- iter+1
  }
  
  # t_prundend <- Sys.time()
  # print(paste("time for pruning process: ", t_prundend-t_prunestart))
  # print("---------------------------")
  
  tau_list[[si]] <- tau
  ALPHA_list[[si]] <- ALPHA_seq
  # save.image(paste('vlmc_alpha_S', si, "_RelSameAlpha005.RData", sep=''))
  
}

time_tree <- Sys.time()

out <- list(tau=tau_list, ALPHA = ALPHA_list)
saveRDS(out, file = "./vlmc_simudata_K=5_R=10_N=5*1e5_VLMC-C.RDS")



# compute stationary distribution to analyze why some context are missing
trans_prob_korder <- matrix(0, nrow = length(X)^K, ncol = length(X)^K)
rownames(trans_prob_korder) <- s_K
colnames(trans_prob_korder) <- s_K
for (k in 1:K) {
  s_k <- substr(rownames(trans_prob_korder), K-k+1, K)
  flag <- as.integer(s_k %in% branches)
  tt1 <- substr(rownames(trans_prob_korder)[flag==1], 2, K)
  tt2 <-  substr(colnames(trans_prob_korder), 1, K-1)
  for (i in 1:sum(flag)){
    trans_prob_korder[which(flag==1)[i], which(tt2==tt1[i])] <- trans_prob[s_k[flag==1][i],]  
  }
}

mu <- eigen(t(trans_prob_korder), symmetric = FALSE)
mu <- mu$vectors[,1]
mu <- Re(mu)
mu <- mu/sum(mu)
names(mu) <- rownames(trans_prob_korder)

P_A <- sum(mu[which(substr(names(mu), 1, 4)=="1033")])
P_XA <- mu[which(substr(names(mu), 1, 4)=="1033")]
P_XA/P_A # this equals to P(X|A)=(0.1856821 0.2473574 0.2426894 0.3242711), not approx (0.273 0.182 0.455 0.091)

# plot results
library(data.tree)

out <- readRDS('./vlmc_simudata_K=5_R=10_N=5*1e5_VLMC-C.RDS')

tt1 <- ToDataFrameTable(out$tau[[1]], "pathString")
tt1 <- strsplit(tt1, "[/]")
vcc.branches <- tt1
for(i in 1:length(tt1)){
  vcc.branches[[i]] <- tt1[[i]][-1]
}

vcc.branches.new <- vcc.branches
for(j in 1:max(lengths(vcc.branches))){
  tt1 <- sapply(vcc.branches, function(x) x[j]) 
  # tt1 <- tt1[!is.na(tt1)] # nodes from j-th levels
  
  if(j>1){
    tt1_parents <- sapply(vcc.branches, function(x) list(x[1:(j-1)]))
    tt1_parents <- sapply(tt1_parents, function(x) paste(x[!is.na(x)], collapse = ""))
    
    tt1_table <- table(tt1)
    tt1_rep_X <- names(tt1_table)[tt1_table>1]
    if(length(tt1_rep_X)==0)
      break()
    
    for(i in 1:length(tt1_rep_X)){
      tt2 <- (tt1==tt1_rep_X[i])
      tt2[is.na(tt2)] <- FALSE
      tt1_index <- seq(length(unique(tt1_parents[tt2])))
      
      for(l in tt1_index){
        tt3 <- tt1_parents==unique(tt1_parents[tt2])[l] & tt1==tt1_rep_X[i]
        tt1[tt3] <- paste(tt1[tt3], "-", l, sep="")
      }
      # tt3 <- which(as.character(tt1)==tt1_rep_X[i])
      # for(l in 1:length(tt3)){
      #   
      #   tt1[tt3[l]] <- paste(tt1[tt3[l]], "2", l, sep="")
      # }
      
    }
  }
  for(i in 1:length(vcc.branches)){
    if(!is.na(tt1[i]))
      vcc.branches.new[[i]][j] <- tt1[i]
  }
  
}

for(i in 1:length(vcc.branches.new)){
  for(j in 1:length(vcc.branches.new[[i]])){
    vcc.branches.new[[i]][j] <- paste("level", j, "-", vcc.branches.new[[i]][j], sep="")  
    
  }
}

library(dplyr)
library(igraph)
# set edge for igraph 
el <- data.frame(from=0, to=0)
cc <- 1
for(i in 1:length(vcc.branches.new)){
  
  for(j in 1:length(vcc.branches.new[[i]])){
    if(j==1){
      el[cc, "from"] <- NA
      el[cc, "to"] <- vcc.branches.new[[i]][j]
      cc <- cc+1
    }else{
      el[cc, "from"] <- vcc.branches.new[[i]][j-1]
      el[cc, "to"] <- vcc.branches.new[[i]][j]
      cc <- cc+1
    }
  }
}
el <- distinct(el)
el[is.na(el)] <- "ax1"
el <- el[order(el[,1]), ]
vcc.label <- c()
# extract labels for each node
for(i in 1:max(lengths(vcc.branches.new))){
  tt1 <- unlist(sapply(vcc.branches.new, function(x) if(!is.na(x[i])) x[i]))
  vcc.label <- c(vcc.label, unique(tt1))
}
tt1 <- vcc.label
tt1 <- c("ax1", tt1)
vcc.label <- strsplit(vcc.label, "[-]")
vcc.label <- sapply(vcc.label, function(x) x[[2]])
vcc.label <- c(NA, vcc.label)
vcc.label[!is.na(vcc.label)] <- sapply(vcc.label[!is.na(vcc.label)], function(x)
  substring(x, first = 1, last=1) )
x <- c("A","C","G","T")
vcc.label <- x[as.integer(vcc.label)+1]

vcc.label <- c(vcc.label, "A", "G", "C", "G")
names(vcc.label) <- c(tt1, "level5-01033", "level5-21033", "level2-10", "level2-20")

#### modify missing context 
el.miss.context <- data.frame(from=c("level4-1033-1", "level4-1033-1", "level1-0", "level1-0"), to=c("level5-01033", "level5-21033", "level2-10", "level2-20")) 
el <- rbind(el, el.miss.context)
el.graph <- graph.edgelist(as.matrix(el))
tt1 <- attr(V(el.graph), "names")
tt2 <- tt1
tt2[which(tt1=="level2-30")] <- "level2-10"
tt2[which(tt1=="level2-10")] <- "level2-20"
tt2[which(tt1=="level2-20")] <- "level2-30"
tt2[which(tt1=="level5-11033")] <- "level5-01033"
tt2[which(tt1=="level5-31033")] <- "level5-11033"
tt2[which(tt1=="level5-01033")] <- "level5-21033"
tt2[which(tt1=="level5-21033")] <- "level5-31033"
V(el.graph)$label <- vcc.label[tt2]

elc.graph <- el.graph

pdf("vcc_context_tree_K=5_N=1e5.pdf", height=7, width=7)
plot(elc.graph, layout=layout_as_tree,
     vertex.label.color='black',
     vertex.color='grey90', vertex.shape="circle", vertex.label.dist=0, vertex.size=13,
     edge.width=0.4, edge.color="black",
     edge.arrow.size=1, edge.arrow.width=0.5,
     # main="vlmc-B"
     # margin=c(-0, -0.3, -0., -0.3)
)
dev.off()

