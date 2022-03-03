rm(list=ls())
# library(VLMC)
library(data.tree)
source("pruning_fun.R")
# S_all <- readRDS('./S_all_simudata.RDS')
S_all <- readRDS("./simulate_data.RDS")
X <- as.character(0:3)
THRESHOLD.GEN <- max(2, round(sqrt(nchar(S_all))/length(X)))


# Generate results of VLMC-C --------------------------------------

# Construct initial tree tau_0
# in each iteration, we can construt a deeper level for the current level, consisting children of 
# nodes in the current level. WIDTH-FIRST!
k <- 1
tau0 <- Node$new("x1")
tt_record_subtree_old <- list()
while (1){
  
  tt_seg <- substring(S_all, 1:(nchar(S_all)-k+1), k:nchar(S_all))
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
                                      obs_seq = S_all, X=X, alpha=ALPHA)
    
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

saveRDS(tau, "./tau_VLMC-C.RDS")
