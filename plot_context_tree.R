rm(llist=ls())
library(dplyr)
library(igraph)

tau <- readRDS('./tau_VLMC-C.RDS')


tt1 <- ToDataFrameTable(tau, "pathString")
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
names(vcc.label) <- tt1

el.graph <- graph.edgelist(as.matrix(el))
tt1 <- attr(V(el.graph), "names")
V(el.graph)$label <- vcc.label[tt1]


pdf("tau_VLMC-C_plot.pdf", height=7, width=7)
plot(el.graph, layout=layout_as_tree,
     vertex.label.color='black',
     vertex.color='grey90', vertex.shape="circle", vertex.label.dist=0, vertex.size=13,
     edge.width=0.4, edge.color="black",
     edge.arrow.size=1, edge.arrow.width=0.5,
     # main="vlmc-B"
     # margin=c(-0, -0.3, -0., -0.3)
)
dev.off()
