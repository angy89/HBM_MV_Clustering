#Given mean and variance return the beta distribution a and b parameters
beta_params = function(mean, variance){
  if(variance > (mean * (1 - mean))){
    stop('Variance cannot be bigger than mean*(1-mean)')
  }
  alpha <- ((1 - mean) / variance - 1 / mean) * mean ^ 2
  beta <- alpha * (1 / mean - 1)
  return(params = list(alpha = alpha, beta = beta))
}

library(igraph)
source("src/sw_cutting.R")
gen_th = 0.3

n_islands = 10
n_nodes = 90
island_size = c(10,5,6,18,4,12,6,10,6,13)
n_weak_edges = 2500
n_edges_between_islands = 3000

cat("Generating network with: \n")
cat("Number of Islands: ",n_islands,"\n")
cat("Islands Size: ",island_size,"\n")
cat("Edges Between Islands: ",n_edges_between_islands,"\n")
cat("Number of Weak Edges: ",n_weak_edges,"\n")

sd_strong_edge = 0.1
sd_weak_edge = 0.1
sd_noise_edge = 0.1

M_se = c()
M_sp = c()
M_acc = c()

for(mean_strong_edge in c(0.7,0.8,0.9)){
  for(mean_weak_edge in c(0.5,0.6,0.7)){
    for(mean_noise_edge in c(0.2,0.3,0.4)){
      cat("MSE: ",mean_strong_edge," MWE: ",mean_weak_edge," MNE: ",mean_noise_edge, "\n")
      
      strong_edge_params = beta_params(mean_strong_edge, sd_strong_edge^2)
      weak_edge_params = beta_params(mean_weak_edge, sd_weak_edge^2)
      noise_edge_params = beta_params(mean_noise_edge, sd_noise_edge^2)
      
      #Each islands of the network is generated as a preferential attachment model
      
      ADJ = matrix(0,nrow = n_nodes,ncol = n_nodes)
      start = 1
      end = island_size[1]
      islands = list()
      
      for(i in 1:n_islands){
        G_pa= sample_pa(island_size[i],directed = TRUE)
        E(G_pa)$weight = rbeta(length(E(G_pa)), strong_edge_params$alpha, 
                               strong_edge_params$beta)
        sub_adj = get.adjacency(G_pa,attr = "weight",type = "both",sparse = FALSE)
        islands[[i]]= start:end
        ADJ[start:end,start:end] = sub_adj
        cat(start, " ",end ,"\n")
        start = end+1
        end = start + island_size[i+1] -1
      }
      
      GG = graph.adjacency(ADJ,mode="directed",weighted = TRUE)
      #plot(GG)
      
      #Some random edges are attacched between the islands;
      # The edges inside the islands and the edges between the islands are added as strong edges
      
      for(i in 1:n_edges_between_islands){
        idx_island_from = sample(x = 1:length(islands),size = 1)
        from = islands[[idx_island_from]][sample(1:length(islands[[idx_island_from]]),1)]
        sample_to=1:length(islands)
        sample_to = sample_to[-idx_island_from]
        idx_island_to = sample(x = sample_to,size = 1)
        to = islands[[idx_island_to]][sample(1:length(islands[[idx_island_to]]),1)]
        ADJ[from, to] = rbeta(1, strong_edge_params$alpha, 
                              strong_edge_params$beta)
      }
      
      #The edges between all the other couple of nodes can be weak edges or noise. A certan number of weak edges are added
      #to the graph respecting the properties that an edge between two nodes i and j is weak if there are less than gen_th stron
      #edges incident on i and j
      
      G_pa = graph.adjacency(ADJ,mode="directed",weighted = TRUE)
      deg = degree(G_pa)
      deg = deg - round(deg * gen_th)
      CandidatesWeakEdges = which(ADJ==0,arr.ind=TRUE)
      pesi = 1/deg[CandidatesWeakEdges[,1]] * 1/deg[CandidatesWeakEdges[,2]]
      CandidatesWeakEdges = cbind(CandidatesWeakEdges,pesi)  
      weak_idx = sample(seq(1, to=nrow(CandidatesWeakEdges)), size = n_weak_edges, replace = F, prob = CandidatesWeakEdges[,3])
      ADJ[CandidatesWeakEdges[weak_idx,1:2]] = rbeta(length(weak_idx), weak_edge_params$alpha, 
                                                     weak_edge_params$beta)
      
      GG = graph.adjacency(ADJ,mode="directed",weighted = TRUE)
      # plot(GG)
      V(GG)$name = paste("V",seq(vcount(GG)),sep="")
      
      #The rest of the edges in the network, in order to have a completelly connected network is noise
      
      noise = matrix(rbeta(vcount(GG)*vcount(GG), noise_edge_params$alpha, 
                           noise_edge_params$beta),
                     nrow=vcount(GG), ncol=vcount(GG))
      true_net = ADJ
      ADJ[ADJ==0] = noise[ADJ==0]

      ths = seq(from=0.1,to=1,by=0.1)
      
    XRes =   mclapply(X = ths,FUN = function(th){
        res = sw_cutting(graph = GG, invert = TRUE, threshold = th,flow = 0)
        res_adj = get.adjacency(res$residualGraph,type = "both",sparse = FALSE,attr = "weight")
        
        true_net[true_net>0] = 1
        res_adj[res_adj>0] = 1
        
        s_pred = res_adj
        s_true = true_net
        
        s_tp = sum(s_pred * s_true)
        s_tn = sum((1 - s_pred) * (1 - s_true))
        s_fp = sum(s_pred * (1 - s_true))
        s_fn = sum(s_true * (1 - s_pred))
        
        se = s_tp / (s_tp + s_fn)
        sp = s_tn / (s_tn + s_fp)
        acc = (s_tp  + s_tn)/ (s_tp + s_fp + s_fn + s_tn)
        
        x = c(se,sp,acc)
        return(x)
        
      },mc.cores = 2)
      
      Mat = do.call(rbind,XRes)
      colnames(Mat) = c("Sensitivity","Specificity","Accuracy")
      # sp = c()
      # se = c()
      # acc = c()
      #for(i in 1:length(ths)){
        # res = sw_cutting(graph = GG, invert = TRUE, threshold = ths[i],flow = 0)
        # res_adj = get.adjacency(res$residualGraph,type = "both",sparse = FALSE,attr = "weight")
        # 
        # true_net[true_net>0] = 1
        # res_adj[res_adj>0] = 1
        # 
        # s_pred = res_adj
        # s_true = true_net
        # 
        # s_tp = sum(s_pred * s_true)
        # s_tn = sum((1 - s_pred) * (1 - s_true))
        # s_fp = sum(s_pred * (1 - s_true))
        # s_fn = sum(s_true * (1 - s_pred))
        # 
        # se = c(se, s_tp / (s_tp + s_fn))
        # sp = c(sp, s_tn / (s_tn + s_fp))
        # acc = c(acc, (s_tp  + s_tn)/ (s_tp + s_fp + s_fn + s_tn))
      #}end for
      
      #
      #cat("MSE: ",mean_strong_edge," MWE: ",mean_weak_edge," MNE: ",mean_noise_edge, "\n")
      
      M_se = rbind(M_se,c(mean_strong_edge,mean_weak_edge,mean_noise_edge,Mat[,1]))
      M_sp = rbind(M_sp,c(mean_strong_edge,mean_weak_edge,mean_noise_edge,Mat[,2]))
      M_acc = rbind(M_acc,c(mean_strong_edge,mean_weak_edge,mean_noise_edge,Mat[,3]))
      
      png(paste("sim_res2/MSE_",mean_strong_edge * 100,"_SSE_",sd_strong_edge * 100,
                "_MWE_",mean_weak_edge * 100,"_SWE_",sd_weak_edge * 100,
                "_MNE_",mean_noise_edge * 100,"_SNE_",sd_noise_edge * 100,
                "_GENTH_",gen_th * 100,
                ".png",sep=""))
      plot(se, col='red', type="l")
      lines(sp, col='blue')
      lines(acc, col='green')
      
      legend(x = "bottomright",legend = c("Specificity","Sensitivity","Accuracy"),
             fill=c("red","blue","green"))
      dev.off()
      save(se,sp,acc, file = paste("sim_res2/MSE_",mean_strong_edge * 100,"_SSE_",sd_strong_edge * 100,
                                   "_MWE_",mean_weak_edge * 100,"_SWE_",sd_weak_edge * 100,
                                   "_MNE_",mean_noise_edge * 100,"_SNE_",sd_noise_edge * 100,
                                   "_GENTH_",gen_th * 100,
                                   ".RData",sep=""))
    }
  }
}

# colnames(M_se) = colnames(M_sp) = c("Mean_Strong","Mean_Weak","Mean_Noise",
#                                     paste("Th:",ths * 100,sep=""))
# 
# png(paste("sensitivity_th",gen_th * 100,n_nodes,"_nodes.png",sep=""))
# rownames(M_se) = paste("MSE:",M_se[,1]*100," MWE:",M_se[,2]*100," MNE:",M_se[,3]*100,sep="")
# heatmap.2(M_se[,4:13],dendrogram="none",main="Sensitivity",
#           margins=c(5,10),cexRow=0.8,cexCol=0.8)
# dev.off()