# library("devtools")
# remove.packages(pkgs = "MVDA")
# remove.packages(pkgs="CorKohonen") 
# install.packages("~/Dropbox/University/Tesi/CorKohonen_1.0.tar.gz", repos = NULL, type = "source")
# install.packages(c('clv', 'amap', 'pvclust', 'fpc', 'infotheo', 'kernlab', 'randomForest', 'ggplot2', 'reshape2'),dependencies = T)
# install.packages("sda")
# install.packages("rminer")
# install.packages("C:/Users/mguida/Dropbox/University/Tesi/MVDA_0.1.tar.gz", repos = NULL, type = "source")

library(clv)
library(amap)
library(pvclust)
library(fpc)
library(infotheo)
library(kernlab)
library(randomForest)
library(ggplot2)
library(reshape2)
library(MVDA)
library(caret)
library(matrixcalc)
library(e1071)
library(rpart)


require(igraph)
source('chromoHBM.r')
source('chromoHBM3C.r')
source('tesi_createListDataset.R')
source('MVDA_clustering.r')
source('tesi_cross_validation.R')
source('tesi_function.R')
source('tesi_cluster_gerarchico.R')
source('tesi_createHBM.R')
source('tesi_createMatrixSimilarity.R')
source('tesi_multi_cross_validation.R')
source('tesi_createHBM3C.R')



library(caret)

lista_nomi_dataset_project<-c("BRCA", "GBM", "KIRK", "OV", "Oxford", "Oxford");

########## creo i nomi dei data set da leggere
lista_nomi_dati_b<-c("variableRNASeq.txt", "variableMirRNASeq.txt");
lista_nomi_dati_g<-c("variable_genes.txt", "variable_miRNA.txt");
lista_nomi_dati_k<-c("variable_genes.txt", "variable_miRNA.txt", "variable_CNV.txt");
lista_nomi_dati_o<-c("variable_genes.txt", "variable_miRNA.txt","variable_protein.txt");
lista_nomi_dati_ox<-c("variable_genes.txt", "variable_miRNA.txt");
lista_nomi_dati_oxp<-c("variable_genes.txt", "variable_miRNA.txt");

link_data_source<-"D:/Users/mguida/Desktop/Progetto_tesi_Maria_domenica_guida"     #path in cui sono inseriti i dati
link_data_result<-paste(link_data_source,"/toSend",sep = "")    #cartella in cui sono contenuti i dataset da testare


lista_nomi_dataset<-list("BRCA"= lista_nomi_dati_b,"GBM"= lista_nomi_dati_g,"KIRK"= lista_nomi_dati_k,"OV"= lista_nomi_dati_o,"Oxford"= lista_nomi_dati_ox, "OxfordPam"= lista_nomi_dati_oxp)

#cartella del data set testato in cui sono contenuti i risultati
listaPath<-c()
link_project<-c()
for(i in 1:length(lista_nomi_dataset)){
  link_project[i]<-paste(link_data_result,lista_nomi_dataset_project[i],sep = "/")   
}
path_data<-c()          #contiene le path di tutti i dataset
for(i in 1:length(lista_nomi_dataset)){
  link_view<-c()
  for(j in 1:length(lista_nomi_dataset[[i]])){
    link_view[j]<-paste(link_project[i],lista_nomi_dataset[[i]][j],sep = "/")
  }
  path_data[[i]]<-link_view
}

path_patient_classes<-c()
path_clust_label <- c()
path_error_label <- c()
for(i in 1:length(lista_nomi_dataset)){
  path_patient_classes[i] <- file.path(paste(link_project[i],"/patient_classes.txt",sep = ""))  
  if(i==5){
    path_clust_label[i] <- file.path(paste(link_project[i],"/matrix_factorization_grave/all_feature/unsupervised/cluster_class.txt",sep = ""))
    path_error_label[i] <- file.path(paste(link_project[i],"/matrix_factorization_grave/all_feature/unsupervised/error.txt",sep = ""))
  }
  if(i == 6){
    path_clust_label[i] <- file.path(paste(link_project[i],"/matrix_factorization_pam50/all_feature/unsupervised/cluster_class.txt",sep = ""))
    path_error_label[i] <- file.path(paste(link_project[i],"/matrix_factorization_pam50/all_feature/unsupervised/error.txt",sep = ""))
  }
  if(i != 5 && i != 6){
    path_clust_label[i] <- file.path(paste(link_project[i],"/matrix_factorization/all_feature/unsupervised/cluster_class.txt",sep = ""))
    path_error_label[i]  <- file.path(paste(link_project[i],"/matrix_factorization/all_feature/unsupervised/error.txt",sep = ""))
    }
  
}
       
####### Ripeto tutto il progetto per ogni dataset
for(proge in 1:length(lista_nomi_dataset)){
 
  mvda_cluster_label<-read.table(path_clust_label[proge])
  mvda_cluster_error<-read.table(path_error_label[proge])
  mvda_cluster_error<-mvda_cluster_error[[1]]

  patient_classes<-read.table(path_patient_classes[proge])
  
  method_similarity <- "pearson"
  indexSim <-0.80       #indice minimo considerato per la matrice di similarit?
  cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
  
  listDataSet<-createListDataset(path_data[[proge]], patient_classes,lista_nomi_dataset_project[proge])
  resultMatrixSimilarity <- createMatrixSimilarity(listDataSet,method_similarity,indexSim)    #contiene sia la matrice di similarit? che quella tagliata
  listMatrixSimilarity<-resultMatrixSimilarity$listSimMatrix
  listMatrixSimilarityTag<-resultMatrixSimilarity$listMatrixT
  size_patients <- ncol(resultMatrixSimilarity$listSimMatrix[[1]])
  
  
  cluster.args <-c()
  
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
  
  
  
  ################## Clusterizzazione gerarchica ##################
  
  #eseguo il codice per calcolare il cluster con il metodo MVDA
  listDissimilarityMatrix<-c()
  for(i in 1:length(listMatrixSimilarity)){
    listDissimilarityMatrix[[i]] <- 1-listMatrixSimilarity[[i]]
  }
  listDissimilarityMatrix[[i+1]]<- result_chromo_hbm_3c
  
  if(proge == 6){
    resultClusterGer <- cluster_gerarchico(listDissimilarityMatrix, "ward.D2", length(mvda_cluster_label$x), colnames(listDissimilarityMatrix[[1]]), patient_classes$label.pazienti.pam50 )
  }
  else{
    resultClusterGer <- cluster_gerarchico(listDissimilarityMatrix, "ward.D2", length(mvda_cluster_label$x), colnames(listDissimilarityMatrix[[1]]), patient_classes$x )
    
  }
  
  
  #cluster_MI<- cluster_mutual_information(pred, mvda_cluster$cluster, mirna)
  
  if(proge == 6){
    lista_nomi_dataset_project[proge] = "Oxford Pam50"
  }
  if(length(lista_nomi_dataset[[proge]]) == 2){
    cat("Risultati clusterizzazione gerarchica \n \n ",lista_nomi_dataset_project[proge]," \n \n \t\t",lista_nomi_dataset[[proge]][1]," \t\t",lista_nomi_dataset[[proge]][2]," \t\t HMB \t\t\t\ MVDA \n \n Errore \t\t",resultClusterGer[[1]]$errore*100,"% \t\t" ,resultClusterGer[[2]]$errore*100,"% \t\t\t" ,resultClusterGer[[3]]$errore*100,"% \t\t", mvda_cluster_error*100,"%\n\n")
     }
  else {
    cat("Risultati clusterizzazione gerarchica \n \n ",lista_nomi_dataset_project[proge]," \n \n \t\t",lista_nomi_dataset[[proge]][1]," \t\t",lista_nomi_dataset[[proge]][2]," \t\t",lista_nomi_dataset[[proge]][3]," \t\t HMB \t\t\t\ MVDA \n \n Errore \t\t",resultClusterGer[[1]]$errore*100,"% \t\t" ,resultClusterGer[[2]]$errore*100,"% \t\t\t\t" ,resultClusterGer[[3]]$errore*100,"% \t\t",resultClusterGer[[4]]$errore*100,"% \t\t", mvda_cluster_error*100,"%\n\n")
    
     }
}




