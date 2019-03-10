#crea la lista dei dataset
createListDataset<-function(listPath, patient_classes,controllo){
  listDataSet<-c()
  for(i in 1:length(listPath)){
    fileRead<-file.path(listPath[i])
    listDataSet[[i]]<- read.table(fileRead)
    # create the patient matrices
    if(controllo == "Oxford"){
      colnames(listDataSet[[i]])<-row.names(patient_classes)
    }
    listDataSet[[i]] <- listDataSet[[i]][row.names(patient_classes)]
  }
  return(listDataSet)
}
