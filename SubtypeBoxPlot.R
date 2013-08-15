SubtypeBoxPlot <- function(master.eset,Wilcox.test,modified.master.eset,output.file){
  
  #######   Indexation of everything   #########
  
  # Indexing everything for proper subtype classification
  index2 <- which(colnames(pData(master.eset))=="subtype")
  subtype.name <- colnames(pData(master.eset))[(index2+1):ncol(pData(master.eset))]
  subtype.comb <- combinations(n=length(subtype.name), r=2)
  
  #Changing clinical data subtype to a more suitable subtype name
  list.index<- list() 
  temp <- as.matrix(pData(master.eset))
  index3 <- which(colnames(temp)=="subtype")
  temp[which(temp[, index3]=="ER+/HER2- High Prolif"),index3] <- "LumB"
  temp[which(temp[, index3]=="ER+/HER2- Low Prolif"),index3] <- "LumA"
  temp[which(temp[, index3]=="ER-/HER2-"),index3] <- "Basal"
  temp[which(temp[, index3]=="HER2+"),index3] <- "Her2"
  
  #Generating List of index for classification
  for (i in 1:length(subtype.name)){    
    if (subtype.name[[i]]=="Global.population"){list.index[["Global"]] <- 1:nrow(pData(master.eset))}
    if (subtype.name[[i]]=="Lums"){list.index[[subtype.name[[i]]]] <- NA}  
    list.index[[subtype.name[[i]]]] <- which(temp[,index3]==subtype.name[i])
  }
  
  #Rearranging List of index
  if (!is.null(list.index$Lums)){list.index$Lums <- c(list.index$LumA , list.index$LumB)} 
  list.index <- list(list.index$Basal,list.index$Her2,list.index$LumB,list.index$LumA,list.index$Lums, list.index$Global)
  names(list.index) <-c("Basal", "Her2","LumB","LumA","Lums","Global")
  
  pdf(output.file)
  for (i in 1:nrow(Wilcox.test)){
    Basal <- exprs(modified.master.eset)[Wilcox.test[i,1],list.index[[1]]]
    Her2 <- exprs(modified.master.eset)[Wilcox.test[i,1],list.index[[2]]]
    LumB <- exprs(modified.master.eset)[Wilcox.test[i,1],list.index[[3]]]
    LumA <- exprs(modified.master.eset)[Wilcox.test[i,1],list.index[[4]]]
    Lums <- exprs(modified.master.eset)[Wilcox.test[i,1],list.index[[5]]]
    Global <- exprs(modified.master.eset)[Wilcox.test[i,1],list.index[[6]]]
  
    boxplot(Basal, Her2,LumB,LumA,Lums,Global, names=c("Basal", "Her2","LumB","LumA","Lums","Global"), xlab="Subtype", ylab="Metagene expression",main=sprintf("Subtype dependance of metagene %s",Wilcox.test[i,1]))
        
  }
  dev.off()
}