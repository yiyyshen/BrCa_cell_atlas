##Code for Signature scores and highest Calls
#Aatish Thennavan Perou Lab

library(Seurat)

# read in scsubtype gene signatures
sigdat <- read.csv("NatGen_Supplementary_table_S4.csv") #define sigdat as the gene list table
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]), #temp_allgenes is a vector that contains 4 vectors of subtyping genes
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""]) #remove the duplicated genes 

#Read in the single cell RDS object as 'Mydata'
# Mydata <- readRDS("path_to_seurat_object.Rdata")
Mydata <- ScaleData(Mydata, features=temp_allgenes) #mydata is scaled and centered based on the genes from the signature gene list???
tocalc<-as.data.frame(Mydata@assays$RNA@scale.data) #extract scaled data???

#calculate mean scsubtype scores
outdat <- matrix(0,
                 nrow=ncol(sigdat), # row is subtype
                 ncol=ncol(tocalc), # column is each individual cell??
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))
for(i in 1:ncol(sigdat)){ # for each subtype
  
  # sigdat[i,!is.na(sigdat[i,])]->module #all rows in a column (all signature genes in a subtype)
  # row <- as.character(unlist(module))
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)}) #Calculate the mean of signature gene counts
  
  outdat[i,]<-as.numeric(temp)

}

final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

##Obtaining the highest call
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames

##Writing out output files (rownames remain the same for both)
write.table(finalm.sweep.t, "Mydata_Scores.txt", sep="\t")
write.table(Finalnames, "Mydata_CALLS.txt", sep="\t")
