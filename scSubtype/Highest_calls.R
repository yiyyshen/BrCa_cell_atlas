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
Mydata <- ScaleData(Mydata, features=temp_allgenes) # The feature genes in Mydata are scaled and centered with mean of 0
tocalc<-as.data.frame(Mydata@assays$RNA@scale.data) # extract scaled data as dataframe named "tocalc"

#calculate mean scsubtype scores
outdat <- matrix(0,
                 nrow=ncol(sigdat), # row is subtype
                 ncol=ncol(tocalc), # column is each individual cell? so row is gene
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))
for(i in 1:ncol(sigdat)){ # for each subtype
  
  # sigdat[i,!is.na(sigdat[i,])]->module #all rows in a sigdat column (all signature genes in a subtype)
  # row <- as.character(unlist(module))
  row <- as.character(sigdat[,i]) #create a variable called row, which contains all genes of a subtype as characters, like row = "gene a" "gene b" "gene c"...
  row<-unique(row[row != ""]) #make sure no empty row and each gene only appears once in each subtype???
  genes<-which(rownames(tocalc) %in% row) #create a new variable called "genes", which includes only the POSITION of the genes that appear in the signature gene list (like 1, 2, 4, ...) 
  
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)}) #Calculate the mean (scaled) gene expression of all signature genes for each cell
  
  outdat[i,]<-as.numeric(temp) #add the mean expression as the subtype score for each cell

}

final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),] #create a variable called "final" and remove the rows with sum=0 from outdat (remove the subtype with 0 total gene expression of all cells)
final<-as.data.frame(final) #make final into a dataframe
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4) # Round final up? what is the 4 doing here??
finalm<-as.matrix(final) #make final into a matrix

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) { #define a function with 2 inputs: data x, and row.w is the default input, and repeat 1 nrow times and divide each by nrow (if there are 4 rows, generate  1 1 1 1, divide each by 4, generate 1/4, 1/4, 1/4, 1/4)
  get_average <- function(v) sum(v * row.w)/sum(row.w) # v is each column in x, to calculate the avergae in each column
  average <- apply(x, 2, get_average) # calculate the average gene expression for each cell (average of all the rows within the specified column) 
  sweep(x, 2, average) # substracting the mean of the column in every row of the column
} #normalizing the subtype score in each row

##Obtaining the highest call 
finalmt<-as.data.frame(t(finalm)) #pick highest number in each column
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames

##Writing out output files (rownames remain the same for both)
write.table(finalm.sweep.t, "Mydata_Scores.txt", sep="\t")
write.table(Finalnames, "Mydata_CALLS.txt", sep="\t")
