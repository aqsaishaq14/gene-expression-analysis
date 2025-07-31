setwd("C:\\Users\\Ar Rehman Laptops\\Downloads\\urba")
getwd()
raw_data_dir <- "C:\\Users\\Ar Rehman Laptops\\Downloads\\urba"
if(!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
}


#1st try
#importing the SDRF file
sdrf_location <- file.path(raw_data_dir,"E-GEOD-23361.sdrf.txt")
SDRF <- read.delim(sdrf_location)
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

#installing oligo to read cel files
options(timeout = 1200)  # set longer download timeout

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("oligo", "pd.huex.1.0.st.v2"))

nBiocManager::install("pd.huex.1.0.st.v2")

library(pd.huex.1.0.st.v2)
library(oligo)
#chatgpt given code bcz my raw data is not found
library(oligo)

# Set to your correct folder where CEL files are stored
setwd("C:\\Users\\Ar Rehman Laptops\\Downloads\\urba")  # change this path!

# Get list of CEL files
cel_files <- list.celfiles(full.names = TRUE)

# Load into raw_data again
raw_data <- read.celfiles(cel_files)

# Now test:
boxplot(raw_data, target = "core")
#####

#2nd Try 
raw_data <- oligo::read.celfiles(filenames=file.path(raw_data_dir,
                                                     SDRF$Array.Data.File),
                                 verbose=FALSE, phenoData=SDRF)
  stopifnot(validObject(raw_data))
  #we made expression set which we called raw data
  #checking col names
  
colnames(exp_raw)
#pData function
head(Biobase::pData(raw_data))
names(Biobase::pData(raw_data))
Biobase::pData(raw_data) <- Biobase::pData(raw_data) [,
                                                      c("Source.Name"  , "Array.Data.File" ,"Factor.Value..DISEASE.STATE." )]

#Quality control of the raw data
Biobase::exprs(raw_data)[1:5, 1:5]
exp_raw <- log2(Biobase::exprs(raw_data))

#boxplot of the raw expression data
oligo::boxplot(raw_data, target="core" ,
               main="Boxplot of log2-intensities for the raw data")

 #Constructing PCA plot of raw expression data
PCA_raw <- prcomp(t(exp_raw), scale= FALSE)
percentvar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentvar[2] / percentvar[1])
dataGG <- data.frame(PC1=PCA_raw$x[,1], PC2=PCA_raw$x[,2],
                     phenotype=pData(raw_data)$Factor.Value..DISEASE.STATE.)
library(ggplot2)
ggplot(dataGG,aes(PC1,PC2))+
  geom_point(aes(shape=phenotype,colour=phenotype))+
  ggtitle("PCA Plot of the log-transformed raw expression data")+
  xlab(paste0("PC1, VarExp: ", percentvar[1] , "%")) +
  ylab(paste0("PC2, VarExp: ", percentvar[2], "%")) +
  theme(plot.title=element_text(hjust=0.5)) +
  coord_fixed(ratio= sd_ratio)+
  scale_shape_manual(values=c(4,15))+
  scale_color_manual(values=c("red" , "purple"))

#chtgpt given bcz it is showing my factor value disease is empty..
nrow(pData(raw_data))
rownames(PCA_raw$x)
rownames(pData(raw_data))
#checking columns
colnames(pData(raw_data))
pData(raw_data)$Factor.Value.disease
grep("disease", colnames(pData(raw_data)), value = TRUE)

