#RMA callibration of the data(Quantile Normalization)
#this code will normalize our data means calulating and adjusting the mean of each column
my_eset_norm <- oligo::rma(raw_data)

#PCA plot of the normalized data
exp_eset <- Biobase::exprs(my_eset_norm)
PCA <- prcomp(t(exp_my_dataset), scale=FALSE)
percentvar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentvar[2] / percentvar[1])
dataGG <- data.frame(PC1=PCA_raw$x[,1], PC2=PCA_raw$x[,2],
                     phenotype=pData(my_eset_norm)$Factor.Value..DISEASE.STATE.)
library(ggplot2)
ggplot(dataGG,aes(PC1,PC2))+
  geom_point(aes(shape=phenotype,colour=phenotype))+
  ggtitle("PCA Plot of the calibrated, summarized data")+
  xlab(paste0("PC1, VarExp: ", percentvar[1] , "%")) +
  ylab(paste0("PC2, VarExp: ", percentvar[2], "%")) +
  theme(plot.title=element_text(hjust=0.5)) +
  coord_fixed(ratio= sd_ratio)+
  scale_shape_manual(values=c(4,15))+
  scale_color_manual(values=c("red" , "purple"))


#boxplot for normalized data
oligo::boxplot(my_eset_norm, target="core",
               las=2,
               main="Boxplot of log2 intensities for the raw data")
