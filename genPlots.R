##########Do not change, scroll below##########################################################
###############################################################################################

#setwd("C:/Users/CIBA-GBU/Desktop/mrfR/new2")
#library('plyr')
library(dplyr)
library(RColorBrewer)
library(gplots)



cds <- read.table("cdsHeatmap.txt", sep="\t", header=T)
cols <- ncol(cds)
rnames <- cds[,1]

cols

cdsmat <- data.matrix(cds[,2:ncol(cds)])
rownames(cdsmat) <- rnames

#cdsmat

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
 numGenomes = 50
} else if (length(args)==1) {
 numGenomes = as.numeric(args[1])
}

numGenomes

hmcol = colorRampPalette(brewer.pal(9, "YlOrRd"))(25)

#png("heatmap3.png", width = 860, height = 550)
png("heatmap_noClustering.png", width = 960, height = 850)
#heatmap.2(cdsdesc25mat, col = hmcol, trace="none", margin=c(30,18))

heatmap.2(cdsmat, col = hmcol, trace="none", margin=c(18,15), cexRow=2.0, cexCol=2.0, Rowv=FALSE,Colv=FALSE)
dev.off()

png("heatmap_clusterOn.png", width = 960, height = 1080)

heatmap.2(cdsmat, col = hmcol, trace="none", margin=c(30,18), cexRow=1.5, cexCol=2.0)

#heatmap.2(cdsmat, col = hmcol, trace="none", margin=c(30,18))
#heatmap.2(cdsmat, col = hmcol, trace="none", margin=c(18,15), cexRow=1.5, cexCol=1.5, Rowv=FALSE,Colv=FALSE)
dev.off()


####All fine above do not change
###working on summary below

summary <- as.data.frame(read.table("summaryAll.txt", sep="\t", header=T))

#detach("package:dplyr")
#library('plyr')
library(ggplot2)
#lostCount <- count(summary, 'Total.missing.CDS.length')
lostCount <- summary %>% dplyr::count(Total_missing_CDS_length)


##works
#ggplot(lostCount, aes(x=freq, y=Total.missing.CDS.length)) + 
#  geom_point(colour="red", aes(size=Total.missing.CDS.length)) ##works but we dont need now
  


##we are going to use plot2 as of now
plot2 <- ggplot(lostCount, aes(x = n, y = Total_missing_CDS_length, colour=Total_missing_CDS_length)) + geom_point(size = 5) + 
    scale_colour_gradientn(colours=rainbow(3)) +
    xlab("No. of genomes") + ylab("Tot. missing CDS length") + theme_classic(base_size = 15)
plot2 

 png("cdsLostvsFreq.png", width = 960, height = 850)
 print(plot2)
 dev.off()

###plot3 also works, just a color variation of plot2, so not using for now
#plot3 <- ggplot(lostCount, aes(x = freq, y = Total.missing.CDS.length, colour=Total.missing.CDS.length)) + geom_point(size = 5) + 
#    scale_colour_gradient(low = "orange", high = "red")
#plot3
  
  
############stacked bar plot working

#summarydf <- as.data.frame(read.table("summaryAll.txt", sep="\t", header=T))

#library('dplyr')
sumdfsmall <- summary %>% select(1, 5:6)

nrow(sumdfsmall)

#detach("package:dplyr")
library('tidyr')
if(nrow(sumdfsmall) <  numGenomes){
#sumdfsmallHd <- head(sumdfsmall, 50)




testtab <- pivot_longer(sumdfsmall, cols=2:3, names_to = "type", values_to = "length")


bar <- ggplot(testtab, aes(x = accession, y = length))+
  geom_col(aes(fill = type), width=1) +  theme_classic(base_size = 15) +
 theme(axis.text.x = element_text(angle = 90, face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14),
legend.title=element_text(size=15), 
    legend.text=element_text(size=14)) 
png("StackedBar.png", width = 1080, height = 850)
 print(bar)
 dev.off() 
} else {

sumdfsmallHd <- head(sumdfsmall, numGenomes)

testtab <- pivot_longer(sumdfsmallHd, cols=2:3, names_to = "type", values_to = "length")


bar <- ggplot(testtab, aes(x = accession, y = length)) +
  geom_col(aes(fill = type), width=1) + theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90,size = 14, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), 
legend.title=element_text(size=15), 
    legend.text=element_text(size=14))
 png("StackedBar.png", width = 1080, height = 850)
 print(bar)
 dev.off()
}

#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



