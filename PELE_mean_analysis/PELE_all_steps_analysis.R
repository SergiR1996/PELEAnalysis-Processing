# Set the current working directory
getwd();setwd("/home/sergi/Documents/MSc_in_Bioinformatics_UAB/Module/Module_6&7/Results/PELE/2HEMHET/WT")
# install.packages("ggExtra") # For density plots added in the plot
library(ggpubr);library(ggcorrplot);library("ggExtra")

# Generate the table of the PELE results. Transform the PELE results with the sed tool
# (sed -i 's/    / /g' *.out) to enable its recognition.
Total <- read.table("report_1.out",header=FALSE,sep=" ")
Total[,length(Total)] <- 1
Total <- Total[5:nrow(Total),2:length(Total)]
for (i in 2:128) {
  table <- read.table(sprintf("report_%s.out",i),header = FALSE,sep=" ")
  table[,length(table)] <- i
  Total <- rbind(Total,table[5:nrow(table),2:length(table)]) # Append the table of one report with 
  # the 4 first elements of the PELE simulation eliminated.
}

# Density plot using kernel trick for PELE metrics

ggplot()+geom_density(data=Total,aes(Total[,8]),
                      bw=0.1,alpha=0.5,fill="#56B2E9")+
  theme_minimal()+ggtitle("Title")+xlab("X label")+
  theme(plot.title = element_text(hjust = 0.5))

# Correlation matrix plotted between different PELE metrics

corr <- round(cor(Total[,5:10]),3)
p.mat <- cor_pmat(Total[,5:10])
ggcorrplot(corr,   ggtheme = ggplot2::theme_gray,
           colors = c("#F8A31B", "white", "#8DD3C7"),lab=TRUE,p.mat = p.mat,title = "Correlation matrix of the PELE quantitative metrics")+
  theme(plot.title = element_text(hjust = 0.5))


# Scatter plot of two metrics with density plots represented in the sides

p <- ggplot(Total,aes(Total[,5],Total[,6]))+geom_point(color="white",size=1,alpha=0.5)+
  stat_density_2d(aes(col = ..level..))+
  theme(panel.background = element_blank(), 
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA),
      legend.position="none")

ggMarginal(p,type = "density",color="darkblue",fill="blue",alpha=0.6)

# As the previous one but used mainly for different subsets in the data.

ggscatterhist(
  table, x = "V8", y = "V9", size = 0.1, alpha = 0.6)
