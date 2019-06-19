# Set the current working directory
getwd();setwd("/home/sergiroda/Documents/Simulations/PETase_Manuel/2HEMHET4/WT/out")
# install.packages("ggExtra") # For density plots added in the plot
library(ggpubr);library(ggcorrplot);library("ggExtra");library(e1071)

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

summary(Total);
sd(Total[,18])
skewness(Total[,9])

# One-sample Wilocoxon test

wilcox.test(Total[,11],mu=3.35,alternative = "less")
wilcox.test(Total[,12],mu=1.8,alternative = "less")

wilcox.test(Total_D238G[,10],Total_D238G[,9],mu=0.2,alternative = "greater")

# Min of a set of columns (distances to substrate)

Total$min <- do.call(pmin, Total[,13:16])

Total <- cbind(Total,(Total$min)/5)
length(Total[,length(Total)][which(Total[,length(Total)]<=1)])/length(Total[,length(Total)])

# Density plot using kernel trick for PELE metrics

ggplot()+geom_density(data=Total,aes(Total$V12),
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

p <- ggplot(Total,aes(Total[,11],Total[,12]))+geom_point(color="white",size=1,alpha=0.5)+
  stat_density_2d(aes(col = ..level..),h=c(0.01,0.01))+
  theme(panel.background = element_blank(), 
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA),
      legend.position="none",
      axis.text = element_text(size =11),
      axis.title = element_text(size =12),
      plot.title = element_text(hjust = 0.5,size = 15))+
  xlim(2.7,2.85)+ylim(1.5,1.65)+
  xlab("Ser-His distance (Å)")+ylab("Asp-His distance (Å)")+ggtitle("2D Density plot of Ser-His distance against Asp-His distance")

ggMarginal(p,type = "density",color="darkblue",fill="blue",alpha=0.6,bw=c(0.01,0.01))

# As the previous one but used mainly for different subsets in the data.

ggscatterhist(
  table, x = "V8", y = "V9", size = 0.1, alpha = 0.6)
