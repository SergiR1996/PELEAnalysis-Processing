# Adding the e1071 library to calculate the skewness of the sets of the data.
#install.packages("e1071")
#install.packages("Peptides", dependencies=TRUE)
#install.packages("ggExtra")
#install.packages("ggridges")
library(e1071);library(Peptides);library(ggpubr);library("ggExtra");library(ggridges)

# The datasets of PELE results are loaded for the different mutants
Results <- read.csv(file="csvfile")

# Comparison of metrics between all PELE simulations using the Kruskal-Wallis test. Also, TukeyHSD test is shown.
?bartlett.test
?kruskal.test
?aov
?TukeyHSD

# Test the correlation between all the sets of the data
#install.packages("ggpubr");install.packages("ggcorrplot")
library(ggcorrplot)
corr <- round(cor(metrics),3)
p.mat <- cor_pmat(metrics)
ggcorrplot(corr,   ggtheme = ggplot2::theme_gray,
           colors = c("#F8A31B", "white", "#8DD3C7"),lab=TRUE,p.mat = p.mat,title = "Correlation matrix of the PELE quantitative metrics")+theme(plot.title = element_text(hjust = 0.5))

# -----------------------------

# Tools to visualize the datasets individually
?summary
?quantile

# Analysis of the normality of the test
?shapiro.test
?qqnorm
?qqline
?hist
?skewness # To check the asimmetry of the real random distribution

# Statistic tests of two samples assuming a normal distribution.
?t.test
?var.test
# Statistic tests of two samples not assuming a normal distribution.
?wilcox.test

# Use this if the PELE report files want to be analyzed step by step
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

# -----------------------------
# MD_simulations_results

# WT data frame
setwd("Working directory for MD simulation")

Results_1 <- read.table(file="Results directory")
Results_2 <- read.table(file="Results directory")
Results <- Results_1
Results_2 <- cbind(Results,Results_2[,2])

WT_MD_results_RMSD <- WT_MD_results_7[,2]
WT_MD_results_RMSD <- cbind(WT_MD_results_RMSD,WT_MD_results_6[,2])
WT_MD_results_RMSD <- cbind(WT_MD_results_RMSD,WT_MD_results_5[,2])
RMSD_WT <- c()
for (i in 1:length(WT_MD_results_5[,2])) {
  RMSD_WT <- c(RMSD_WT,mean(WT_MD_results_RMSD[i,1:3]))
}
for (i in (length(WT_MD_results_5[,2])+1):length(WT_MD_results_6[,2])){
  RMSD_WT <- c(RMSD_WT,mean(WT_MD_results_RMSD[i,1:2]))
}
for (i in (length(WT_MD_results_6[,2])+1):length(WT_MD_results_7[,2])){
  RMSD_WT <- c(RMSD_WT,mean(WT_MD_results_RMSD[i,1]))
}
RMSD_WT <- as.data.frame(RMSD_WT)
rmsd <- cbind(WT_MD_results_7[,1],RMSD_WT)
write.table(rmsd, "rmsd_tot.xvg", sep="\t",append = FALSE)

names(WT_MD_results) <-c("names of columns")

?plot(density)

ggplot() +geom_density(data=Results,aes(ColumnName,fill="ColumnName"),
                      bw=0.1,alpha=0.5)+
  geom_density(data=Results2,aes(ColumNName2,fill="ColumnName2"),
               bw=0.1,alpha=0.5)+theme_minimal()+
  guides(fill=guide_legend(title="Metric"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))

figure <- ggarrange(plot1, plot2,
                    labels = c(1,2),
                    ncol = 2, nrow = 1)

annotate_figure(figure,
                top = text_grob("Title of figure", size = 15))

p <- ggplot(WT_MD_results,aes(WT_MD_results[,3],WT_MD_results[,4]))+geom_point(color="black",size=1,alpha=0.5)+
  geom_density_2d(color="blue")+theme_minimal()

ggMarginal(p,type = "density",color="blue",fill="blue",alpha=0.5)

# As the previous one but used mainly for different subsets in the data.

ggscatterhist(
  table, x = "V8", y = "V9", size = 0.1, alpha = 0.6)

# Find the minimum around a set of metrics in a data frame
ComplexMin <- do.call(pmin, metrics)

# Ridgeline plots

ggplot()+
  geom_density_ridges(data=DF1,aes(DF1[,column],y="DF1",fill="DF1",scale=0.1))+
  geom_density_ridges(data=DF2,aes(DF2[,column],y="DF2",fill="DF2",scale=0.1))+
  geom_density_ridges(data=DF3,aes(DF3[,column],y="DF3",fill="DF3",scale=0.1))+
  theme_minimal()+ggtitle("Title of plot")+xlab("X axis title")+ylab("Y axis title")+
  guides(fill=guide_legend(title="Y axis title"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"))

