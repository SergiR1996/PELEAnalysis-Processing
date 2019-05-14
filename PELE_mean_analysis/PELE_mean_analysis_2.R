# Adding the e1071 library to calculate the skewness of the sets of the data.
install.packages("e1071")
library(e1071)

# The datasets of PELE results are loaded for the different mutants
WT_results <- read.csv(file="Results2/Results_WT_2HEMHET4.csv")
D238SL272D_results <- read.csv(file="Results2/Results_T124D_2HEMHET4.csv")
D238SW269D_results <- read.csv(file="Results2/Results_Q137E_2HEMHET4.csv")
D238S278K_results <- read.csv(file="Results2/Results_DE_2HEMHET4.csv")

summary(WT_results);summary(D238SL272D_results);summary(D238S278K_results)

ServsHis_data <- data.frame(
  ServsHis=c(WT_results[,10],D238SL272D_results[,7],D238SW269D_results[,7],D238S278K_results[,7]),
  mutants =factor(rep(c("WT", "A", "B","C"), times=c(length(WT_results[,10]),length(D238SL272D_results[,7]),
                                                                         length(D238SW269D_results[,7]),length(D238S278K_results[,7])))))
ANOVA_ServsHis <- aov(ServsHis~mutants,data = ServsHis_data)
TukeyHSD(ANOVA_ServsHis);plot(TukeyHSD(ANOVA_ServsHis),las=2,cex=0.5,col="darkblue")+mtext("Serine-Histidine distance of main active site",line=0.5)

HisvsAsp_data <- data.frame(
  HisvsAsp=c(WT_results[,11],D238SL272D_results[,8],D238SW269D_results[,8],D238S278K_results[,8]),
  mutants =factor(rep(c("WT", "A", "B","C"), times=c(length(WT_results[,11]),length(D238SL272D_results[,8]),
                                                         length(D238SW269D_results[,8]),length(D238S278K_results[,8])))))
ANOVA_HisvsAsp <- aov(HisvsAsp~mutants,data = HisvsAsp_data)
TukeyHSD(ANOVA_HisvsAsp);plot(TukeyHSD(ANOVA_HisvsAsp),las=2,cex=0.5,col="darkblue")+mtext("Histidine-Aspartate distance of main active site",line=0.5)

Ser238vsHis_data <- data.frame(
  Ser238vsHis=c(WT_results[,10],D238SL272D_results[,13],D238SW269D_results[,13],D238S278K_results[,13]),
  mutants =factor(rep(c("WT", "A", "B","C"), times=c(length(WT_results[,10]),length(D238SL272D_results[,13]),
                                                         length(D238SW269D_results[,13]),length(D238S278K_results[,13])))))
ANOVA_Ser238vsHis <- aov(Ser238vsHis~mutants,data = Ser238vsHis_data)
TukeyHSD(ANOVA_Ser238vsHis);plot(TukeyHSD(ANOVA_Ser238vsHis),las=2,cex=0.5,col="darkblue")+mtext("Serine133-Histidine distance of potential active site",line=0.5)

His238vsAsp_data <- data.frame(
  His238vsAsp=c(WT_results[,11],D238SL272D_results[,14],D238SW269D_results[,14],D238S278K_results[,16],D238S278K_results[,17]),
  mutants =factor(rep(c("WT", "A", "B","C","C2"), times=c(length(WT_results[,11]),length(D238SL272D_results[,14]),
                                                         length(D238SW269D_results[,14]),length(D238S278K_results[,16]),length(D238S278K_results[,17])))))
ANOVA_His238vsAsp <- aov(His238vsAsp~mutants,data = His238vsAsp_data)
TukeyHSD(ANOVA_His238vsAsp);plot(TukeyHSD(ANOVA_His238vsAsp),las=2,cex=0.5,col="darkblue")+mtext("Histidine-Acid distance of potential active site",line=0.5)

SASASer238vsHis_data <- data.frame(
  SASASer238vsHis=c(WT_results[,13],D238SL272D_results[,16],D238SW269D_results[,16],D238S278K_results[,18]),
  mutants =factor(rep(c("WT", "A", "B","C"), times=c(length(WT_results[,13]),length(D238SL272D_results[,16]),
                                                         length(D238SW269D_results[,16]),length(D238S278K_results[,18])))))
ANOVA_SASASer238vsHis <- aov(SASASer238vsHis~mutants,data = SASASer238vsHis_data)
TukeyHSD(ANOVA_SASASer238vsHis);plot(TukeyHSD(ANOVA_SASASer238vsHis),las=2,cex=0.5,col="darkblue")+mtext("SASA of nucleophile O of Ser of potential active site",line=0.5)


# Test the correlation between all the sets of the data
install.packages("ggpubr");install.packages("ggcorrplot")
library(ggpubr);library(ggcorrplot)
corr <- round(cor(WT_results[,6:10]),3)
p.mat <- cor_pmat(WT_results[,6:10])
ggcorrplot(corr,   ggtheme = ggplot2::theme_gray,
           colors = c("#F8A31B", "white", "#8DD3C7"),lab=TRUE,p.mat = p.mat,title = "Correlation matrix of the PELE quantitative metrics")+theme(plot.title = element_text(hjust = 0.5))

# -----------------------------
# MD_simulations_results

# WT data frame

WT_MD_results_1 <- read.table(file="Results/MD_WT/ServsHis.dat")
WT_MD_results_2 <- read.table(file="Results/MD_WT/AspOD1vsHis.dat")
WT_MD_results_3 <- read.table(file="Results/MD_WT/AspOD2vsHis.dat")
WT_MD_results <- WT_MD_results_1
WT_MD_results <- cbind(WT_MD_results,WT_MD_results_2[,2])
WT_MD_results <- cbind(WT_MD_results,WT_MD_results_3[,2])
names(WT_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist")

# T124D data frame

L272D_MD_results_1 <- read.table(file="Results2/MD_T124D/ServsHis.dat")
L272D_MD_results_2 <- read.table(file="Results2/MD_T124D/AspOD1vsHis.dat")
L272D_MD_results_3 <- read.table(file="Results2/MD_T124D/AspOD2vsHis.dat")
L272D_MD_results_4 <- read.table(file="Results2/MD_T124D/Ser133vsHis.dat")
L272D_MD_results_5 <- read.table(file="Results2/MD_T124D/Asp124OD1vsHis.dat")
L272D_MD_results_6 <- read.table(file="Results2/MD_T124D/Asp124OD2vsHis.dat")
L272D_MD_results <- L272D_MD_results_1
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_2[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_3[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_4[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_5[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_6[,2])
names(L272D_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser133vsHis_dist","Asp124OD1vsHis_dist","Asp124OD2vsHis_dist")

summary(WT_MD_results)
summary(L272D_MD_results)

ServsHisdiffL272D <- wilcox.test(L272D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=4.25)$p.value
AspOD1vsHisdiffL272D <- wilcox.test(L272D_MD_results[,4],WT_MD_results[,4], alternative="greater",mu=2)$p.value
Ser238vsHisdiffL272D <- wilcox.test(L272D_MD_results[,5],L272D_MD_results[,2], alternative="less",mu=-2.5)$p.value
Asp272vsHisdiffL272D <- wilcox.test(L272D_MD_results[,6],L272D_MD_results[,4], alternative="greater",mu=2.8)$p.value

# Q137E data frame

L272D_MD_results_1 <- read.table(file="Results2/MD_Q137E/ServsHis.dat")
L272D_MD_results_2 <- read.table(file="Results2/MD_Q137E/AspOD1vsHis.dat")
L272D_MD_results_3 <- read.table(file="Results2/MD_Q137E/AspOD2vsHis.dat")
L272D_MD_results_4 <- read.table(file="Results2/MD_Q137E/Ser133vsHis.dat")
L272D_MD_results_5 <- read.table(file="Results2/MD_Q137E/GluOE1vsHis.dat")
L272D_MD_results_6 <- read.table(file="Results2/MD_Q137E/GluOE2vsHis.dat")
L272D_MD_results <- L272D_MD_results_1
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_2[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_3[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_4[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_5[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_6[,2])
names(L272D_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","GluOE1vsHis_dist","GluOE2vsHis_dist")

ServsHisdiffL272D <- wilcox.test(L272D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.7)$p.value
AspOD1vsHisdiffL272D <- wilcox.test(L272D_MD_results[,4],WT_MD_results[,4], alternative="less")$p.value
Ser238vsHisdiffL272D <- wilcox.test(L272D_MD_results[,5],L272D_MD_results[,2], alternative="greater",mu=6.9)$p.value
Asp272vsHisdiffL272D <- wilcox.test(L272D_MD_results[,6],L272D_MD_results[,4], alternative="greater",mu=4.6)$p.value

# DE data frame

L272D_MD_results_1 <- read.table(file="Results2/MD_DE/ServsHis.dat")
L272D_MD_results_2 <- read.table(file="Results2/MD_DE/AspOD1vsHis.dat")
L272D_MD_results_3 <- read.table(file="Results2/MD_DE/AspOD2vsHis.dat")
L272D_MD_results_4 <- read.table(file="Results2/MD_DE/Ser133vsHis.dat")
L272D_MD_results_5 <- read.table(file="Results2/MD_DE/Asp124OD1vsHis.dat")
L272D_MD_results_6 <- read.table(file="Results2/MD_DE/Asp124OD2vsHis.dat")
L272D_MD_results_7 <- read.table(file="Results2/MD_DE/GluOE1vsHis.dat")
L272D_MD_results_8 <- read.table(file="Results2/MD_DE/GluOE2vsHis.dat")
L272D_MD_results <- L272D_MD_results_1
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_2[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_3[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_4[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_5[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_6[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_7[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_8[,2])
names(L272D_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser133vsHis_dist","Asp124OD1vsHis_dist","Asp124OD2vsHis_dist","GluOE1vsHis_dist","GluOE2vsHis_dist")

ServsHisdiffL272D <- wilcox.test(L272D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.15)$p.value
AspOD1vsHisdiffL272D <- wilcox.test(L272D_MD_results[,4],WT_MD_results[,4], alternative="greater",mu=0.01)$p.value
Ser238vsHisdiffL272D <- wilcox.test(L272D_MD_results[,5],L272D_MD_results[,2], alternative="greater",mu=1.4)$p.value
Asp272vsHisdiffL272D <- wilcox.test(L272D_MD_results[,6],L272D_MD_results[,4], alternative="greater",mu=3.4)$p.value
GluvsHisdiffL272D <- wilcox.test(L272D_MD_results[,9],L272D_MD_results[,4], alternative="greater",mu=0.7)$p.value
