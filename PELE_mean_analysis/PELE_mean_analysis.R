# Adding the e1071 library to calculate the skewness of the sets of the data.
install.packages("e1071")
library(e1071)

# The datasets of PELE results are loaded for the different mutants
WT_results <- read.csv(file="Results/Results_WT_2HEMHET.csv")
D238SL272D_results <- read.csv(file="Results/Results_D238SL272D_2HEMHET.csv")
D238SW269D_results <- read.csv(file="Results/Results_D238SW269D_2HEMHET.csv")
D238S278K_results <- read.csv(file="Results/Results_D238S278K_2HEMHET.csv")
D238S278R_results <- read.csv(file="Results/Results_D238S278R_2HEMHET.csv")

# Comparison of catalytic distances between all PELE simulations using the Kruskal-Wallis test. Also, TukeyHSD test is shown.
bartlett.test(list(WT_results[,6],D238SL272D_results[,8],D238SW269D_results[,8],D238S278K_results[,8],D238S278R_results[,8]))
kruskal.test(list(WT_results[,6],D238SL272D_results[,8],D238SW269D_results[,8],D238S278K_results[,8],D238S278R_results[,8]))
?kruskal.test

ServsHis_data <- data.frame(
  ServsHis=c(WT_results[,6],D238SL272D_results[,8],D238SW269D_results[,8],D238S278K_results[,8],D238S278R_results[,8]),
  mutants =factor(rep(c("WT", "A", "B","C","D"), times=c(length(WT_results[,6]),length(D238SL272D_results[,8]),
                                                                         length(D238SW269D_results[,8]),length(D238S278K_results[,8]),length(D238S278R_results[,8])))))
ANOVA_ServsHis <- aov(ServsHis~mutants,data = ServsHis_data)
TukeyHSD(ANOVA_ServsHis);plot(TukeyHSD(ANOVA_ServsHis),las=2,cex=0.5,col="darkblue")+mtext("Serine-Histidine distance of main active site",line=0.5)

HisvsAsp_data <- data.frame(
  HisvsAsp=c(WT_results[,7],D238SL272D_results[,11],D238SW269D_results[,11],D238S278K_results[,11],D238S278R_results[,11]),
  mutants =factor(rep(c("WT", "A", "B","C","D"), times=c(length(WT_results[,7]),length(D238SL272D_results[,11]),
                                                         length(D238SW269D_results[,11]),length(D238S278K_results[,11]),length(D238S278R_results[,11])))))
ANOVA_HisvsAsp <- aov(HisvsAsp~mutants,data = HisvsAsp_data)
TukeyHSD(ANOVA_HisvsAsp);plot(TukeyHSD(ANOVA_HisvsAsp),las=2,cex=0.5,col="darkblue")+mtext("Histidine-Aspartate distance of main active site",line=0.5)

Ser238vsHis_data <- data.frame(
  Ser238vsHis=c(WT_results[,6],D238SL272D_results[,9],D238SW269D_results[,9],D238S278K_results[,9],D238S278R_results[,9]),
  mutants =factor(rep(c("WT", "A", "B","C","D"), times=c(length(WT_results[,6]),length(D238SL272D_results[,9]),
                                                         length(D238SW269D_results[,9]),length(D238S278K_results[,9]),length(D238S278R_results[,9])))))
ANOVA_Ser238vsHis <- aov(Ser238vsHis~mutants,data = Ser238vsHis_data)
TukeyHSD(ANOVA_Ser238vsHis);plot(TukeyHSD(ANOVA_Ser238vsHis),las=2,cex=0.5,col="darkblue")+mtext("Serine238-Histidine distance of potential active site",line=0.5)

His238vsAsp_data <- data.frame(
  His238vsAsp=c(WT_results[,7],D238SL272D_results[,14],D238SW269D_results[,14],D238S278K_results[,14],D238S278R_results[,14]),
  mutants =factor(rep(c("WT", "A", "B","C","D"), times=c(length(WT_results[,7]),length(D238SL272D_results[,14]),
                                                         length(D238SW269D_results[,14]),length(D238S278K_results[,14]),length(D238S278R_results[,14])))))
ANOVA_His238vsAsp <- aov(His238vsAsp~mutants,data = His238vsAsp_data)
TukeyHSD(ANOVA_His238vsAsp);plot(TukeyHSD(ANOVA_His238vsAsp),las=2,cex=0.5,col="darkblue")+mtext("Histidine-Aspartate distance of potential active site",line=0.5)

SASASer238vsHis_data <- data.frame(
  SASASer238vsHis=c(WT_results[,9],D238SL272D_results[,16],D238SW269D_results[,16],D238S278K_results[,16],D238S278R_results[,16]),
  mutants =factor(rep(c("WT", "A", "B","C","D"), times=c(length(WT_results[,9]),length(D238SL272D_results[,16]),
                                                         length(D238SW269D_results[,16]),length(D238S278K_results[,16]),length(D238S278R_results[,16])))))
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

# K281L data frame

K281L_MD_results_1 <- read.table(file="Results/MD_K281L/ServsHis.dat")
K281L_MD_results_2 <- read.table(file="Results/MD_K281L/AspOD1vsHis.dat")
K281L_MD_results_3 <- read.table(file="Results/MD_K281L/AspOD2vsHis.dat")
K281L_MD_results_4 <- read.table(file="Results/MD_K281L/Ser238vsHis.dat")
K281L_MD_results_5 <- read.table(file="Results/MD_K281L/Asp272OD1vsHis.dat")
K281L_MD_results_6 <- read.table(file="Results/MD_K281L/Asp272OD2vsHis.dat")
K281L_MD_results <- K281L_MD_results_1
K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_2[,2]);K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_3[,2])
K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_4[,2]);K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_5[,2])
K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_6[,2])
names(K281L_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp272OD1vsHis_dist","Asp272OD2vsHis_dist")

ServsHisdiffK281L <- wilcox.test(K281L_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.11)$p.value
AspOD1vsHisdiffK281L <- wilcox.test(K281L_MD_results[,3],WT_MD_results[,4], alternative="greater",mu=0.06)$p.value
Ser238vsHisdiffK281L <- wilcox.test(K281L_MD_results[,5],K281L_MD_results[,2], alternative="less",mu=-0.75)$p.value
Asp272vsHisdiffK281L <- wilcox.test(K281L_MD_results[,6],K281L_MD_results[,3], alternative="greater",mu=2.56)$p.value

# -----------------------------
  
# Tools to visualize the datasets individually
summary(D238SL272D_results)
quantile(WT_results[,10])
ServsHis <- D238SL272D_results[,8]
Ser238vsHis <- D238SL272D_results[,9]
boxplot(PELE_results[,8],PELE_results[,9])
pairs(PELE_results[,10:13])

# Analysis of the normality of the test
shapiro.test(PELE_results[,8])
qqnorm(PELE_results[,8])
qqline(PELE_results[,8])
hist(PELE_results[,9])
skewness(PELE_results[,9]) # To check the asimmetry of the real random distribution

# Statistic tests of two samples of both assuming a normal distribution or not.
t.test(PELE_results[,8], PELE_results[,9], alternative="greater", var.equal=FALSE)
t.test(PELE_results[,13], PELE_results[,11], alternative="greater", var.equal=FALSE,mu=1.5)
# Wrong since they are not normal. We therefore use Wilcox test
wilcox.test(PELE_results[,9], PELE_results[,8], alternative="greater",mu=0.055)
boxplot(PELE_results[,9],PELE_results[,8])

# Statistic tests of more than two samples of not assuming a normal distribution.
bartlett.test(list(PELE_results[,4],PELE_results[,5]))
aov(PELE_results[,9]~PELE_results[,8])
kruskal.test(list(ServsHis,Ser238vsHis))
