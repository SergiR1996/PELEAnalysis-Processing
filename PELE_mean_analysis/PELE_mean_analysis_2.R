# Adding the e1071 library to calculate the skewness of the sets of the data.
# install.packages("e1071")
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
# install.packages("ggpubr");install.packages("ggcorrplot")
library(ggpubr);library(ggcorrplot)
corr <- round(cor(WT_results[,6:10]),3)
p.mat <- cor_pmat(WT_results[,6:10])
ggcorrplot(corr,   ggtheme = ggplot2::theme_gray,
           colors = c("#F8A31B", "white", "#8DD3C7"),lab=TRUE,p.mat = p.mat,title = "Correlation matrix of the PELE quantitative metrics")+theme(plot.title = element_text(hjust = 0.5))

# -----------------------------
# MD_simulations_results

setwd("/home/sergiroda/repos/PELEAnalysis-Processing/PELE_mean_analysis")

# WT data frame

WT_MD_results_1 <- read.table(file="Results/MD_WT/ServsHis_tot.dat")
WT_MD_results_2 <- read.table(file="Results/MD_WT/AspOD1vsHis_tot.dat")
WT_MD_results_3 <- read.table(file="Results/MD_WT/AspOD2vsHis_tot.dat")
WT_MD_results_4 <- read.table(file="Results/MD_WT/localrmsd_tot.xvg")
WT_MD_results_5 <- read.table(file="Results/MD_WT/rmsd_tot.xvg")
WT_MD_results <- WT_MD_results_1
WT_MD_results <- cbind(WT_MD_results,WT_MD_results_2[,2])
WT_MD_results <- cbind(WT_MD_results,WT_MD_results_3[,2])
names(WT_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist")

# T124D data frame

T124D_MD_results_1 <- read.table(file="Results2/MD_T124D/ServsHis.dat")
T124D_MD_results_2 <- read.table(file="Results2/MD_T124D/AspOD1vsHis.dat")
T124D_MD_results_3 <- read.table(file="Results2/MD_T124D/AspOD2vsHis.dat")
T124D_MD_results_4 <- read.table(file="Results2/MD_T124D/Ser133vsHis.dat")
T124D_MD_results_5 <- read.table(file="Results2/MD_T124D/Asp124OD1vsHis.dat")
T124D_MD_results_6 <- read.table(file="Results2/MD_T124D/Asp124OD2vsHis.dat")
T124D_MD_results_7 <- read.table(file="Results2/MD_T124D/localrmsd.xvg")
T124D_MD_results_8 <- read.table(file="Results2/MD_T124D/localrmsd_potential.xvg")
T124D_MD_results_9 <- read.table(file="Results2/MD_T124D/rmsd.xvg")
T124D_MD_results <- T124D_MD_results_1
T124D_MD_results <- cbind(T124D_MD_results,T124D_MD_results_2[,2]);T124D_MD_results <- cbind(T124D_MD_results,T124D_MD_results_3[,2])
T124D_MD_results <- cbind(T124D_MD_results,T124D_MD_results_4[,2]);T124D_MD_results <- cbind(T124D_MD_results,T124D_MD_results_5[,2])
T124D_MD_results <- cbind(T124D_MD_results,T124D_MD_results_6[,2])
names(T124D_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser133vsHis_dist","Asp124OD1vsHis_dist","Asp124OD2vsHis_dist")

ServsHisdiffL272D <- wilcox.test(L272D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=4.25)$p.value
AspOD1vsHisdiffL272D <- wilcox.test(L272D_MD_results[,4],WT_MD_results[,4], alternative="greater",mu=2)$p.value
Ser238vsHisdiffL272D <- wilcox.test(L272D_MD_results[,5],L272D_MD_results[,2], alternative="less",mu=-2.5)$p.value
Asp272vsHisdiffL272D <- wilcox.test(L272D_MD_results[,6],L272D_MD_results[,4], alternative="greater",mu=2.8)$p.value

# Q137E data frame

Q137E_MD_results_1 <- read.table(file="Results2/MD_Q137E/ServsHis.dat")
Q137E_MD_results_2 <- read.table(file="Results2/MD_Q137E/AspOD1vsHis.dat")
Q137E_MD_results_3 <- read.table(file="Results2/MD_Q137E/AspOD2vsHis.dat")
Q137E_MD_results_4 <- read.table(file="Results2/MD_Q137E/Ser133vsHis.dat")
Q137E_MD_results_5 <- read.table(file="Results2/MD_Q137E/GluOE1vsHis.dat")
Q137E_MD_results_6 <- read.table(file="Results2/MD_Q137E/GluOE2vsHis.dat")
Q137E_MD_results_7 <- read.table(file="Results2/MD_Q137E/localrmsd.xvg")
Q137E_MD_results_8 <- read.table(file="Results2/MD_Q137E/localrmsd_potential.xvg")
Q137E_MD_results_9 <- read.table(file="Results2/MD_Q137E/rmsd.xvg")
Q137E_MD_results <- Q137E_MD_results_1
Q137E_MD_results <- cbind(Q137E_MD_results,Q137E_MD_results_2[,2]);Q137E_MD_results <- cbind(Q137E_MD_results,Q137E_MD_results_3[,2])
Q137E_MD_results <- cbind(Q137E_MD_results,Q137E_MD_results_4[,2]);Q137E_MD_results <- cbind(Q137E_MD_results,Q137E_MD_results_5[,2])
Q137E_MD_results <- cbind(Q137E_MD_results,Q137E_MD_results_6[,2])
names(Q137E_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","GluOE1vsHis_dist","GluOE2vsHis_dist")

ServsHisdiffL272D <- wilcox.test(L272D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.7)$p.value
AspOD1vsHisdiffL272D <- wilcox.test(L272D_MD_results[,4],WT_MD_results[,4], alternative="less")$p.value
Ser238vsHisdiffL272D <- wilcox.test(L272D_MD_results[,5],L272D_MD_results[,2], alternative="greater",mu=6.9)$p.value
Asp272vsHisdiffL272D <- wilcox.test(L272D_MD_results[,6],L272D_MD_results[,4], alternative="greater",mu=4.6)$p.value

# Q137K data frame

Q137K_MD_results_1 <- read.table(file="Results2/MD_Q137K/ServsHis.dat")
Q137K_MD_results_2 <- read.table(file="Results2/MD_Q137K/AspOD1vsHis.dat")
Q137K_MD_results_3 <- read.table(file="Results2/MD_Q137K/AspOD2vsHis.dat")
Q137K_MD_results_4 <- read.table(file="Results2/MD_Q137K/Ser133vsHis.dat")
Q137K_MD_results_5 <- read.table(file="Results2/MD_Q137K/Asp124OD1vsHis.dat")
Q137K_MD_results_6 <- read.table(file="Results2/MD_Q137K/Asp124OD2vsHis.dat")
Q137K_MD_results_7 <- read.table(file="Results2/MD_Q137K/localrmsd.xvg")
Q137K_MD_results_8 <- read.table(file="Results2/MD_Q137K/localrmsd_potential.xvg")
Q137K_MD_results_9 <- read.table(file="Results2/MD_Q137K/rmsd.xvg")
Q137K_MD_results <- Q137K_MD_results_1
Q137K_MD_results <- cbind(Q137K_MD_results,Q137K_MD_results_2[,2]);Q137K_MD_results <- cbind(Q137K_MD_results,Q137K_MD_results_3[,2])
Q137K_MD_results <- cbind(Q137K_MD_results,Q137K_MD_results_4[,2]);Q137K_MD_results <- cbind(Q137K_MD_results,Q137K_MD_results_5[,2])
Q137K_MD_results <- cbind(Q137K_MD_results,Q137K_MD_results_6[,2])
names(Q137K_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp124OD1vsHis_dist","Asp124OD2vsHis_dist")

# DE data frame

DE_MD_results_1 <- read.table(file="Results2/MD_DE/ServsHis.dat")
DE_MD_results_2 <- read.table(file="Results2/MD_DE/AspOD1vsHis.dat")
DE_MD_results_3 <- read.table(file="Results2/MD_DE/AspOD2vsHis.dat")
DE_MD_results_4 <- read.table(file="Results2/MD_DE/Ser133vsHis.dat")
DE_MD_results_5 <- read.table(file="Results2/MD_DE/Asp124OD1vsHis.dat")
DE_MD_results_6 <- read.table(file="Results2/MD_DE/Asp124OD2vsHis.dat")
DE_MD_results_7 <- read.table(file="Results2/MD_DE/GluOE1vsHis.dat")
DE_MD_results_8 <- read.table(file="Results2/MD_DE/GluOE2vsHis.dat")
DE_MD_results_9 <- read.table(file="Results2/MD_DE/localrmsd.xvg")
DE_MD_results_10 <- read.table(file="Results2/MD_DE/localrmsd_potential.xvg")
DE_MD_results_11 <- read.table(file="Results2/MD_DE/rmsd.xvg")
DE_MD_results <- DE_MD_results_1
DE_MD_results <- cbind(DE_MD_results,DE_MD_results_2[,2]);DE_MD_results <- cbind(DE_MD_results,DE_MD_results_3[,2])
DE_MD_results <- cbind(DE_MD_results,DE_MD_results_4[,2]);DE_MD_results <- cbind(DE_MD_results,DE_MD_results_5[,2])
DE_MD_results <- cbind(DE_MD_results,DE_MD_results_6[,2]);DE_MD_results <- cbind(DE_MD_results,DE_MD_results_7[,2])
DE_MD_results <- cbind(DE_MD_results,DE_MD_results_8[,2])
names(DE_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser133vsHis_dist","Asp124OD1vsHis_dist","Asp124OD2vsHis_dist","GluOE1vsHis_dist","GluOE2vsHis_dist")

ServsHisdiffL272D <- wilcox.test(L272D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.15)$p.value
AspOD1vsHisdiffL272D <- wilcox.test(L272D_MD_results[,4],WT_MD_results[,4], alternative="greater",mu=0.01)$p.value
Ser238vsHisdiffL272D <- wilcox.test(L272D_MD_results[,5],L272D_MD_results[,2], alternative="greater",mu=1.4)$p.value
Asp272vsHisdiffL272D <- wilcox.test(L272D_MD_results[,6],L272D_MD_results[,4], alternative="greater",mu=3.4)$p.value
GluvsHisdiffL272D <- wilcox.test(L272D_MD_results[,9],L272D_MD_results[,4], alternative="greater",mu=0.7)$p.value

# T124L

mean(T124L_MD_results_6[,2]);sd(T124L_MD_results_6[,2])

T124L_MD_results_1 <- read.table(file="Results2/MD_T124L/ServsHis.dat")
T124L_MD_results_2 <- read.table(file="Results2/MD_T124L/AspOD1vsHis.dat")
T124L_MD_results_3 <- read.table(file="Results2/MD_T124L/AspOD2vsHis.dat")
T124L_MD_results_4 <- read.table(file="Results2/MD_T124L/Ser133vsHis.dat")
T124L_MD_results_5 <- read.table(file="Results2/MD_T124L/GluOE1vsHis.dat")
T124L_MD_results_6 <- read.table(file="Results2/MD_T124L/GluOE2vsHis.dat")
T124L_MD_results_7 <- read.table(file="Results2/MD_T124L/localrmsd.xvg")
T124L_MD_results_8 <- read.table(file="Results2/MD_T124L/localrmsd_potential.xvg")
T124L_MD_results_9 <- read.table(file="Results2/MD_T124L/rmsd.xvg")
T124L_MD_results <- T124L_MD_results_1
T124L_MD_results <- cbind(T124L_MD_results,T124L_MD_results_2[,2]);T124L_MD_results <- cbind(T124L_MD_results,T124L_MD_results_3[,2])
T124L_MD_results <- cbind(T124L_MD_results,T124L_MD_results_4[,2]);T124L_MD_results <- cbind(T124L_MD_results,T124L_MD_results_5[,2])
T124L_MD_results <- cbind(T124L_MD_results,T124L_MD_results_6[,2])
names(T124L_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser133vsHis_dist","Asp124OD1vsHis_dist","Asp124OD2vsHis_dist")

mean(T124L_MD_results_4[,2])

plot(density(WT_MD_results[,2]))
plot(density(L272D_MD_results_3[,2]),col=3)
lines(density(L272D_MD_results_5[,2]),col=2)

# T124L total

T124LT_MD_results_1 <- read.table(file="Results2/MD_T124L/ServsHis_tot.dat")
T124LT_MD_results_2 <- read.table(file="Results2/MD_T124L/AspOD1vsHis_tot.dat")
T124LT_MD_results_3 <- read.table(file="Results2/MD_T124L/AspOD2vsHis_tot.dat")
T124LT_MD_results_4 <- read.table(file="Results2/MD_T124L/Ser133vsHis_tot.dat")
T124LT_MD_results_5 <- read.table(file="Results2/MD_T124L/GluOE1vsHis_tot.dat")
T124LT_MD_results_6 <- read.table(file="Results2/MD_T124L/GluOE2vsHis_tot.dat")
T124LT_MD_results_7 <- read.table(file="Results2/MD_T124L/localrmsd_tot.xvg")
T124LT_MD_results_8 <- read.table(file="Results2/MD_T124L/localrmsd_potential_tot.xvg")
T124LT_MD_results_9 <- read.table(file="Results2/MD_T124L/rmsd_tot.xvg")
T124LT_MD_results_10 <- read.table(file="Results2/MD_T124L/rmsd_2.xvg")
T124LT_MD_results_11 <- read.table(file="Results2/MD_T124L/rmsd_3.xvg")
T124LT_MD_results_RMSD <- T124LT_MD_results_10[,2]
T124LT_MD_results_RMSD <- cbind(T124LT_MD_results_RMSD,T124LT_MD_results_11[,2])
T124LT_MD_results_RMSD <- cbind(T124LT_MD_results_RMSD,T124L_MD_results_9[,2])
RMSD_WT <- c()
for (i in 1:length(T124L_MD_results_9[,2])) {
   RMSD_WT <- c(RMSD_WT,mean(T124LT_MD_results_RMSD[i,1:3]))
 }
 for (i in (length(T124L_MD_results_9[,2])+1):length(T124LT_MD_results_11[,2])){
   RMSD_WT <- c(RMSD_WT,mean(T124LT_MD_results_RMSD[i,1:2]))
}
for (i in (length(T124LT_MD_results_11[,2])+1):length(T124LT_MD_results_10[,2])){
  RMSD_WT <- c(RMSD_WT,mean(T124LT_MD_results_RMSD[i,1]))
}
RMSD_WT <- as.data.frame(RMSD_WT)
rmsd <- cbind(T124LT_MD_results_10[,1],RMSD_WT)
write.table(rmsd, "rmsd_tot.xvg", sep="\t",append = FALSE)

# Density plots with the ggplot package
ggplot()+geom_density(data=L272D_MD_results,aes(ServsHis_dist,fill="ServsHis"),
                      bw=0.1,alpha=0.5)+
  geom_density(data=L272D_MD_results,aes(AspOD2vsHis_dist,fill="AspvsHis"),
               bw=0.1,alpha=0.5)+theme_minimal()+
  guides(fill=guide_legend(title="Metric"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))


df <- as.data.frame(L272D_MD_results_4)
data <- stack(df)
names(data) <- c("values","Metric")

ggplot(data,aes(x=values))+geom_density(aes(group=Metric,colour=Metric,fill=Metric),
                      bw=0.1,alpha=0.5)


#------------------------------------------------------

# RMSD of the simulation

ggplot(data=WT_MD_results_8,aes(x=WT_MD_results_8[,1],y=WT_MD_results_8[,2],
                                color="WT"))+geom_line(alpha=0.7)+
  geom_line(data=T124D_MD_results_9,
            aes(x=T124D_MD_results_9[,1],
                y=T124D_MD_results_9[,2],color="D129HT124DQ137N"),alpha=0.7)+
  geom_line(data=DE_MD_results_11,
            aes(x=DE_MD_results_11[,1],
                y=DE_MD_results_11[,2],color="D129HT124DQ137E"),alpha=0.7)+
  geom_line(data=T124L_MD_results_9,
            aes(x=T124L_MD_results_9[,1],
                y=T124L_MD_results_9[,2],color="D129HT124LQ137E"),alpha=0.7)+
  geom_line(data=Q137K_MD_results_9,
            aes(x=Q137K_MD_results_9[,1],
                y=Q137K_MD_results_9[,2],color="D129HT124DQ137K"),alpha=0.7)+
  geom_line(data=Q137E_MD_results_9,
            aes(x=Q137E_MD_results_9[,1],
                y=Q137E_MD_results_9[,2],color="D129HQ137E"),alpha=0.7)+
  theme_minimal()+ggtitle("Global RMSD (least-squares fit to backbone)")+xlab("Time (ps)")+ylab("RMSD (nm)")+
  guides(color=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

# LOCAL RMSD MAIN ACTIVE SITE

ggplot()+geom_density(data=WT_MD_results_4,aes(WT_MD_results_4[,2],fill="WT"),
                      bw=0.01,alpha=0.5)+
  geom_density(data=T124D_MD_results_7,aes(T124D_MD_results_7[,2],fill="D129HT124DQ137N"),
               bw=0.01,alpha=0.5)+
  geom_density(data=Q137E_MD_results_7,aes(Q137E_MD_results_7[,2],fill="D129HQ137E"),
               bw=0.01,alpha=0.5)+
  geom_density(data=Q137K_MD_results_7,aes(Q137K_MD_results_7[,2],fill="D129HT124DQ137K"),
               bw=0.01,alpha=0.5)+
  geom_density(data=DE_MD_results_9,aes(DE_MD_results_9[,2],fill="D129HT124DQ137E"),
               bw=0.01,alpha=0.5)+
  geom_density(data=T124L_MD_results_7,aes(T124L_MD_results_7[,2],fill="D129HT124LQ137E"),
               bw=0.01,alpha=0.5)+
  theme_minimal()+ggtitle("Local RMSD for the main catalytic triad")+xlab("Local RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

# LOCAL RMSD POTENTIAl ACTIVE SITE

ggplot()+
  geom_density(data=T124D_MD_results_8,aes(T124D_MD_results_8[,2],fill="D129HT124DQ137N"),
               bw=0.01,alpha=0.5)+
  geom_density(data=Q137E_MD_results_8,aes(Q137E_MD_results_8[,2],fill="D129HQ137E"),
               bw=0.01,alpha=0.5)+
  geom_density(data=Q137K_MD_results_8,aes(Q137K_MD_results_8[,2],fill="D129HT124DQ137K"),
               bw=0.01,alpha=0.5)+
  geom_density(data=DE_MD_results_10,aes(DE_MD_results_10[,2],fill="D129HT124DQ137E"),
               bw=0.01,alpha=0.5)+
  geom_density(data=T124L_MD_results_8,aes(T124L_MD_results_8[,2],fill="D129HT124LQ137E"),
               bw=0.01,alpha=0.5)+
  theme_minimal()+ggtitle("Local RMSD for the potential catalytic triad")+xlab("Local RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

# Catalytic distances MAIN

ggplot()+
  geom_density(data=WT_MD_results_1,aes(WT_MD_results_1[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124D_MD_results_1,aes(T124D_MD_results_1[,2],fill="D129HT124DQ137N"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137E_MD_results_1,aes(Q137E_MD_results_1[,2],fill="D129HQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137K_MD_results_1,aes(Q137K_MD_results_1[,2],fill="D129HT124DQ137K"),
               bw=0.1,alpha=0.5)+
  geom_density(data=DE_MD_results_1,aes(DE_MD_results_1[,2],fill="D129HT124DQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_1,aes(T124L_MD_results_1[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His distance for the main catalytic triad")+xlab("Ser-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

ggplot()+
  geom_density(data=WT_MD_results_3,aes(WT_MD_results_3[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124D_MD_results_3,aes(T124D_MD_results_3[,2],fill="D129HT124DQ137N"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137E_MD_results_3,aes(Q137E_MD_results_3[,2],fill="D129HQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137K_MD_results_3,aes(Q137K_MD_results_3[,2],fill="D129HT124DQ137K"),
               bw=0.1,alpha=0.5)+
  geom_density(data=DE_MD_results_3,aes(DE_MD_results_3[,2],fill="D129HT124DQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_3,aes(T124L_MD_results_3[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Asp-His distance for the main catalytic triad")+xlab("Asp-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

# Catalytic distances POTENTIAL

ggplot()+
  geom_density(data=WT_MD_results_1,aes(WT_MD_results_1[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124D_MD_results_4,aes(T124D_MD_results_4[,2],fill="D129HT124DQ137N"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137E_MD_results_4,aes(Q137E_MD_results_4[,2],fill="D129HQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137K_MD_results_4,aes(Q137K_MD_results_4[,2],fill="D129HT124DQ137K"),
               bw=0.1,alpha=0.5)+
  geom_density(data=DE_MD_results_4,aes(DE_MD_results_4[,2],fill="D129HT124DQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_4,aes(T124L_MD_results_4[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His distance for the potential catalytic triad")+xlab("Ser-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

summary(Q137K_MD_results_6[,2])
sd(Q137K_MD_results_6[,2])

ggplot()+
  geom_density(data=WT_MD_results_3,aes(WT_MD_results_3[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124D_MD_results_5,aes(T124D_MD_results_5[,2],fill="D129HT124DQ137N"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137E_MD_results_6,aes(Q137E_MD_results_6[,2],fill="D129HQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Q137K_MD_results_6,aes(Q137K_MD_results_6[,2],fill="D129HT124DQ137K"),
               bw=0.1,alpha=0.5)+
  geom_density(data=DE_MD_results_8,aes(DE_MD_results_8[,2],fill="D129HT124DQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_5,aes(T124L_MD_results_5[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Acid-His distance for the potential catalytic triad")+xlab("Acid-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

summary(WT_MD_results)
summary(DE_MD_results)
sd(DE_MD_results$Ser133vsHis_dist);sd(DE_MD_results$Asp124OD1vsHis_dist)
sd(DE_MD_results$Asp124OD2vsHis_dist)
sd(DE_MD_results$GluOE1vsHis_dist);sd(DE_MD_results$GluOE2vsHis_dist)
sd(T124L_MD_results_5[,2])
summary(T124LT_MD_results_8);sd(T124LT_MD_results_8[,2])
summary(Q137K_MD_results_8[,2])
summary(DE_MD_results_8)

var.test(WT_MD_results[,2],DE_MD_results_4[,2])
skewness(DE_MD_results_9[,2])

t.test(DE_MD_results_4[,2],WT_MD_results[,2]
       ,var.equal = FALSE,alternative = "greater")

var.test(WT_MD_results[,4],DE_MD_results_7[,2])

t.test(DE_MD_results_7[,2],WT_MD_results[,4]
       ,var.equal = FALSE,alternative = "greater")

#------------------------------------------------------  
  
ComplexL_MD_results_1 <- read.table(file="Results2/MD_Complex/ServsHis.dat")
ComplexL_MD_results_2 <- read.table(file="Results2/MD_Complex/GluOE1vsHis.dat")
ComplexL_MD_results_3 <- read.table(file="Results2/MD_Complex/GluOE2vsHis.dat")
ComplexL_MD_results_4 <- read.table(file="Results2/MD_Complex/ServsC1.dat")
ComplexL_MD_results_5 <- read.table(file="Results2/MD_Complex/ServsC8.dat")
ComplexL_MD_results_6 <- read.table(file="Results2/MD_Complex/ServsC11.dat")
ComplexL_MD_results_7 <- read.table(file="Results2/MD_Complex/ServsC23.dat")
ComplexL_MD_results_8 <- read.table(file="Results2/MD_Complex/ServsC30.dat")
ComplexL_MD_results_9 <- read.table(file="Results2/MD_Complex/ServsC33.dat")
ComplexL_MD_results_10 <- read.table(file="Results2/MD_Complex/localrmsd.xvg")
ComplexL_MD_results_11 <- read.table(file="Results2/MD_Complex/localrmsd_potential.xvg")

ComplexL_MD_results <- ComplexL_MD_results_4
ComplexL_MD_results <- cbind(ComplexL_MD_results,Complex_MD_results_5[,2]);ComplexL_MD_results <- cbind(ComplexL_MD_results,ComplexL_MD_results_6[,2])
ComplexL_MD_results <- cbind(ComplexL_MD_results,Complex_MD_results_7[,2]);ComplexL_MD_results <- cbind(ComplexL_MD_results,ComplexL_MD_results_8[,2])
ComplexL_MD_results <- cbind(ComplexL_MD_results,Complex_MD_results_9[,2])

ComplexL_MD_results$min <- do.call(pmin, ComplexL_MD_results[,2:7])

# Total

ComplexTL_MD_results_1 <- read.table(file="Results2/MD_Complex/ServsHis_tot.dat")
ComplexTL_MD_results_2 <- read.table(file="Results2/MD_Complex/GluOE1vsHis_tot.dat")
ComplexTL_MD_results_3 <- read.table(file="Results2/MD_Complex/GluOE2vsHis_tot.dat")
ComplexTL_MD_results_4 <- read.table(file="Results2/MD_Complex/ServsC1_tot.dat")
ComplexTL_MD_results_5 <- read.table(file="Results2/MD_Complex/ServsC8_tot.dat")
ComplexTL_MD_results_6 <- read.table(file="Results2/MD_Complex/ServsC11_tot.dat")
ComplexTL_MD_results_7 <- read.table(file="Results2/MD_Complex/ServsC23_tot.dat")
ComplexTL_MD_results_8 <- read.table(file="Results2/MD_Complex/ServsC30_tot.dat")
ComplexTL_MD_results_9 <- read.table(file="Results2/MD_Complex/ServsC33_tot.dat")
ComplexTL_MD_results_10 <- read.table(file="Results2/MD_Complex/localrmsd_tot.xvg")
ComplexTL_MD_results_11 <- read.table(file="Results2/MD_Complex/localrmsd_potential_tot.xvg")

ComplexTL_MD_results <- ComplexTL_MD_results_4
ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_5[,2]);ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_6[,2])
ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_7[,2]);ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_8[,2])
ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_9[,2])

ComplexTL_MD_results$min <- do.call(pmin, ComplexTL_MD_results[,2:7])

ggplot()+geom_density(data=Complex_MD_results,aes(min),
                      bw=0.1,alpha=0.5,fill="#56B2E9")+
  theme_minimal()+ggtitle("Ser-Substrate for the potential catalytic triad")+xlab("Ser-Substrate distance (Ang)")+
  theme(plot.title = element_text(hjust = 0.5))

# Catalytic distances in complex

ggplot()+
  geom_density(data=Complex_MD_results_1,aes(Complex_MD_results_1[,2],fill="T124LD129HQ137E_with_ligand"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_4,aes(T124L_MD_results_4[,2],fill="T124LD129HQ137E"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His distance for the potential catalytic triad")+xlab("Ser-His distance (Ang)")+
  guides(fill=guide_legend(title="MD simulation"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

ggplot()+
  geom_density(data=Complex_MD_results_2,aes(Complex_MD_results_2[,2],fill="T124LD129HQ137E_with_ligand"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_5,aes(T124L_MD_results_5[,2],fill="T124LD129HQ137E"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("His-Asp distance for the potential catalytic triad")+xlab("His-Asp distance (Ang)")+
  guides(fill=guide_legend(title="MD simulation"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

# RMSD in complex

ggplot()+geom_density(data=Complex_MD_results_11,aes(Complex_MD_results_11[,2]),
                      bw=0.01,alpha=0.5,fill="#56B2E9")+
  theme_minimal()+ggtitle("Local RMSD for the potential catalytic triad")+xlab("Local RMSD (nm)")+
  theme(plot.title = element_text(hjust = 0.5))
