# Adding the e1071 library to calculate the skewness of the sets of the data.
#install.packages("e1071")
#install.packages("Peptides", dependencies=TRUE)
library(e1071);library(Peptides);library(ggpubr)

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
  mutants =factor(rep(c("WT", "A", "B","C","D","E"), times=c(length(WT_results[,9]),length(D238SL272D_results[,16]),
                                                         length(D238SW269D_results[,16]),length(D238S278K_results[,16]),length(D238S278R_results[,16])))))
ANOVA_SASASer238vsHis <- aov(SASASer238vsHis~mutants,data = SASASer238vsHis_data)
TukeyHSD(ANOVA_SASASer238vsHis);plot(TukeyHSD(ANOVA_SASASer238vsHis),las=2,cex=0.5,col="darkblue")+mtext("SASA of nucleophile O of Ser of potential active site",line=0.5)

# Test the correlation between all the sets of the data
#install.packages("ggpubr");install.packages("ggcorrplot")
library(ggcorrplot)
corr <- round(cor(WT_results[,6:10]),3)
p.mat <- cor_pmat(WT_results[,6:10])
ggcorrplot(corr,   ggtheme = ggplot2::theme_gray,
           colors = c("#F8A31B", "white", "#8DD3C7"),lab=TRUE,p.mat = p.mat,title = "Correlation matrix of the PELE quantitative metrics")+theme(plot.title = element_text(hjust = 0.5))

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

# -----------------------------
# MD_simulations_results

# WT data frame
setwd("/home/sergiroda/repos/PELEAnalysis-Processing/PELE_mean_analysis")

WT_MD_results_1 <- read.table(file="Results/MD_WT/ServsHis_tot.dat")
WT_MD_results_2 <- read.table(file="Results/MD_WT/AspOD1vsHis_tot.dat")
WT_MD_results_3 <- read.table(file="Results/MD_WT/AspOD2vsHis_tot.dat")
WT_MD_results_4 <- read.table(file="Results/MD_WT/localrmsd_tot.xvg")
WT_MD_results_5 <- read.table(file="Results/MD_WT/rmsd.xvg")
WT_MD_results_6 <- read.table(file="Results/MD_WT/rmsd_2.xvg")
WT_MD_results_7 <- read.table(file="Results/MD_WT/rmsd_3.xvg")
WT_MD_results_8 <- read.table(file="Results/MD_WT/rmsd_tot.xvg")
WT_MD_results <- WT_MD_results_1
WT_MD_results <- cbind(WT_MD_results,WT_MD_results_2[,2])
WT_MD_results <- cbind(WT_MD_results,WT_MD_results_3[,2])
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

names(WT_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist")

# D238SL272D data frame

L272D_MD_results_1 <- read.table(file="Results/MD_L272D/ServsHis.dat")
L272D_MD_results_2 <- read.table(file="Results/MD_L272D/AspOD1vsHis.dat")
L272D_MD_results_3 <- read.table(file="Results/MD_L272D/AspOD2vsHis.dat")
L272D_MD_results_4 <- read.table(file="Results/MD_L272D/Ser238vsHis.dat")
L272D_MD_results_5 <- read.table(file="Results/MD_L272D/Asp272OD1vsHis.dat")
L272D_MD_results_6 <- read.table(file="Results/MD_L272D/Asp272OD2vsHis.dat")
L272D_MD_results_7 <- read.table(file="Results/MD_L272D/localrmsd.xvg")
L272D_MD_results_8 <- read.table(file="Results/MD_L272D/localrmsd_potential.xvg")
L272D_MD_results_9 <- read.table(file="Results/MD_L272D/rmsd.xvg")
L272D_MD_results <- L272D_MD_results_1
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_2[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_3[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_4[,2]);L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_5[,2])
L272D_MD_results <- cbind(L272D_MD_results,L272D_MD_results_6[,2])
names(L272D_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp272OD1vsHis_dist","Asp272OD2vsHis_dist")

xlim <- range(W269D_MD_results[,2],L272D_MD_results[,5])
plot(density(W269D_MD_results[,2]),xlim=xlim)
summary(W269D_MD_results)
lines(density(W269D_MD_results[,5]),col=3)

ServsHisdiffL272D <- wilcox.test(L272D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.06)$p.value
AspOD1vsHisdiffL272D <- wilcox.test(L272D_MD_results[,4],WT_MD_results[,4], alternative="less",mu=-0.02)$p.value
Ser238vsHisdiffL272D <- wilcox.test(L272D_MD_results[,5],L272D_MD_results[,2], alternative="greater",mu=0.9)$p.value
Asp272vsHisdiffL272D <- wilcox.test(L272D_MD_results[,6],L272D_MD_results[,4], alternative="greater",mu=3.4)$p.value

# D238SW269D data frame

W269D_MD_results_1 <- read.table(file="Results/MD_W269D/ServsHis.dat")
W269D_MD_results_2 <- read.table(file="Results/MD_W269D/AspOD1vsHis.dat")
W269D_MD_results_3 <- read.table(file="Results/MD_W269D/AspOD2vsHis.dat")
W269D_MD_results_4 <- read.table(file="Results/MD_W269D/Ser238vsHis.dat")
W269D_MD_results_5 <- read.table(file="Results/MD_W269D/Asp269OD1vsHis.dat")
W269D_MD_results_6 <- read.table(file="Results/MD_W269D/Asp269OD2vsHis.dat")
W269D_MD_results_7 <- read.table(file="Results/MD_W269D/localrmsd.xvg")
W269D_MD_results_8 <- read.table(file="Results/MD_W269D/localrmsd_potential.xvg")
W269D_MD_results_9 <- read.table(file="Results/MD_W269D/rmsd.xvg")
W269D_MD_results <- W269D_MD_results_1
W269D_MD_results <- cbind(W269D_MD_results,W269D_MD_results_2[,2]);W269D_MD_results <- cbind(W269D_MD_results,W269D_MD_results_3[,2])
W269D_MD_results <- cbind(W269D_MD_results,W269D_MD_results_4[,2]);W269D_MD_results <- cbind(W269D_MD_results,W269D_MD_results_5[,2])
W269D_MD_results <- cbind(W269D_MD_results,W269D_MD_results_6[,2])
names(W269D_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp269OD1vsHis_dist","Asp269OD2vsHis_dist")
summary(W269D_MD_results)
ServsHisdiffW269D <- wilcox.test(W269D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.08)$p.value
AspOD1vsHisdiffW269D <- wilcox.test(W269D_MD_results[,4],WT_MD_results[,4], alternative="greater",mu=0.04)$p.value
Ser238vsHisdiffW269D <- wilcox.test(W269D_MD_results[,5],W269D_MD_results[,2], alternative="less",mu=-0.6)$p.value
Asp269vsHisdiffW269D <- wilcox.test(W269D_MD_results[,6],W269D_MD_results[,4], alternative="greater",mu=0.125)$p.value
summary(W269D_MD_results)
sd(W269D_MD_results[,2]);sd(W269D_MD_results[,5])
sd(W269D_MD_results[,4]);sd(W269D_MD_results[,6])

ggplot() +geom_density(data=W269D_MD_results,aes(ServsHis_dist,fill="ServsHis"),
                      bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results,aes(Ser238vsHis_dist,fill="Ser238vsHis"),
               bw=0.1,alpha=0.5)+theme_minimal()+
  guides(fill=guide_legend(title="Metric"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))

ggplot()+geom_density(data=W269D_MD_results,aes(AspOD2vsHis_dist,fill="AspvsHis"),
                      bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results,aes(Asp269OD1vsHis_dist,fill="Asp269vsHis"),
               bw=0.1,alpha=0.5)+theme_minimal()+
  guides(fill=guide_legend(title="Metric"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))

# W269D total

W269DT_MD_results_1 <- read.table(file="Results/MD_W269D/ServsHis_tot.dat")
W269DT_MD_results_2 <- read.table(file="Results/MD_W269D/AspOD1vsHis_tot.dat")
W269DT_MD_results_3 <- read.table(file="Results/MD_W269D/AspOD2vsHis_tot.dat")
W269DT_MD_results_4 <- read.table(file="Results/MD_W269D/Ser238vsHis_tot.dat")
W269DT_MD_results_5 <- read.table(file="Results/MD_W269D/Asp269OD1vsHis_tot.dat")
W269DT_MD_results_6 <- read.table(file="Results/MD_W269D/Asp269OD2vsHis_tot.dat")
W269DT_MD_results_7 <- read.table(file="Results/MD_W269D/localrmsd_tot.xvg")
W269DT_MD_results_8 <- read.table(file="Results/MD_W269D/localrmsd_potential_tot.xvg")
W269DT_MD_results_9 <- read.table(file="Results/MD_W269D/rmsd_tot.xvg")

# W269DT_MD_results_10 <- read.table(file="Results/MD_W269D/rmsd_2.xvg")
# W269DT_MD_results_11 <- read.table(file="Results/MD_W269D/rmsd_3.xvg")
# W269D_MD_results_RMSD <- W269D_MD_results_9[,2]
# W269D_MD_results_RMSD <- cbind(W269D_MD_results_RMSD,W269D_MD_results_11[,2])
# W269D_MD_results_RMSD <- cbind(W269D_MD_results_RMSD,W269D_MD_results_10[,2])
# RMSD_WT <- c()
# for (i in 1:length(W269DT_MD_results_10[,2])) {
#   RMSD_WT <- c(RMSD_WT,mean(W269D_MD_results_RMSD[i,1:3]))
# }
# for (i in (length(W269DT_MD_results_10[,2])+1):length(W269DT_MD_results_11[,2])){
#   RMSD_WT <- c(RMSD_WT,mean(W269D_MD_results_RMSD[i,1:2]))
# }
# for (i in (length(W269DT_MD_results_11[,2])+1):length(W269D_MD_results_9[,2])){
#   RMSD_WT <- c(RMSD_WT,mean(W269D_MD_results_RMSD[i,1]))
# }
# RMSD_WT <- as.data.frame(RMSD_WT)
# rmsd <- cbind(W269D_MD_results_9[,1],RMSD_WT)
# write.table(rmsd, "rmsd_tot.xvg", sep="\t",append = FALSE)

summary(W269D_MD_results_1);summary(W269DT_MD_results_1)

# Global RMSD

ggplot(data=W269D_MD_results_9,aes(x=W269D_MD_results_9[,1],y=W269D_MD_results_9[,2],
                                color="D238SW269D"))+geom_line(alpha=0.7)+
  geom_line(data=W269DT_MD_results_9,
            aes(x=W269DT_MD_results_9[,1],
                y=W269DT_MD_results_9[,2],color="D238SW269D: Average"),alpha=0.7)+
  geom_line(data=T124L_MD_results_9,
            aes(x=T124L_MD_results_9[,1],
                y=T124L_MD_results_9[,2],color="D129HT124LQ137E"),alpha=0.7)+
  geom_line(data=T124LT_MD_results_9,
            aes(x=T124LT_MD_results_9[,1],
                y=T124LT_MD_results_9[,2],color="D129HT124LQ137E: Average"),alpha=0.7)+
  theme_minimal()+ggtitle("Global RMSD (least-squares fit to backbone)")+xlab("Time (ps)")+ylab("RMSD (nm)")+
  guides(color=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

# Local RMSD main

LRM <- ggplot()+
  geom_density(data=WT_MD_results_4,aes(WT_MD_results_4[,2],fill="WT"),
                bw=0.01,alpha=0.5)+
  geom_density(data=W269DT_MD_results_7,aes(W269DT_MD_results_7[,2],fill="D238SW269D: Average"),
                bw=0.01,alpha=0.5)+
  geom_density(data=W269D_MD_results_7,aes(W269D_MD_results_7[,2],fill="D238SW269D"),
                bw=0.01,alpha=0.5)+
  geom_density(data=T124L_MD_results_7,aes(T124L_MD_results_7[,2],fill="D129HT124LQ137E"),
                bw=0.01,alpha=0.5)+
  geom_density(data=T124LT_MD_results_7,aes(T124LT_MD_results_7[,2],fill="D129HT124LQ137E: Average"),
                bw=0.01,alpha=0.5)+
  theme_minimal()+xlab("RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

# Local RMSD potential

LRP <- ggplot()+
  geom_density(data=W269DT_MD_results_8,aes(W269DT_MD_results_8[,2],fill="D238SW269D: Average"),
               bw=0.01,alpha=0.5)+
  geom_density(data=W269D_MD_results_8,aes(W269D_MD_results_8[,2],fill="D238SW269D"),
               bw=0.01,alpha=0.5)+
  geom_density(data=T124L_MD_results_8,aes(T124L_MD_results_8[,2],fill="D129HT124LQ137E"),
               bw=0.01,alpha=0.5)+
  geom_density(data=T124LT_MD_results_8,aes(T124LT_MD_results_8[,2],fill="D129HT124LQ137E: Average"),
               bw=0.01,alpha=0.5)+
  theme_minimal()+xlab("RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))

figure <- ggarrange(LRM, LRP,
                    labels = c(1,2),
                    ncol = 2, nrow = 1)

annotate_figure(figure,
                top = text_grob("Local RMSD", size = 15))

# Catalytic distances

SER <- ggplot()+
  geom_density(data=WT_MD_results_1,aes(WT_MD_results_1[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269DT_MD_results_4,aes(W269DT_MD_results_4[,2],fill="D238SW269D: Average"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_4,aes(W269D_MD_results_4[,2],fill="D238SW269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_4,aes(T124L_MD_results_4[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124LT_MD_results_4,aes(T124LT_MD_results_4[,2],fill="D129HT124LQ137E: Average"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His")+xlab("Ser-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))

ASP <- ggplot()+
  geom_density(data=WT_MD_results_3,aes(WT_MD_results_3[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269DT_MD_results_5,aes(W269DT_MD_results_5[,2],fill="D238SW269D: Average"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_5,aes(W269D_MD_results_5[,2],fill="D238SW269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124L_MD_results_5,aes(T124L_MD_results_5[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=T124LT_MD_results_5,aes(T124LT_MD_results_5[,2],fill="D129HT124LQ137E: Average"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Asp-His")+xlab("Asp-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.position="none",
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))

figure <- ggarrange(SER, ASP,
                    labels = c(1,2),
                    ncol = 2, nrow = 1)

annotate_figure(figure,
                top = text_grob("Catalytic distances", size = 15))

# S278K data frame

S278K_MD_results_1 <- read.table(file="Results/MD_S278K/ServsHis.dat")
S278K_MD_results_2 <- read.table(file="Results/MD_S278K/AspOD1vsHis.dat")
S278K_MD_results_3 <- read.table(file="Results/MD_S278K/AspOD2vsHis.dat")
S278K_MD_results_4 <- read.table(file="Results/MD_S278K/Ser238vsHis.dat")
S278K_MD_results_5 <- read.table(file="Results/MD_S278K/Asp272OD1vsHis.dat")
S278K_MD_results_6 <- read.table(file="Results/MD_S278K/Asp272OD2vsHis.dat")
S278K_MD_results_7 <- read.table(file="Results/MD_S278K/localrmsd.xvg")
S278K_MD_results_8 <- read.table(file="Results/MD_S278K/localrmsd_potential.xvg")
S278K_MD_results_9 <- read.table(file="Results/MD_S278K/rmsd.xvg")
S278K_MD_results <- S278K_MD_results_1
S278K_MD_results <- cbind(S278K_MD_results,S278K_MD_results_2[,2]);S278K_MD_results <- cbind(S278K_MD_results,S278K_MD_results_3[,2])
S278K_MD_results <- cbind(S278K_MD_results,S278K_MD_results_4[,2]);S278K_MD_results <- cbind(S278K_MD_results,S278K_MD_results_5[,2])
S278K_MD_results <- cbind(S278K_MD_results,S278K_MD_results_6[,2])
names(S278K_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp272OD1vsHis_dist","Asp272OD2vsHis_dist")

ServsHisdiffS278K <- wilcox.test(S278K_MD_results[,2],WT_MD_results[,2], alternative="less",mu=-0.06)$p.value
AspOD1vsHisdiffS278K <- wilcox.test(S278K_MD_results[,4],WT_MD_results[,4], alternative="less",mu=-0.03)$p.value
Ser238vsHisdiffS278K <- wilcox.test(S278K_MD_results[,5],S278K_MD_results[,2], alternative="greater",mu=0.55)$p.value
Asp272vsHisdiffS278K <- wilcox.test(S278K_MD_results[,7],S278K_MD_results[,4], alternative="greater",mu=2.65)$p.value

# S278D data frame

S278D_MD_results_1 <- read.table(file="Results/MD_S278D/ServsHis.dat")
S278D_MD_results_2 <- read.table(file="Results/MD_S278D/AspOD1vsHis.dat")
S278D_MD_results_3 <- read.table(file="Results/MD_S278D/AspOD2vsHis.dat")
S278D_MD_results_4 <- read.table(file="Results/MD_S278D/Ser238vsHis.dat")
S278D_MD_results_5 <- read.table(file="Results/MD_S278D/Asp272OD1vsHis.dat")
S278D_MD_results_6 <- read.table(file="Results/MD_S278D/Asp272OD2vsHis.dat")
S278D_MD_results_7 <- read.table(file="Results/MD_S278D/localrmsd.xvg")
S278D_MD_results_8 <- read.table(file="Results/MD_S278D/localrmsd_potential.xvg")
S278D_MD_results_9 <- read.table(file="Results/MD_S278D/rmsd.xvg")
S278D_MD_results <- S278D_MD_results_1
S278D_MD_results <- cbind(S278D_MD_results,S278D_MD_results_2[,2]);S278D_MD_results <- cbind(S278D_MD_results,S278D_MD_results_3[,2])
S278D_MD_results <- cbind(S278D_MD_results,S278D_MD_results_4[,2]);S278D_MD_results <- cbind(S278D_MD_results,S278D_MD_results_5[,2])
S278D_MD_results <- cbind(S278D_MD_results,S278D_MD_results_6[,2])
names(S278D_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp272OD1vsHis_dist","Asp272OD2vsHis_dist")

# S278R data frame

S278R_MD_results_1 <- read.table(file="Results/MD_S278R/ServsHis.dat")
S278R_MD_results_2 <- read.table(file="Results/MD_S278R/AspOD1vsHis.dat")
S278R_MD_results_3 <- read.table(file="Results/MD_S278R/AspOD2vsHis.dat")
S278R_MD_results_4 <- read.table(file="Results/MD_S278R/Ser238vsHis.dat")
S278R_MD_results_5 <- read.table(file="Results/MD_S278R/Asp272OD1vsHis.dat")
S278R_MD_results_6 <- read.table(file="Results/MD_S278R/Asp272OD2vsHis.dat")
S278R_MD_results_7 <- read.table(file="Results/MD_S278R/localrmsd.xvg")
S278R_MD_results_8 <- read.table(file="Results/MD_S278R/localrmsd_potential.xvg")
S278R_MD_results_9 <- read.table(file="Results/MD_S278R/rmsd.xvg")
S278R_MD_results <- S278R_MD_results_1
S278R_MD_results <- cbind(S278R_MD_results,S278R_MD_results_2[,2]);S278R_MD_results <- cbind(S278R_MD_results,S278R_MD_results_3[,2])
S278R_MD_results <- cbind(S278R_MD_results,S278R_MD_results_4[,2]);S278R_MD_results <- cbind(S278R_MD_results,S278R_MD_results_5[,2])
S278R_MD_results <- cbind(S278R_MD_results,S278R_MD_results_6[,2])
names(S278R_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp272OD1vsHis_dist","Asp272OD2vsHis_dist")

ServsHisdiffS278R <- wilcox.test(S278R_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.08)$p.value
AspOD1vsHisdiffS278R <- wilcox.test(S278R_MD_results[,4],WT_MD_results[,4], alternative="greater",mu=0.16)$p.value
Ser238vsHisdiffS278R <- wilcox.test(S278R_MD_results[,5],S278R_MD_results[,2], alternative="less",mu=-0.25)$p.value
Asp272vsHisdiffS278R <- wilcox.test(S278R_MD_results[,6],S278R_MD_results[,3], alternative="greater",mu=2.5)$p.value

summary(S278R_MD_results)

# K281L data frame

K281L_MD_results_1 <- read.table(file="Results/MD_K281L/ServsHis.dat")
K281L_MD_results_2 <- read.table(file="Results/MD_K281L/AspOD1vsHis.dat")
K281L_MD_results_3 <- read.table(file="Results/MD_K281L/AspOD2vsHis.dat")
K281L_MD_results_4 <- read.table(file="Results/MD_K281L/Ser238vsHis.dat")
K281L_MD_results_5 <- read.table(file="Results/MD_K281L/Asp272OD1vsHis.dat")
K281L_MD_results_6 <- read.table(file="Results/MD_K281L/Asp272OD2vsHis.dat")
K281L_MD_results_7 <- read.table(file="Results/MD_K281L/localrmsd.xvg")
K281L_MD_results_8 <- read.table(file="Results/MD_K281L/localrmsd_potential.xvg")
K281L_MD_results_9 <- read.table(file="Results/MD_K281L/rmsd.xvg")
K281L_MD_results <- K281L_MD_results_1
K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_2[,2]);K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_3[,2])
K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_4[,2]);K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_5[,2])
K281L_MD_results <- cbind(K281L_MD_results,K281L_MD_results_6[,2])
names(K281L_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp272OD1vsHis_dist","Asp272OD2vsHis_dist")

ServsHisdiffK281L <- wilcox.test(K281L_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.12)$p.value
AspOD1vsHisdiffK281L <- wilcox.test(K281L_MD_results[,3],WT_MD_results[,4], alternative="greater",mu=0.06)$p.value
Ser238vsHisdiffK281L <- wilcox.test(K281L_MD_results[,5],K281L_MD_results[,2], alternative="less",mu=-0.75)$p.value
Asp272vsHisdiffK281L <- wilcox.test(K281L_MD_results[,6],K281L_MD_results[,3], alternative="greater",mu=2.56)$p.value
summary(K281L_MD_results)

# K281G data frame

K281G_MD_results_1 <- read.table(file="Results/MD_K281G/ServsHis.dat")
K281G_MD_results_2 <- read.table(file="Results/MD_K281G/AspOD1vsHis.dat")
K281G_MD_results_3 <- read.table(file="Results/MD_K281G/AspOD2vsHis.dat")
K281G_MD_results_4 <- read.table(file="Results/MD_K281G/Ser238vsHis.dat")
K281G_MD_results_5 <- read.table(file="Results/MD_K281G/Asp272OD1vsHis.dat")
K281G_MD_results_6 <- read.table(file="Results/MD_K281G/Asp272OD2vsHis.dat")
K281G_MD_results_7 <- read.table(file="Results/MD_K281G/localrmsd.xvg")
K281G_MD_results_8 <- read.table(file="Results/MD_K281G/localrmsd_potential.xvg")
K281G_MD_results_9 <- read.table(file="Results/MD_K281G/rmsd.xvg")
K281G_MD_results <- K281G_MD_results_1
K281G_MD_results <- cbind(K281G_MD_results,K281G_MD_results_2[,2]);K281G_MD_results <- cbind(K281G_MD_results,K281G_MD_results_3[,2])
K281G_MD_results <- cbind(K281G_MD_results,K281G_MD_results_4[,2]);K281G_MD_results <- cbind(K281G_MD_results,K281G_MD_results_5[,2])
K281G_MD_results <- cbind(K281G_MD_results,K281G_MD_results_6[,2])
names(K281G_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp272OD1vsHis_dist","Asp272OD2vsHis_dist")

# W269D_A240L

A240L_MD_results_1 <- read.table(file="Results/MD_A240L/ServsHis.dat")
A240L_MD_results_2 <- read.table(file="Results/MD_A240L/AspOD1vsHis.dat")
A240L_MD_results_3 <- read.table(file="Results/MD_A240L/AspOD2vsHis.dat")
A240L_MD_results_4 <- read.table(file="Results/MD_A240L/Ser238vsHis.dat")
A240L_MD_results_5 <- read.table(file="Results/MD_A240L/Asp269OD1vsHis.dat")
A240L_MD_results_6 <- read.table(file="Results/MD_A240L/Asp269OD2vsHis.dat")
A240L_MD_results_7 <- read.table(file="Results/MD_A240L/localrmsd.xvg")
A240L_MD_results_8 <- read.table(file="Results/MD_A240L/localrmsd_potential.xvg")
A240L_MD_results_9 <- read.table(file="Results/MD_A240L/rmsd.xvg")
A240L_MD_results <- A240L_MD_results_1
A240L_MD_results <- cbind(A240L_MD_results,A240L_MD_results_2[,2]);A240L_MD_results <- cbind(A240L_MD_results,A240L_MD_results_3[,2])
A240L_MD_results <- cbind(A240L_MD_results,A240L_MD_results_4[,2]);A240L_MD_results <- cbind(A240L_MD_results,A240L_MD_results_5[,2])
A240L_MD_results <- cbind(A240L_MD_results,A240L_MD_results_6[,2])
names(A240L_MD_results) <-c("step","ServsHis_dist","AspOD1vsHis_dist","AspOD2vsHis_dist","Ser238vsHis_dist","Asp269OD1vsHis_dist","Asp269OD2vsHis_dist")
summary(A240L_MD_results)
ServsHisdiffW269D <- wilcox.test(W269D_MD_results[,2],WT_MD_results[,2], alternative="greater",mu=0.08)$p.value
AspOD1vsHisdiffW269D <- wilcox.test(W269D_MD_results[,4],WT_MD_results[,4], alternative="greater",mu=0.04)$p.value
Ser238vsHisdiffW269D <- wilcox.test(W269D_MD_results[,5],W269D_MD_results[,2], alternative="less",mu=-0.6)$p.value
Asp269vsHisdiffW269D <- wilcox.test(W269D_MD_results[,6],W269D_MD_results[,4], alternative="greater",mu=0.125)$p.value
summary(W269D_MD_results)
sd(W269D_MD_results[,2]);sd(W269D_MD_results[,5])
sd(W269D_MD_results[,4]);sd(W269D_MD_results[,6])

ggplot()+geom_density(data=W269D_MD_results,aes(ServsHis_dist,fill="ServsHis"),
                      bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results,aes(Ser238vsHis_dist,fill="Ser238vsHis"),
               bw=0.1,alpha=0.5)+theme_minimal()+
  guides(fill=guide_legend(title="Metric"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))

ggplot()+geom_density(data=W269D_MD_results,aes(AspOD2vsHis_dist,fill="AspvsHis"),
                      bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results,aes(Asp269OD2vsHis_dist,fill="Asp269vsHis"),
               bw=0.1,alpha=0.5)+theme_minimal()+
  guides(fill=guide_legend(title="Metric"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#56B2E9"))

# W269D_D239S_Complex

Complex_MD_results_1 <- read.table(file="Results/MD_Complex/ServsHis.dat")
Complex_MD_results_2 <- read.table(file="Results/MD_Complex/AspOD1vsHis.dat")
Complex_MD_results_3 <- read.table(file="Results/MD_Complex/AspOD2vsHis.dat")
Complex_MD_results_4 <- read.table(file="Results/MD_Complex/ServsC1.dat")
Complex_MD_results_5 <- read.table(file="Results/MD_Complex/ServsC8.dat")
Complex_MD_results_6 <- read.table(file="Results/MD_Complex/ServsC11.dat")
Complex_MD_results_7 <- read.table(file="Results/MD_Complex/ServsC23.dat")
Complex_MD_results_8 <- read.table(file="Results/MD_Complex/ServsC30.dat")
Complex_MD_results_9 <- read.table(file="Results/MD_Complex/ServsC33.dat")
Complex_MD_results_10 <- read.table(file="Results/MD_Complex/localrmsd.xvg")
Complex_MD_results_11 <- read.table(file="Results/MD_Complex/localrmsd_potential.xvg")
Complex_MD_results_12 <- read.table(file="Results/MD_Complex/rmsd.xvg")

Complex_MD_results <- Complex_MD_results_4
Complex_MD_results <- cbind(Complex_MD_results,Complex_MD_results_5[,2]);Complex_MD_results <- cbind(Complex_MD_results,Complex_MD_results_6[,2])
Complex_MD_results <- cbind(Complex_MD_results,Complex_MD_results_7[,2]);Complex_MD_results <- cbind(Complex_MD_results,Complex_MD_results_8[,2])
Complex_MD_results <- cbind(Complex_MD_results,Complex_MD_results_9[,2])

Complex_MD_results$min <- do.call(pmin, Complex_MD_results[,2:7])

summary(Complex_MD_results$min);sd(Complex_MD_results$min)

# Total

ComplexT_MD_results_1 <- read.table(file="Results/MD_Complex/ServsHis_tot.dat")
ComplexT_MD_results_2 <- read.table(file="Results/MD_Complex/AspOD1vsHis_tot.dat")
ComplexT_MD_results_3 <- read.table(file="Results/MD_Complex/AspOD2vsHis_tot.dat")
ComplexT_MD_results_4 <- read.table(file="Results/MD_Complex/ServsC1_tot.dat")
ComplexT_MD_results_5 <- read.table(file="Results/MD_Complex/ServsC8_tot.dat")
ComplexT_MD_results_6 <- read.table(file="Results/MD_Complex/ServsC11_tot.dat")
ComplexT_MD_results_7 <- read.table(file="Results/MD_Complex/ServsC23_tot.dat")
ComplexT_MD_results_8 <- read.table(file="Results/MD_Complex/ServsC30_tot.dat")
ComplexT_MD_results_9 <- read.table(file="Results/MD_Complex/ServsC33_tot.dat")
ComplexT_MD_results_10 <- read.table(file="Results/MD_Complex/localrmsd_tot.xvg")
ComplexT_MD_results_11 <- read.table(file="Results/MD_Complex/localrmsd_potential_tot.xvg")
ComplexT_MD_results_12 <- read.table(file="Results/MD_Complex/rmsd_2.xvg")
ComplexT_MD_results_13 <- read.table(file="Results/MD_Complex/rmsd_3.xvg")
ComplexT_MD_results_14 <- read.table(file="Results/MD_Complex/rmsd_tot.xvg")
ComplexT_MD_results_RMSD  <- ComplexT_MD_results_13[,2]
ComplexT_MD_results_RMSD <- cbind(ComplexT_MD_results_RMSD,ComplexT_MD_results_12[,2])
ComplexT_MD_results_RMSD <- cbind(ComplexT_MD_results_RMSD,Complex_MD_results_12[,2])
RMSD_WT <- c()
for (i in 1:length(Complex_MD_results_12[,2])) {
   RMSD_WT <- c(RMSD_WT,mean(ComplexT_MD_results_RMSD[i,1:3]))
}
for (i in (length(Complex_MD_results_12[,2])+1):length(ComplexT_MD_results_12[,2])){
   RMSD_WT <- c(RMSD_WT,mean(ComplexT_MD_results_RMSD[i,1:2]))
}
for (i in (length(ComplexT_MD_results_12[,2])+1):length(ComplexT_MD_results_13[,2])){
   RMSD_WT <- c(RMSD_WT,mean(ComplexT_MD_results_RMSD[i,1])) 
}
RMSD_WT <- as.data.frame(RMSD_WT)
rmsd <- cbind(ComplexT_MD_results_13[,1],RMSD_WT)
write.table(rmsd, "rmsd_tot.xvg", sep="\t",append = FALSE)


ComplexT_MD_results <- ComplexT_MD_results_4
ComplexT_MD_results <- cbind(ComplexT_MD_results,ComplexT_MD_results_5[,2]);ComplexT_MD_results <- cbind(ComplexT_MD_results,ComplexT_MD_results_6[,2])
ComplexT_MD_results <- cbind(ComplexT_MD_results,ComplexT_MD_results_7[,2]);ComplexT_MD_results <- cbind(ComplexT_MD_results,ComplexT_MD_results_8[,2])
ComplexT_MD_results <- cbind(ComplexT_MD_results,ComplexT_MD_results_9[,2])

ComplexT_MD_results$min <- do.call(pmin, ComplexT_MD_results[,2:7])

ggplot()+geom_density(data=Complex_MD_results,aes(min),
                      bw=0.1,alpha=0.5,fill="#56B2E9")+
  theme_minimal()+ggtitle("Ser-Substrate for the potential catalytic triad")+xlab("Ser-Substrate distance (Ang)")+
  theme(plot.title = element_text(hjust = 0.5))

# Catalytic distances in complex

ggplot()+
  geom_density(data=Complex_MD_results_1,aes(Complex_MD_results_1[,2],fill="D238SW269D_with_ligand"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_4,aes(W269D_MD_results_4[,2],fill="D238SW269D"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His distance for the potential catalytic triad")+xlab("Ser-His distance (Ang)")+
  guides(fill=guide_legend(title="MD simulation"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

ggplot()+
  geom_density(data=Complex_MD_results_2,aes(Complex_MD_results_2[,2],fill="D238SW269D_with_ligand"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_5,aes(W269D_MD_results_5[,2],fill="D238SW269D"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("His-Asp distance for the potential catalytic triad")+xlab("His-Asp distance (Ang)")+
  guides(fill=guide_legend(title="MD simulation"))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00","#082D8E","#ff001a"))

# RMSD in complex

ggplot()+geom_density(data=ComplexT_MD_results_11,aes(ComplexT_MD_results_11[,2]),
                      bw=0.01,alpha=0.5,fill="#56B2E9")+
  theme_minimal()+ggtitle("Local RMSD for the potential catalytic triad")+xlab("Local RMSD (nm)")+
  theme(plot.title = element_text(hjust = 0.5))

# RMSD of the simulation

ggplot(data=WT_MD_results_8,aes(x=WT_MD_results_8[,1],y=WT_MD_results_8[,2],
                                color="WT"))+geom_line(alpha=0.7)+
  geom_line(data=K281G_MD_results_9,
            aes(x=K281G_MD_results_9[,1],
                y=K281G_MD_results_9[,2],color="D238SL272DK281G"),alpha=0.7)+
  geom_line(data=K281L_MD_results_9,
            aes(x=K281L_MD_results_9[,1],
                y=K281L_MD_results_9[,2],color="D238SL272DK281L"),alpha=0.7)+
  geom_line(data=S278K_MD_results_9,
            aes(x=S278K_MD_results_9[,1],
                y=S278K_MD_results_9[,2],color="D238SL272DS278KK281L"),alpha=0.7)+
  geom_line(data=S278R_MD_results_9,
            aes(x=S278R_MD_results_9[,1],
                y=S278R_MD_results_9[,2],color="D238SL272DS278RK281L"),alpha=0.7)+
  geom_line(data=S278D_MD_results_9,
            aes(x=S278D_MD_results_9[,1],
                y=S278D_MD_results_9[,2],color="D238SL272DS278D"),alpha=0.7)+
  theme_minimal()+ggtitle("Global RMSD (least-squares fit to backbone)")+xlab("Time (ps)")+ylab("RMSD (nm)")+
  guides(color=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

ggplot(data=WT_MD_results_8,aes(x=WT_MD_results_8[,1],y=WT_MD_results_8[,2],
                                color="WT"))+geom_line(alpha=0.7)+
  geom_line(data=L272D_MD_results_9,
            aes(x=L272D_MD_results_9[,1],
                y=L272D_MD_results_9[,2],color="D238SL272D"),alpha=0.7)+
  geom_line(data=W269D_MD_results_9,
            aes(x=W269D_MD_results_9[,1],
                y=W269D_MD_results_9[,2],color="D238SW269D"),alpha=0.7)+
  geom_line(data=A240L_MD_results_9,
            aes(x=A240L_MD_results_9[,1],
                y=A240L_MD_results_9[,2],color="D238SW269DA240L"),alpha=0.7)+
  theme_minimal()+ggtitle("Global RMSD (least-squares fit to backbone)")+xlab("Time (ps)")+ylab("RMSD (nm)")+
  guides(color=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

# LOCAL RMSD MAIN ACTIVE SITE

ggplot()+geom_density(data=WT_MD_results_4,aes(WT_MD_results_4[,2],fill="WT"),
                      bw=0.01,alpha=0.5)+
  geom_density(data=K281L_MD_results_7,aes(K281L_MD_results_7[,2],fill="D238SL272DK281L"),
               bw=0.01,alpha=0.5)+
  geom_density(data=K281G_MD_results_7,aes(K281G_MD_results_7[,2],fill="D238SL272DK281G"),
               bw=0.01,alpha=0.5)+
  geom_density(data=S278K_MD_results_7,aes(S278K_MD_results_7[,2],fill="D238SL272DS278KK281L"),
               bw=0.01,alpha=0.5)+
  geom_density(data=S278R_MD_results_7,aes(S278R_MD_results_7[,2],fill="D238SL272DS278RK281L"),
               bw=0.01,alpha=0.5)+
  geom_density(data=S278D_MD_results_7,aes(S278D_MD_results_7[,2],fill="D238SL272DS278D"),
               bw=0.01,alpha=0.5)+
  theme_minimal()+ggtitle("Local RMSD for the main catalytic triad")+xlab("Local RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

ggplot()+geom_density(data=WT_MD_results_4,aes(WT_MD_results_4[,2],fill="WT"),
                      bw=0.01,alpha=0.5)+
  geom_density(data=L272D_MD_results_7,aes(L272D_MD_results_7[,2],fill="D238L272D"),
               bw=0.01,alpha=0.5)+
  geom_density(data=W269D_MD_results_7,aes(W269D_MD_results_7[,2],fill="D238W269D"),
               bw=0.01,alpha=0.5)+
  geom_density(data=A240L_MD_results_7,aes(A240L_MD_results_7[,2],fill="D238SW269DA240L"),
               bw=0.01,alpha=0.5)+
  theme_minimal()+ggtitle("Local RMSD for the main catalytic triad")+xlab("Local RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00"))

# LOCAL RMSD POTENTIAl ACTIVE SITE

ggplot()+
  geom_density(data=K281L_MD_results_8,aes(K281L_MD_results_8[,2],fill="D238SL272DK281L"),
               bw=0.01,alpha=0.5)+
  geom_density(data=K281G_MD_results_8,aes(K281G_MD_results_8[,2],fill="D238SL272DK281G"),
               bw=0.01,alpha=0.5)+
  geom_density(data=S278K_MD_results_8,aes(S278K_MD_results_8[,2],fill="D238SL272DS278KK281L"),
               bw=0.01,alpha=0.5)+
  geom_density(data=S278R_MD_results_8,aes(S278R_MD_results_8[,2],fill="D238SL272DS278RK281L"),
               bw=0.01,alpha=0.5)+
  geom_density(data=S278D_MD_results_8,aes(S278D_MD_results_8[,2],fill="D238SL272DS278D"),
               bw=0.01,alpha=0.5)+
  theme_minimal()+ggtitle("Local RMSD for the potential catalytic triad")+xlab("Local RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#ffdab9","#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

ggplot()+
  geom_density(data=L272D_MD_results_8,aes(L272D_MD_results_8[,2],fill="D238L272D"),
               bw=0.01,alpha=0.5)+
  geom_density(data=W269D_MD_results_8,aes(W269D_MD_results_8[,2],fill="D238W269D"),
               bw=0.01,alpha=0.5)+
  geom_density(data=A240L_MD_results_8,aes(A240L_MD_results_8[,2],fill="D238SW269DA240L"),
               bw=0.01,alpha=0.5)+
  theme_minimal()+ggtitle("Local RMSD for the potential catalytic triad")+xlab("Local RMSD (nm)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a","#ffdab9","#a3ff00"))

# Catalytic distances MAIN

ggplot()+
  geom_density(data=WT_MD_results_1,aes(WT_MD_results_1[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=L272D_MD_results_1,aes(L272D_MD_results_1[,2],fill="D238L272D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281L_MD_results_1,aes(K281L_MD_results_1[,2],fill="D238SL272DK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281G_MD_results_1,aes(K281G_MD_results_1[,2],fill="D238SL272DK281G"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278K_MD_results_1,aes(S278K_MD_results_1[,2],fill="D238SL272DS278KK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278R_MD_results_1,aes(S278R_MD_results_1[,2],fill="D238SL272DS278RK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278D_MD_results_1,aes(S278D_MD_results_1[,2],fill="D238SL272DS278D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_1,aes(W269D_MD_results_1[,2],fill="D238W269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=A240L_MD_results_1,aes(A240L_MD_results_1[,2],fill="D238SW269DA240L"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His distance for the main catalytic triad")+xlab("Ser-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#ffdab9","#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

ggplot()+
  geom_density(data=WT_MD_results_3,aes(WT_MD_results_3[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=L272D_MD_results_3,aes(L272D_MD_results_3[,2],fill="D238L272D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281L_MD_results_2,aes(K281L_MD_results_2[,2],fill="D238SL272DK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281G_MD_results_3,aes(K281G_MD_results_3[,2],fill="D238SL272DK281G"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278K_MD_results_3,aes(S278K_MD_results_3[,2],fill="D238SL272DS278KK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278R_MD_results_3,aes(S278R_MD_results_3[,2],fill="D238SL272DS278RK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278D_MD_results_3,aes(S278D_MD_results_3[,2],fill="D238SL272DS278D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_2,aes(W269D_MD_results_2[,2],fill="D238W269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=A240L_MD_results_3,aes(A240L_MD_results_3[,2],fill="D238SW269DA240L"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Asp-His distance for the main catalytic triad")+xlab("Asp-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#ffdab9","#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

# Catalytic distances POTENTIAL

ggplot()+
  geom_density(data=WT_MD_results_1,aes(WT_MD_results_1[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=L272D_MD_results_4,aes(L272D_MD_results_4[,2],fill="D238L272D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_4,aes(W269D_MD_results_4[,2],fill="D238W269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=A240L_MD_results_4,aes(A240L_MD_results_4[,2],fill="D238SW269DA240L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281L_MD_results_4,aes(K281L_MD_results_4[,2],fill="D238SL272DK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281G_MD_results_4,aes(K281G_MD_results_4[,2],fill="D238SL272DK281G"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278K_MD_results_4,aes(S278K_MD_results_4[,2],fill="D238SL272DS278KK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278R_MD_results_4,aes(S278R_MD_results_4[,2],fill="D238SL272DS278RK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278D_MD_results_4,aes(S278D_MD_results_4[,2],fill="D238SL272DS278D"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His distance for the potential catalytic triad")+xlab("Ser-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#ffdab9","#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

ggplot()+
  geom_density(data=K281L_MD_results_4,aes(K281L_MD_results_4[,2],fill="D238SL272DK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281G_MD_results_4,aes(K281G_MD_results_4[,2],fill="D238SL272DK281G"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278K_MD_results_4,aes(S278K_MD_results_4[,2],fill="D238SL272DS278KK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278R_MD_results_4,aes(S278R_MD_results_4[,2],fill="D238SL272DS278RK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278D_MD_results_4,aes(S278D_MD_results_4[,2],fill="D238SL272DS278D"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Ser-His distance for the potential catalytic triad")+xlab("Ser-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#ffdab9","#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))

ggplot()+
  geom_density(data=WT_MD_results_3,aes(WT_MD_results_3[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=L272D_MD_results_5,aes(L272D_MD_results_5[,2],fill="D238L272D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281L_MD_results_5,aes(K281L_MD_results_5[,2],fill="D238SL272DK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=K281G_MD_results_5,aes(K281G_MD_results_5[,2],fill="D238SL272DK281G"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278K_MD_results_6,aes(S278K_MD_results_6[,2],fill="D238SL272DS278KK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278R_MD_results_5,aes(S278R_MD_results_5[,2],fill="D238SL272DS278RK281L"),
               bw=0.1,alpha=0.5)+
  geom_density(data=S278D_MD_results_5,aes(S278D_MD_results_5[,2],fill="D238SL272DS278D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=W269D_MD_results_5,aes(W269D_MD_results_5[,2],fill="D238W269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=A240L_MD_results_6,aes(A240L_MD_results_6[,2],fill="D238SW269DA240L"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Asp-His distance for the potential catalytic triad")+xlab("Asp-His distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#ffdab9","#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))+
  scale_color_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                              ,"#ffdab9","#a3ff00","#afeeee",
                              "#9910f9","#136d14","#f66909"))


summary(WT_MD_results);sd(WT_MD_results$ServsHis_dist);
sd(WT_MD_results$AspOD1vsHis_dist);sd(WT_MD_results$AspOD2vsHis_dist)
summary(L272D_MD_results)
sd(L272D_MD_results$Ser238vsHis_dist);sd(L272D_MD_results$Asp272OD1vsHis_dist)
sd(L272D_MD_results$Asp272OD2vsHis_dist)
summary(K281L_MD_results);sd(K281L_MD_results$Ser238vsHis_dist);sd(K281L_MD_results$Asp272OD1vsHis_dist)
sd(K281L_MD_results$Asp272OD2vsHis_dist)
summary(K281G_MD_results);sd(K281G_MD_results$Ser238vsHis_dist);sd(K281G_MD_results$Asp272OD1vsHis_dist)
sd(K281G_MD_results$Asp272OD2vsHis_dist)
summary(S278D_MD_results);sd(S278D_MD_results$Ser238vsHis_dist);sd(S278D_MD_results$Asp272OD1vsHis_dist)
sd(S278D_MD_results$Asp272OD2vsHis_dist)
summary(S278K_MD_results);sd(S278K_MD_results$Ser238vsHis_dist);sd(S278K_MD_results$Asp272OD1vsHis_dist)
sd(S278K_MD_results$Asp272OD2vsHis_dist)
summary(S278R_MD_results);sd(S278R_MD_results$Ser238vsHis_dist);sd(S278R_MD_results$Asp272OD1vsHis_dist)
sd(S278R_MD_results$Asp272OD2vsHis_dist)
summary(W269D_MD_results);sd(W269D_MD_results$Ser238vsHis_dist);sd(W269D_MD_results$Asp269OD1vsHis_dist)
sd(W269D_MD_results$Asp269OD2vsHis_dist)
summary(A240L_MD_results);sd(A240L_MD_results$Ser238vsHis_dist);sd(A240L_MD_results$Asp269OD1vsHis_dist)
sd(A240L_MD_results$Asp269OD2vsHis_dist)

skewness(WT_MD_results$ServsHis_dist);
skewness(WT_MD_results$AspOD1vsHis_dist);skewness(WT_MD_results$AspOD2vsHis_dist)
skewness(L272D_MD_results$Ser238vsHis_dist);skewness(L272D_MD_results$Asp272OD1vsHis_dist)
skewness(L272D_MD_results$Asp272OD2vsHis_dist)
skewness(K281L_MD_results$Ser238vsHis_dist);skewness(K281L_MD_results$Asp272OD1vsHis_dist)
skewness(K281L_MD_results$Asp272OD2vsHis_dist)
skewness(K281G_MD_results$Ser238vsHis_dist);skewness(K281G_MD_results$Asp272OD1vsHis_dist)
skewness(K281G_MD_results$Asp272OD2vsHis_dist)
skewness(S278D_MD_results$Ser238vsHis_dist);skewness(S278D_MD_results$Asp272OD1vsHis_dist)
skewness(S278D_MD_results$Asp272OD2vsHis_dist)
skewness(S278K_MD_results$Ser238vsHis_dist);skewness(S278K_MD_results$Asp272OD1vsHis_dist)
skewness(S278K_MD_results$Asp272OD2vsHis_dist)
skewness(S278R_MD_results$Ser238vsHis_dist);skewness(S278R_MD_results$Asp272OD1vsHis_dist)
skewness(S278R_MD_results$Asp272OD2vsHis_dist)
skewness(W269D_MD_results$Ser238vsHis_dist);skewness(W269D_MD_results$Asp269OD1vsHis_dist)
skewness(W269D_MD_results$Asp269OD2vsHis_dist)
skewness(A240L_MD_results$Ser238vsHis_dist);skewness(A240L_MD_results$Asp269OD1vsHis_dist)
skewness(A240L_MD_results$Asp269OD2vsHis_dist)

p <- ggplot(WT_MD_results,aes(WT_MD_results[,3],WT_MD_results[,4]))+geom_point(color="black",size=1,alpha=0.5)+
  geom_density_2d(color="blue")+theme_minimal()

ggMarginal(p,type = "density",color="blue",fill="blue",alpha=0.5)

# Inference of MD simulations

var.test(c(K281L_MD_results_8[,2],K281G_MD_results_8[,2],
              S278K_MD_results_8[,2],S278R_MD_results_8[,2],
              S278D_MD_results_8[,2],L272D_MD_results_8[,2]),
         c(A240L_MD_results_8[,2],W269D_MD_results_8[,2]))

t.test(c(K281L_MD_results_8[,2],K281G_MD_results_8[,2],
           S278K_MD_results_8[,2],S278R_MD_results_8[,2],
           S278D_MD_results_8[,2],L272D_MD_results_8[,2]),
         c(A240L_MD_results_8[,2],W269D_MD_results_8[,2]),var.equal = FALSE,alternative = "greater")

sd(c(K281L_MD_results_8[,2],K281G_MD_results_8[,2],
  S278K_MD_results_8[,2],S278R_MD_results_8[,2],
  S278D_MD_results_8[,2],L272D_MD_results_8[,2]))
sd(c(A240L_MD_results_8[,2],W269D_MD_results_8[,2]))

var.test(WT_MD_results[,2],A240L_MD_results_4[,2])

t.test(A240L_MD_results_4[,2],WT_MD_results[,2]
       ,var.equal = FALSE,alternative = "less")

var.test(WT_MD_results[,4],A240L_MD_results_6[,2])

t.test(A240L_MD_results_6[,2],WT_MD_results[,4]
       ,var.equal = FALSE,alternative = "greater")
