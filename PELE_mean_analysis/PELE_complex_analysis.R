# W269D_D238S_Complex

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

ComplexMin <- do.call(pmin, Complex_MD_results[,2:7])

Complex_MD_results$min <- do.call(pmin, Complex_MD_results[,2:7])

Complex_Asp <- Complex_MD_results_2
Complex_Asp <- cbind(Complex_Asp,Complex_MD_results_3[,2])

Complex_Asp$min <- do.call(pmin, Complex_Asp[,2:3])

summary(Complex_MD_results$min);sd(Complex_MD_results$min)

# W269D_D238S_Complex_2

Complex_MD_results2_1 <- read.table(file="Results/MD_Complex2/ServsHis.dat")
Complex_MD_results2_2 <- read.table(file="Results/MD_Complex2/AspOD1vsHis.dat")
Complex_MD_results2_3 <- read.table(file="Results/MD_Complex2/AspOD2vsHis.dat")
Complex_MD_results2_4 <- read.table(file="Results/MD_Complex2/ServsC1.dat")
Complex_MD_results2_5 <- read.table(file="Results/MD_Complex2/ServsC8.dat")
Complex_MD_results2_6 <- read.table(file="Results/MD_Complex2/ServsC11.dat")
Complex_MD_results2_7 <- read.table(file="Results/MD_Complex2/ServsC23.dat")
Complex_MD_results2_8 <- read.table(file="Results/MD_Complex2/ServsC30.dat")
Complex_MD_results2_9 <- read.table(file="Results/MD_Complex2/ServsC33.dat")
Complex_MD_results2_10 <- read.table(file="Results/MD_Complex2/localrmsd.xvg")
Complex_MD_results2_11 <- read.table(file="Results/MD_Complex2/localrmsd_potential.xvg")
Complex_MD_results2_12 <- read.table(file="Results/MD_Complex2/rmsd_2.xvg")

Complex_MD_results2 <- Complex_MD_results2_4
Complex_MD_results2 <- cbind(Complex_MD_results2,Complex_MD_results2_5[,2]);Complex_MD_results2 <- cbind(Complex_MD_results2,Complex_MD_results2_6[,2])
Complex_MD_results2 <- cbind(Complex_MD_results2,Complex_MD_results2_7[,2]);Complex_MD_results2 <- cbind(Complex_MD_results2,Complex_MD_results2_8[,2])
Complex_MD_results2 <- cbind(Complex_MD_results2,Complex_MD_results2_9[,2])

Complex_MD_results2$min <- do.call(pmin, Complex_MD_results2[,2:7])

Complex2_Asp <- Complex_MD_results2_2
Complex2_Asp <- cbind(Complex2_Asp,Complex_MD_results2_3[,2])

Complex2_Asp$min <- do.call(pmin, Complex2_Asp[,2:3])

summary(Complex_MD_results2$min);sd(Complex_MD_results2$min)

# W269D_D238S_Complex_3

Complex_MD_results3_1 <- read.table(file="Results/MD_Complex3/ServsHis.dat")
Complex_MD_results3_2 <- read.table(file="Results/MD_Complex3/AspOD1vsHis.dat")
Complex_MD_results3_3 <- read.table(file="Results/MD_Complex3/AspOD2vsHis.dat")
Complex_MD_results3_4 <- read.table(file="Results/MD_Complex3/ServsC1.dat")
Complex_MD_results3_5 <- read.table(file="Results/MD_Complex3/ServsC8.dat")
Complex_MD_results3_6 <- read.table(file="Results/MD_Complex3/ServsC11.dat")
Complex_MD_results3_7 <- read.table(file="Results/MD_Complex3/ServsC23.dat")
Complex_MD_results3_8 <- read.table(file="Results/MD_Complex3/ServsC30.dat")
Complex_MD_results3_9 <- read.table(file="Results/MD_Complex3/ServsC33.dat")
Complex_MD_results3_10 <- read.table(file="Results/MD_Complex3/localrmsd.xvg")
Complex_MD_results3_11 <- read.table(file="Results/MD_Complex3/localrmsd_potential.xvg")
Complex_MD_results3_12 <- read.table(file="Results/MD_Complex3/rmsd_3.xvg")

Complex_MD_results3 <- Complex_MD_results3_4
Complex_MD_results3 <- cbind(Complex_MD_results3,Complex_MD_results3_5[,2]);Complex_MD_results3 <- cbind(Complex_MD_results3,Complex_MD_results3_6[,2])
Complex_MD_results3 <- cbind(Complex_MD_results3,Complex_MD_results3_7[,2]);Complex_MD_results3 <- cbind(Complex_MD_results3,Complex_MD_results3_8[,2])
Complex_MD_results3 <- cbind(Complex_MD_results3,Complex_MD_results3_9[,2])

ComplexMin3 <- do.call(pmin, Complex_MD_results3[,2:7])

Complex_MD_results3$min <- do.call(pmin, Complex_MD_results3[,2:7])

Complex3_Asp <- Complex_MD_results3_2
Complex3_Asp <- cbind(Complex3_Asp,Complex_MD_results3_3[,2])

Complex3_Asp$min <- do.call(pmin, Complex3_Asp[,2:3])

summary(c(Complex_MD_results3$min,Complex_MD_results$min));sd(c(Complex_MD_results3$min,Complex_MD_results$min))

summary(c(Complex_MD_results3_1[,2],Complex_MD_results_1[,2]));sd(c(Complex_MD_results3_1[,2],Complex_MD_results_1[,2]))
skewness(c(Complex_MD_results3_1[,2],Complex_MD_results_1[,2]))
summary(Complex_MD_results2_1[,2]);sd(Complex_MD_results2_1[,2])
skewness(Complex_MD_results2_1[,2])

var.test(c(Complex_MD_results3_1[,2],Complex_MD_results_1[,2]),
         Complex_MD_results2_1[,2])
t.test(c(Complex_MD_results3_1[,2],Complex_MD_results_1[,2]),
       Complex_MD_results2_1[,2],
       var.equal = FALSE,alternative = "less")

summary(c(Complex_Asp$min,Complex3_Asp$min));sd(c(Complex_Asp$min,Complex3_Asp$min))
skewness(c(Complex_Asp$min,Complex3_Asp$min))
summary(ComplexT_Asp$min);sd(ComplexT_Asp$min)
skewness(W269DT_MD_results_4[,2])

var.test(c(Complex_Asp$min,Complex3_Asp$min),
         Complex2_Asp$min)
t.test(c(Complex_Asp$min,Complex3_Asp$min),
       Complex2_Asp$min,
       var.equal = FALSE,alternative = "less")

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

ComplexT_Asp <- ComplexT_MD_results_2
ComplexT_Asp <- cbind(ComplexT_Asp,ComplexT_MD_results_3[,2])

ComplexT_Asp$min <- do.call(pmin, ComplexT_Asp[,2:3])

# D129HT124LQ137E

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
ComplexL_MD_results_12 <- read.table(file="Results2/MD_Complex/rmsd.xvg")

ComplexL_MD_results <- ComplexL_MD_results_4
ComplexL_MD_results <- cbind(ComplexL_MD_results,ComplexL_MD_results_5[,2]);ComplexL_MD_results <- cbind(ComplexL_MD_results,ComplexL_MD_results_6[,2])
ComplexL_MD_results <- cbind(ComplexL_MD_results,ComplexL_MD_results_7[,2]);ComplexL_MD_results <- cbind(ComplexL_MD_results,ComplexL_MD_results_8[,2])
ComplexL_MD_results <- cbind(ComplexL_MD_results,ComplexL_MD_results_9[,2])

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
ComplexTL_MD_results_12 <- read.table(file="Results2/MD_Complex/rmsd_2.xvg")
ComplexTL_MD_results_13 <- read.table(file="Results2/MD_Complex/rmsd_3.xvg")
ComplexTL_MD_results_14 <- read.table(file="Results2/MD_Complex/rmsd_tot.xvg")
ComplexTL_MD_results_RMSD  <- ComplexTL_MD_results_13[,2]
ComplexTL_MD_results_RMSD <- cbind(ComplexTL_MD_results_RMSD,ComplexTL_MD_results_12[,2])
ComplexTL_MD_results_RMSD <- cbind(ComplexTL_MD_results_RMSD,ComplexL_MD_results_12[,2])
RMSD_WT <- c()
for (i in 1:length(ComplexL_MD_results_12[,2])) {
  RMSD_WT <- c(RMSD_WT,mean(ComplexTL_MD_results_RMSD[i,1:3]))
}
for (i in (length(ComplexL_MD_results_12[,2])+1):length(ComplexTL_MD_results_12[,2])){
  RMSD_WT <- c(RMSD_WT,mean(ComplexTL_MD_results_RMSD[i,1:2]))
}
for (i in (length(ComplexTL_MD_results_12[,2])+1):length(ComplexTL_MD_results_13[,2])){
  RMSD_WT <- c(RMSD_WT,mean(ComplexTL_MD_results_RMSD[i,1])) 
}
RMSD_WT <- as.data.frame(RMSD_WT)
rmsd <- cbind(ComplexTL_MD_results_13[,1],RMSD_WT)
write.table(rmsd, "rmsd_tot.xvg", sep="\t",append = FALSE)

ComplexTL_MD_results <- ComplexTL_MD_results_4
ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_5[,2]);ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_6[,2])
ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_7[,2]);ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_8[,2])
ComplexTL_MD_results <- cbind(ComplexTL_MD_results,ComplexTL_MD_results_9[,2])

ComplexTL_MD_results$min <- do.call(pmin, ComplexTL_MD_results[,2:7])

ggplot(data=WT_MD_results_8,aes(x=WT_MD_results_8[,1],y=WT_MD_results_8[,2],
                                color="WT"))+geom_line(alpha=0.7)+
  geom_line(data=Complex_MD_results_12,
            aes(x=Complex_MD_results_12[,1],
                y=Complex_MD_results_12[,2],color="D238SW269D"),alpha=0.7)+
  geom_line(data=ComplexT_MD_results_14,
            aes(x=ComplexT_MD_results_14[,1],
                y=ComplexT_MD_results_14[,2],color="D238SW269D: Average"),alpha=0.7)+
  geom_line(data=ComplexL_MD_results_12,
            aes(x=ComplexL_MD_results_12[,1],
                y=ComplexL_MD_results_12[,2],color="D129HT124LQ137E"),alpha=0.7)+
  geom_line(data=ComplexTL_MD_results_14,
            aes(x=ComplexTL_MD_results_14[,1],
                y=ComplexTL_MD_results_14[,2],color="D129HT124LQ137E: Average"),alpha=0.7)+
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
  geom_density(data=ComplexT_MD_results_10,aes(ComplexT_MD_results_10[,2],fill="D238SW269D: Average"),
               bw=0.01,alpha=0.5)+
  geom_density(data=Complex_MD_results_10,aes(Complex_MD_results_10[,2],fill="D238SW269D"),
               bw=0.01,alpha=0.5)+
  geom_density(data=ComplexL_MD_results_10,aes(ComplexL_MD_results_10[,2],fill="D129HT124LQ137E"),
               bw=0.01,alpha=0.5)+
  geom_density(data=ComplexTL_MD_results_10,aes(ComplexTL_MD_results_10[,2],fill="D129HT124LQ137E: Average"),
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
  geom_density(data=ComplexT_MD_results_11,aes(ComplexT_MD_results_11[,2],fill="D238SW269D: Average"),
               bw=0.01,alpha=0.5)+
  geom_density(data=Complex_MD_results_11,aes(Complex_MD_results_11[,2],fill="D238SW269D"),
               bw=0.01,alpha=0.5)+
  geom_density(data=ComplexL_MD_results_11,aes(ComplexL_MD_results_11[,2],fill="D129HT124LQ137E"),
               bw=0.01,alpha=0.5)+
  geom_density(data=ComplexTL_MD_results_11,aes(ComplexTL_MD_results_11[,2],fill="D129HT124LQ137E: Average"),
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
  geom_density(data=ComplexT_MD_results_1,aes(ComplexT_MD_results_1[,2],fill="D238SW269D: Average"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Complex_MD_results_1,aes(Complex_MD_results_1[,2],fill="D238SW269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=ComplexL_MD_results_1,aes(ComplexL_MD_results_1[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=ComplexTL_MD_results_1,aes(ComplexTL_MD_results_1[,2],fill="D129HT124LQ137E: Average"),
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

summary(ComplexT_MD_results_2[,2]);sd(ComplexT_MD_results_2[,2])

ASP <- ggplot()+
  geom_density(data=WT_MD_results_3,aes(WT_MD_results_3[,2],fill="WT"),
               bw=0.1,alpha=0.5)+
  geom_density(data=ComplexT_MD_results_2,aes(ComplexT_MD_results_2[,2],fill="D238SW269D: Average"),
               bw=0.1,alpha=0.5)+
  geom_density(data=Complex_MD_results_2,aes(Complex_MD_results_2[,2],fill="D238SW269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=ComplexL_MD_results_2,aes(ComplexL_MD_results_2[,2],fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=ComplexTL_MD_results_2,aes(ComplexTL_MD_results_2[,2],fill="D129HT124LQ137E: Average"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Acid-His")+xlab("Acid-His distance (Å)")+
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

# ------------

ggplot()+
  geom_density(data=ComplexT_MD_results,aes(min,fill="D238SW269D: Average"),
                       bw=0.1,alpha=0.5)+
  geom_density(data=Complex_MD_results,aes(min,fill="D238SW269D"),
               bw=0.1,alpha=0.5)+
  geom_density(data=ComplexL_MD_results,aes(min,fill="D129HT124LQ137E"),
               bw=0.1,alpha=0.5)+
  geom_density(data=ComplexTL_MD_results,aes(min,fill="D129HT124LQ137E: Average"),
               bw=0.1,alpha=0.5)+
  theme_minimal()+ggtitle("Serine-Substrate distance for the potential catalytic triad")+xlab("Serine-Substrate distance (Å)")+
  guides(fill=guide_legend(title="Variant"))+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_fill_manual(values=c("#56B2E9", "#D69F00", "#7e082a"
                             ,"#a3ff00","#afeeee",
                             "#9910f9","#136d14","#f66909"))

summary(ComplexT_MD_results_10[,2]);skewness(ComplexT_MD_results_10[,2])

wilcox.test(ComplexT_MD_results_10[,2],Complex_MD_results_10[,2],mu=0.0025,alternative="greater")

# -------------------------
# Oxyanion distances

names(Complex_MD_results) <-c("step","ServsC1","ServsC8","ServsC11",
                              "ServsC23","ServsC30","ServsC33")

Complex_MD_results_Asn1 <- read.table(file="Results/MD_Complex/Asn301HD1vsC1.dat")
Complex_MD_results_Asn2 <- read.table(file="Results/MD_Complex/Asn301HD2vsC1.dat")
Complex_MD_results_Asn3 <- read.table(file="Results/MD_Complex/Asn301HD1vsC8.dat")
Complex_MD_results_Asn4 <- read.table(file="Results/MD_Complex/Asn301HD2vsC8.dat")
Complex_MD_results_Asn5 <- read.table(file="Results/MD_Complex/Asn301HD1vsC11.dat")
Complex_MD_results_Asn6 <- read.table(file="Results/MD_Complex/Asn301HD2vsC11.dat")
Complex_MD_results_Asn7 <- read.table(file="Results/MD_Complex/Asn301HD1vsC23.dat")
Complex_MD_results_Asn8 <- read.table(file="Results/MD_Complex/Asn301HD2vsC23.dat")
Complex_MD_results_Asn9 <- read.table(file="Results/MD_Complex/Asn301HD1vsC30.dat")
Complex_MD_results_Asn10 <- read.table(file="Results/MD_Complex/Asn301HD2vsC30.dat")
Complex_MD_results_Asn11 <- read.table(file="Results/MD_Complex/Asn301HD1vsC33.dat")
Complex_MD_results_Asn12 <- read.table(file="Results/MD_Complex/Asn301HD2vsC33.dat")

Complex_MD_results_Asn <- Complex_MD_results_Asn1[,2]
Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn2[,2]);Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn3[,2])
Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn4[,2]);Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn5[,2])
Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn6[,2]);Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn7[,2])
Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn8[,2]);Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn9[,2])
Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn10[,2]);Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn11[,2])
Complex_MD_results_Asn <- cbind(Complex_MD_results_Asn,Complex_MD_results_Asn12[,2])

Complex_Oxyanion <- c()
for (i in 1:length(ComplexMin)){
  a <- which(Complex_MD_results==ComplexMin[i],arr.ind=T)# Get column index
  if (a[,2]==2){
  Complex_Oxyanion <- c(Complex_Oxyanion,min(Complex_MD_results_Asn[a[,1],1:2]))}
  else if (a[,2]==3){
    Complex_Oxyanion <- c(Complex_Oxyanion,min(Complex_MD_results_Asn[a[,1],3:4]))}
  else if (a[,2]==4){
    Complex_Oxyanion <- c(Complex_Oxyanion,min(Complex_MD_results_Asn[a[,1],5:6]))}
  else if (a[,2]==5){
    Complex_Oxyanion <- c(Complex_Oxyanion,min(Complex_MD_results_Asn[a[,1],7:8]))}
  else if (a[,2]==6){
    Complex_Oxyanion <- c(Complex_Oxyanion,min(Complex_MD_results_Asn[a[,1],9:10]))}
  else if (a[,2]==7){
    Complex_Oxyanion <- c(Complex_Oxyanion,min(Complex_MD_results_Asn[a[,1],11:12]))}
}

names(Complex_MD_results3) <-c("step","ServsC1","ServsC8","ServsC11",
                              "ServsC23","ServsC30","ServsC33")

Complex_MD_results3_Asn1 <- read.table(file="Results/MD_Complex3/Asn301HD1vsC1.dat")
Complex_MD_results3_Asn2 <- read.table(file="Results/MD_Complex3/Asn301HD2vsC1.dat")
Complex_MD_results3_Asn3 <- read.table(file="Results/MD_Complex3/Asn301HD1vsC8.dat")
Complex_MD_results3_Asn4 <- read.table(file="Results/MD_Complex3/Asn301HD2vsC8.dat")
Complex_MD_results3_Asn5 <- read.table(file="Results/MD_Complex3/Asn301HD1vsC11.dat")
Complex_MD_results3_Asn6 <- read.table(file="Results/MD_Complex3/Asn301HD2vsC11.dat")
Complex_MD_results3_Asn7 <- read.table(file="Results/MD_Complex3/Asn301HD1vsC23.dat")
Complex_MD_results3_Asn8 <- read.table(file="Results/MD_Complex3/Asn301HD2vsC23.dat")
Complex_MD_results3_Asn9 <- read.table(file="Results/MD_Complex3/Asn301HD1vsC30.dat")
Complex_MD_results3_Asn10 <- read.table(file="Results/MD_Complex3/Asn301HD2vsC30.dat")
Complex_MD_results3_Asn11 <- read.table(file="Results/MD_Complex3/Asn301HD1vsC33.dat")
Complex_MD_results3_Asn12 <- read.table(file="Results/MD_Complex3/Asn301HD2vsC33.dat")

Complex_MD_results3_Asn <- Complex_MD_results3_Asn1[,2]
Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn2[,2]);Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn3[,2])
Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn4[,2]);Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn5[,2])
Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn6[,2]);Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn7[,2])
Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn8[,2]);Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn9[,2])
Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn10[,2]);Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn11[,2])
Complex_MD_results3_Asn <- cbind(Complex_MD_results3_Asn,Complex_MD_results3_Asn12[,2])

Complex_Oxyanion3 <- c()
for (i in 1:length(ComplexMin3)){
  a <- which(Complex_MD_results3==ComplexMin3[i],arr.ind=T) # Get column index
if (a[,2]==2){
  Complex_Oxyanion3 <- c(Complex_Oxyanion3,min(Complex_MD_results3_Asn[a[,1],1:2]))}
else if (a[,2]==3){
  Complex_Oxyanion3 <- c(Complex_Oxyanion3,min(Complex_MD_results3_Asn[a[,1],3:4]))}
else if (a[,2]==4){
  Complex_Oxyanion3 <- c(Complex_Oxyanion3,min(Complex_MD_results3_Asn[a[,1],5:6]))}
else if (a[,2]==5){
  Complex_Oxyanion3 <- c(Complex_Oxyanion3,min(Complex_MD_results3_Asn[a[,1],7:8]))}
else if (a[,2]==6){
  Complex_Oxyanion3 <- c(Complex_Oxyanion3,min(Complex_MD_results3_Asn[a[,1],9:10]))}
else if (a[,2]==7){
  Complex_Oxyanion3 <- c(Complex_Oxyanion3,min(Complex_MD_results3_Asn[a[,1],11:12]))}
}

summary(Complex_Oxyanion3)
summary(Complex_Oxyanion)
CO <- c(Complex_Oxyanion,Complex_Oxyanion3)
CO <- as.data.frame(CO)
summary(CO);sd(CO$CO)
skewness(CO$CO)

co <- as.data.frame(Complex_Oxyanion)

ggplot()+geom_density(data=data.frame(CO),aes(CO),
                      bw=0.2,alpha=0.5,fill="#56B2E9")+
  theme_minimal()+ggtitle("Oxyanion hole distance")+xlab("Asn301 - carbonyl oxygen distance (Å)")+
  theme(axis.text = element_text(size =11),
        axis.title = element_text(size =13),
        legend.title = element_text(size =13),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5,size = 15))
