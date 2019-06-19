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