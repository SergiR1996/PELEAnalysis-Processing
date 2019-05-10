install.packages("e1071")
library(e1071)

PELE_results <- read.csv(file="PETase_Manuel_2HEMHET_D238SL272D.csv")
summary(PELE_results);sd(PELE_results[,5])
quantile(PELE_results[,10])
ServsHis <- PELE_results[,8]
Ser238vsHis <- PELE_results[,9]
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
