# Set the current working directory
getwd();setwd("/home/sergiroda/Documents/Simulations/PETase_Manuel/Multiple_active_sites")

# Generate the table of the PELE results. Transform the PELE results with the sed tool
# (sed -i 's/    / /g' *.out) to enable its recognition.
Total <- read.table("out/report_1.out",header=FALSE,sep=" ")
Total[,26] <- 1
Total <- Total[5:nrow(Total),]
for (i in 2:128) {
  table <- read.table(sprintf("out/report_%s.out",i),header = FALSE,sep=" ")
  table[,26] <- i
  Total <- rbind(Total,table[5:nrow(table),]) # Append the table of one report with 
  # the 4 first elements of the PELE simulation eliminated.
}

# Density plot of the different active sites
install.packages("sm");library(sm)

# First active site
xlim <- range(Total[,5:10])
D1 <- density(Total[,5],from = xlim[1], to = xlim[2]);D2 <- density(Total[,6],from = xlim[1], to = xlim[2])
D3 <- density(Total[,7],from = xlim[1], to = xlim[2]);D4 <- density(Total[,8],from = xlim[1], to = xlim[2])
D5 <- density(Total[,9],from = xlim[1], to = xlim[2]);D6 <- density(Total[,10],from = xlim[1], to = xlim[2])

ylim <- range(D1$y, D2$y, D3$y, D4$y, D5$y, D6$y)
plot(D1, col = 1, lwd = 2, type = "l",xlim = xlim, ylim = ylim)
lines(D2, col = 2, lwd = 2);lines(D3, col = 3, lwd = 2)
lines(D4, col = 4, lwd = 2);lines(D5, col = 5, lwd = 2);lines(D6, col = 6, lwd = 2)

# Second active site
xlim <- range(Total[,11:16])
D1 <- density(Total[,11],from = xlim[1], to = xlim[2]);D2 <- density(Total[,12],from = xlim[1], to = xlim[2])
D3 <- density(Total[,13],from = xlim[1], to = xlim[2]);D4 <- density(Total[,14],from = xlim[1], to = xlim[2])
D5 <- density(Total[,15],from = xlim[1], to = xlim[2]);D6 <- density(Total[,16],from = xlim[1], to = xlim[2])

ylim <- range(D1$y, D2$y, D3$y, D4$y, D5$y, D6$y)
plot(D1, col = 1, lwd = 2, type = "l",xlim = xlim, ylim = ylim)
lines(D2, col = 2, lwd = 2);lines(D3, col = 3, lwd = 2)
lines(D4, col = 4, lwd = 2);lines(D5, col = 5, lwd = 2);lines(D6, col = 6, lwd = 2)

# Third active site
xlim <- range(Total[,17:22])
D1 <- density(Total[,17],from = xlim[1], to = xlim[2]);D2 <- density(Total[,18],from = xlim[1], to = xlim[2])
D3 <- density(Total[,19],from = xlim[1], to = xlim[2]);D4 <- density(Total[,20],from = xlim[1], to = xlim[2])
D5 <- density(Total[,21],from = xlim[1], to = xlim[2]);D6 <- density(Total[,22],from = xlim[1], to = xlim[2])

ylim <- range(D1$y, D2$y, D3$y, D4$y, D5$y, D6$y)
plot(D1, col = 1, lwd = 2, type = "l",xlim = xlim, ylim = ylim)
lines(D2, col = 2, lwd = 2);lines(D3, col = 3, lwd = 2)
lines(D4, col = 4, lwd = 2);lines(D5, col = 5, lwd = 2);lines(D6, col = 6, lwd = 2)

summary(Total[,11:16])
summary(Total[,17:22])

Total[which.min(Total$V22),]
Total[which.min(Total$V13),]