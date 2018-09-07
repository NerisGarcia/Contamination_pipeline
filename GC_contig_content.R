library(seqinr)
library(ggplot2)
args<-commandArgs(TRUE)
# arg[1] = name of the sample
# arg[2] = path to sample contigs (.fasta)
# arg[3] = name of the output
name= args[1]
contigs = read.fasta(file=args[2])

names <- getAnnot(contigs)
dfContigs <- matrix(ncol = 3, nrow = length(names))
colnames(dfContigs) <- (c("Contig", "Length", "GC"))


for (i in 1:length(names)){
  res = summary.SeqFastadna(contigs[[i]])
  #SUMM <- cbind((c(names[[i]], res$length, res$GC)))
  dfContigs[i, 1] <- c(names[[i]])
  dfContigs[i, 2] <- c(res$length)
  dfContigs[i, 3] <- c(res$GC)
}
(res$GC)
df <- as.data.frame(dfContigs)
#Data frame with large contigs >10.000
dfLC <- df[as.numeric(as.character(df$Length))>=10000,]

#Data frame with short contigs
dfSC <-df[as.numeric(as.character(df$Length))<10000,]
#LC mean

Promedio <- sum(as.numeric(as.character(dfLC$Length))*as.numeric(as.character(dfLC$GC)))/sum(as.numeric(as.character(dfLC$Length)))
#Desviacion tipica
dfLC$Promedio = Promedio
dfLC$DV = "NA"
dfLC$DVp = "NA"
dfLC$DVn = "NA"
dfLC$Select = "NA"

dfSC$Promedio = Promedio
dfSC$DV = 3*sqrt(as.numeric(as.character(dfSC$Length))*Promedio*(1-Promedio))

#df$DV = as.numeric(as.character(df$Length))*df$Promedio*(1-df$Promedio)
#valor intervalo positivo
dfSC$DVp = Promedio + dfSC$DV/ as.numeric(as.character(dfSC$Length))
#valor intervalo negativo
dfSC$DVn = Promedio - dfSC$DV/(as.numeric(as.character(dfSC$Length)))

#True --> Los que nos quedamos
#False --> los que hay que eliminar
dfSC$Select = (as.numeric(as.character(dfSC$GC)) < as.numeric(as.character(dfSC$DVp))) & (as.numeric(as.character(dfSC$GC)) > as.numeric(as.character(dfSC$DVn)))

final <- rbind(dfLC, dfSC)

write.csv(final, paste(args[3])) 
#stats
fTrue <- final[final$Select==TRUE,]
#write.csv(final, paste(args[3], "_TRUE.csv")) 

fFalse <- final[final$Select==FALSE,]
#write.csv(final, paste(args[3], "_FALSE.csv")) 



plot <- final[, c(2,3)]

interval00.05 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.0) & (as.numeric(as.character(plot$GC))<=0.05) ,][,1])))
interval05.10 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.05) & (as.numeric(as.character(plot$GC))<=0.10) ,][,1])))

interval10.15 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.10) & (as.numeric(as.character(plot$GC))<=0.15) ,][,1])))
interval15.20 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.15) & (as.numeric(as.character(plot$GC))<=0.20) ,][,1])))

interval20.25 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.20) & (as.numeric(as.character(plot$GC))<=0.25) ,][,1])))
interval25.30 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.25) & (as.numeric(as.character(plot$GC))<=0.30) ,][,1])))

interval30.35 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.30) & (as.numeric(as.character(plot$GC))<=0.35) ,][,1])))
interval35.40 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.35) & (as.numeric(as.character(plot$GC))<=0.40) ,][,1])))

interval40.45 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.40) & (as.numeric(as.character(plot$GC))<=0.45) ,][,1])))
interval45.50 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.45) & (as.numeric(as.character(plot$GC))<=0.50) ,][,1])))

interval50.55 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.50) & (as.numeric(as.character(plot$GC))<=0.55) ,][,1])))
interval55.60 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.55) & (as.numeric(as.character(plot$GC))<=0.60) ,][,1])))

interval60.65 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.60) & (as.numeric(as.character(plot$GC))<=0.65) ,][,1])))
interval65.70 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.65) & (as.numeric(as.character(plot$GC))<=0.70) ,][,1])))

interval70.75 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.70) & (as.numeric(as.character(plot$GC))<=0.75) ,][,1])))
interval75.80 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.75) & (as.numeric(as.character(plot$GC))<=0.80) ,][,1])))

interval80.85 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.80) & (as.numeric(as.character(plot$GC))<=0.85) ,][,1])))
interval85.90 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.85) & (as.numeric(as.character(plot$GC))<=0.90) ,][,1])))

interval90.95 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.90) & (as.numeric(as.character(plot$GC))<=0.95) ,][,1])))
interval95.100 <- sum(as.numeric(as.character(plot[(as.numeric(as.character(plot$GC))>0.95) & (as.numeric(as.character(plot$GC))<=0.100) ,][,1])))


fplot <- as.data.frame(rbind(c("00-05", interval00.05),c("05-10", interval05.10), c("10-15", interval10.15),c("15-20", interval15.20),c("20-25", interval20.25),c("25-30", interval25.30),c("30-35", interval30.35),c("35-40", interval35.40),c("40-45", interval40.45),c("45-50", interval45.50),c("50-55", interval50.55),c("55-60", interval55.60), c("60-65", interval60.65), c("65-70", interval65.70),c("70-75", interval70.75),c("75-80", interval75.80),c("80-85", interval80.85),c("85-90", interval85.90),c("90-95", interval90.95),c("95-100", interval95.100)))

colnames(fplot)<-(c("GC interval", "Number of bases"))  

L = sum(as.numeric(as.character(fplot$`Number of bases`)))
fplot$percentaje <- (as.numeric(as.character(fplot$`Number of bases`))*100)/L


#save the plot in a variable image to be able to export to svg
finalplot <-ggplot(data=fplot, aes(x=fplot$`GC interval`, y =(as.numeric(as.character(fplot$`Number of bases`))), fill=((fplot$`GC interval`)))) + geom_bar(stat="identity" )  + ylab("Number of bases") + xlab("GC content") +  geom_text(aes(label=round(fplot$percentaje, signif(1))), vjust=-0.3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)  + ggtitle(args[1]) + theme(plot.title = element_text(hjust = 0.5))

#This actually save the plot in a image
fileplotname = paste(args[3], "plot.png")
ggsave(fileplotname, width = 16, height = 9, dpi = 100)
