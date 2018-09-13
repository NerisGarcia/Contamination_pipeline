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
