#Load the librairies
library(miic)
library(qgraph)
library(corrplot)
library(Hmisc)

#Load the data
Hemato=hematoData

#Explore data
head(Hemato)
View(Hemato)
dim(Hemato)
summary(Hemato)
colnames(Hemato)
sum(is.na(Hemato)) #There is no missing value in our data

#Convert our data to numeric
sapply(Hemato, class)
data <- apply(Hemato,2,function(x) as.numeric(as.character(x)))
View(data)                    
sapply(data, class)

#TF Correlation network
res_cor <- cor(data) 
corrplot(res_cor, type = "upper", order = "hclust",
         col = c("black", "white"), bg = "lightblue", tl.col="black", tl.srt=80)
View(round(res_cor, 2))

#Choice of several thresholds according to the mean of r
data <- as.matrix(hematoData)
result_corr <- rcorr(data, type="spearman")
result_corr
mean(result_corr$r) #0,13

#Plot the correlation network while trying different threshold to filter the edges
names = colnames(data)
Graph_cor <- qgraph(res_cor, graph = "cor", layout = "spring", threshold = 0.05,
                    alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                    legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                    negCol = "red",details = TRUE)
Graph_cor_01 <- qgraph(res_cor, graph = "cor", layout = "spring", threshold = 0.04,
                    alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                    legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                    negCol = "red",details = TRUE)
Graph_cor_02 <- qgraph(res_cor, graph = "cor", layout = "spring", threshold = 0.02,
                    alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                    legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                    negCol = "red",details = TRUE)
Graph_cor_03 <- qgraph(res_cor, graph = "cor", layout = "spring", threshold = 0.01,
                    alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                    legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                    negCol = "red",details = TRUE)
Graph_cor_05 <- qgraph(res_cor, graph = "cor", layout = "spring", threshold = 0.3,
                    alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                    legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                    negCol = "red",details = TRUE)

#TF Partial Correlation network
#Convert data to numeric
data <- apply(Hemato,2,function(x) as.numeric(as.character(x)))
#Covariance matrix
res_cov = cov(data)
#View(round(res_cov, 2))
diag(res_cov)
det(res_cov)

#Regularization by adding a small value lambda to the covariance matrix diagonal
#lambda = 1e-5
for (i in seq(1,33))
{
  res_cov[i,i]=res_cov[i,i]+1e-5
  
}
#Invert the covariance matrix
inv_cov = solve(res_cov)
inv_cov

#Plot the partial correlation network while trying different threshold to filter the edges
Graph_pcor <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.05,
                    alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                    legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                    negCol = "red",details = TRUE)
Graph_pcor_01 <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.04,
                       alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                       legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                       negCol = "red",details = TRUE)
Graph_pcor_02 <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.02,
                       alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                       legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                       negCol = "red",details = TRUE)
Graph_pcor_03 <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.01,
                       alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                       legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                       negCol = "red",details = TRUE)

#Trying with lambda = 1e-10
for (i in seq(1,33))
{
  res_cov[i,i]=res_cov[i,i]+1e-10
}

inv_cov = solve(res_cov)

Graph_pcor <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.05,
                     alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                     legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                     negCol = "red",details = TRUE)
Graph_pcor_01 <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.04,
                        alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                        legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                        negCol = "red",details = TRUE)
Graph_pcor_02 <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.02,
                        alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                        legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                        negCol = "red",details = TRUE)
Graph_pcor_03 <- qgraph(inv_cov, graph = "pcor", layout = "spring", threshold = 0.01,
                        alpha = 0.05, nodeNames = names,sampleSize = nrow(data),
                        legend.cex = 0.3,vsize = 3.5,label.cex=3,posCol = "black",
                        negCol = "red",details = TRUE)
