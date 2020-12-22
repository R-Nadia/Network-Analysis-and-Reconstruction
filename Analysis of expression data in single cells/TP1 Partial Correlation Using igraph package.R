## Part 2 - Partial Correlation
detach("package:miic", unload=TRUE)
library(miic)

data("hematoData")
hematoData
is.factor(hematoData[,1])

hematoData_int = hematoData
for(iCol in c(1:ncol(hematoData))){
  hematoData_int[,iCol] <- as.numeric(as.character(hematoData[,iCol]))
}

is.factor(hematoData_int[,1])
is.numeric(hematoData_int[,1])



#Definir un parametre lambda pour calculer la matrice inverse de la matrice covariance
res_cov = cov(hematoData_int)
#View(round(res_cov, 2))
diag(res_cov)
det(res_cov)

plot_graph_p <- function(graph, thresh, lambda)
{
  for (i in seq(1,33))
  {
    res_cov[i,i]=res_cov[i,i]+lambda
    
  }
  res_corr_filter <- -solve(res_cov)
  res_corr_filter[abs(res_corr_filter) < thresh] <- 0
  graph <- graph_from_adjacency_matrix(res_corr_filter, diag = FALSE, weight = TRUE, mode = "undirected")
  V(graph)$color <- ifelse(attr(V(graph), "names") %in% Hematoporetic, "red", 
                    ifelse(attr(V(graph), "names") %in% Endothelial , "violet", 
                    ifelse(attr(V(graph), "names") %in% Hox , "gray", "deepskyblue")))
  E(graph)$color <- ifelse(E(graph)$weight > 0, "blue","red")
  #E(graph)$width <- abs(E(graph)$weight) * 3
  plot(graph, layout=layout_with_dh)
  
}

#Jouer sur les parametres seuil et lambda
plot_graph_p(res_corr_filter, 0.01, 1e-5)






