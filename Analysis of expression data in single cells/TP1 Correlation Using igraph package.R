#Les bibliotheques dont on a besoin
library(miic)
library(corrplot)
library(igraph)
library(Hmisc)

#telecharger le jeu de donnees
data(hematoData)
#initialiser une liste des differents genes
Endothelial <- c("Erg", "Sox17", "HoxB4", "Sox7", "Tbx20", "Tbx3", "Notch1")
Hematoporetic <- c("Cbfa2t3h", "Sfpi1","Gfi1", "Nfe2", "Gfi1b", "Gata1", "Mitf", "Myb", "Ikaros", "Runx1")
Hox <- c("HoxD8", "HoxB2")

#Analyse descriptive
data <- hematoData
head(data)
ncol(data) #33 genes
nrow(data) #3934 cells
colnames(data) #genes' names
sum(is.na(data)) #verifier si il y a des valeurs manquantes

##### CORRELATION 
#Calculer les coefficients de correlation
data <- as.matrix(data) #transformer la data frame en une matrice
res_corr <- rcorr(data, type="spearman") #utiliser la méthode spearman car on ne sait pas si il la relation entre les donnees et linaire ou pas
res_corr 

#Plot le reseau de correlation
mean(res_corr$r) #calculer la moyenne des coefficients pour pouvoir determiner un seuil

#On cree une fonction qui permet de tracer les graphes avec differents seuils
plot_graph <- function(graph, thresh)
{
  res_corr_filter[abs(res_corr$r) < thresh] <- 0
  graph <- graph_from_adjacency_matrix(res_corr_filter, diag = FALSE, weight = TRUE, mode = "undirected")
  V(graph)$color <- ifelse(attr(V(graph), "names") %in% Hematoporetic, "red", 
                    ifelse(attr(V(graph), "names") %in% Endothelial , "violet", 
                    ifelse(attr(V(graph), "names") %in% Hox , "gray", "deepskyblue")))
  E(graph)$color <- ifelse(E(graph)$weight > 0, "blue","red")
  E(graph)$width <- abs(E(graph)$weight) * 3
  plot(graph, layout=layout_with_dh)
  
}
res_corr_filter <- res_corr$r 

#changer la valeur du seuil pour en tester plusieurs
plot_graph(res_corr_filter, 0.01)
plot_graph(res_corr_filter, 0.02)
plot_graph(res_corr_filter, 0.04)
plot_graph(res_corr_filter, 0.05)
plot_graph(res_corr_filter, 0.3)


