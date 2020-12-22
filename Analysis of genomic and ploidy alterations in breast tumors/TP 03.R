###Load libraries
library("miic")
library("bnlearn")
library("janitor")
library("igraph")
library("qgraph")
library("pcalg")

###Functions

#Plot graph
plot_graph <- function(graph){
  graph = delete.vertices(graph, which(degree(graph)==0))
  e = get.edgelist(graph, names=FALSE)
  l = qgraph.layout.fruchtermanreingold(e, vcount=vcount(graph),
                                        area=8*(vcount(graph)^2), repulse.rad=(vcount(graph)^3.1))
  # color nodes
  colors = rep("orange", length(V(graph)))
  colors[names(V(graph)) == tolower(names(V(graph)))] = "yellow"
  colors[names(V(graph)) == "Ploidy"] = "green"
  
  plot(graph, vertex.color=colors,layout=l, vertex.size=9,
       vertex.label.cex=1, edge.arrow.size=.05)
}

#Idendify the mutated genes that are significantly related to gene expression
lc_most_related_to_uc = function(graph) {
  l=list()
  for(node in names(V(graph))) {
    if (node == tolower(node)) { l[node]=0 }
  }
  for(i in 1:nrow(get.edgelist(graph))) {
    r <- get.edgelist(graph)[i,]
    if (r[1] == toupper(r[1]) & r[2] == tolower(r[2])) { l[r[2]]=l[[r[2]]]+1}
    if (r[2] == toupper(r[2]) & r[1] == tolower(r[1])) { l[r[1]]=l[[r[1]]]+1}
  }
  return(sort(unlist(l), decreasing=T))
}

###Load & explore data
Data=cosmicCancer
View(cosmicCancer)
head(cosmicCancer)
dim(cosmicCancer)
colnames(cosmicCancer)

###1. Network reconstruction with the hill-climbing approach
hc<-hc(cosmicCancer)

##Cleaning data
#Missing value
sum(is.na(cosmicCancer)) # 8 na
summary(cosmicCancer)
dataCleaned<- cosmicCancer[complete.cases(cosmicCancer), ]
dim(dataCleaned)
#Constant variables
same <- sapply(dataCleaned, function(.col){ all(.col[1L] == .col)})
dataCleaned <- dataCleaned[!same]
dim(dataCleaned)
#'Ploidy' variable
dataCleaned[,"Ploidy"]=as.factor(dataCleaned[,"Ploidy"])
class(dataCleaned[,"Ploidy"])
levels(dataCleaned[,"Ploidy"])

##Apply Hill-Climbing
#hc model
hc<-hc(dataCleaned)
hc

#Graph network
adj_mat <- amat(hc)
graph_hc = graph_from_adjacency_matrix(adj_mat)
plot_graph(graph_hc)

#Idendify the lower case nodes most related to upper case nodes
lc_most_related_to_uc(graph_hc)[1:10]

#Identify the hubs
sort(hub_score(graph_hc)$vector, decreasing = TRUE)[1:10]

#Top 10 nodes and edges interms of betweenness centrality measure.
#Nodes
sort(betweenness(graph_hc), decreasing = TRUE)[1:10]
#Edges
ed=edge_betweenness(graph_hc)
i=order(ed, decreasing=T)[1:10]
edges=get.edgelist(graph_hc)[i,]
colnames(edges) = c('From', 'To')
View(edges)

###2. Network reconstruction with the PC approach
##alpha=0.01
#Graph network
data_pc=data.matrix(dataCleaned)
data_pc=data_pc-1
nlevels=apply(data_pc,2,function(x) length(attr(as.factor(x), "levels")))
suffStat=list(dm = data_pc, nlev=nlevels, adaptDF = FALSE)
pc_model <- pc(suffStat, indepTest=disCItest, alpha=0.01, labels=colnames(data_pc))
pc_model
pc<-as.bn(pc_model, check.cycles = FALSE)
adj_pc<-amat(pc)
graph_pc = graph_from_adjacency_matrix(adj_pc)
plot_graph(graph_pc)

#Idendify the lower case nodes most related to upper case nodes
lc_most_related_to_uc(graph_pc)[1:6]

#Identify the hubs
sort(hub_score(graph_pc)$vector, decreasing = TRUE)[1:10]

#Top 10 nodes and edges interms of betweenness centrality measure.
#Nodes
sort(betweenness(graph_pc), decreasing = TRUE)[1:10]
#Edges
ed=edge_betweenness(graph_pc)
i=order(ed, decreasing=T)[1:10]
edges=get.edgelist(graph_pc)[i,]
colnames(edges) = c('From', 'To')
View(edges)

##alpha=0.06
#Graph network
data_pc=data.matrix(dataCleaned)
data_pc=data_pc-1
nlevels=apply(data_pc,2,function(x) length(attr(as.factor(x), "levels")))
suffStat=list(dm = data_pc, nlev=nlevels, adaptDF = FALSE)
pc_model <- pc(suffStat, indepTest=disCItest, alpha=0.06, labels=colnames(data_pc))
pc_model
pc<-as.bn(pc_model, check.cycles = FALSE)
adj_pc<-amat(pc)
graph_pc = graph_from_adjacency_matrix(adj_pc)
plot_graph(graph_pc)

#Idendify the lower case nodes most related to upper case nodes
lc_most_related_to_uc(graph_pc)[1:6]

#Identify the hubs
sort(hub_score(graph_pc)$vector, decreasing = TRUE)[1:10]

#Top 10 nodes and edges interms of betweenness centrality measure.
#Nodes
sort(betweenness(graph_pc), decreasing = TRUE)[1:10]
#Edges
ed=edge_betweenness(graph_pc)
i=order(ed, decreasing=T)[1:10]
edges=get.edgelist(graph_pc)[i,]
colnames(edges) = c('From', 'To')
View(edges)

###3. Network reconstruction with the MIIC approach
data(cosmicCancer_stateOrder)

##Execute MIIC (reconstruct graph)
#n_shuffles = 100 & conf_threshold = 0.001
miic.res <- miic(
  input_data = cosmicCancer, state_order = cosmicCancer_stateOrder, latent = "yes",
  n_shuffles = 100, conf_threshold = 0.001)

graph_miic = graph_from_adjacency_matrix(miic.res$adj_matrix)
plot_graph(graph_miic)

#n_shuffles = 200 & conf_threshold = 0.01
miic.res <- miic(
  input_data = cosmicCancer, state_order = cosmicCancer_stateOrder, latent = "yes",
  n_shuffles = 200, conf_threshold = 0.01)


graph_miic = graph_from_adjacency_matrix(miic.res$adj_matrix)
plot_graph(graph_miic)

#n_shuffles = 500 & conf_threshold = 0.001
miic.res <- miic(
  input_data = cosmicCancer, state_order = cosmicCancer_stateOrder, latent = "yes",
  n_shuffles = 500, conf_threshold = 0.001)


graph_miic = graph_from_adjacency_matrix(miic.res$adj_matrix)
plot_graph(graph_miic)

#Idendify the lower case nodes most related to upper case nodes
lc_most_related_to_uc(graph_miic)[1:6]

#Identify the hubs
sort(hub_score(graph_miic)$vector, decreasing = TRUE)[1:10]

#Top 10 nodes and edges interms of betweenness centrality measure.
#Nodes
sort(betweenness(graph_miic), decreasing = TRUE)[1:10]
#Edges
ed=edge_betweenness(graph_miic)
i=order(ed, decreasing=T)[1:10]
edges=get.edgelist(graph_miic)[i,]
colnames(edges) = c('From', 'To')
View(edges)
