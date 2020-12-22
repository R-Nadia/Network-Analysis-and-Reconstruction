library(bnlearn)
library(igraph)
library(pcalg)
install.packages("BiocManager")
BiocManager::install("Rgraphviz")
library("grid")

## Function to plot graphs
plot_graph <- function(adj_mat){
  graph <- graph_from_adjacency_matrix(adj_mat)
  plot(graph, edge.arrow.size=0.1, vertex.size=10, vertex.label.cex=0.8)
}

## Function to calculate scores
scores <- function(real_network, model) {
  cm = bnlearn::compare(bnlearn::skeleton(real_network), bnlearn::skeleton(model))
  precision = cm$tp / (cm$tp + cm$fp)
  recall = cm$tp /(cm$tp + cm$fn)
  fscore = (2 * precision * recall) /(precision + recall)
  return(c(cm, "precision" = precision, "recall" = recall, "fscore" = fscore))
}
#########
##  1
########
### a
data = insurance
head(data)
View(data)
dim(data)  #20000 * 27
colnames(data) 

### b
modelstring = paste0("[Age][Mileage][SocioEcon|Age][GoodStudent|Age:SocioEcon]",
                     "[RiskAversion|Age:SocioEcon][OtherCar|SocioEcon][VehicleYear|SocioEcon:RiskAversion]",
                     "[MakeModel|SocioEcon:RiskAversion][SeniorTrain|Age:RiskAversion]",
                     "[HomeBase|SocioEcon:RiskAversion][AntiTheft|SocioEcon:RiskAversion]",
                     "[RuggedAuto|VehicleYear:MakeModel][Antilock|VehicleYear:MakeModel]",
                     "[DrivingSkill|Age:SeniorTrain][CarValue|VehicleYear:MakeModel:Mileage]",
                     "[Airbag|VehicleYear:MakeModel][DrivQuality|RiskAversion:DrivingSkill]",
                     "[Theft|CarValue:HomeBase:AntiTheft][Cushioning|RuggedAuto:Airbag]",
                     "[DrivHist|RiskAversion:DrivingSkill][Accident|DrivQuality:Mileage:Antilock]",
                     "[ThisCarDam|RuggedAuto:Accident][OtherCarCost|RuggedAuto:Accident]",
                     "[MedCost|Age:Accident:Cushioning][ILiCost|Accident]",
                     "[ThisCarCost|ThisCarDam:Theft:CarValue][PropCost|ThisCarCost:OtherCarCost]")
dag = model2network(modelstring)

### c
class(dag)
dag

### d
adj_mat = bnlearn::amat(dag) 

### e 
plot_graph(adj_mat)




#########
##  2
########
### b
data_hc = bnlearn::hc(data)
class(data_hc)
data_hc

### c
adj_matHC = bnlearn::amat(data_hc)

### d
plot_graph(adj_matHC)

### e
scores(dag, data_hc)

### f 
diff.args = list(show.first=FALSE)
graphviz.compare(bnlearn::skeleton(dag), bnlearn::skeleton(data_hc),
                 diff.args=diff.args, shape = "ellipse", 
                 sub=c("", paste("Prediction of real graph using HC with: \n - red",
                                 "for false positive edges \n - blue" ,
                                 "for false negative edges \n PS:",
                                 "orientations aren't taken into account")))
### g
graphviz.compare((dag), (data_hc),
                 diff.args=diff.args, shape = "ellipse", 
                 sub=c("", paste("Prediction of real graph using HC with: \n - red",
                                 "for false positive edges \n - blue" ,
                                 "for false negative edges \n PS:",
                                 "orientations are taken into account")))

#########
##  3
########
### a
#convert to data matrix
data_pc = data.matrix(data) 
#Make the categories start from 0
data_pc = data_pc - 1 
#Compute the number of levels for each variable
nlevels = apply(data_pc, 2, function(x) length(unique(x))) 
#Prepare the suffStat object
suffStat = list(dm=data_pc, nlev=nlevels, adaptDF=FALSE)
pcl <- pc(suffStat, disCItest, labels = colnames(data_pc), alpha=0.05,  verbose = FALSE)
pcl
### b
pcl_bn <- as.bn(pcl)
adj_matPC <- bnlearn::amat(pcl_bn)

### c 
plot_graph(adj_matPC)

### d
scores(dag, pcl_bn)

### e
graphviz.compare(bnlearn::skeleton(dag), bnlearn::skeleton(pcl_bn),
                 diff.args=diff.args, shape = "ellipse", 
                 sub=c("", paste("Prediction of real graph using PC with: \n - red",
                                 "for false positive edges \n - blue" ,
                                 "for false negative edges \n PS:",
                                 "orientations aren't taken into account")))


#########
##  4
########
### a
data_Aracne <- bnlearn::aracne(data)
class(data_Aracne)
data_Aracne
### b
Adj_matAR <- bnlearn::amat(data_Aracne)

### c
plot_graph(Adj_matAR)

### d
scores(dag, data_Aracne)

### e
graphviz.compare(bnlearn::skeleton(dag), bnlearn::skeleton(data_Aracne),
                 diff.args=diff.args, shape="ellipse",
                 sub=c("", paste("Prediction of real graph using ARACNE with: \n - red",
                                 "for false positive edges \n - blue" ,
                                 "for false negative edges \n PS:",
                                 "orientations aren't taken into account")))

















