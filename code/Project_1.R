
library(tidyverse) ; library(igraph) ; library(network); library(sna)
edges = read.csv("dolphin_edges.csv")
nodes = read.csv("dolphin_nodes.csv")

dol<-matrix( rep(NA, dim(nodes)[1]*dim(nodes)[1]), nrow=dim(nodes)[1], ncol=dim(nodes)[1])
for (i in 1:dim(edges)[1]){
  dol[edges$x[i], edges$y[i]]=1; 
  dol[edges$y[i], edges$x[i]]=1
}
dol[is.na(dol)]=0
dol

dol_graph = graph_from_adjacency_matrix(dol)

ecount(dol_graph); vcount(dol_graph)

par(mfrow=c(1,2))
hist(degree(dol), col="lightblue", xlim=c(0,50),
     xlab="Vertex Degree", ylab="Frequency", main="")

hist(graph.strength(dol_graph), col="pink", xlim=c(0,50),
     xlab="Vertex Strength", ylab="Frequency", main="")
dev.off()

dyad.census(dol)

transitivity(dol_graph) ; reciprocity(dol_graph, mode = 'default')
average.path.length(dol_graph); diameter(dol_graph)

nodes.viz = mutate(nodes, id = row_number()) %>% select(id, everything())
dol.viz<-as.network(edges,matrix.type="edgelist", directed=FALSE, vertex.attr=nodes.viz)
dcolors = c("red", "blue", "black")[match(nodes$sex, c("F", "M", "U"))]
plot(dol.viz, vertex.col = dcolors )
plot(dol.viz)

# hergm
library("hergm")
dolnet<-as.network(dol, directed = FALSE, matrix.type="adjacency" )
set.seed(0)
dol.local <- hergm(dolnet~ edges_ij + triangle_ijk, max_number = 9, sample_size=50000, 
                   parametic=TRUE, posterior.thinning = 2)

summary(dol.local)

gof(dol.local)

plot(dol.local)

# SBM_bernoulli
library(blockmodels)
set.seed(0)
dol.sbm= BM_bernoulli("SBM_sym",dol)
dol.sbm$estimate()
which.max(dol.sbm$ICL)
dol.sbm$plot_obs_pred(2)
dol.sbm$memberships[[2]]$Z
# vizualization
sbm.Z<-dol.sbm$memberships[[2]]$Z
sbZ<-rep(2, 62)
for (i in 1:62){
  if (sbm.Z[i, 1]>sbm.Z[i, 2]){sbZ[i]=1}
}
dcolors= c("red", "blue")[match(sbZ, c(1,2))]
plot(dol.viz, vertex.col = dcolors)
sbZ

# DCSBM_poisson
library(randnet)
BHMC.estimate(dol, 4)
LRBIC(dol, 4, lambda = NULL, model = "both")
ssc <- reg.SSP(dol,K=2,lap=FALSE)
NMI(ssc$cluster,sbZ)
ssc$cluster
ssc$cluster==sbZ
dcolors2= c("red", "blue")[match(ssc$cluster, c(1,2))]
plot(dol.viz, vertex.col=dcolors2 )
