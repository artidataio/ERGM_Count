library(intergraph)
library(data.table)
library(network)
library(igraph)
library(ergm)
library(ergm.count)
library(rprojroot)

root <- find_root(is_rstudio_project)

data_LS <- readRDS(paste(root,"/processed data/data_LS",sep=""))

fit_LS <- list()
coef_LS <- list()
  
#Edgewise Models

#model 1 : Poisson
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["poisson_LS"]] <- LS
coef_LS[["poisson_DT"]] <- DT 

#model 2: nonzero
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1) + nonzero,
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["nonzero_LS"]] <- LS
coef_LS[["nonzero_DT"]] <- DT 

#model 3: CMP
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1) + CMP,
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["cmp_LS"]] <- LS
coef_LS[["cmp_DT"]] <- DT 

#model 4: nodeocov
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~nodeocov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["nodeocov_LS"]] <- LS
coef_LS[["nodeocov_DT"]] <- DT 

#model 5: nodeicov
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["nodeicov_LS"]] <- LS
coef_LS[["nodeicov_DT"]] <- DT 

#model 6: edgecov
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ edgecov(home_mat),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["edgecov_LS"]] <- LS
coef_LS[["edgecov_DT"]] <- DT 

#model 7: exco
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_7_LS"]] <- LS
coef_LS[["model_7_DT"]] <- DT 

#Model 8: exco +sum
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ sum + edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_8_LS"]] <- LS
coef_LS[["model_8_DT"]] <- DT 

#model 9: exco +nonzero
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ nonzero + edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_9_LS"]] <- LS
coef_LS[["model_9_DT"]] <- DT 

#model 10: exco + cmp
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ CMP + edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_10_LS"]] <- LS
coef_LS[["model_10_DT"]] <- DT 

#model 11: exco + cmp + nonzero
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ nonzero + CMP + edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_11_LS"]] <- LS
coef_LS[["model_11_DT"]] <- DT 


#model 12: nabsdiff
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1) + mutual("nabsdiff"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["nabsdiff_LS"]] <- LS
coef_LS[["nabsdiff_DT"]] <- DT 

#model 13: minimum
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1) + mutual("min"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["min_LS"]] <- LS
coef_LS[["min_DT"]] <- DT 

#model 14: geomean
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1) + mutual("geometric"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["geomean_LS"]] <- LS
coef_LS[["geomean_DT"]] <- DT 

#model 15: exco+mutual 
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ mutual("nabsdiff")+edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_15_LS"]] <- LS
coef_LS[["model_15_DT"]] <- DT 

#model 16: exco+mutual + nonzero
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ nonzero+mutual("nabsdiff")+edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_16_LS"]] <- LS
coef_LS[["model_16_DT"]] <- DT 

#model 17: exco +mutual +cmp
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ CMP +mutual("nabsdiff")+edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_17_LS"]] <- LS
coef_LS[["model_17_DT"]] <- DT 

#model 18: exco + mutual +nonzero +cmp
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  home_mat <- as.matrix(as_adjacency_matrix(data_LS[[i]]$igraph,attr="home"))
  model <- ergm(training~ CMP +nonzero +mutual("nabsdiff")+edgecov(home_mat) + nodeocov("market_value") + nodeicov("market_value"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["model_18_LS"]] <- LS
coef_LS[["model_18_DT"]] <- DT 

#model 19: transitivity
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1) + transitiveweights("min","max","min"),
                response = "weight",reference = ~Poisson)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["transitivity_LS"]] <- LS
coef_LS[["transitivity_DT"]] <- DT 

#model 20: geometric
set.seed(240193)
LS <- list()
DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1),
                response = "weight",reference = ~Geometric)
  LS[[i]] <- model 
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  setnames(coef_DT, "rn", "parameter")
  coef_DT <- cbind(coef_DT,year = data_LS[[i]]$name, season = data_LS[[i]]$season)
  DT <- rbind(DT,coef_DT)
}
fit_LS[["geometric_LS"]] <- LS
coef_LS[["geometric_DT"]] <- DT 


saveRDS(fit_LS,paste(root,"/processed data/fit_LS",sep=""))
saveRDS(coef_LS,paste(root,"/processed data/coef_LS",sep=""))
