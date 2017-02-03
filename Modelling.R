library(intergraph)
library(data.table)
library(network)
library(igraph)
library(ergm)
library(ergm.count)
library(xtable)

data_LS <- readRDS(paste(getwd(),"/processed data/data_LS",sep=""))

fit_LS <- list()
coef_LS <- list()
  
#Edgewise Models

#model 1 : Poisson
#fitting Poisson
set.seed(240193)
poisson_LS <- list()
poisson_DT <- data.table()
for(i in 1:6){
  training <- asNetwork(data_LS[[i]]$igraph)
  model <- ergm(training~sum(pow=1),
                response = "weight",reference = ~Poisson)
  coef_DT <- setDT(summary(model)$coef,keep.rownames=T)
  
}
fit_LS[["poisson_LS"]]<-poisson_LS

#fitting Geometric
set.seed(240193)
geom_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  geom_LS[[i]] <- ergm(training_i~sum(pow=1),
                          response = "weight",reference = ~Geometric)
}
fit_LS[["geom_LS"]]<-geom_LS


#fitting CMP
set.seed(240193)
cmp_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  cmp_LS[[i]] <- ergm(training_i~
                        sum(pow=1) + 
                        CMP,
                          response = "weight",reference = ~Poisson)
}
fit_LS[["cmp_LS"]] <- cmp_LS

#fitting nonzero
set.seed(240193)
nonzero_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  nonzero_LS[[i]] <- ergm(training_i~
                            sum(pow=1) + 
                            nonzero,
                          response = "weight",reference = ~Poisson)
}
fit_LS[["nonzero_LS"]] <- nonzero_LS

#fitting mutual
#geomean
set.seed(240193)
geomean_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  geomean_LS[[i]] <- ergm(training_i~
                              sum(pow=1) + 
                              mutual("geometric"),
                            response = "weight",reference = ~Poisson)
}
fit_LS[["geomean_LS"]] <- geomean_LS

#minimum
set.seed(240193)
minimum_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  minimum_LS[[i]] <- ergm(training_i~
                            sum(pow=1) + 
                            mutual("min"),
                          response = "weight",reference = ~Poisson)
}
fit_LS[["minimum_LS"]] <- minimum_LS

#nabsdiff
set.seed(240193)
nabsdiff_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  nabsdiff_LS[[i]] <- ergm(training_i~
                             sum(pow=1) + 
                             mutual("nabsdiff"),
                           response = "weight",reference = ~Poisson)
}
fit_LS[["nabsdiff_LS"]] <- nabsdiff_LS

#fitting transitivity
set.seed(240193)
transitivity_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  transitivity_LS[[i]] <- ergm(training_i~
                                 sum(pow=1) + 
                                 transitiveweights("geomean","sum","geomean"),
                               response = "weight",reference = ~Poisson)
}
fit_LS[["transitivity_LS"]] <- transitivity_LS



#Exogenous Covariate
exco_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  exco_LS[[i]] <- ergm(training_i~
                         sum(pow=1) +
                         nodeocov("market_value")+
                         nodeicov("market_value"),
                       response = "weight",reference = ~Poisson)
}
fit_LS[["exco_LS"]] <- exco_LS


#Predictive Modelling
#Best edgewise
set.seed(240193)
edgewise_LS <- list()
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  edgewise_LS[[i]] <- ergm(training_i~
                         sum(pow=1) +
                         nodeocov("market_value")+
                         nodeicov("market_value"),
                       response = "weight",reference = ~Poisson)
}
fit_LS[["edgewise_LS"]] <- edgewise_LS

#Best edgewise
set.seed(240193)
pairwise_LS <- list() 
for(i in 1:6){
  training_i <- asNetwork(data_LS[[i]]$igraph)
  pairwise_LS[[i]] <- ergm(training_i~
                             sum(pow=1) +
                             nodeocov("market_value")+
                             nodeicov("market_value")+
                             mutual("nabsdiff"),
                           response = "weight",reference = ~Poisson)
}
fit_LS[["pairwise_LS"]] <- pairwise_LS 


saveRDS(fit_LS,"fit_LS")






estimate <- unlist(lapply(poisson_LS, function(x){x$coef}),use.names = F)
standard_error <- unlist(lapply(poisson_LS, function(x){sqrt(x$covar[1,1])}),use.names = F)
parameter <- unlist(lapply(poisson_LS, function(x){names(x$coef)}),use.names = F)
model <- paste(rep("fit",6),1:6)
type <- c(rep("estimate",6),rep("standard error",6))
table_DT <- data.table(model,parameter,type,value = c(estimate,standard_error))

