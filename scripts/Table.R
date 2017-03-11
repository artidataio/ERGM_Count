library(data.table)
library(ztable)
library(stringi)
library(rprojroot)

#set up root
root <- find_root(is_rstudio_project)


coef_LS <- readRDS(paste(root,"/processed data/coef_LS",sep=""))
table_LS <- list()

for(names in names(coef_LS)){
  DT <- coef_LS[[names]]
  rgroup <- (1:6)
  n.rgroup <- rep(nrow(DT)/6,6)

  z <- ztable(DT[,c(1,2,3,5)],digits = 4)
  z <- addRowColor(z,DT[`p-value`<0.1,which=T]+1,color="platinum")
  z <- addrgroup(z,rgroup=rgroup,n.rgroup=n.rgroup,cspan.rgroup=NULL)
  
  table_LS[[stri_sub(names,1,-4)]] <- z
}

poisson_DT <- coef_LS$poisson_DT
geometric_DT <- coef_LS$geometric_DT
pois_geom_table <- ztable(cbind(poisson_DT[,c(1,2,3,5)],geometric_DT[,c(2,3,5)]))
pois_geom_table <- addCellColor(pois_geom_table,
                                rows = rep(poisson_DT[`p-value`<0.1,which=T]+1,3),
                                cols = rep(c(3,4,5),each = length(poisson_DT[`p-value`<0.1,which=T])),
                                color ="platinum")
pois_geom_table <- addCellColor(pois_geom_table,
                                rows = rep(geometric_DT[`p-value`<0.1,which=T]+1,3),
                                cols = rep(c(6,7,8),each = length(geometric_DT[`p-value`<0.1,which=T])),
                                color ="platinum")

pois_geom_table <- addcgroup(pois_geom_table,
                             cgroup = c("", "Poisson","Geometric"),
                             n.cgroup = c(1,3,3))
pois_geom_table <- vlines(pois_geom_table,
                          add=c(3,6))

attr(pois_geom_table$x,"footer")<-"Highlighted parameters are significant at 0.1 level"
saveRDS(table_LS,paste(root,"/processed data/table_LS",sep=""))

