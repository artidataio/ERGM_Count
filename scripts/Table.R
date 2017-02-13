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
  rgroup <- as.roman(1:6)
  n.rgroup <- rep(nrow(DT)/6,6)

  z <- ztable(DT[,c(1,2,3,5)],digits = 4)
  z <- addRowColor(z,DT[`p-value`<0.1,which=T]+1,color="platinum")
  z <- addrgroup(z,rgroup=rgroup,n.rgroup=n.rgroup,cspan.rgroup=NULL)
  
  table_LS[[stri_sub(names,1,-4)]] <- z
}

saveRDS(table_LS,paste(root,"/processed data/table_LS",sep=""))

