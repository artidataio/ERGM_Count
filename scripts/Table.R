library(data.table)
library(ztable)
library(stringr)
library(rprojroot)

#set up root
root <- find_root(is_rstudio_project)
options(ztable.type = "latex")

coef_LS <- readRDS(paste(root,"/processed data/coef_LS",sep=""))
table_LS <- list()

for(names in names(coef_LS)){
  DT <- coef_LS[[names]]
  DT$parameter <- str_replace(DT$parameter, "\\..+", "")
  rgroup <- c("I",
              "II",
              "III",
              "IV",
              "V",
              "VI")
  
  n.rgroup <- rep(nrow(DT)/6,6)

  z <- ztable(DT[,c(1,2,3,5)],digits = 4)
  z <- addRowColor(z,DT[`p-value`<0.1,which=T]+1,color="platinum")
  z <- addrgroup(z,rgroup=rgroup,n.rgroup=n.rgroup,cspan.rgroup=1)
  z <- update_ztable(z,size = 2, tablewidth = 0.4,caption.placement = "bottom",placement="H")
  
  attr(z$x,"footer")<-"Highlighted parameters are significant at 0.1 level"
  table_LS[[str_sub(names,1,-4)]] <- z
}

saveRDS(table_LS,paste(root,"/processed data/table_LS",sep=""))

