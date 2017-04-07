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
  setnames(DT,"parameter", "Statistic")
  setnames(DT,"Estimate", "Parameter Estimate")
  setnames(DT,"Std. Error", "Standard Error")
  
  
  n.rgroup <- rep(nrow(DT)/6,6)

  z <- ztable(DT[,c(1,2,3,5)],digits = 4)
  z <- addRowColor(z,DT[`p-value`<0.05,which=T]+1,color="platinum")
  z <- addrgroup(z,rgroup=rgroup,n.rgroup=n.rgroup,cspan.rgroup=1)
  z <- update_ztable(z,
                     size = 2, 
                     tablewidth = 0.3,
                     caption.placement = "bottom",
                     placement="H",
                     include.rownames = F,
                     align="cccc")
  
  attr(z$x,"footer")<-"Highlighted parameters are significant at 0.05 level"
  table_LS[[str_sub(names,1,-4)]] <- z
}

saveRDS(table_LS,paste(root,"/processed data/table_LS",sep=""))

####################
#Summary Statistics#
####################

data_LS <- readRDS(paste(root,"/processed data/data_LS",sep=""))
DT <- data.table()
for(i in 1:6){
  foo_DT <- data.table(data_LS[[i]]$edges_DT,
                       season = data_LS[[i]]$season,
                       time = data_LS[[i]]$time)
  DT <- rbind(DT,foo_DT)
}

summary_DT <- DT[,.(mean = mean(weight),
                    sd = sd(weight)),by=.(season,time)]

summary_table <- data.frame(t(summary_DT[,c(3,4)]))
colnames(summary_table) <- summary_DT$time
cgroup <- unique(summary_DT$season)
summary_table <- ztable(summary_table, include.rownames=T) 
summary_table <- addcgroup(summary_table, 
                           cgroup= cgroup,
                           n.cgroup = c(2,2,2))

saveRDS(summary_table,paste(root,"/processed data/summary_table",sep=""))



lambda_table <- cbind(summary_DT,V2 = table_LS$poisson$x$`Parameter Estimate`)
lambda_table[,V2:=exp(V2),]
setnames(lambda_table,"V2","$\\hat{\\lambda}$")
lambda_table <- t(lambda_table[,c(5,3,4)])
colnames(lambda_table)<- paste0("\\textbf{",c("I","II","III","IV","V","VI") ,"}")
lambda_table <- ztable(lambda_table)

saveRDS(lambda_table,paste(root,"/processed data/lambda_table",sep=""))
