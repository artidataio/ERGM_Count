#Best model choice for each training data
edgewise_aic <- NULL
model_name <- NULL
for(name in names(fit_LS)){
  coef_name <- gsub("\\..*$", "",names(fit_LS[[name]][[1]]$coef))
  model <- paste(coef_name,collapse=" + ")
  model_name <- c(model_name,model)
  aic <- lapply(fit_LS[[name]],AIC)
  edgewise_aic <- rbind(edgewise_aic,unlist(aic)) 
}
edgewise_aic_DT <- data.table(edgewise_aic)
colnames(edgewise_aic_DT) <- paste("data",1:6,sep="_")
edgewise_aic_DT <- cbind(model_name,edgewise_aic_DT)

best_model <- apply(edgewise_aic, 2, which.min)
best_edgewise_LS <- list()
for(i in 1:6){
  model <- fit_LS$edgewise_LS[[best_model[[i]]]][[i]]
  best_edgewise_LS[[i]]<- model
}
