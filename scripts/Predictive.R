library(data.table)
library(rprojroot)
library(stringr)
library(ztable)

root <- find_root(is_rstudio_project)

data_LS <- readRDS(paste(root,"/processed data/data_LS",sep=""))
coef_LS <- readRDS(paste(root,"/processed data/coef_LS",sep=""))

coef_LS$transitivity_DT <- NULL
coef_LS$geometric_DT <- NULL
coef_LS$min_DT <- NULL
coef_LS$geomean_DT <- NULL

DT <- data.table()
for(i in 1:6){
  DT_i <- data_LS[[i]]$data
  DT_i[,":="(name = data_LS[[i]]$name,
             season = data_LS[[i]]$season),]
  foo_DT_i <- data_LS[[i]]$vertices_DT #contain amv
  
  setnames(foo_DT_i,c("HomeTeam","home_amv"))
  DT_i <- merge(DT_i,foo_DT_i,by= "HomeTeam",all.x=T)
  
  setnames(foo_DT_i,c("AwayTeam","away_amv"))
  DT_i <- merge(DT_i,foo_DT_i,by= "AwayTeam",all.x=T)
  
  DT <- rbind(DT,DT_i)
}

make_p_MAT <- function(SUM,
                       nonzero,
                       CMP,
                       nodeocov,
                       nodeicov,
                       edgecov,
                       mutual,
                       home_amv,
                       away_amv){
  max_goal <- 10
  p_MAT <- matrix(numeric((max_goal+1)^2),ncol= max_goal+1)
  #i is home
  for(i in 0:max_goal){
    for(j in 0:max_goal){
      p_MAT[i+1,j+1] <- exp(SUM*(i+j)+
                              nonzero*((i>0)+(j>0))+
                              nodeocov*( home_amv*i+ away_amv*j) +
                              nodeicov*( away_amv*i+ home_amv*j) +
                              edgecov*(i)- 
                              mutual*(abs(i-j)))*(factorial(i)*factorial(j))^(CMP-1)
    }
  }
  p_MAT <- p_MAT/sum(p_MAT)
  return(p_MAT)
}


get_logarithmic_score<- function(p_MAT,obs_i,obs_j){
  #since the the first row and column store p(0)
  return(-log(p_MAT[obs_i+1,obs_j+1]))
}
get_quadratic_score<- function(p_MAT,obs_i,obs_j){
  #since the the first row and column store p(0)
  return(-2*p_MAT[obs_i+1,obs_j+1]+sum(p_MAT^2))
} 
get_spherical_score<- function(p_MAT,obs_i,obs_j){
  #since the the first row and column store p(0)
  return(-p_MAT[obs_i+1,obs_j+1]/sum(p_MAT^2))
}

logarithmic_MAT <- matrix(nrow = length(coef_LS),ncol=6)
quadratic_MAT <- matrix(nrow = length(coef_LS),ncol=6)
spherical_MAT <- matrix(nrow = length(coef_LS),ncol=6)

for(row in 1:length(coef_LS)){
  models <- coef_LS[[row]]
  for(col  in 1:6){ #Because of 6 networks  that we have
    print(str_c(row,col,sep=" "))
    num_param <- nrow(models)/6
    model <- models[seq(num_param*(col-1)+1,num_param*col),]
    model$parameter <- str_replace(model$parameter, "\\..+", "")
    
    #getting the parameter
    SUM <- ifelse("sum" %in% model$parameter,
                  model[parameter=="sum",Estimate],0)
    
    nonzero <- ifelse("nonzero" %in% model$parameter,
                      model[parameter=="nonzero",Estimate],0)
    
    CMP <- ifelse("CMP" %in% model$parameter,
                  model[parameter=="CMP",Estimate],0)
    
    
    nodeocov <- ifelse("nodeocov" %in% model$parameter,
                       model[parameter=="nodeocov",Estimate],0)
    
    edgecov <- ifelse("edgecov" %in% model$parameter,
                      model[parameter=="edgecov",Estimate],0)
    
    nodeicov <- ifelse("nodeicov" %in% model$parameter,
                       model[parameter=="nodeicov",Estimate],0)
    
    mutual <- ifelse("mutual" %in% model$parameter,
                     model[parameter=="mutual",Estimate],0)
    
    print(str_c(SUM,nonzero,CMP,edgecov,nodeocov,nodeicov,sep=" "))
    logarithmic <- numeric(nrow(DT))
    quadratic <- numeric(nrow(DT))
    spherical <- numeric(nrow(DT))
    
    for(match in 1:nrow(DT)){
      p_MAT <- make_p_MAT(SUM,
                          nonzero,
                          CMP, 
                          nodeocov, 
                          nodeicov,
                          edgecov,
                          mutual,
                          DT[match,]$home_amv,
                          DT[match,]$away_amv)
      
      
      obs_i <- DT[match,FTHG]
      obs_j <- DT[match,FTAG]
      
      p_pred <- p_MAT[obs_i+1,obs_j+1]
      l2norm <- sum(p_MAT^2)
    
      logarithmic[match] <- -log(p_pred) #NOTE on the scoring method
      quadratic[match]<- -2*p_pred+l2norm
      spherical[match]<- -p_pred/l2norm
    }
    
    logarithmic_MAT[row,col]<- mean(logarithmic) 
    quadratic_MAT[row,col]<- mean(quadratic)
    spherical_MAT[row,col]<- mean(spherical)
  }
}

model <- NULL
for(names in names(coef_LS)){
  coef_name <- str_replace(coef_LS[[names]][,unique(parameter)],"\\..*$", "")
  coef_name <- str_c(coef_name,collapse = " + ")
  model <- c(model,coef_name)
}
colnames(pearson_DT) <- c("model","I","II","III","IV","V","VI")
logarithmic_DT <- data.table(model,logarithmic_MAT)
colnames(logarithmic_DT) <- c("model","I","II","III","IV","V","VI")

quadratic_DT <- data.table(model,quadratic_MAT)
colnames(quadratic_DT) <- c("model","I","II","III","IV","V","VI")

spherical_DT <- data.table(model,spherical_MAT)
colnames(spherical_DT) <- c("model","I","II","III","IV","V","VI")


###################
#prediction part 2#
###################
logarithmic2_MAT <- matrix(nrow = length(coef_LS),ncol=6)
quadratic2_MAT <- matrix(nrow = length(coef_LS),ncol=6)
spherical2_MAT <- matrix(nrow = length(coef_LS),ncol=6)

for(row in 1:length(coef_LS)){
  models <- coef_LS[[row]]
  for(col  in 1:6){ #Because of 6 networks  that we have
    print(str_c(row,col,sep=" "))
    num_param <- nrow(models)/6
    model <- models[seq(num_param*(col-1)+1,num_param*col),]
    model$parameter <- str_replace(model$parameter, "\\..+", "")
    
    #getting the parameter
    SUM <- ifelse("sum" %in% model$parameter,
                  model[parameter=="sum",Estimate],0)
    
    nonzero <- ifelse("nonzero" %in% model$parameter,
                      model[parameter=="nonzero",Estimate],0)
    
    CMP <- ifelse("CMP" %in% model$parameter,
                  model[parameter=="CMP",Estimate],0)
    
    
    nodeocov <- ifelse("nodeocov" %in% model$parameter,
                       model[parameter=="nodeocov",Estimate],0)
    
    edgecov <- ifelse("edgecov" %in% model$parameter,
                      model[parameter=="edgecov",Estimate],0)
    
    nodeicov <- ifelse("nodeicov" %in% model$parameter,
                       model[parameter=="nodeicov",Estimate],0)
    
    mutual <- ifelse("mutual" %in% model$parameter,
                     model[parameter=="mutual",Estimate],0)
    
    print(str_c(SUM,nonzero,CMP,edgecov,nodeocov,nodeicov,sep=" "))
    logarithmic2 <- numeric(nrow(DT))
    quadratic2 <- numeric(nrow(DT))
    spherical2 <- numeric(nrow(DT))
    
    for(match in 1:nrow(DT)){
      p_MAT <- make_p_MAT(SUM,
                          nonzero,
                          CMP, 
                          nodeocov, 
                          nodeicov,
                          edgecov,
                          mutual,
                          DT[match,]$home_amv,
                          DT[match,]$away_amv)
      
      p_home<- sum(p_MAT[lower.tri(p_MAT)])
      p_draw<- sum(diag(p_MAT))
      p_away<- sum(p_MAT[upper.tri(p_MAT)])
      
      p <- c(p_home,p_draw,p_away) 
      l2norm <- sum(p^2)
        
      obs_i <- DT[match,FTHG]
      obs_j <- DT[match,FTAG]
      
      p_x <- p[ifelse(obs_i>obs_j,1,ifelse(obs_j>obs_i,3,2))]
      
      logarithmic2[match] <- -log(p_x) #NOTE on the scoring method
      quadratic2[match]<- -2*p_x+l2norm
      spherical2[match]<- -p_x/l2norm
    }
    
    logarithmic2_MAT[row,col]<- mean(logarithmic2) 
    quadratic2_MAT[row,col]<- mean(quadratic2)
    spherical2_MAT[row,col]<- mean(spherical2)
  }
}

model <- NULL
for(names in names(coef_LS)){
  coef_name <- str_replace(coef_LS[[names]][,unique(parameter)],"\\..*$", "")
  coef_name <- str_c(coef_name,collapse = " + ")
  model <- c(model,coef_name)
}

logarithmic2_DT <- data.table(model,logarithmic2_MAT)
colnames(logarithmic2_DT) <- c("model","I","II","III","IV","V","VI")

quadratic2_DT <- data.table(model,quadratic2_MAT)
colnames(quadratic2_DT) <- c("model","I","II","III","IV","V","VI")

spherical2_DT <- data.table(model,spherical2_MAT)
colnames(spherical2_DT) <- c("model","I","II","III","IV","V","VI")

assessment_LS <- list()
assessment_LS[["logarithmic_DT"]]<- logarithmic_DT
assessment_LS[["quadratic_DT"]]<- quadratic_DT
assessment_LS[["spherical_DT"]]<- spherical_DT

assessment_LS[["logarithmic2_DT"]]<- logarithmic2_DT
assessment_LS[["quadratic2_DT"]]<- quadratic2_DT
assessment_LS[["spherical2_DT"]]<- spherical2_DT


saveRDS(assessment_LS,paste(root,"/processed data/assessment_LS",sep=""))


