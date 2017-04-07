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

DT <- fread(paste(root,"/dataset/prediction/EPL 2016-17.txt",sep=""))
DT <- DT[,c(3,4,5,6)]
(team <- DT[,unique(HomeTeam)])
#[1] "Burnley"        "Crystal Palace" "Everton"        "Hull"          
#[5] "Man City"       "Middlesbrough"  "Southampton"    "Arsenal"       
#[9] "Bournemouth"    "Chelsea"        "Man United"     "Leicester"     
#[13] "Stoke"          "Swansea"        "Tottenham"      "Watford"       
#[17] "West Brom"      "Sunderland"     "West Ham"       "Liverpool" 
market_value <- c(3.70, 6.67, 9.96,3.89,
                  21.89, 3.98, 7.93, 18.22,
                  4.08, 20.60, 19.42, 8.21,
                  6.73, 4.72, 17.88,4.89,
                  4.70, 3.64, 7.76, 14.28)

vertices_DT <- data.table(team,market_value)
setnames(vertices_DT,c("HomeTeam","home_amv"))
DT <- merge(DT,vertices_DT,by= "HomeTeam",all.x=T)
setnames(vertices_DT,c("AwayTeam","away_amv"))
DT <- merge(DT,vertices_DT,by= "AwayTeam",all.x=T)


#Important Functions
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
logarithmic_DT <- data.table(model,logarithmic_MAT)
colnames(logarithmic_DT) <- c("Sufficient Statistics","I","II","III","IV","V","VI")

quadratic_DT <- data.table(model,quadratic_MAT)
colnames(quadratic_DT) <- c("Sufficient Statistics","I","II","III","IV","V","VI")

spherical_DT <- data.table(model,spherical_MAT)
colnames(spherical_DT) <- c("Sufficient Statistics","I","II","III","IV","V","VI")


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
colnames(logarithmic2_DT) <- c("Sufficient Statistics","I","II","III","IV","V","VI")

quadratic2_DT <- data.table(model,quadratic2_MAT)
colnames(quadratic2_DT) <- c("Sufficient Statistics","I","II","III","IV","V","VI")

spherical2_DT <- data.table(model,spherical2_MAT)
colnames(spherical2_DT) <- c("Sufficient Statistics","I","II","III","IV","V","VI")

assessment_LS <- list()
assessment_LS[["logarithmic_DT"]]<- logarithmic_DT
assessment_LS[["quadratic_DT"]]<- quadratic_DT
assessment_LS[["spherical_DT"]]<- spherical_DT

assessment_LS[["logarithmic2_DT"]]<- logarithmic2_DT
assessment_LS[["quadratic2_DT"]]<- quadratic2_DT
assessment_LS[["spherical2_DT"]]<- spherical2_DT


saveRDS(assessment_LS,paste(root,"/processed data/assessment_LS",sep=""))




