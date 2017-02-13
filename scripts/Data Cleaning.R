library(data.table)
library(RColorBrewer)
library(igraph)

#insert 3 different simiar files
temp <- list.files(paste(getwd(),"/dataset/",sep=""),pattern = "*.txt")
data_LS <- list()
for(i in 1:length(temp)){
  dir <- paste(getwd(),"/dataset/",sep="")
  DT <- fread(paste(dir,temp[i],sep=""))
  DT <- DT[,.(Date,HomeTeam,AwayTeam,FTHG,FTAG)]
  
  #Here we guess that half season denoted by the middle row 
  DT_1 <- DT[1:(nrow(DT)/2),,]
  DT_2 <- DT[(nrow(DT)/2+1):nrow(DT),,]
  
  test_1<- c(DT_1$HomeTeam,DT_1$AwayTeam)
  #Test whether each team has played the same amount of time
  #   Mode    TRUE    NA's 
  #logical      20       0 
  if(sum(table(test_1)==19) != 20){
    print("fark")
    break
  }
  print(head(test_1))
  
  foo1 <- paste(DT_1$HomeTeam,DT_1$AwayTeam)
  foo2 <- paste(DT_1$AwayTeam,DT_1$HomeTeam)
  
  test_2 <- c(foo1,foo2)
  if(sum(table(test_2)==1)!=380){
    print("fark_2")
    break
  }
  
  DT_LS_1 <- list(data = DT_1,name = substr(temp[i],1,13), season = "first")
  DT_LS_2 <- list(data = DT_2,name = substr(temp[i],1,13), season = "second")
  print(head(test_2))
  data_LS<- append(data_LS,list(DT_LS_1))
  data_LS<- append(data_LS,list(DT_LS_2))
}

(team <- sort(unique(data_LS[[2]]$data[,HomeTeam,])))
#[1] "Arsenal"        "Aston Villa"    "Cardiff"        "Chelsea"       
#[5] "Crystal Palace" "Everton"        "Fulham"         "Hull"          
#[9] "Liverpool"      "Man City"       "Man United"     "Newcastle"     
#[13] "Norwich"        "Southampton"    "Stoke"          "Sunderland"    
#[17] "Swansea"        "Tottenham"      "West Brom"      "West Ham"
data_LS[[2]]$name
#[1] "EPL 2013-2014"
market_value <- c(8.99, 3.42, 2.37, 12.15,
                  1.57, 7.25, 3.37, 2.12,
                  8.79, 17.74, 15.05, 5.72,
                  3.30, 3.31, 4.34, 4.09,
                  3.83, 8.22, 3.07, 3.05)

data_LS[[1]]$vertices_DT <- data.table(team,market_value) 
data_LS[[2]]$vertices_DT <- data.table(team,market_value)

(team <- sort(unique(data_LS[[3]]$data[,HomeTeam,])))
#[1] "Arsenal"        "Aston Villa"    "Burnley"        "Chelsea"       
#[5] "Crystal Palace" "Everton"        "Hull"           "Leicester"     
#[9] "Liverpool"      "Man City"       "Man United"     "Newcastle"     
#[13] "QPR"            "Southampton"    "Stoke"          "Sunderland"    
#[17] "Swansea"        "Tottenham"      "West Brom"      "West Ham"
data_LS[[3]]$name
#[1] "EPL 2014-2015"
market_value <- c(9.63, 3.34, 1.23, 13.84,
                  2.13, 4.44, 3.49, 1.51,
                  8.14, 13.28, 12.15, 3.33,
                  2.81, 4.06, 4.30, 3.39,
                  2.99, 7.10, 2.96, 3.42)

data_LS[[3]]$vertices_DT <- data.table(team,market_value) 
data_LS[[4]]$vertices_DT <- data.table(team,market_value)

(team <- sort(unique(data_LS[[5]]$data[,HomeTeam,])))
#[1] "Arsenal"        "Aston Villa"    "Bournemouth"    "Chelsea"       
#[5] "Crystal Palace" "Everton"        "Leicester"      "Liverpool"     
#[9] "Man City"       "Man United"     "Newcastle"      "Norwich"       
#[13] "Southampton"    "Stoke"          "Sunderland"     "Swansea"       
#[17] "Tottenham"      "Watford"        "West Brom"      "West Ham" 
data_LS[[5]]$name
#[1] "EPL 2015-2016"
market_value <- c(9.97, 2.55, 1.39, 14.50,
                  3.18, 5.23, 2.61, 7.39,
                  11.91, 8.50, 4.78, 2.53,
                  5.36, 4.10, 3.05, 3.84,
                  8.19, 2.82, 3.34, 3.70)

data_LS[[5]]$vertices_DT <- data.table(team,market_value) 
data_LS[[6]]$vertices_DT <- data.table(team,market_value)

#Setting up color DT
max_goal <- 8
color <- c("white",brewer.pal(max_goal,"Greens"),"white",brewer.pal(max_goal,"Reds"))
weight <- rep(0:max_goal,2)
home <- c(rep(1,max_goal+1),rep(0,max_goal+1))
color_DT <- data.table(color,weight,home)

for(i in 1:6){
  DT <- data_LS[[i]]$data
  edges_DT <- rbind(DT[,.(from = HomeTeam,
                          to = AwayTeam, 
                          weight = FTHG,
                          home=1),],
                    DT[,.(from = AwayTeam,
                          to = HomeTeam, 
                          weight = FTAG,
                          home=0),])
  
  edges_DT <- merge(edges_DT,color_DT,by.x=c("weight","home"),by.y=c("weight","home"))
  setcolorder(edges_DT,c("from","to","weight","home","color"))
  data_LS[[i]]$edges_DT <- edges_DT
  data_LS[[i]]$igraph <- graph_from_data_frame(d = edges_DT,
                                               directed = T,
                                               vertices = data_LS[[i]]$vertices_DT)
}

saveRDS(data_LS,paste(getwd(),"/processed data/data_LS",sep=""))


