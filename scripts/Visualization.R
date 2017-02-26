library(data.table)
library(ggthemes)
library(ggplot2)
library(latex2exp)
library(rprojroot)
library(igraph)
library(RColorBrewer)

library(gridExtra)
#setup
root <- find_root(is_rstudio_project)


#Geometric and Poisson
n <- 10
probability <- c(dgeom(0:n,2/5),dpois(0:n,3/2))
type <- c(rep("Geometric",n+1),rep("Poisson",n+1))
y_ij <- c(rep(0:n,2))
graph_DT <- data.table(probability,type,y_ij)
geom_pois <- ggplot(data=graph_DT)+
  geom_bar(aes(x=as.factor(y_ij),y=probability),stat="identity",width=1,color="white")+
  geom_rangeframe(aes(x=as.factor(y_ij),y=probability))+
  facet_grid(.~type,scales = "free")+
  theme_tufte(base_size = 15)+
  labs(x=expression(y["ij"]),y=expression(paste("Pr (",Y["ij"]," = ",y["ij"]," | ",mu," = ","3/2",")" )))

saveRDS(geom_pois, paste(root,"/plots/geom_pois",sep=""))

#Zero Modified Poisson
n <- 10
x_grid <- 3
y_grid <- 3
y_ij <- rep(rep(0:n,x_grid),y_grid)
theta_1 <- NULL
for(i in seq(0,1,length.out =y_grid)){theta_1 <- c(theta_1,rep(i,(n+1)*x_grid))}
theta_2 <- NULL
for(i in seq(-1,1,length.out =x_grid)){theta_2 <- c(theta_2,rep(i,(n+1)))}
theta_2 <- rep(theta_2,y_grid)

graph_DT <- data.table(y_ij,theta_1,theta_2)
graph_DT[,probability:=exp(theta_1*y_ij+theta_2*(y_ij>0))/factorial(y_ij)]
graph_DT[,probability:=probability/sum(probability),by=.(theta_1,theta_2)]
graph_DT[,":="(y_ij=as.factor(y_ij),
               theta_1=as.factor(theta_1),
               theta_2=as.factor(theta_2)),]

for(i in 1:length(levels(graph_DT$theta_1))){
  levels(graph_DT$theta_1)[i] <- TeX(paste("$\\theta_{1} = ",levels(graph_DT$theta_1)[i],"$"))
}
for(i in 1:length(levels(graph_DT$theta_2))){
  levels(graph_DT$theta_2)[i] <- TeX(paste("$\\theta_{2} = ",levels(graph_DT$theta_2)[i],"$"))
}

zero_modified <- ggplot(data=graph_DT)+
  geom_bar(aes(x=y_ij,y=probability),stat="identity",width=1,color="white")+
  geom_bar(data=graph_DT[y_ij==0],aes(x=y_ij,y=probability),
           stat="identity",width=1,color="white",fill="red")+
  geom_rangeframe(data=graph_DT[probability>=0.001],aes(x=y_ij,y=probability))+
  facet_grid(theta_1~theta_2,labeller = label_parsed)+
  labs(x=TeX("$y_{ij}$"), y= TeX("$Pr(Y_{ij}=y_{ij}|\\theta_{1},\\theta_{2})$"))+
  theme_tufte()

saveRDS(zero_modified,paste(root,"/plots/zero_modified",sep=""))


# Conway-Maxwell-Poisson
n <- 20
x_grid <- 3
y_grid <- 3
y_ij <- rep(rep(0:n,x_grid),y_grid)
theta_1 <- NULL
for(i in seq(0,1,length.out =y_grid)){theta_1 <- c(theta_1,rep(i,(n+1)*x_grid))}
theta_2 <- NULL
for(i in seq(-0.5,0.5,length.out =x_grid)){theta_2 <- c(theta_2,rep(i,(n+1)))}
theta_2 <- rep(theta_2,y_grid)

graph_DT <- data.table(y_ij,theta_1,theta_2)
graph_DT[,probability:=exp(theta_1*y_ij)*(factorial(y_ij))^(theta_2-1)]
graph_DT[,probability:=probability/sum(probability),by=.(theta_1,theta_2)]
graph_DT[,":="(theta_1=as.factor(theta_1),
               theta_2=as.factor(theta_2)),]

for(i in 1:length(levels(graph_DT$theta_1))){
  levels(graph_DT$theta_1)[i] <- TeX(paste("$\\theta_{1} = ",levels(graph_DT$theta_1)[i],"$"))
}
for(i in 1:length(levels(graph_DT$theta_2))){
  levels(graph_DT$theta_2)[i] <- TeX(paste("$\\theta_{2} = ",levels(graph_DT$theta_2)[i],"$"))
}

cmp <- ggplot(data=graph_DT)+
  geom_bar(aes(x=y_ij,y=probability),stat="identity",width=1,color="white")+
  geom_rangeframe(data=graph_DT[probability>=0.001],aes(x=y_ij,y=probability))+
  facet_grid(theta_1~theta_2,labeller = label_parsed)+
  labs(x=TeX("$y_{ij}$"), y= TeX("$Pr(Y_{ij}=y_{ij}|\\theta_{1},\\theta_{2})$"))+
  theme_tufte()

saveRDS(cmp,paste(root,"/plots/cmp",sep=""))

#comparison of mutuality
n <- 20
x_grid <- 3
y_ij <- rep(rep(0:n,n+1),x_grid)
y_ji <- NULL
for(i in 0:n){y_ji <- c(y_ji,rep(i,n+1))}

y_ji <- rep(y_ji,x_grid)
theta_1 <- rep(1, (n+1)^2*x_grid)
theta_2 <- NULL
for(i in seq(-1,1,length.out = x_grid)){theta_2 <- c(theta_2,rep(i,(n+1)^2))}

graph_DT <- data.table(y_ij,y_ji,theta_1,theta_2)
graph_DT[,minimum:= ifelse(y_ij <= y_ji,y_ij,y_ji),]

statistics <- NULL
for(i in c("geomean","min","nabsdiff")){statistics <- c(statistics,rep(i,(n+1)^2*x_grid)) }

graph_DT <- data.table(graph_DT,statistics)

graph_DT[statistics=="min",probability:= exp(theta_1*(y_ij+y_ji)+theta_2*minimum) /
           (factorial(y_ij)*factorial(y_ji)),]
graph_DT[statistics=="min",probability:= probability/sum(probability),by = .(theta_2)]

graph_DT[statistics=="geomean",probability:= exp(theta_1*(y_ij+y_ji)+theta_2*sqrt(y_ij)*sqrt(y_ji)) /
           (factorial(y_ij)*factorial(y_ji)),]
graph_DT[statistics=="geomean",probability:= probability/sum(probability),by=.(theta_2)]

graph_DT[statistics=="nabsdiff",probability:= exp(theta_1*(y_ij+y_ji)-theta_2*abs(y_ij-y_ji)) /
           (factorial(y_ij)*factorial(y_ji)),]
graph_DT[statistics=="nabsdiff",probability:= probability/sum(probability),by=.(theta_2)]

graph_DT[,":="(theta_2=as.factor(theta_2)),]

for(i in 1:length(levels(graph_DT$theta_2))){
  levels(graph_DT$theta_2)[i] <- TeX(paste("$\\theta_{2} = ",levels(graph_DT$theta_2)[i],"$"))
}

mutual<- ggplot(data=graph_DT)+
    geom_tile(aes(x=y_ij,y=y_ji,fill=probability),colour="white") +
    geom_rangeframe(aes(x=y_ij,y=y_ji))+
    theme_tufte(base_size = 15)+
    facet_grid(statistics~theta_2,labeller = label_parsed)+
    labs(x=TeX("$y_{ij}$"),y=TeX("$y_{ji}$"),fill=TeX("$Pr(Y_{ij},Y_{ji}|\\theta_{1}=1,\\theta_{2})$"))+
    scale_fill_gradientn(colours=c("white","black"))+
    coord_equal()+
    guides(fill = guide_legend(
      title.position = "left",
      title.theme = element_text(angle = 90)))

saveRDS(mutual,paste(root,"/plots/mutual",sep=""))

data_LS <- readRDS(paste(root,"/processed data/data_LS",sep=""))

#1st plot
i<-1
graph_DT <- rbind(data.table(data_LS[[i*2]]$edges_DT),data.table(data_LS[[i*2-1]]$edges_DT))
graph_DT[,home:= ifelse(home==1,"home","away")]
graph_DT[,weight:=as.character(weight)]
graph_DT[,color:=NULL,]
  
amv_DT <- data.table(data_LS[[i*2]]$vertices_DT)
  
na_DT <- data.table(from = amv_DT$team, 
                    to = amv_DT$team, 
                    weight ="na" , 
                    home = rep(c("home","away"),each = nrow(amv_DT)))
graph_DT <- rbind(graph_DT,na_DT)
  
setkey(amv_DT,"market_value")
  
graph_DT$home <- factor(graph_DT$home,levels=c("home","away"))
graph_DT$from <-  factor(graph_DT$from, levels = amv_DT$team)
graph_DT$to <- factor(graph_DT$to,levels = amv_DT$team)
amv_DT$team <- factor(amv_DT$team,levels = amv_DT$team)
  
  
heatmap_1_gg <- ggplot(graph_DT)+
  geom_tile(aes(y=from,x=to,fill=weight),color="black")+
  xlim(rev(levels(graph_DT$to)))+
  scale_fill_manual(values = c("na" = "black", "0" = "white",
                               "1" = '#f7fbff',"2" = '#deebf7',"3" = '#c6dbef',"4" = '#9ecae1',
                               "5" = '#6baed6',"6" = '#4292c6',"7" = '#2171b5',"8" = '#084594'))+
  theme_tufte()+
  theme(axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  facet_grid(.~home)
  
barplot_1_gg <- ggplot(amv_DT)+
  geom_col(aes(x = team, y = market_value),fill = "black",color = "white")+
  geom_text(data = amv_DT[market_value>=3],
            aes(x = team, y= market_value-1.5,label = market_value),color="white",size=3)+
  geom_text(data = amv_DT[market_value<3],
            aes(x = team, y= market_value+1.5,label = market_value),color="black",size=3)+
  coord_fixed()+
  coord_flip()+
  scale_y_continuous(position = "top",limits = c(0,20))+
  theme_tufte()+
  theme(axis.ticks = element_blank(),
        axis.text =  element_blank(),
        axis.title.y = element_blank())

#2nd plot
i<-2
graph_DT <- rbind(data.table(data_LS[[i*2]]$edges_DT),data.table(data_LS[[i*2-1]]$edges_DT))
graph_DT[,home:= ifelse(home==1,"home","away")]
graph_DT[,weight:=as.character(weight)]
graph_DT[,color:=NULL,]

amv_DT <- data.table(data_LS[[i*2]]$vertices_DT)

na_DT <- data.table(from = amv_DT$team, 
                    to = amv_DT$team, 
                    weight ="na" , 
                    home = rep(c("home","away"),each = nrow(amv_DT)))
graph_DT <- rbind(graph_DT,na_DT)

setkey(amv_DT,"market_value")

graph_DT$home <- factor(graph_DT$home,levels=c("home","away"))
graph_DT$from <-  factor(graph_DT$from, levels = amv_DT$team)
graph_DT$to <- factor(graph_DT$to,levels = amv_DT$team)
amv_DT$team <- factor(amv_DT$team,levels = amv_DT$team)


heatmap_2_gg <- ggplot(graph_DT)+
  geom_tile(aes(y=from,x=to,fill=weight),color="black")+
  xlim(rev(levels(graph_DT$to)))+
  scale_fill_manual(values = c("na"= "black", "0"= "white",
                               "1"= '#f7fcf5', "2"= '#e5f5e0', "3"= '#c7e9c0',"4"= '#a1d99b',
                               "5"= '#74c476',"6"='#41ab5d',"7"= '#238b45',"8"='#005a32'))+
  theme_tufte()+
  theme(axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  facet_grid(.~home)

barplot_2_gg <- ggplot(amv_DT)+
  geom_col(aes(x = team, y = market_value),fill = "black",color = "white")+
  geom_text(data = amv_DT[market_value>=3],
            aes(x = team, y= market_value-1.5,label = market_value),color="white",size=3)+
  geom_text(data = amv_DT[market_value<3],
            aes(x = team, y= market_value+1.5,label = market_value),color="black",size=3)+
  coord_fixed()+
  coord_flip()+
  scale_y_continuous(position = "top",limits = c(0,20))+
  theme_tufte()+
  theme(axis.ticks = element_blank(),
        axis.text =  element_blank(),
        axis.title.y = element_blank())

#3rd plot
i<-3
graph_DT <- rbind(data.table(data_LS[[i*2]]$edges_DT),data.table(data_LS[[i*2-1]]$edges_DT))
graph_DT[,home:= ifelse(home==1,"home","away")]
graph_DT[,weight:=as.character(weight)]
graph_DT[,color:=NULL,]

amv_DT <- data.table(data_LS[[i*2]]$vertices_DT)

na_DT <- data.table(from = amv_DT$team, 
                    to = amv_DT$team, 
                    weight ="na" , 
                    home = rep(c("home","away"),each = nrow(amv_DT)))
graph_DT <- rbind(graph_DT,na_DT)

setkey(amv_DT,"market_value")

graph_DT$home <- factor(graph_DT$home,levels=c("home","away"))
graph_DT$from <-  factor(graph_DT$from, levels = amv_DT$team)
graph_DT$to <- factor(graph_DT$to,levels = amv_DT$team)
amv_DT$team <- factor(amv_DT$team,levels = amv_DT$team)


heatmap_3_gg <- ggplot(graph_DT)+
  geom_tile(aes(y=from,x=to,fill=weight),color="black")+
  xlim(rev(levels(graph_DT$to)))+
  scale_fill_manual(values = c("na"= "black", "0"= "white",
                               "1"= '#fff5f0',"2"= '#fee0d2',"3"= '#fcbba1',"4"= '#fc9272',
                               "5"= '#fb6a4a',"6"= '#ef3b2c',"7"= '#cb181d',"8" ='#99000d'))+
  theme_tufte()+
  theme(axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  facet_grid(.~home)

barplot_3_gg <- ggplot(amv_DT)+
  geom_col(aes(x = team, y = market_value),fill = "black",color = "white")+
  geom_text(data = amv_DT[market_value>=3],
            aes(x = team, y= market_value-1.5,label = market_value),color="white",size=3)+
  geom_text(data = amv_DT[market_value<3],
            aes(x = team, y= market_value+1.5,label = market_value),color="black",size=3)+
  coord_fixed()+
  coord_flip()+
  scale_y_continuous(position = "top",limits = c(0,20))+
  theme_tufte()+
  theme(axis.ticks = element_blank(),
        axis.text =  element_blank(),
        axis.title.y = element_blank())

  
grid.arrange(heatmap_1_gg, barplot_1_gg, 
             heatmap_2_gg, barplot_2_gg,
             heatmap_3_gg, barplot_3_gg,
             layout_matrix = rbind(c(1,1,1,2),c(3,3,3,4),c(5,5,5,6)))
  