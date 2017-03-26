library(ggthemes)
library(ggplot2)
library(latex2exp)
library(rprojroot)
library(igraph)
library(RColorBrewer)
library(gridExtra)
library(data.table)

#setup
root <- find_root(is_rstudio_project)
data_LS <- readRDS(paste(root,"/processed data/data_LS",sep=""))
coef_LS <- readRDS(paste(root,"/processed data/coef_LS",sep=""))

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

#plot of EPL 2013-14
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
amv_DT$team <- factor(amv_DT$team,levels = rev(amv_DT$team))
  
colormap <- c("na" = "black", "0" = "white",
              "1" = '#f7fbff',"2" = '#deebf7',"3" = '#c6dbef',"4" = '#9ecae1',
              "5" = '#6baed6',"6" = '#4292c6',"7" = '#2171b5',"8" = '#084594')

legend_DT <- data.table(from = "foo", to=names(colormap))

legend_1_gg <- ggplot(legend_DT)+
  geom_tile(aes(y=from,x=to,fill=to),color="black")+
  scale_fill_manual(values=colormap)+
  coord_equal()+
  theme_tufte()+
  theme(legend.position = "none",
        axis.title= element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=40))

heatmap_1_gg <- ggplot(graph_DT)+
  geom_tile(aes(y=from,x=to,fill=weight),color="black")+
  xlim(rev(levels(graph_DT$to)))+
  scale_fill_manual(values = colormap)+
  theme_tufte()+
  coord_equal()+
  theme(axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  facet_grid(.~home,switch="x")


barplot_1_gg <- ggplot(amv_DT)+
  geom_col(aes(x = team, y = market_value),fill = "black",color = "white")+
  scale_y_continuous(limits = c(0,20))+
  geom_hline(data =data.table(y=seq(0,20,5)) ,aes(yintercept=y),col="white")+
  theme_tufte(ticks=F)+
  theme(axis.ticks = element_blank(),
        axis.text.x =  element_blank(),
        axis.title = element_blank(),
        text = element_text(size=20))

#Plot of EPL 2014-15
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
amv_DT$team <- factor(amv_DT$team,levels = rev(amv_DT$team))

colormap <- c("na"= "black", "0"= "white",
              "1"= '#f7fcf5', "2"= '#e5f5e0', "3"= '#c7e9c0',"4"= '#a1d99b',
              "5"= '#74c476',"6"='#41ab5d',"7"= '#238b45',"8"='#005a32')

legend_DT <- data.table(from = "foo", to=names(colormap))

legend_2_gg <- ggplot(legend_DT)+
  geom_tile(aes(y=from,x=to,fill=to),color="black")+
  scale_fill_manual(values=colormap)+
  coord_equal()+
  theme_tufte()+
  theme(legend.position = "none",
        axis.title= element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=40))

heatmap_2_gg <- ggplot(graph_DT)+
  geom_tile(aes(y=from,x=to,fill=weight),color="black")+
  xlim(rev(levels(graph_DT$to)))+
  scale_fill_manual(values = colormap)+
  theme_tufte()+
  coord_equal()+
  theme(axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  facet_grid(.~home,switch="x")


barplot_2_gg <- ggplot(amv_DT)+
  geom_col(aes(x = team, y = market_value),fill = "black",color = "white")+
  scale_y_continuous(limits = c(0,20))+
  geom_hline(data =data.table(y=seq(0,20,5)) ,aes(yintercept=y),col="white")+
  theme_tufte(ticks=F)+
  theme(axis.ticks = element_blank(),
        axis.text.x =  element_blank(),
        axis.title = element_blank(),
        text = element_text(size=20))

#Plot of EPL 2015-16
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
amv_DT$team <- factor(amv_DT$team,levels = rev(amv_DT$team))

colormap <- c("na"= "black", "0"= "white",
              "1"= '#fff5f0',"2"= '#fee0d2',"3"= '#fcbba1',"4"= '#fc9272',
              "5"= '#fb6a4a',"6"= '#ef3b2c',"7"= '#cb181d',"8" ='#99000d')

legend_DT <- data.table(from = "foo", to=names(colormap))

legend_3_gg <- ggplot(legend_DT)+
  geom_tile(aes(y=from,x=to,fill=to),color="black")+
  scale_fill_manual(values=colormap)+
  coord_equal()+
  theme_tufte()+
  theme(legend.position = "none",
        axis.title= element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=40))

heatmap_3_gg <- ggplot(graph_DT)+
  geom_tile(aes(y=from,x=to,fill=weight),color="black")+
  xlim(rev(levels(graph_DT$to)))+
  scale_fill_manual(values = colormap)+
  theme_tufte()+
  coord_equal()+
  theme(axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  facet_grid(.~home,switch="x")


barplot_3_gg <- ggplot(amv_DT)+
  geom_col(aes(x = team, y = market_value),fill = "black",color = "white")+
  scale_y_continuous(limits = c(0,20))+
  geom_hline(data =data.table(y=seq(0,20,5)) ,aes(yintercept=y),col="white")+
  theme_tufte(ticks=F)+
  theme(axis.ticks = element_blank(),
        axis.text.x =  element_blank(),
        axis.title = element_blank(),
        text = element_text(size=20))


#Network representation
library(igraph)
library(scales)
library(RColorBrewer)

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
par(mfrow = c(1,1),mar=c(3, 3, 3, 3))
data_igraph <- data_LS[[2]]$igraph

vertex_order <- order(V(data_igraph)$market_value,decreasing=T)
lab_locs <- radian.rescale(order(vertex_order), direction=-1, start=0)

plot(data_igraph, 
     layout=layout_in_circle(data_igraph, order = vertex_order),
     vertex.label.dist=0.8,
     vertex.label.degree=lab_locs,
     vertex.size=sqrt(V(data_igraph)$market_value)*15/max(sqrt(V(data_igraph)$market_value)),
     vertex.frame.color="NA",
     vertex.color="#081d58",
     edge.arrow.size=.5,
     edge.curved = 0.1,
     edge.lty =1,
     edge.width=0.5*E(data_igraph)$weight)

legend(x=1.2, y=1.3, legend=0:8,lwd=0.5*(0:8), cex =1,
       col = c("white",brewer.pal(8,"Greens")),
       y.intersp = .5,bty = "n",title="Goals(Home)")

legend(x=1.2, y=1.3, legend=0:8,lwd=0.5*(0:8), cex = 1,
       col = c("white",brewer.pal(8,"Reds")),
       y.intersp = .5,bty ="n",title= "Goals(Away)")
#making the legend of the circle is hard



par(mfrow = c(1,1),mar=c(3, 3, 3, 3))
data_igraph <- data_LS[[3]]$igraph

vertex_order <- order(V(data_igraph)$market_value,decreasing=T)
lab_locs <- radian.rescale(order(vertex_order), direction=-1, start=0)

plot(data_igraph, 
     layout=layout_in_circle(data_igraph, order = vertex_order),
     vertex.label.dist=0.8,
     vertex.label.degree=lab_locs,
     vertex.size=sqrt(V(data_igraph)$market_value)*15/max(sqrt(V(data_igraph)$market_value)),
     vertex.frame.color="NA",
     vertex.color="#081d58",
     edge.arrow.size=.5,
     edge.curved = 0.1,
     edge.lty =1,
     edge.width=0.5*E(data_igraph)$weight)


###Independent edgewise Visualization

max_goal <- 10

#poisson
get_p_geometric<- function(SUM){
  p<- exp(SUM*(0:max_goal))
  return(p/sum(p))
}

get_p_poisson <- function(SUM){
  p <- exp(SUM*(0:max_goal))/factorial(0:10)
  return(p/sum(p))
}

get_p_nonzero <- function(SUM,nonzero){
  p <- exp(SUM*(0:max_goal)+nonzero*(0:10 >0))/factorial(0:10)
  return(p/sum(p))
}

get_p_CMP <- function(SUM, CMP){
  p <- exp(SUM*(0:max_goal))*factorial(0:10)^(CMP-1)
  return(p/sum(p))
}

graph_DT <-data.table()

#poisson
#model #data #goal #p #goal #coef1 #coef2
for(i in 1:6){
  model <- "geometric"
  data <- i
  coef_1 <- coef_LS$geometric_DT$Estimate[i]
  coef_2 <- NA
  goal <- 0:max_goal
  p <- get_p_geometric(coef_1)
  graph_DT <- rbind(graph_DT,data.table(model,data,coef_1,coef_2,goal,p))
}
for(i in 1:6){
  model <- "poisson"
  data <- i
  coef_1 <- coef_LS$poisson_DT$Estimate[i]
  coef_2 <- NA
  goal <- 0:max_goal
  p <- get_p_poisson(coef_1)
  graph_DT <- rbind(graph_DT,data.table(model,data,coef_1,coef_2,goal,p))
}
for(i in 1:6){
  model <- "nonzero"
  data <- i
  coef_1 <- coef_LS$nonzero_DT$Estimate[i*2-1]
  coef_2 <- coef_LS$nonzero_DT$Estimate[i*2]
  goal <- 0:max_goal
  p <- get_p_nonzero(coef_1,coef_2)
  graph_DT <- rbind(graph_DT,data.table(model,data,coef_1,coef_2,goal,p))
}
for(i in 1:6){
  model <- "cmp"
  data <- i
  coef_1 <- coef_LS$cmp_DT$Estimate[i*2-1]
  coef_2 <- coef_LS$cmp_DT$Estimate[i*2]
  goal <- 0:max_goal
  p <- get_p_CMP(coef_1,coef_2)
  graph_DT <- rbind(graph_DT,data.table(model,data,coef_1,coef_2,goal,p))
}

col_DT <- data.table()
for(i in 1:length(data_LS)){
  foo_DT <- data.table(data = i,data_LS[[i]]$edges_DT)
  col_DT <- rbind(col_DT,foo_DT)
}

col_DT <- col_DT[,.(count=.N),by=.(data,weight)]
col_DT[,total:=sum(count),by=data]
col_DT[,p:=count/total ,]
graph_DT[,goal:=goal+0.5]

edgewise <- ggplot()+
  geom_col(data=col_DT,aes(x=as.factor(weight),y=p),width=1,color="white")+
  geom_step(data=graph_DT,aes(x=goal,y=p,group = model,color = model),size=1.2)+
  geom_rangeframe(data=graph_DT,aes(x=goal,y=p))+
  facet_wrap(~data,ncol=2)+
  theme_tufte()+
  theme(legend.position="top")+
  labs(x=TeX("$y_{ij}$"), y= TeX("$Pr(Y_{ij}=y_{ij})$"))
  
saveRDS(edgewise, paste(root,"/plots/edgewise",sep=""))

