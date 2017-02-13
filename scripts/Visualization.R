library(data.table)
library(ggthemes)
library(ggplot2)
library(latex2exp)
library(rprojroot)

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
  theme_tufte()+
  labs(x=expression(y["ij"]),y=expression(paste("Pr (",Y["ij"]," = ",y["ij"]," | ",mu," = ","3/2",")" )))

saveRDS(geom_pois, paste(root,"/plots/geom_pois",sep=""))

#Zero Modified Poisson
n <- 10
x_grid <- 5
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
x_grid <- 5
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
n <- 30
x_grid <- 5
y_ij <- rep(rep(0:n,n+1),x_grid)
y_ji <- NULL
for(i in 0:n){y_ji <- c(y_ji,rep(i,n+1))}

y_ji <- rep(y_ji,x_grid)
theta_1 <- rep(1, (n+1)^2*x_grid)
theta_2 <- NULL
for(i in seq(-2,2,length.out = x_grid)){theta_2 <- c(theta_2,rep(i,(n+1)^2))}

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
    geom_raster(aes(x=y_ij,y=y_ji,fill=probability)) +
    geom_rangeframe(aes(x=y_ij,y=y_ji))+
    theme_tufte()+
    facet_grid(statistics~theta_2,labeller = label_parsed)+
    labs(x=TeX("$y_{ij}$"),y=TeX("$y_{ji}$"),fill=TeX("$Pr(Y_{ij},Y_{ji}|\\theta_{1}=1,\\theta_{2})$"))+
    scale_fill_gradientn(colours=c("white","black"))+
    guides(fill = guide_legend(
      title.position = "left",
      title.theme = element_text(angle = 90)))

saveRDS(mutual,paste(root,"/plots/mutual",sep=""))
