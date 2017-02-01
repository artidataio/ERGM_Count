library(data.table)
library(ggthemes)
library(ggplot2)
library(latex2exp)
#comparison of mutuality

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

saveRDS(geom_pois, paste(getwd(),"/plots/geom_pois",sep=""))
