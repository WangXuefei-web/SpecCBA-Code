# SpecCBA: Real data application
library(rTensor)
library(ggplot2)
library(ggstar)
library(reshape2)
library(Matrix)
library(tidyverse)
library(plotly)
library(clue)
#library(RSpectra)
source('SpecCovCBA_fun.R')
library(maps)
library(geosphere)
library(ggpubr)

load(file="data/FAO_data.Rdata")
load(file="FAO_SpecCBA_result.Rdata")
set.seed(123)
### Application 1: 94 countries and 19 kinds of processed food ###
# Apply SpecCBA
N <- processed.arrT@modes[1]
K <- 4
lb.re <- SpecCBA.app(processed.arrT,alpha.seq=seq(-1,5,0.2),K)
M.mat <- SpecCBA.A(processed.arrT,lb.re$alpha,K)$M #symmetric

# Draw map
world_map <- map_data("world")
# trade_map <- country_map %>% mutate(Community=as.character(lb.re$label.hat))
trade_map <- processed.map #saved result

# Connect map
processed.S3 <- modeSum(processed.arrT,3,drop=TRUE)@data
S.scale <- processed.S3/processed.arrT@modes[3]
S.scale[lower.tri(S.scale,diag=TRUE)] <- 0
S.dfup <- as.data.frame(S.scale)
colnames(S.dfup) <- c(1:N)
#med_val <- median(S.scale[upper.tri(S.scale,diag=FALSE)])
S.df <- S.dfup %>% mutate(row_ID=c(1:N)) %>% gather(key="col_ID",value="value", -row_ID) %>% filter(value>0)
S.df$col_ID <- as.numeric(S.df$col_ID)


# M.scale <- (M.mat-min(M.mat)+1e-6)/(max(M.mat)-min(M.mat))
# M.scale[lower.tri(M.scale,diag=TRUE)] <- 0
# M.dfup <- as.data.frame(M.scale)
# colnames(M.dfup) <- c(1:N)
# med_val <- median(M.scale[upper.tri(M.scale,diag=FALSE)])
# M.df <- M.dfup %>% mutate(row_ID=c(1:N)) %>% gather(key="col_ID",value="value", -row_ID) %>% filter(value>med_val)
# M.df$col_ID <- as.numeric(M.df$col_ID)

trade.df <- trade_map %>% select(long,lat) %>% mutate(row_ID=c(1:N),col_ID=c(1:N))

connect_map <- S.df %>% left_join(trade.df,by="row_ID") %>% select(value,long1=long,lat1=lat,col_ID=col_ID.x) %>% 
  left_join(trade.df,by="col_ID") %>% mutate(alpha=value/5) %>% select(alpha,value,long1,lat1,long2=long,lat2=lat)

# c.exchange <- function(x){
#   if(x==1)
#   {return(4)}
#   else if(x==4)
#   {return(1)}
#   else
#   {return(x)}
# }
# 
# trade_map <- trade_map %>% mutate(C1=map_int(Community,function(x){as.numeric(x)})) %>% 
#   mutate(Community=as.character(map_int(C1,c.exchange))) %>% select(-C1)
#processed.map <- trade_map

png("result/real_processed.png", height=1000,width=2000, res=168)

ggplot() + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),col="#f2f2f2",bg="white",lwd=0.05) +
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())+
  geom_curve(data=connect_map,aes(x=long1,y=lat1,xend=long2,yend=lat2,alpha=alpha),color="lightblue",linewidth=0.001,show.legend = FALSE)+#curvature=0
  geom_point(data=trade_map,aes(x = long, y = lat, group=Community, color=Community,shape = Community),size=2)+
  scale_shape_manual(values = c(15:18,25))+
  scale_color_brewer(palette = 'Set1')+
  geom_text(data=trade_map,aes(x = long, y = lat, label=region),color="black",vjust=1,size=2)+
  theme(legend.position = c(0.1,0.3))

dev.off()




### Application 2: 94 countries and 11 kinds of unprocessed food ###
# Apply SpecCBA
N <- unprocessed.arrT@modes[1]
K <- 4
lb.re <- SpecCBA.app(unprocessed.arrT,alpha.seq=seq(-1,5,0.2),K)
M.mat <- SpecCBA.A(unprocessed.arrT,lb.re$alpha,K)$M #symmetric

# Draw map
world_map <- map_data("world")
# trade_map <- country_map %>% mutate(Community=as.character(lb.re$label.hat))
trade_map <- unprocessed.map #saved result

# Connect map
unprocessed.S3 <- modeSum(unprocessed.arrT,3,drop=TRUE)@data

S.scale <- unprocessed.S3/(unprocessed.arrT@modes[3])
S.scale[lower.tri(S.scale,diag=TRUE)] <- 0
S.dfup <- as.data.frame(S.scale)
colnames(S.dfup) <- c(1:N)
#med_val <- median(S.scale[upper.tri(S.scale,diag=FALSE)])
S.df <- S.dfup %>% mutate(row_ID=c(1:N)) %>% gather(key="col_ID",value="value", -row_ID) %>% filter(value>0)
S.df$col_ID <- as.numeric(S.df$col_ID)

trade.df <- trade_map %>% select(long,lat) %>% mutate(row_ID=c(1:N),col_ID=c(1:N))

connect_map <- S.df %>% left_join(trade.df,by="row_ID") %>% select(value,long1=long,lat1=lat,col_ID=col_ID.x) %>% 
  left_join(trade.df,by="col_ID") %>% mutate(alpha=value/20) %>% select(alpha,value,long1,lat1,long2=long,lat2=lat)

# c.exchange <- function(x){
#   if(x==3)
#     {return(1)}
#   else if(x==1)
#   {return(4)}
#   else if(x==4)
#   {return(3)}
#   else
#     {return(x)}
# }
# 
# trade_map <- trade_map %>% mutate(C1=map_int(Community,function(x){as.numeric(x)})) %>%
#   mutate(Community=as.character(map_int(C1,c.exchange))) %>% select(-C1)
# unprocessed.map <- trade_map

png("result/real_unprocessed.png", height=1000,width=2000, res=168)

ggplot() + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),col="#f2f2f2",bg="white",lwd=0.05) +
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())+
  geom_curve(data=connect_map,aes(x=long1,y=lat1,xend=long2,yend=lat2,alpha=alpha),color="lightblue",linewidth=0.001,show.legend = FALSE)+#curvature=0
  geom_point(data=trade_map,aes(x = long, y = lat, group=Community, color=Community,shape = Community),size=2)+
  scale_shape_manual(values = c(15:18,25))+
  scale_color_brewer(palette = 'Set1')+
  geom_text(data=trade_map,aes(x = long, y = lat, label=region),color="black",vjust=1,size=2)+
  theme(legend.position = c(0.1,0.3))

dev.off()

save(processed.map,unprocessed.map,core.map,food_diff,file="FAO_SpecCBA_result.Rdata")



### Analysis ###
food_diff <- processed.map %>% select(nodeID, proC=Community) %>% left_join(unprocessed.map,by="nodeID") %>%
  select(nodeID,region,proC,unproC=Community) %>% mutate(change=map2_int(proC,unproC,function(x,y){as.numeric(x!=y)}))

draw_trade_heat <- function(trade_map,S.scale,showlegend=FALSE){
  N <- dim(S.scale)[1]
  S.scale[lower.tri(S.scale,diag=TRUE)] <- 0
  S.dfup <- as.data.frame(S.scale)
  colnames(S.dfup) <- c(1:N)
  #med_val <- median(S.scale[upper.tri(S.scale,diag=FALSE)])
  S.df <- S.dfup %>% mutate(row_ID=c(1:N)) %>% gather(key="col_ID",value="value", -row_ID) %>% filter(value>0)
  S.df$col_ID <- as.numeric(S.df$col_ID)
  
  
  # reorder nodes by label
  lb.df <- data.frame(label = as.numeric(trade_map$Community)) %>% mutate(row_ID=c(1:N),degree=rowMeans(S.scale)) %>%
    arrange(label,desc(degree)) %>% mutate(reID=c(1:N))
  
  sumAdj.map <- S.df %>% rename(col_ID=row_ID,row_ID=col_ID) %>% bind_rows(S.df) %>% 
    left_join(lb.df,by="row_ID") %>% rename(row_re = reID) %>% select(-label,-degree) %>% 
    left_join(lb.df,join_by(col_ID == row_ID)) %>% rename(col_re = reID) %>% select(-label,-degree)
  
  # add empty row to satisfies low="white"
  if(min(sumAdj.map$value)==max(sumAdj.map$value)){
    emptyrow <- sumAdj.map[1,] %>% mutate(col_ID=row_ID,col_re=row_re,value=0)
    sumAdj.map <- sumAdj.map %>% bind_rows(emptyrow)
  }
  
  # axis intercept
  count <- as.numeric(table(trade_map$Community))
  cut <- (cumsum(count)+0.5)[-length(count)]
  
  p <- ggplot(sumAdj.map,aes(x=row_re,y=-col_re,fill=value))+
    geom_tile(show.legend = showlegend)+
    scale_fill_gradient2(low="white",high="red")+
    coord_equal()+
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    theme_void()+
    theme(legend.direction="horizontal",legend.position = c(0,1),legend.key.size = unit(0.5,"cm"))+
    geom_vline(xintercept = cut,linewidth=0.8)+
    geom_hline(yintercept = -cut,linewidth=0.8)
  return(p)
}

## Heatmap for processed food 19
# S.scale = meanAdj
trade_map <- processed.map
processed.S3 <- modeSum(processed.arrT,3,drop=TRUE)@data
S.scale <- processed.S3/(processed.arrT@modes[3])
#filestr <- "processed_layers/MeanAdj"
#filename <- paste("result/",filestr, ".png", sep = "", collapse = NULL)

# draw
p <- draw_trade_heat(trade_map,S.scale,TRUE)
p.mean <- p+
     ggtitle("MeanAdj")+
     theme(plot.title = element_text(hjust = 0.5))
# print(p.mean)
# png(filename, height=600,width=700, res=168)
# pt <- p+
#   ggtitle("Mean of layer-wise adjacency matrices")+
#   theme(plot.title = element_text(hjust = 0.5))
# print(pt)
# dev.off()
plist <- list()
# S.scale = l-th layer of Adj
for (l in (1:processed.arrT@modes[3])){#
  
  S.scale <- processed.arrT@data[,,l]
  # filestr=paste("processed_layers/layer",as.character(l), sep = "", collapse = NULL)
  # filename <- paste("result/",filestr, ".png", sep = "", collapse = NULL)
  # draw
  p <- draw_trade_heat(trade_map,S.scale)
  pt <- p+
    ggtitle(as.character(processed_food$layerLabel[l]))+
    theme(plot.title = element_text(hjust = 0.5)) 
  plist[[l]] <-pt
  # png(filename, height=650,width=600, res=168)
  # print(pt)
  # dev.off()
}

parr <- ggarrange(p.mean,plist[[1]],plist[[2]],plist[[3]],
          plist[[4]],plist[[5]],plist[[6]],plist[[7]],
          plist[[8]],plist[[9]],plist[[10]],plist[[11]],
          plist[[12]],plist[[13]],plist[[14]],plist[[15]],
          plist[[16]],plist[[17]],plist[[18]],plist[[19]],
          ncol = 4, nrow = 5,common.legend = T)

png("result/layers_heatmap_processed.png", height=1500,width=1000, res=100)
print(parr)
dev.off()

## Heatmap for unprocessed food
trade_map <- unprocessed.map
unprocessed.S3 <- modeSum(unprocessed.arrT,3,drop=TRUE)@data
S.scale <- unprocessed.S3/(unprocessed.arrT@modes[3])
#filestr <- "unprocessed_layers/MeanAdj"
#filename <- paste("result/",filestr, ".png", sep = "", collapse = NULL)

# draw
p <- draw_trade_heat(trade_map,S.scale,TRUE)
p.mean <- p+
  ggtitle("MeanAdj")+
  theme(plot.title = element_text(hjust = 0.5))
# print(p.mean)
# png(filename, height=600,width=700, res=168)
# pt <- p+
#   ggtitle("Mean of layer-wise adjacency matrices")+
#   theme(plot.title = element_text(hjust = 0.5))
# print(pt)
# dev.off()
plist <- list()
# S.scale = l-th layer of Adj
for (l in (1:unprocessed.arrT@modes[3])){#
  
  S.scale <- unprocessed.arrT@data[,,l]
  # filestr=paste("unprocessed_layers/layer",as.character(l), sep = "", collapse = NULL)
  # filename <- paste("result/",filestr, ".png", sep = "", collapse = NULL)
  # draw
  p <- draw_trade_heat(trade_map,S.scale)
  pt <- p+
    ggtitle(as.character(unprocessed_food$layerLabel[l]))+
    theme(plot.title = element_text(hjust = 0.5)) 
  plist[[l]] <-pt
  # png(filename, height=650,width=600, res=168)
  # print(pt)
  # dev.off()
}

parr <- ggarrange(p.mean,plist[[1]],plist[[2]],plist[[3]],
                  plist[[4]],plist[[5]],plist[[6]],plist[[7]],
                  plist[[8]],plist[[9]],plist[[10]],plist[[11]],
                  ncol = 4, nrow = 3,common.legend = T)

png("result/layers_heatmap_unprocessed.png", height=900,width=1000, res=100)
print(parr)
dev.off()


### Whole 30 foods ###
load(file="data/FAO_data.Rdata")
load(file="FAO_SpecCBA_result.Rdata")
set.seed(123)

# Apply SpecCBA
N <- core.arrT@modes[1]
K <- 4
lb.re <- SpecCBA.app(core.arrT,alpha.seq=seq(-1,5,0.2),K)
M.mat <- SpecCBA.A(core.arrT,lb.re$alpha,K)$M #symmetric

# Draw map
world_map <- map_data("world")
# trade_map <- country_map %>% mutate(Community=as.character(lb.re$label.hat))
trade_map <- core.map #saved result

# Connect map
core.S3 <- modeSum(core.arrT,3,drop=TRUE)@data
S.scale <- core.S3/core.arrT@modes[3]
S.scale[lower.tri(S.scale,diag=TRUE)] <- 0
S.dfup <- as.data.frame(S.scale)
colnames(S.dfup) <- c(1:N)
#med_val <- median(S.scale[upper.tri(S.scale,diag=FALSE)])
S.df <- S.dfup %>% mutate(row_ID=c(1:N)) %>% gather(key="col_ID",value="value", -row_ID) %>% filter(value>0)
S.df$col_ID <- as.numeric(S.df$col_ID)


# M.scale <- (M.mat-min(M.mat)+1e-6)/(max(M.mat)-min(M.mat))
# M.scale[lower.tri(M.scale,diag=TRUE)] <- 0
# M.dfup <- as.data.frame(M.scale)
# colnames(M.dfup) <- c(1:N)
# med_val <- median(M.scale[upper.tri(M.scale,diag=FALSE)])
# M.df <- M.dfup %>% mutate(row_ID=c(1:N)) %>% gather(key="col_ID",value="value", -row_ID) %>% filter(value>med_val)
# M.df$col_ID <- as.numeric(M.df$col_ID)

trade.df <- trade_map %>% select(long,lat) %>% mutate(row_ID=c(1:N),col_ID=c(1:N))

connect_map <- S.df %>% left_join(trade.df,by="row_ID") %>% select(value,long1=long,lat1=lat,col_ID=col_ID.x) %>% 
  left_join(trade.df,by="col_ID") %>% mutate(alpha=value/10) %>% select(alpha,value,long1,lat1,long2=long,lat2=lat)

# c.exchange <- function(x){
#   if(x==3)
#   {return(2)}
#   else if(x==3)
#   {return(2)}
#   else
#   {return(x)}
# }
# 
# trade_map <- trade_map %>% mutate(C1=map_int(Community,function(x){as.numeric(x)})) %>%
#   mutate(Community=as.character(map_int(C1,c.exchange))) %>% select(-C1)
# core.map <- trade_map

png("result/real_core.png", height=1000,width=2000, res=168)

ggplot() + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),col="#f2f2f2",bg="white",lwd=0.05) +
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())+
  geom_curve(data=connect_map,aes(x=long1,y=lat1,xend=long2,yend=lat2,alpha=alpha),color="lightblue",linewidth=0.05,show.legend = FALSE)+#curvature=0
  geom_star(data=trade_map,aes(x = long, y = lat, group=Community, fill= Community, color= Community, starshape = Community),size=2)+
  scale_starshape_manual(values = c(1,11,13,15,28))+
  scale_color_brewer(palette = 'Set1')+
  scale_fill_brewer(palette = 'Set1')+
  geom_text(data=trade_map,aes(x = long, y = lat, label=region),color="black",vjust=1,size=3)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme(legend.position = c(0.1,0.3))

dev.off()

## Heatmap for core food 30
trade_map <- core.map
core.S3 <- modeSum(core.arrT,3,drop=TRUE)@data
S.scale <- core.S3/(core.arrT@modes[3])
#filestr <- "processed_layers/MeanAdj"
#filename <- paste("result/",filestr, ".png", sep = "", collapse = NULL)

# draw
p <- draw_trade_heat(trade_map,S.scale,TRUE)
p.mean <- p+
  ggtitle("MeanAdj")+
  theme(plot.title = element_text(hjust = 0.5))
# print(p.mean)
# png(filename, height=600,width=700, res=168)
# pt <- p+
#   ggtitle("Mean of layer-wise adjacency matrices")+
#   theme(plot.title = element_text(hjust = 0.5))
# print(pt)
# dev.off()
plist <- list()
# S.scale = l-th layer of Adj
for (l in (1:core.arrT@modes[3])){#
  
  S.scale <- core.arrT@data[,,l]
  # draw
  p <- draw_trade_heat(trade_map,S.scale)
  pt <- p+
    ggtitle(as.character(core_food$layerLabel[l]))+
    theme(plot.title = element_text(hjust = 0.5)) 
  plist[[l]] <-pt
  # png(filename, height=650,width=600, res=168)
  # print(pt)
  # dev.off()
}

parr <- ggarrange(plist[[1]],plist[[2]],plist[[3]],
                  plist[[4]],plist[[5]],plist[[6]],plist[[7]],
                  plist[[8]],plist[[9]],plist[[10]],plist[[11]],
                  plist[[12]],plist[[13]],plist[[14]],plist[[15]],
                  plist[[16]],plist[[17]],plist[[18]],plist[[19]],
                  plist[[20]],plist[[21]],plist[[22]],plist[[23]],
                  plist[[24]],plist[[25]],plist[[26]],plist[[27]],
                  plist[[28]],plist[[29]],plist[[30]],
                  ncol = 6, nrow = 5,common.legend = T)

png("result/layers_heatmap_core.png", height=1200,width=1400, res=100)
print(parr)
dev.off()

print((core.map%>% filter(Community=="1") %>% select(region))$region)
region.degree <- modeSum(modeSum(core.arrT,m=2,drop=TRUE),m=2,drop=TRUE)@data
C.prep <- core.map %>% arrange(nodeID) %>% mutate(degree=region.degree) %>% arrange(desc(degree)) %>% select(region,Community) 
C1 <- C.prep %>% filter(Community=="1") %>% select(region)
C2 <- C.prep %>% filter(Community=="2") %>% select(region)
C3 <- C.prep %>% filter(Community=="3") %>% select(region)
C4 <- C.prep %>% filter(Community=="4") %>% select(region)
### Data preprocessing ###
# colnames(country_care) <- c("nodeID","region")
# country_care=country_map[,1:2]
# country_map <- country_care %>% left_join(world_map, by = "region") %>% group_by(nodeID,region) %>% 
#   summarise(long=mean(long),lat=mean(lat))
# write.csv(country_map, file="data/country_map.csv ", row.names= FALSE )
# 
# SELECTREGION1 <- world_map %>% filter(region %in% c("South Korea","Russia","Saudi Arabia","UK","USA","United Arab Emirates",
#                                                     "Bosnia and Herzegovina","Cape Verde","Republic of Congo","Costa Rica",
#                                                     "Czech Republic","Moldova","North Macedonia","Venezuela","El Salvador",
#                                                     "New Zealand","South Africa","Sri Lanka","Tanzania","Trinidad","Tobago",
#                                                     "Syria")) %>% group_by(region) %>% summarise(long=mean(long),lat=mean(lat))
# SELECTREGION2 <- world_map %>% filter(region %in% c("Trinidad","Tobago")) %>% summarise(long=mean(long),lat=mean(lat))
# # deleted layer 7,8,9,44,68
# country_map <- read.csv(file="data/country_map.csv",header = TRUE)
# processed_food
# processed.arrT
# save(country_map,processed_food,processed.arrT,file="data/FAO_99_19.Rdata")
# 
# country_order <- country_map %>% mutate(orderID=order(nodeID)) %>% select(nodeID,orderID)
# 
# # unprocessed arrT
# unprocessed_food_order <- unprocessed_food %>%  mutate(orderID=order(layerID)) %>% select(layerID,orderID)
# 
# unprocessed.edges <- unprocessed_food %>% left_join(full_data,by="layerID") %>% filter(from %in% country_map$nodeID) %>% 
#   filter(to %in% country_map$nodeID)  %>% filter(weight>8)%>% filter(from!=to) %>% select(-layerLabel,-weight) %>% 
#   mutate(min=pmap_int(list(from,to),min),max=pmap_int(list(from,to),max)) %>% select(layerID,min,max) %>% 
#   distinct() 
# 
# unprocessed.ind <-unprocessed.edges %>% rename(nodeID=min) %>% left_join(country_order,by="nodeID") %>%
#   rename(x=orderID) %>% select(-nodeID) %>% rename(nodeID=max) %>% left_join(country_order,by="nodeID") %>%
#   rename(y=orderID) %>% select(-nodeID) %>% left_join(unprocessed_food_order,by="layerID") %>% 
#   rename(z=orderID) %>% select(-layerID)
# 
# 
# unprocessed.arrT <- array(0,dim=c(94,94,11))
# unprocessed.arrT[cbind(unprocessed.ind$x,unprocessed.ind$y,unprocessed.ind$z)] <- 1
# unprocessed.arrT[cbind(unprocessed.ind$y,unprocessed.ind$x,unprocessed.ind$z)] <- 1
# unprocessed.arrT <- as.tensor(unprocessed.arrT)
# # processed arrT
# processed_food_order <- processed_food %>%  mutate(orderID=order(layerID)) %>% select(layerID,orderID)
# 
# processed.edges <- processed_food %>% left_join(full_data,by="layerID") %>% filter(from %in% country_map$nodeID) %>% 
#   filter(to %in% country_map$nodeID)  %>% filter(weight>8)%>% filter(from!=to) %>% select(-layerLabel,-weight) %>% 
#   mutate(min=pmap_int(list(from,to),min),max=pmap_int(list(from,to),max)) %>% select(layerID,min,max) %>% 
#   distinct() 
# 
# processed.ind <-processed.edges %>% rename(nodeID=min) %>% left_join(country_order,by="nodeID") %>%
#   rename(x=orderID) %>% select(-nodeID) %>% rename(nodeID=max) %>% left_join(country_order,by="nodeID") %>%
#   rename(y=orderID) %>% select(-nodeID) %>% left_join(processed_food_order,by="layerID") %>% 
#   rename(z=orderID) %>% select(-layerID)
# 
# processed.arrT <- array(0,dim=c(94,94,19))
# processed.arrT[cbind(processed.ind$x,processed.ind$y,processed.ind$z)] <- 1
# processed.arrT[cbind(processed.ind$y,processed.ind$x,processed.ind$z)] <- 1
# processed.arrT <- as.tensor(processed.arrT)
# 
# # core arrT
core_food <- processed_food %>% bind_rows(unprocessed_food) %>% arrange(layerID)
corefood.edges <- processed.edges %>% bind_rows(unprocessed.edges) %>% arrange(layerID)
country_order <- country_map %>% mutate(orderID=order(nodeID)) %>% select(nodeID,orderID)

core_food_order <- core_food %>%  mutate(orderID=order(layerID)) %>% select(layerID,orderID)

core.edges <- core_food %>% left_join(full_data,by="layerID") %>% filter(from %in% country_map$nodeID) %>%
  filter(to %in% country_map$nodeID)  %>% filter(weight>8)%>% filter(from!=to) %>% select(-layerLabel,-weight) %>%
  mutate(min=pmap_int(list(from,to),min),max=pmap_int(list(from,to),max)) %>% select(layerID,min,max) %>%
  distinct()

core.ind <-core.edges %>% rename(nodeID=min) %>% left_join(country_order,by="nodeID") %>%
  rename(x=orderID) %>% select(-nodeID) %>% rename(nodeID=max) %>% left_join(country_order,by="nodeID") %>%
  rename(y=orderID) %>% select(-nodeID) %>% left_join(core_food_order,by="layerID") %>%
  rename(z=orderID) %>% select(-layerID)


core.arrT <- array(0,dim=c(94,94,30))
core.arrT[cbind(core.ind$x,core.ind$y,core.ind$z)] <- 1
core.arrT[cbind(core.ind$y,core.ind$x,core.ind$z)] <- 1
core.arrT <- as.tensor(core.arrT)
save(full_data,country_map,food_type,processed_food,unprocessed_food, core_food,
     unprocessed.edges,processed.edges, corefood.edges, unprocessed.arrT,
     processed.arrT,core.arrT,file="data/FAO_data.Rdata")




