#####################################
                                    #
      #Network Comparison           #
      #MDS and GCM-11               #
      #by Korryn Bodner             #
                                    #
#####################################

#############################################
#MDS & 3D PLOTS
#############################################
#packages for mds and 3dplots
library(car)
library(rgl)

#read in gcd11 file created using code from Yaveroglu et al. 2014 
#instructions and code for GCD-11 found here: http://www0.cs.ucl.ac.uk/staff/natasa/GCD/
gcd.11<-read.csv("gcd11.csv",row.names=1,header=T)

#MDS
d <- dist(gcd.11) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim

#create dataframe of MDS in 3 dimensions
mds.3d<-data.frame(fit$points[,1:3])
mds.3d$names<-row.names(mds.3d)

#name variable so all networks within the same group are categorized the same
mds.3d$final.names<-NA
mds.3d$final.names[which(grepl("rewire_fifty*", mds.3d$names))]<-"Rewired stage"
mds.3d$final.names[which(grepl("rewire_adult*", mds.3d$names))]<-"Rewired adult"
mds.3d$final.names[which(grepl("adult_network", mds.3d$names))]<-"Adult"
mds.3d$final.names[which(grepl("fifty_network", mds.3d$names))]<-"Stage"

mds.3d$final.names<-as.factor(mds.3d$final.names)

#plot interactive 3d graph
scatter3d(x=mds.3d$X1,y=mds.3d$X2,z=mds.3d$X3, sphere.size = 2.5, 
          groups=mds.3d$final.names,bty=f,
          surface.col=c("darkgoldenrod4", "gold1","pink2", "magenta4"),surface=F,data=mds.3d,
          xlab="",ylab="",zlab="",axis.col = c("black","black","black"),axis.scales = F)
legend3d("topright",pch=19, legend = c("Stage",
                                       "Adult", "Rewired stage",
                                       "Rewired adult"),
         col = c("magenta4","darkgoldenrod4","pink3","gold1"), cex=1.3, inset=c(0.17))

################################################
#GCM-11s & GCD-11 Comparisons
################################################
##Packages

library(gdata) #GCD-11 comparisons
library(corrplot) #GCM-11 matrices
#source("http://www.sthda.com/upload/rquery_cormat.r")

####Calculate GCD-11s for each network group

#gcd_11 with only rewired adult networks (100x100)
random_adult_index<-grepl("rewire_adult",colnames(gcd.11))
random_adult<-subset(gcd.11,subset=random_adult_index,select=which(random_adult_index))

#gcd_11 with only rewired stage-structured networks (100x100)
random_stages_index<-grepl("rewire_fifty",colnames(gcd.11))
random_stages<-subset(gcd.11,subset=random_stages_index,select=which(random_stages_index))

#gcd_11 with all rewired networks, rewired adult x rewired stage (100x100)
all_randoms<-subset(gcd.11,subset=random_adult_index,select=which(random_stages_index))

#adult network with all networks (202x1)
adult_index<-grepl("adult_network", colnames(gcd.11))
adult<-data.frame(subset(gcd.11,select=which(adult_index)))
  
#stage-structured network with all networks (202x1)
stages_index<-grepl("fifty_network", colnames(gcd.11)) 
stages<-data.frame(subset(gcd.11,select=which(stages_index)))

#upper matrix for random networks
adult.avg.up.triangle<-upperTriangle(random_adult)
stage.avg.up.triangle<-upperTriangle(random_stages)

#Avg. within-group GCD for random networks
random.adult<-sum(adult.avg.up.triangle)/length(adult.avg.up.triangle)
random.stage<-sum(stage.avg.up.triangle)/length(adult.avg.up.triangle)

#Avg between-group GCD for random networks
random.between.group<-sum(all_randoms)/(100*100) #matrix is 100x100

#GCD scores for adults
adult.with.stage<-as.numeric(subset(adult,row.names(adult)=="fifty_network"))
adult_rewireadult_subset<-subset(adult,grepl("rewire_adult",row.names(adult)))  
adult.with.avg.adult<-sum(adult_rewireadult_subset$adult_network)/length(adult_rewireadult_subset$adult_network)

#GCD scores for stage 
stage_rewirestage_subset<-subset(stages,grepl("rewire_fifty",row.names(stages)))
stage.with.avg.stage<-sum(stage_rewirestage_subset$fifty_network)/length(stage_rewirestage_subset$fifty_network)


#create table
gcd_scores<-c(adult.with.stage,stage.with.avg.stage,
  adult.with.avg.adult,random.between.group,
  random.stage,random.adult)

pairwise_comparisons<-c("Stage - Adult", "Stage - Avg. rewired staged",
        "Adult - Avg. rewired adult", "Avg. rewired adult - Avg. rewired stage",
        "Random rewired stage", "Radom rewired adult")
gcd_table<-as.table(cbind(pairwise_comparisons,round(gcd_scores,2)))
colnames(gcd_table)<-c("Network pairwise comparisons", "GCD-11 scores")

#view table
gcd_table

####Plot GCM-11s

#read in adult and stage-structured .ndump2 files 
# .ndump2 files created using GCD-11 code by Yaveroglu et al. 2014

filelist.adult = list.files(pattern = "rewire_adult.*.ndump2")
files.adult<-lapply(filelist.adult, FUN=read.table, header=F) #need to set it up to read files in git
filelist.stage = list.files(pattern = "rewire_fifty.*.ndump2")
files.stage<-lapply(filelist.stage, FUN=read.table, header=F) #need to set it up to read files in git

#read in stage-structured and adult inferred networks
adult.orbits<-read.csv("adult.gcd.csv",header=T)
stage.orbits<-read.csv("stages.gcd.csv",header=T)

#Step 1: Create average node positions for rewired networks

#function for node averages
average_node_pos<-function(list_of_files){

  n=length(list_of_files) #length of files
  m=nrow(list_of_files[[1]]) #number of nodes - 31 for adult, for stage
 
  #create table for heatmap
  table.gcd<-matrix(0,m,10)
  
  for(i in 1:n){
    #subset from V1-V15 (inclusive)
    table<-list_of_files[[i]][,1:15]
    
    #get rid of V4 (orbit 3) and > V12
    temp.table.gcd<-table[c(1:2,4:12)]
    
    table.gcd=table.gcd+temp.table.gcd
  }
  #average summed total of each node
  average.table.gcd<-table.gcd/100 #100 rewired networks
  return(average.table.gcd)
}

avg.adult.orbits<-average_node_pos(files.adult)
colnames(avg.adult.orbits)<-colnames(adult.orbits)
avg.stage.orbits<-average_node_pos(files.stage)
colnames(avg.stage.orbits)<-colnames(stage.orbits)

#Step 2: Reorder orbits to optomize patterns
#reorder matrix
heatmap.order<-function(dataset){
  temp.df<-dataset
  temp.order<-data.frame()
  temp.order<-cbind(temp.df$orbit_2, temp.df$orbit_11,temp.df$orbit_0,
                    temp.df$orbit_5, temp.df$orbit_7, temp.df$orbit_4,
                    temp.df$orbit_8,temp.df$orbit_10,temp.df$orbit_1,
                    temp.df$orbit_9, temp.df$orbit_6)
  colnames(temp.order)<-c("orbit 2","orbit 11", "orbit 0",
                          "orbit 5", "orbit 7", "orbit 4",
                          "orbit 8", "orbit 10", "orbit 1", 
                          "orbit 9","orbit 6")
  return(temp.order)
}

stage.orbits.ordered<-heatmap.order(stage.orbits)
adult.orbits.ordered<-heatmap.order(adult.orbits)
avg.adult.orbits.ordered<-heatmap.order(avg.adult.orbits)
avg.stage.orbits.ordered<-heatmap.order(avg.stage.orbits)


#Step 2: Modify corrplot function for spearman correlation

#code provided by (http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...,method="spearman")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#Step 3: Plot GCM-11s

#color palette
col=c("#000080", "#0000EE", "#0080FF", "#87CEFF", "#E0F3F8", "#FFEC8B", "#FFA500",
       "#FF4500",
       "#D73027", "#AB0000")

#GCM-11 for stage-structured network
p.stages <- cor.mtest(stage.orbits.ordered) #warnings are just for for spearman correlation
stage.only <-cor(stage.orbits.ordered,method="spearman")
corrplot(stage.only, type="full", order="original",
         tl.col="black",
         col=col,sig.level=0.05,
         p.mat=p.stages)

#GCM-11 for adult network
p.adult <- cor.mtest(adult.orbits.ordered)
adult.only <-cor(adult.orbits.ordered,method="spearman")
corrplot(adult.only, type="full", order="original",
         tl.col="black", 
         col=col,sig.level=0.05,
         p.mat=p.adult)

#GCM-11 for average rewired stage-structured network
p.stages.rewired <- cor.mtest(avg.stage.orbits.ordered)
stage.rewired.only <-cor(avg.stage.orbits.ordered,method="spearman")
corrplot(stage.rewired.only, type="full", order="original",
         tl.col="black", 
         col=col,sig.level=0.05,
         p.mat=p.stages.rewired)

#GCM-11 for average adult network
p.adult.rewired <- cor.mtest(avg.adult.orbits.ordered)
adult.rewired.only <-cor(avg.adult.orbits.ordered,method="spearman")
corrplot(adult.rewired.only, type="full", order="original",
         tl.col="black", 
         col=col,sig.level=0.05,
         p.mat=p.adult.rewired)


