                                #
  #Dissimilarity Analysis       #
      #by Korryn Bodner         #
#################################

library(ggplot2)
library(ggrepel)


df.ind_fish<-read.csv("NEON_fish_inds.csv", header=T)

#read in stage-structure network's edgelist
#duplicate 1st and 2nd column and switch col #s to allow for easy counting of all undirected interactions

stages<-read.table("stage.edgelist.txt")
stages.copy<-stages[c(2,1)]
colnames(stages.copy)<-colnames(stages)

#FUNCTION avg.size.metrics
#calculates avg adult-juv size diff. and avg. adult length per species
avg.size.metrics<-function(dataset){
  
  fish.names<-levels(as.factor(dataset$scientificName))
  final.df<-data.frame(character(length(fish.names)),numeric(length(fish.names)),numeric(length(fish.names)),stringsAsFactors=FALSE)
  colnames(final.df)<-c("scientificName","size.diff","avg.adult.length")
  for (i in 1:length(fish.names)){
    final.df$scientificName[i]<-as.character(fish.names[i])
    final.df$size.diff[i]<-NA
    temp<-subset(dataset,scientificName==fish.names[i])
    adult.length<-round(mean(temp$fishTotalLength[temp$fishLifeStage=="adult"],na.rm=T))
    final.df$avg.adult.length[i]<-adult.length
    if(length(unique(temp$fishLifeStage))>1){
      #print(summary(aov(fishTotalLength~fishLifeStage,data=temp)))
      juv.length<-mean(temp$fishTotalLength[temp$fishLifeStage=="juvenile"],na.rm=T)
      diff<-adult.length - juv.length
      final.df$size.diff[i]<-diff
    }
    i=i+1
  }
  return(final.df)
}

#create size and size differences dataset
fish.char.dataset<-avg.size.metrics(df.ind_fish)

#need to append all connections
stages.full<-rbind(stages,stages.copy)
colnames(stages.full)<-c("Species", "Connections")

#copy down species name, life stage
stages.full$species.name<-unlist(lapply(strsplit(stages.full$Species, split="_"), '[', 1)) #takes first element
stages.full$life.stage<-unlist(lapply(strsplit(stages.full$Species, split="_"), '[', 2)) #takes second element

#the stage-structured fish (i.e. the largest 50)
largestFish.50 = c("A.melas","A.natalis","H.etowanum","I.gagei","L.aepyptera",
                   "L.chrysocephalus","L.cyanellus","L.macrochirus","L.megalotis",
                   "M.carinatum","N.funebris","N.leptacanthus","N.leptocephalus",
                   "S.atromaculatus","S.fontinalis","S.trutta")
#piscivore fish
piscFish = c("A.natalis","C.carolinae","L.cyanellus","L.megalotis",
             "S.trutta", "S.fontinalis","S.atromaculatus","M.dolomieu")

#indicate if in pisc
stages.full$pisc<-stages.full$species.name %in% piscFish

#indicate if in 50% largest
stages.full$fifty.largest<-stages.full$species.name %in% largestFish.50


#######Jaccard's Index Prep###
#subset only 50 largest fish
df.largest.50<-subset(stages.full,!is.na(stages.full$life.stage))
species.names.50<-unique(df.largest.50$species.name)

#create empty dataframe
similarity.df<-data.frame(species.names.50)
similarity.df$jaccard<-NA
similarity.df$pisc<-NA

#calculate jaccard score for each larger species (and record if piscivore)
for(i in 1:length(species.names.50)){
  df.temp<-subset(df.largest.50,df.largest.50$species.name == species.names.50[i])
  adult<-subset(df.temp,life.stage=="adult")
  juv<-subset(df.temp,life.stage=="juvenile")
  
  #prep for similarity
  a=sum(adult$Connections %in% juv$Connections) #in both
  b=nrow(adult) - a #in adult but not juv
  c=nrow(juv) - a #in juv but not adult
  
  jaccard.index= 1-(a/(a+b+c))
  
  #add values to dataset
  similarity.df$jaccard[i]<-jaccard.index
  similarity.df$pisc[i]<-df.temp$pisc[1]
}

#connect average length/weight and difference (between juv and adult) with similiarity metrics?
new.data<-merge(similarity.df,fish.char.dataset, by.x = "species.names.50", 
      by.y = "scientificName",all.x=T,all.y=F)
new.data.noNA<-subset(new.data,!is.na(size.diff))

#linear regressions
#Jaccard index ~ adult-juv length difference (mm)
size.diff<-lm(jaccard~size.diff,data=new.data.noNA)
summary(size.diff)

#regression with L.megalotis removed
data.subset<-subset(new.data.noNA,species.names.50!="L.megalotis")
summary(lm(jaccard~size.diff,data=data.subset))


#plot Jaccard dissimilarity against adult-juv size-difference
jaccard.plot<-ggplot(new.data.noNA, aes(x=size.diff, y=jaccard),alpha=0.5) + 
  geom_smooth(method=lm, se=F,col="grey69") + theme_classic() +
  geom_point(aes(size=avg.adult.length,col=pisc)) +
  scale_color_manual(values = c("cyan4", "darkred"),labels=c("Non-piscivore","Piscivore")) +
  scale_size(breaks=c(100,125,150,175)) + labs(colour="Feeding behaviour", size="Avg. adult. length (mm)")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2), name="Jaccard dissimilarity index")+
  scale_x_continuous(limits=c(10,80),name="Avg. length difference between adult and juvenile (mm)") +
  geom_text_repel(aes(label = species.names.50, color =pisc), 
                  size = 3,segment.alpha = 0,show.legend = F) 


jaccard.plot
