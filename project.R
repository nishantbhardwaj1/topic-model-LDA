#Loading necessary libraries, set working directory, and reading in the data. 
library(bibliometrix)
library(tm)
library(topicmodels)
library(slam)
library(lme4)
getwd()


#reading the files
k <- convert2df(file=c("literature.bib"), dbsource = "isi", format = "bibtex")

#finding the dimension of the matrix
dim(k)
#sort the vector or factor into ascending or descending
sort(names(k))
#removing articles without abstracts, first finding out how many articles are without abstracts
length(which(is.na(k$AB)))
k <- k[-which(is.na(k$AB)),]
files <- k$AB
summary(files)


#creating a corpus
docs <- Corpus(VectorSource(files))
length(docs)
########################### tidying up the corpus
# Tranform to lower case
docs <-tm_map(docs,content_transformer(tolower))
writeLines(as.character(docs[[1]]))

## Remove ELSEVIER LTD. ALL RIGHTS RESERVED
pubtripe <- c("(c)", "elsevier", "ltd", "all rights reserved")
docs <- tm_map(docs, removeWords, pubtripe)

# Remove stop words
sort(stopwords("english"))
stop_add<-c("also", "although","can","will","one","two","use","whether", "however","may","within")
stop_words <- c(stopwords('english'), stop_add)
docs <- tm_map(docs, removeWords, stop_words)


# Change hyphen and slash to space
toSpace <- content_transformer(function(x, pattern) { return (gsub(pattern, " ", x))})
docs <- tm_map(docs, toSpace, "-")
docs <- tm_map(docs, toSpace, "/")

# Remove punctuation
docs <- tm_map(docs, removePunctuation)

# Strip digits
docs <- tm_map(docs, removeNumbers)

# Remove whitespace
docs <- tm_map(docs, stripWhitespace)

writeLines(as.character(docs[[1]]))
######stemming the document
docs <- tm_map(docs,stemDocument) 
writeLines(as.character(docs[[1]]))
######creating a document term matrix
dtm <- DocumentTermMatrix(docs)
dim(dtm)
# Convert rownames to article IDs
rownames(dtm) <- k$UT
##removing words which occur in just N documents
## Create binary matrix to count number of articles that each word appears in
dtmbin <- dtm
dtmbin[which(apply(dtmbin, 2, function(x) x>=1))] <- 1

## Have a look at the proportion of words contained in 5 or fewer articles...
length(which(col_sums(dtmbin)<=1))/ncol(dtmbin)
length(which(col_sums(dtmbin)<=2))/ncol(dtmbin)
length(which(col_sums(dtmbin)<=5))/ncol(dtmbin) 
length(which(col_sums(dtmbin)<=10))/ncol(dtmbin)
## Remove words occuring in only 5 or fewer documents
dtm <- dtm[,-which(col_sums(dtmbin)<=5)]  
dim(dtm)


##################################LDAtuning for finding out suitable number of topics
install.packages("ldatuning")
library(ldatuning)
result <- FindTopicsNumber(
  dtm,
  topics = seq(from = 5, to = 75, by = 5),
  metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
  method = "Gibbs",
  control = list(seed = 77),
  mc.cores = 4L,
  verbose = TRUE
)
FindTopicsNumber_plot(result)
############################################Topic modelling
# fixing number of toipcs and iterations and using gibbs sampling and 2000 iterations for ensuring convergence
ntopics <- 30
tmod <- LDA(dtm, k=ntopics, method="Gibbs", control=list(burnin=1000, thin=100, iter=2000, best=T))
save(tmod, file="tmod.RData")

#inspecting 20 highest weighted words
Tt<-terms(tmod, 20)
###################Wrting a csv file for the terms
write.csv(Tt,'tmod3.csv')
topics(tmod)






##############################################################################################
#visualising the fitted lda model on web based interface
topicmodels_json_ldavis <- function(tmod, docs, dtm){
  # Required packages
  library(topicmodels)
  library(dplyr)
  library(stringi)
  library(tm)
  library(LDAvis)
library(servr)
library(slam)
library(gistr)
  install.packages("devtools")
  library(devtools)
  devtools::install_github("cpsievert/LDAvis")
  

  # Find required quantities
  phi <- posterior(tmod)$terms %>% as.matrix
  theta <- posterior(tmod)$topics %>% as.matrix
  vocab <- colnames(phi)
  term_freq <- slam::col_sums(dtm)
  # Convert to json
  json_lda4 <- LDAvis::createJSON(phi = phi, theta = theta,
                                 vocab = vocab,
                                 doc.length = as.vector(table(dtm$i)),
                                 term.frequency = term_freq, reorder.topics = FALSE)
  
  return(json_lda4)
}

serVis(json_lda4, out.dir = 'vism',as.gist=TRUE,open.browser = TRUE)
  
  
######################################################################################################




### Extract the weight of each word within each topic and seeing how word weight decresd after 20-25 words so 20 was good choice
lda.post <- posterior(tmod)
M1 <- lda.post$terms

# Unnecessarily long code for compiling dataframe
nwo <- seq(1,500,1)
weightsum <- data.frame()
for(t in 1:dim(M1)[1]){ 
  topr <- as.numeric(sort(M1[t,], decreasing=T))
  for(w in 1:length(nwo)){
    out <- data.frame(topic=t, nwords=nwo[w], sumweight=sum(topr[1:nwo[w]]), singlewordweight=topr[nwo[w]])
    weightsum <- rbind(weightsum, out)
  }  }

# Plot word weight against word rank 
plot(singlewordweight ~ nwords, weightsum,col="royalblue",
     xlab="Word rank", ylab="Word weight", 
     cex.lab=1.4, cex=0.5,col.lab="plum4")
#plotting word weight and rank as boxplot and limitng the x axis to 40 words 
library(ggplot2)
ggplot(data=weightsum, aes(group=nwords,x=nwords, y=singlewordweight))+geom_boxplot()+xlim(0,40)+scale_y_continuous(breaks=seq(0,0.30, by =0.05))+xlab("Word rank")+ylab("Word weight")
ggsave("plot1.png", last_plot())

#########################################plotting top 5 words of every topic
library(tidyverse)
library(tidytext)
congress_lda_td <- tidy(tmod)
congress_lda_td
#combining two dataframes so that we can have topic names in congress_lda_td
congress_lda_td<- 
  merge(congress_lda_td, topic.names, by="topic")
top_terms <- congress_lda_td %>%
  group_by(topic) %>%
  top_n(5, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)
top_terms
top_terms %>%
  mutate(topic = factor(topic),
         term = reorder_within(term, beta, topic)) %>%
transform(top_terms,
            label=factor(label,levels=c("Site comparison",
                                        "Seasonality", 
                                        "Mustelid poulation dynamics",
                                        "Reproduction", 
                                        "Habitat selection",
                                        "Wildlife Management","Home range","Distribution","Deer","Sampling","Population dynamics","Population genetics","Movement ecology","Habitat","Dormice nest","Studies interpretation",
                                        "Epidemiology","Hare population trends","Urban ecology","Squirrel interactions","Small mammals","Fox ecology and control",
                                        "Badger","Modelling&data","Temporal pattern","Prey selection","Vole cycling","Foraging & diet","Bat ecology","Rabbit impacts" 
            ))) %>% ggplot(aes(term, beta))+
  geom_bar(alpha = 0.8, stat = "identity", show.legend = FALSE) +
  scale_x_reordered() +
  facet_wrap(~ label, scales = "free", ncol = 4) +
  coord_flip()
ggsave("topic.png", last_plot())


######################################frequency of topics within the corpus
##Each article has a probability of being ‘about’ each topic (for each article, these probabilities sum to one across all topics). Articles can therefore be assigned to a topic by identifying the topic with the highest probability.
# Plot probability of each topic for first article in dataset
plot(lda.post$topics[1,]~seq(1,ntopics),
     xlab="Topic number", ylab="Probability", pch=16,
     main=paste("Article ID =", row.names(lda.post$topics)[1]))
# Check the topic assigned to that article
topics(tmod)[1]

# Extract which topic each article has been assigned to
mtops <- topics(tmod)

# Summarise data for a barplot, how many articles have been assigned to each topic
topicfreq <- tapply(mtops, mtops, length)
topicfreq <- sort(topicfreq, decreasing=T)
barplot(topicfreq,col="steelblue",las=2,cex.names = 0.8,
        xlab="Topic number", ylab="Number of articles")



#######################Topic labelling
# Create topic labels
topic.names <- tibble(
  topic = 1 : 30,
  label = c(
    "Site comparison",
    "Seasonality", 
    "Mustelid poulation dynamics",
    "Reproduction", 
    "Habitat selection",
    "Wildlife Management","Home range","Distribution","Deer","Sampling","Population dynamics","Population genetics","Movement ecology","Habitat","Dormice nest","Studies interpretation",
    "Epidemiology","Hare population trends","Urban ecology","Squirrel interactions","Small mammals","Fox ecology and control",
    "Badger","Modelling&data","Temporal pattern","Prey selection","Vole cycling","Foraging & diet","Bat ecology","Rabbit impacts" 
    
  )
)
#adding theme info to topic.names dataframe
library(tidyverse)
topic.names1 <- topic.names %>%
  mutate(theme = case_when(
    topic %in% Generic ~ "Generic",
    topic %in% Mammals ~ "Mammals",
    topic %in% Population ~ "Population",
    topic %in% Basicecology ~ "Basicecology",
    topic %in% Management ~ "Management",
    topic %in% Spatialpattern ~ "Spatial"
  ))



#################################################Topic similarity
#installing the library ape
library(ape)
##########M1 is matrix of weight of each word within each topic, I log 10 transformed the matrix M1
MS<-log10(M1)
#calculating euclidean distance
ms.dist1 <- hclust(dist(MS, method = "euclidean")) #D1
ms.dist.p1 <-as.phylo(ms.dist1)
plot(ms.dist.p1)

#tiplabels
ms.dist.p1$tip.label<-as.character(topic.names1[[2]])
my_colours <- c("tomato1", "goldenrod","springgreen3","turquoise3","steelblue2","mediumorchid1")
plot(ms.dist.p1, use.edge.length=F, label.offset=4,tip.color =my_colours[as.factor(topic.names1$theme)])
tiplabels(topic.names1$high.topic,adj=-.2, frame = "none",cex=1, font=2) #0.0005 0.1
add.scale.bar(cex = 0.7, font = 2, col = "black")
legend("bottomleft",title="Theme",legend = c("Basic ecology","Generic","Mammals","Management","Population","Spatial"),fill = my_colours,xpd = T, cex = 0.9)


##############################################Topic co-occurence/research gaps
library(ape)

D1 <- dist(M1, method = "euclidean")
D2 <- dist(M2, method = "euclidean")
# scale D1 and D2 to values b/w 0 and 1

D1.scale <- (D1-min(D1))/(max(D1)-min(D1))  
D2.scale <- (D2-min(D2))/(max(D2)-min(D2))   

#calc product of D1 and D2
D.p <- D1.scale*D2.scale

#plotting
dst <- data.matrix(D.p)
dim <- ncol(dst)
#heatmap

par(mar=c(9,9,1,1))
image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="", col = rev(heat.colors(10)))
axis(1, 1:dim, topic.names[,2], cex.axis = 0.7, las=3)
axis(2, 1:dim, topic.names[,2], cex.axis = 0.7, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=0.6)















##############################################Topic generality vs specificity

M2 <- lda.post$topics
M2 <- t(M2)
M2bin <- matrix(0, nrow=dim(M2)[1], ncol=dim(M2)[2])
M2bin[which(apply(M2, 2, function(x) x == max(x)))] <- 1 
M2select <- matrix(0, nrow=dim(M2)[1], ncol=dim(M2)[2])
M2select[which(M2bin==1)] <- M2[which(M2bin==1)]
M2unsel <- matrix(0, nrow=dim(M2)[1], ncol=dim(M2)[2])
M2unsel[which(M2bin==0)] <- M2[which(M2bin==0)]
meansel <- apply(M2select, 1, function(x) sum(x)/length(which(x!=0)))
meanunsel <- apply(M2unsel, 1, function(x) sum(x)/length(which(x!=0)))
generality <- data.frame(meansel=meansel, meanunsel=meanunsel, TopicN=seq(1,ntopics))

## Plotting
{plot(meansel ~ meanunsel, generality, 
      ylab="Mean weight (selected articles)", xlab="Mean weight (unselected articles)",
      col="white",col.axis="black", las=2,col.lab="plum4", cex.lab=1.3, cex.axis=1.0,srt=45)
  text(generality$meanunsel, generality$meansel, 
       labels=generality$TopicN, font=2, cex=1.1,col="black")
  title(main="Specific topics", adj=0, font.main=2, cex.main=1,col.main="blue")
  title(sub="General topics", adj=1, font.sub=2, cex.main=1,col.sub="red")}


# putting topics in theme
Generic <- c("2","1", "24", "16", "25","10")

Mammals <- c("9", "15","20", "21", "22", "23", "29","30")

Population <- c("3", "11","12", "18", "27")

Basicecology <- c("4", "26", "28")

Management <- c("6", "17", "19")

Spatialpattern <- c("5", "7", "8", "13", "14")

#adding theme info to generality dataframe to colour the topic numbers according to themes
generality2 <- generality %>%
  mutate(theme = case_when(
    TopicN %in% Generic ~ "Generic",
    TopicN %in% Mammals ~ "Mammals",
    TopicN %in% Population ~ "Population",
    TopicN %in% Basicecology ~ "Basic ecology",
    TopicN %in% Management ~ "Management",
    TopicN %in% Spatialpattern ~ "Spatial"
  ))

#ggggggggplotttttttt
ggplot(generality2, aes(x=meanunsel, y=meansel, colour=theme, label=generality2$TopicN), alpha=0.8)+geom_point(alpha=0.1)+
  geom_text(size=6)+labs(x="Mean weight (unselected articles)",y= "Mean weight (selected articles)")







##################################################topic popularity
topartyear <- data.frame(TopicN=topics(tmod), UT=k$UT, PY=k$PY)
nartyr <- aggregate(k$PY, by=list(k$PY), length)
names(nartyr) <- c("PY","nart")
topartyear$decade<-floor(topartyear$PY/10)*10
nartyr$decade<-floor(nartyr$PY/10)*10
nartdec<-aggregate(nartyr$nart,by=list(decade=nartyr$decade),sum)
names(nartdec) <- c("decade","nart")

topicyr <- aggregate(topartyear$TopicN, by=list(topartyear$TopicN, topartyear$decade), length)
names(topicyr) <- c("TopicN","decade","narticles")
topicyr <- merge(topicyr, nartyr, by="decade")
topicyr$propart <- topicyr$narticles / topicyr$nart
topicyr$topicF <- as.factor(topicyr$TopicN)
topicyr$PYscaled <- scale(topicyr$decade)

## Model
library(lme4)
popmod <- glmer(narticles ~ PYscaled + (1 + PYscaled|topicF), data=topicyr, 
                family=poisson(link = "log")) 
summary(popmod)
table(topicyr$topicF,topicyr$PY)

random.df <- ranef(popmod, condVar=TRUE)$topicF
names(random.df) <- c("Intercept","Slope")
random.df$TopicN <- seq(1,ntopics)

#adding the label of topic to random.df dataframe to allow the labeling of topic numbers in the plot
library(tidyverse)
random.df3 <- random.df %>%
  mutate(Topic = case_when(
    TopicN %in%Site_comparison  ~ "Site comparison",
    TopicN %in% Seasonality ~ "Seasonality",
    TopicN %in% Mustelid_population_dynamics ~ "Mustelid population dynamics",
    TopicN %in% Reproduction ~ "Reproduction",
    TopicN %in% Habitat_selection ~ "Habitat selection",
    TopicN %in% Wildife_management ~ "Wildlife management",
    TopicN %in% Home_range~"Home range",
    TopicN %in% Distribution~"Distribution",
    TopicN %in% Deer~"Deer",
    TopicN %in% sampling~"Sampling",
    TopicN %in% Population_dynamics~"Population dynamics",
    TopicN %in% Population_genetics~"Population genetics",
    TopicN %in% Movement_ecology~"Movement ecology",
    TopicN %in% Habitat~"Habitat",
    TopicN %in% Dormice_nest~"Dormice nest",
    TopicN %in% Studies_interpretation~"Studies interpretation",
    TopicN %in% Epidemiology~"Epidemiology",
    TopicN %in% Hare_population_trends~"Hare population trends",
    TopicN %in% Urban_ecology~"Urban ecology",
    TopicN %in% Squirrel_interactions~"Squirrel interactions",
    TopicN %in% Small_mammals~"Small mammals",
    TopicN %in% Fox_ecology_and_control~"Fox ecology and control",
    TopicN %in% Badger~"Badger",
    TopicN %in% Modelling_data~"Modelling & data",
    TopicN %in% Temporal_pattern~"Temporal pattern",
    TopicN %in% Prey_selection~"Prey selection",
    TopicN %in% Vole_cycling~ "Vole cycling",
    TopicN %in% Foraging_diet~"Foraging & diet",
    TopicN %in% Bat_ecology~"Bat ecology",
    TopicN %in% Rabbit_impacts~"Rabbit impacts"
    
  ))



## Plot
{plot(Slope ~ Intercept, random.df, 
      ylab="Slope(growth rate)",xlab="Intercept(number of publications)", 
      col="black", pch=20,cex.lab=1.0, cex.axis=1.2
     # xlim=c(-0.16,0.16), ylim=c(-0.05,0.05)
      )
  
  abline(h=0,v=0,lty=3)
  text(random.df$Intercept, random.df$Slope, labels=random.df3$Topic,pos=1,
       font=2, cex=0.7)
  text(0.4,-0.5, "Large number but declining", col="cadetblue",font=3, cex =1.1)
  text(0.4,0.5, "Large number and increasing",col="firebrick3", font=3,cex =1.1)
  text(-0.3,-0.5, "Small number and declining",col="lightslateblue", font=3,cex =1.1)
  text(-0.4,0.5, "Small number but increasing",col="peachpuff4", font=3,cex =1.1)}


#plotting topic populrarity over the years, no of articles over different years
g2<- ggplot(topicyr,aes(x=PY,y=narticles,colour=topicF))+geom_line(size=0.8)+facet_wrap(topicF~.)+
  
  labs(title="topic popularity", x="publication year", y="number of articles")

g2


####adding a row of topic names to topicyr, all the fuss to just add topic names to topicyr dataframe
Site_comparison <- ("1")
Seasonality <- ("2")
Mustelid_population_dynamics <- ("3")
Reproduction <-("4")
Habitat_selection <- ("5")
Wildife_management<- ("6")
Home_range<- ("7")
Distribution<- ("8")
Deer<- ("9")
sampling<- ("10")
Population_dynamics<- ("11")
Population_genetics<- ("12")
Movement_ecology<- ("13")
Habitat<- ("14")
Dormice_nest<- ("15")
Studies_interpretation<-("16")
Epidemiology<- ("17")
Hare_population_trends<-("18")
Urban_ecology<-("19")
Squirrel_interactions<-("20")
Small_mammals<-("21")
Fox_ecology_and_control<-("22")
Badger<-("23")
Modelling_data<-("24")
Temporal_pattern<-("25")
Prey_selection<-("26")
Vole_cycling<-("27")
Foraging_diet<-("28")
Bat_ecology<-("29")
Rabbit_impacts<-("30")

topart5 <- topicyr %>%
  mutate(Topic = case_when(
    topicF %in%Site_comparison  ~ "Site comparison",
    topicF %in% Seasonality ~ "Seasonality",
    topicF %in% Mustelid_population_dynamics ~ "Mustelid population dynamics",
    topicF %in% Reproduction ~ "Reproduction",
    topicF %in% Habitat_selection ~ "Habitat selection",
    topicF %in% Wildife_management ~ "Wildlife management",
    topicF %in% Home_range~"Home range",
    topicF %in% Distribution~"Distribution",
    topicF %in% Deer~"Deer",
    topicF %in% sampling~"Sampling",
    topicF %in% Population_dynamics~"Population dynamics",
    topicF %in% Population_genetics~"Population genetics",
    topicF %in% Movement_ecology~"Movement ecology",
    topicF %in% Habitat~"Habitat",
    topicF %in% Dormice_nest~"Dormice nest",
    topicF %in% Studies_interpretation~"Studies interpretation",
    topicF %in% Epidemiology~"Epidemiology",
    topicF %in% Hare_population_trends~"Hare population trends",
    topicF %in% Urban_ecology~"Urban ecology",
    topicF %in% Squirrel_interactions~"Squirrel interactions",
    topicF %in% Small_mammals~"Small mammals",
    topicF %in% Fox_ecology_and_control~"Fox ecology and control",
    topicF %in% Badger~"Badger",
    topicF %in% Modelling_data~"Modelling & data",
    topicF %in% Temporal_pattern~"Temporal pattern",
    topicF %in% Prey_selection~"Prey selection",
    topicF %in% Vole_cycling~ "Vole cycling",
    topicF %in% Foraging_diet~"Foraging & diet",
    topicF %in% Bat_ecology~"Bat ecology",
    topicF %in% Rabbit_impacts~"Rabbit impacts"
    
  ))
g9<- ggplot(topart5,aes(x=PY,y=narticles))+geom_line(size=0.6, colour="darkslategray4")+facet_wrap(Topic~.)+
  
  labs(title="topic popularity", x="publication year", y="number of articles")

g9


#plotting few selected topics from the list using subset

gg<- ggplot(subset(topart5, topicF %in% c(17, 19, 20, 21,23,29)), aes( x=PY,y=narticles))+geom_line(size=0.6, colour="steelblue")+facet_wrap(. ~ Topic) +
  labs(title="topic popularity", x="publication year", y="number of articles")# change topic number to plot other topics. 
gg






#####################plotting themes over time/years
topart2 <- topicyr %>%
  mutate(theme = case_when(
    topicF %in% Generic ~ "Generic",
    topicF %in% Mammals ~ "Mammals",
    topicF %in% Population ~ "Population",
    topicF %in% Basicecology ~ "Basic ecology",
    topicF %in% Management ~ "Management",
    topicF %in% Spatialpattern ~ "Spatial"
  ))
g7<- ggplot(topart2,aes(x=PY,y=narticles))+geom_line(colour='steelblue')+facet_wrap(theme~.)+
  
  labs( x="publication year", y="number of articles")

g7

#####putting years into decades and then plotting  themes over years because it was quite messy and confusing earlier in g7
# topics by decades
dec1970 <- c("1971","1972", "1973", "1974", "1975", "1976","1977", "1978", "1979", "1980")

dec1980 <- c("1981","1982", "1983", "1984", "1985", "1986","1987", "1988", "1989", "1990")

dec1990 <- c("1991","1992", "1993", "1994", "1995", "1996","1997", "1998", "1999", "2000")

dec2000 <- c("2001", "2002", "2003", "2004","2005", "2006", "2007","2007", "2008", "2009","2010")

dec2010 <- c("2011", "2012", "2013", "2014","2015", "2016", "2017","2017", "2018", "2019","2020")
topart3 <- topart2 %>%
  mutate(decade = case_when(
    PY %in% dec1970 ~ "1971-80",
    PY %in% dec1980 ~ "1981-90",
    PY %in% dec1990 ~ "1991-00",
    PY %in% dec2000 ~ "2001-10",
    PY %in% dec2010 ~ "2011-20"
    
  ))
g6<- ggplot(topart3,aes(x=decade, y=narticles,fill=decade))+geom_boxplot(size=0.7, alpha=0.5)+facet_wrap(theme~.)+
  
  labs( x="publication decade", y="number of articles")

g6
library(ggplot2)

topart4 <- topart2 %>%
  mutate(decade = case_when(
    PY %in% dec1970 ~ "1971-80",
    PY %in% dec1980 ~ "1981-90",
    PY %in% dec1990 ~ "1991-00",
    PY %in% dec2000 ~ "2001-10",
    PY %in% dec2010 ~ "2011-20"
    
  ))
g3<- ggplot(topart4,aes(x=decade, y=narticles))+geom_col(size=0.5)+scale_y_continuous(breaks=seq(0,150, by =50))+facet_wrap(theme~.)+
  
  labs( x="publication decade", y="number of articles")

g3

ggsave("theme.png", last_plot())
ggsave('theme2.png',last_plot())
topics(tmod)[51]

library (tidyverse)

