# Predicting the birth of a word
# Brandon C. Roy, Michael C. Frank, Philip DeCamp, Matthew Miller, Deb Roy
#
# Supplementary materials
# 
# Supporting functions for constructing the word panels in Figure 3
#

library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)

#### Load all the data required for this figure ----
md.temp <- local({
  d.temp <- read.csv(file="data/HSP_temporal_dists.tsv",sep="\t",header=TRUE,stringsAsFactors=FALSE)
  names(d.temp) <- gsub(names(d.temp),pattern="^X",replacement="")
  md.temp <- melt(d.temp,id.vars=c("word","total"),variable.name="hour",value.name="count")
  md.temp$hour <- as.numeric(levels(md.temp$hour))[md.temp$hour]
  return(md.temp)
})

topic.topwords <- read.table(file="data/HSP_topic_topwords.tsv",stringsAsFactors=FALSE,header=TRUE)

m.topic.dists <- local({
  topic.dists <- read.table(file="data/HSP_topic_dists.tsv",sep="\t",header=TRUE,stringsAsFactors=FALSE)
  m.topic.dists <- melt(topic.dists,id.vars=c("word","count.ep","count.tok"))
  return(m.topic.dists)
})

# This is important - so that all topic plots are ordered by the background in descending background probability
m.topic.dists$variable <- factor(m.topic.dists$variable,
                                 levels=levels(with(subset(m.topic.dists,word=="BACKGROUND"),
                                                    reorder(variable,-value,ordered=TRUE))))


#### Prepare the plotting and utility functions ----

# The main function that actually draws the whole panel for the target word
draw_panel <- function(wrd) {
  stopifnot(require(grid))
  
  # Prepare the data we'll need
  img <- load_spatial_image(wrd)
  
  #
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 4, widths=c(.09,.7,.8,1.1))))
  
  # Write the word text
  pushViewport(viewport(layout.pos.col = 1))
  word.title <- wrd   
  grid.text(word.title,x=.6,y=.5,rot=90,gp=gpar(fontsize=9,lineheight=.8))
  popViewport()
  
  # draw the spatial distribution image
  pushViewport(viewport(layout.pos.col = 2))
  grid.raster(img[55:805,100:1300,])
  popViewport()
  
  # Plot the temporal distribution
  pushViewport(viewport(layout.pos.col = 3))
  plt <- build_temporalplot(wrd)
  print(plt,newpage=FALSE)
  popViewport()
  
  # Draw the topic distribution
  pushViewport(viewport(layout.pos.col = 4))
  topic_fig <- build_topicfig(wrd)
  grid.draw(topic_fig)
  popViewport()
  
  popViewport()  
}

# Utility function that computes the KL-divergence components
kl_components <- function(p,q) {
  z <- p*(log(p) - log(q))
  return(z)
}

# Function that loads the spatial image
load_spatial_image <- function(wrd) {
  stopifnot(require(png))
  return(readPNG(paste("data/spatial_images/sd_projected_",wrd,".png",sep="")))
}

# Builds the temporal distinctiveness figure element
build_temporalplot <- function(wrd,color.fill="#FA6A03") {
  plt <- ggplot(data=md.temp, aes(x=hour,y=count/total,group=word)) + 
    geom_area(data=subset(md.temp,word=="BACKGROUND"),stat="identity",fill="gray80") +
    geom_bar(data=subset(md.temp,word==wrd), stat="identity",fill=color.fill, alpha=.9, aes(width=.6)) +
    scale_x_continuous(limits=c(6,22), breaks=seq(8,20,4)) + 
    scale_y_continuous(breaks=c(0,.15,.3),limits=c(0,.33)) +
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=8),
          axis.text.y = element_blank(), axis.ticks.y=element_blank(), # WITHOUT ticks
          axis.text.x = element_text(hjust=-.2,vjust=0),
          plot.margin = unit(c(0,.15,-.2,-.25),"cm"), # top, right, bottom, left
          axis.ticks.margin= unit(-.1,"cm"))
  return(plt)
}

# Gets the top words for the top topics for the given word
get_topwords <- function(wrd,mt.dists,topic.topwords) {
  stopifnot(require(reshape2))
  stopifnot(all(c("variable","value") %in% colnames(mt.dists)))
  
  k = 3
  pq <- as.data.frame(t(data.frame(row.names="variable",
                                   dcast(m.topic.dists, variable ~ word, 
                                         subset=.(word %in% c(wrd,"BACKGROUND"))))))
  z <- kl_components(pq[wrd,],pq["BACKGROUND",])
  foo <- order(z,decreasing=TRUE)
  topic.names <- colnames(z)[foo[1:k]]
  
  return(topic.topwords[1:5,tolower(topic.names)])
}

# Builds the linguistic topic distinctiveness figure element
build_topicfig <- function(wrd) {
  t.words <- get_topwords(wrd,m.topic.dists,topic.topwords)
  bar.3 <- build_topicbars(wrd)
  
  # Do NOT create a new page since this is to be embedded
  pushViewport(viewport(layout=grid.layout(2,1,heights=c(.75,1),widths=c(1,1))))
  
  # Draw the graph
  pushViewport(viewport(layout.pos.row = 2))
  print(bar.3,newpage=FALSE)
  popViewport()
  
  # Draw the top words
  pushViewport(viewport(layout.pos.row = 1))
  grid.draw(build_topicwords(wrd))
  popViewport()
  
  # Pop the outermost grid viewport
  popViewport()
  return (grid.grab(wrap=TRUE))
}

# Builds the ling. topic bar graph
build_topicbars <- function(wrd,color.fill="#D8489D") {
  # Build the topic bar graph
  
  # This bit about making this thing a data frame is important so that the topic ids are carried through the kl_components function
  pq <- as.data.frame(t(data.frame(row.names="variable",
                                   dcast(m.topic.dists, variable ~ word, subset=.(word %in% c(wrd,"BACKGROUND"))))))
  kl.c <- kl_components(pq[wrd,],pq["BACKGROUND",])
  kl.c <- melt(kl.c,measure.vars=c(1:ncol(kl.c)),variable.name="topic")
  kl.c$rank <- order(order(kl.c$value,decreasing=TRUE))
  
  top.3 <- kl.c$topic[order(kl.c$rank)[1:3]]
  
  bar.3 <- ggplot(data=m.topic.dists,aes(x=variable,y=value)) + 
    geom_bar(subset = .(word=="BACKGROUND"), stat="identity", width=1, fill="grey") +
    geom_bar(subset = .(word==wrd), stat="identity", width=.6, alpha=.7, fill=color.fill) +
    geom_bar(data=subset(m.topic.dists,word==wrd), subset = .(variable %in% top.3),
             stat="identity", width=.6, alpha=1, fill=color.fill) +
    scale_x_discrete(breaks=top.3, labels=sub(pattern="^T",replacement="",top.3)) +
    scale_y_continuous(limits=c(0,.34), breaks=c(0,.15,.3)) +
    theme(panel.background=element_blank(), panel.grid=element_blank(), 
          legend.position="none",
          axis.title=element_blank(), axis.text=element_text(size=8),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_text(hjust=-.15,vjust=0),
          axis.ticks.margin= unit(-.1,"cm"),
          plot.margin = unit(c(0,.1,-.2,-.2),"cm"))
  
  return(bar.3)
}

# Builds the graphical representation of the ling. topic top words
build_topicwords <- function(wrd,font.size=9) {
  t <- gTree(name="Top.topic.words")
  t.words <- get_topwords(wrd,mt.dists=m.topic.dists,topic.topwords=topic.topwords)

  # Build the grobs
  text.x = unit(0,"npc")
  text.x2 <- grobWidth(textGrob("T25: ",hjust=0,gp=gpar(fontface="bold",fontsize=font.size)))
  text.y = unit(1,"npc")
  
  for (i in 1:ncol(t.words)) {
    titleGrob <- textGrob(toupper(colnames(t.words)[i]),
                          x=text.x, y=text.y, just=c("left","top"),
                          gp=gpar(fontface="bold",fontsize=font.size))
    
    wordsGrob <- textGrob(paste(t.words[,i],collapse=", "),
                          x=text.x2, y=text.y, just=c("left","top"),
                          gp=gpar(fontsize=font.size))
    
    text.y <- text.y - unit(.9,"lines")    
    
    t <- addGrob(t,titleGrob)
    t <- addGrob(t,wordsGrob)
  }
  return(t)
}
