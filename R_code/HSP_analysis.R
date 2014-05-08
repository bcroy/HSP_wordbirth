# Predicting the birth of a word
# Brandon C. Roy, Michael C. Frank, Philip DeCamp, Matthew Miller, Deb Roy
#
# Supplementary materials
# Reproduce analyses and generate figures

library(reshape2)
library(ggplot2)
library(plyr)
library(Hmisc)
library(xtable)

rm(list=ls())

#### Load the base data and prepare the derived predictors ----
d <- read.csv("data/HSP_predictors_base.csv",stringsAsFactors=FALSE)

# Prepare the frequency, phoneme and MLU data
d$sln.freq.pre <- with(d, scale(ifelse(is.infinite(log(n.freq.pre)), NA, log(n.freq.pre))))
d$s.cmu.phon <- with(d,scale(cmu.phon))
d$s.uttlen.pre <- with(d,scale(uttlen.pre))

# Prepare the spatial, temporal and linguistic topic distinctiveness predictors
d$srl.sp.KL <- scale(resid(lm(log(sp.KL) ~ log(sp.count), data=d, na.action=na.exclude)))
d$srl.temp.KL <- scale(resid(lm(log(temp.KL) ~ log(temp.count), data=d, na.action=na.exclude)))
d$srl.topic.KL <- scale(resid(lm(log(topic.KL) ~ log(topic.count.ep),data=d, na.action=na.exclude)))

# Prepare the control predictors: imageability, concreteness, familiarity
d$s.mrc.imag <- with(d, scale(mrc.imag))
d$s.mrc.conc <- with(d, scale(mrc.conc))
d$s.mrc.fam <- with(d, scale(mrc.fam))

#### Direct predictor correlations with AoA ----
cor.test(~ s.cmu.phon + aoa, data=d) # 0.249
cor.test(~ s.uttlen.pre + aoa, data=d) # 0.192
cor.test(~ sln.freq.pre + aoa, data=d) # -0.177, same as sqrt(summary(lm.f)$r.squared) for regression model below

cor.test(~ srl.sp.KL + aoa, data=d) # -0.405
cor.test(~ srl.temp.KL + aoa, data=d) # -0.335
cor.test(~ srl.topic.KL + aoa, data=d) # -0.277

#### SI section 4.4: Control analyses ----

# Correlations between imageability, concreteness, familiarity, frequency and distinctiveness predictors
(mrc.corrs <- rcorr(as.matrix(d[,c("s.mrc.imag","s.mrc.conc","s.mrc.fam",
                                   "sln.freq.pre","srl.sp.KL","srl.temp.KL","srl.topic.KL")])))

# Number of common words when including imageability, concreteness and familiarity: 430
cc <- complete.cases(d[,c("s.mrc.imag","s.mrc.conc","s.mrc.fam","sln.freq.pre",
                          "srl.sp.KL","srl.temp.KL","srl.topic.KL","s.cmu.phon","s.uttlen.pre")])
summary(cc)

# Impact of imageability and concreteness on baseline + spatial model
(classic.lm.space <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.sp.KL,data=subset(d,cc)))
(classic.lm.space.imag <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.sp.KL+s.mrc.imag,data=subset(d,cc)))
(classic.lm.space.conc <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.sp.KL+s.mrc.conc,data=subset(d,cc)))

# Impact of imageability and concreteness on baseline + temporal model
(classic.lm.time <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.temp.KL,data=subset(d,cc)))
(classic.lm.time.imag <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.temp.KL+s.mrc.imag,data=subset(d,cc)))
(classic.lm.time.conc <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.temp.KL+s.mrc.conc,data=subset(d,cc)))

# Impact of imageability and concreteness on baseline + linguistic model
(classic.lm.ling <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.topic.KL,data=subset(d,cc)))
(classic.lm.ling.imag <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.topic.KL+s.mrc.imag,data=subset(d,cc)))
(classic.lm.ling.conc <- lm(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.topic.KL+s.mrc.conc,data=subset(d,cc)))

#### Build regression models, compute predicted AoA, and construct tables for SI section 4.3 ----
lm.f <- lm(aoa ~ sln.freq.pre, data=d, na.action=na.exclude)
lm.f.p <- lm(aoa ~ sln.freq.pre + s.cmu.phon, data=d, na.action=na.exclude)
lm.f.p.m <- lm(aoa ~ sln.freq.pre + s.cmu.phon + s.uttlen.pre, data=d, na.action=na.exclude)
lm.f.p.m.sp <- lm(aoa ~ sln.freq.pre + s.cmu.phon + s.uttlen.pre + srl.sp.KL, data=d, na.action=na.exclude)
lm.f.p.m.topic <- lm(aoa ~ sln.freq.pre + s.cmu.phon + s.uttlen.pre + srl.topic.KL, data=d, na.action=na.exclude)
lm.f.p.m.temp <- lm(aoa ~ sln.freq.pre + s.cmu.phon + s.uttlen.pre + srl.temp.KL, data=d, na.action=na.exclude)

# only for SI
lm.f.p.m.all <- lm(aoa ~ sln.freq.pre + s.cmu.phon + s.uttlen.pre + 
                     srl.temp.KL + srl.sp.KL + srl.topic.KL, data=d, na.action=na.exclude)

# generate SI latex - Tables 1-5
xtable(lm.f.p.m)
xtable(lm.f.p.m.sp)
xtable(lm.f.p.m.temp)
xtable(lm.f.p.m.topic)
xtable(lm.f.p.m.all)

# Add to dataset and reshape
d$Spatial <- fitted(lm.f.p.m.sp)
d$Temporal <- fitted(lm.f.p.m.temp)
d$Topic <- fitted(lm.f.p.m.topic)
d$Frequency <- fitted(lm.f)
d$Phonemes <- fitted(lm.f.p)
d$MLU <- fitted(lm.f.p.m)

d$base <- fitted(lm.f.p.m)

md <- melt(d, id.vars=c("word","aoa","small.cat","base"),
           measure.vars=c("Frequency","Phonemes","MLU","Spatial","Temporal","Topic"),
           variable.name="model",
           value.name="pred.aoa")

md$delta[md$model %in% c("Frequency")] <- d$Frequency - mean(d$aoa)
md$delta[md$model %in% c("Phonemes")] <- d$Phonemes - d$Frequency
md$delta[md$model %in% c("MLU")] <- d$MLU - d$Phonemes
md$delta[md$model %in% c("Spatial","Temporal","Topic")] <- 
  md$pred.aoa[md$model %in% c("Spatial","Temporal","Topic")] - md$base[md$model %in% c("Spatial","Temporal","Topic")]

#### Figure 1A and 1B: Model coefficient plot ----
datas.names <- c("All","Nouns","Predicates","Closed Class Words")
formulas.names <- c("None","Spatial","Temporal","Linguistic")
datas <- list(d,
              subset(d,small.cat=="Nouns"),
              subset(d,small.cat=="Adjectives"|small.cat=="Verbs"),
              subset(d,small.cat=="Closed class"))
formulas <- list(aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre,
                 aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.sp.KL,
                 aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.temp.KL,
                 aoa ~ s.cmu.phon+sln.freq.pre+s.uttlen.pre+srl.topic.KL)
counts <- c(nrow(d),
            nrow(subset(d,small.cat=="Nouns")),
            nrow(subset(d,small.cat=="Verbs" | small.cat=="Adjectives")),
            nrow(subset(d,small.cat=="Closed class")))

coefs <- local({
  coefs <- data.frame()
  for (i in 1:length(datas)) {
    for (f in 1:length(formulas)) {
      print(paste("******",datas.names[i],formulas.names[f],"*******"))
      l <- lm(formulas[[f]],data=datas[[i]])
      cs <- round(coef(l),2)
      ps <- round(summary(l)$coef[,4],3)
      new.coefs <- data.frame(dataset=datas.names[i],
                              predictor=formulas.names[f],
                              pred=names(coef(l)),
                              coef=coef(l),
                              se=coef(summary(l))[, 2])
      
      if (f==1) {
        new.coefs <- rbind(new.coefs, data.frame(dataset=datas.names[i],
                                                 predictor=formulas.names[f],
                                                 pred="None",
                                                 coef=0,
                                                 se=0))
      }
      coefs <- rbind(coefs, new.coefs)
    }
  }
  return(coefs)
})

# Set up common color scale and labels
color.scale <- scale_fill_manual(name="", values=c("lightgray","gray","black","white","#ff2f2f","#fa6a03","#d8489d"),
                                 labels=c("# Phonemes","MLU","Frequency","",
                                          "Spatial distinctiveness",
                                          "Temporal distinctiveness","Linguistic distinctiveness"))

# Plot figure 1A
quartz(width=6,height=2.75)
qplot(predictor, coef, 
      position="dodge", fill=pred, facets=~dataset, geom="bar", stat="identity",
      data=subset(coefs,dataset == "All" & pred!="(Intercept)")) + 
  geom_linerange(aes(ymin=coef-se,ymax=coef+se), position=position_dodge(width=.95)) + 
  ylim(c(-25,25)) + 
  xlab("Distinctiveness Model") + ylab("Coefficient (Days AoFP/SD)") +
  color.scale + theme_bw() + theme(legend.key = element_rect(color="white",fill="white"))

# Plot figure 1B
quartz(width=10.1,height=2.75)
qplot(predictor, coef, 
      position="dodge", fill=pred, facets=~dataset, geom="bar", stat="identity",
      data=subset(coefs,dataset != "All" & pred != "(Intercept)")) + 
  geom_linerange(aes(ymin=coef-se,ymax=coef+se),position=position_dodge(width=.95)) + 
  xlab("Distinctiveness Model") + ylab("Coefficient (Days AoFP/SD)") + 
  color.scale + theme_bw() + theme(legend.position = "none")

#### Figure 2: Word scatter plot for each model ----

words.fig2 <- c("bath","beautiful","breakfast","bye","car",
                "cat","color","come","cow","fish","hi","if",
                "kick","moon","motorcycle","no","read","scared",
                "something","that","wheel","with","yes","you")

md$model <- revalue(md$model, c("Frequency"="Frequency",
                                "Phonemes"="Frequency + Phonemes",
                                "MLU"="Frequency + Phonemes + MLU (Base)",
                                "Spatial"="Base + Spatial Distinctiveness",
                                "Temporal"="Base + Temporal Distinctiveness",
                                "Topic"="Base + Linguistic Distinctiveness"))

# Convert all dates from days to months (month \approx 30.3 days for display purposes)
md$aoa <- md$aoa / 30.3
md$pred.aoa <- md$pred.aoa / 30.3
md$delta <- md$delta / 30.3
# Identify words to call out visually in scatter plot
md <- ddply(md,.(model), transform, selected = word %in% words.fig2)

quartz(width=9,height=7)
ggplot(md,aes(x=aoa, y=pred.aoa)) + 
  geom_point(color="gray80",alpha=.3) + 
  geom_abline(slope=1,intercept=0,colour="gray",lty=2) +
  geom_smooth(method="lm",colour="gray",fill="lightgray",lty=3) + 
  facet_wrap(~model) +
  geom_linerange(aes(ymax=pred.aoa-delta,ymin=pred.aoa,colour=small.cat),alpha=.5,data=subset(md,selected)) + 
  geom_point(aes(colour=small.cat),alpha=.9,data=subset(md,selected)) +
  xlab("True Age of First Production (Months)") + 
  ylab("Predicted Age of First Production (Months)") + 
  xlim(9,25.5) + ylim(14,23) +
  theme_bw() + scale_color_discrete("") +
  theme(legend.position="top") +
  geom_text(subset=.(selected), aes(label=word,colour=small.cat), alpha=.75,
            show_guide=FALSE,angle=0,size=3,hjust=0,vjust=0) + 
  theme(panel.grid=element_blank())

#### Figure 3: Word distinctiveness plots ----

# Load the functions used for creating this figure into an environment
env.fig3 <- new.env()
sys.source(file="R_code/HSP_plotfunctions.R",env.fig3)

attach(env.fig3)

words.fig3 <- c("beautiful","breakfast","bye","car","cow","fish","hi","kick","moon","with")
for (wrd in words.fig3) {
  quartz(width=6,height=.9)
  draw_panel(wrd)
}

detach(env.fig3)
rm(env.fig3)

#### SI Figure 9: Pairwise predictor correlation matrix ----

vars <- c("s.cmu.phon","sln.freq.pre","s.uttlen.pre","srl.sp.KL","srl.temp.KL","srl.topic.KL")
var.names <- c("phonemes","log(freq)","utt length","spatial","temporal","topical")

# Compute the pairwise predictor correlations
cs <- rcorr(as.matrix(d[,vars]))
corrs <- cs$r
corrs[upper.tri(corrs,diag=TRUE)] <- NA
corrs <- melt(corrs,varnames=c("var1","var2"),value.name="corr")
corrs <- na.omit(corrs)
corrs$var1 <- mapvalues(corrs$var1, vars, var.names)
corrs$var2 <- mapvalues(corrs$var2, vars, var.names)
corrs$rounded <- round(corrs$corr,digits=2)

# Draw the figure
quartz(width=6,height=5)
qplot(var1,var2,fill=abs(corr),geom="tile",data=corrs) +
  geom_text(aes(label=rounded),colour="white",size=4) +
  scale_fill_gradient(limits=c(0,.6),name="correlation") +
  xlab("") + ylab("") + theme_bw() + 
  theme(legend.position="none", panel.background=element_blank(),
        panel.border=element_blank(), axis.ticks=element_blank(),        
        plot.background=element_blank(),panel.grid=element_blank())
