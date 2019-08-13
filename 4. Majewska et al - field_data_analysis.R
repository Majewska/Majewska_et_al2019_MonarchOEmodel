## 
## Code to accompany: Multiple transmission routes sustain high prevalence of a virulent parasite in a butterfly host 
## (Majewska, Sims, Schneider, Altizer, Hall 2019, Proceedings B;  DOI: 10.1098/rspb.2019.1630) 
## Code was written by AAM and contains field data analyses
##

####
####last updated 7.26.19


# Part 4 Field Data Analysis

rm(list=ls()) 
graphics.off()


require(doBy)
require(lme4)
require(visreg)
require(dplyr)
require(GGally)
require(MuMIn)
require(multcomp)


###########################################################################

# Infection data 

larv<-read.csv("data/Larva _data.csv", na.strings=c(""," ","NA"))
adults<-read.csv("data/Adult_data.csv", na.strings=c(""," ","NA"))

# select the needed columns from each 
l<-dplyr::select(larv, month, week, Sex, OE.Y.N, stage, Animal.ID)
a<-dplyr::select(adults, month, week, Sex, OE.Y.N, stage, Animal.ID )


# Join adult and larva data
df<-rbind(l,a)


# Infection data analysis

# Is OE infection status status influenced by sex, month of collection and collection stage?
summary(model1.2<-glm(OE.Y.N ~ Sex + month * stage, family = binomial, data = df, na.action = na.omit))


# Calculte proportion of infected males and females 

sumpart<-df %>% 
  na.omit() %>% 
  group_by(Sex,OE.Y.N) %>% 
summarize(totalsampled = n()) 

tots<-aggregate(totalsampled~Sex, data=sumpart, FUN=sum) 


tot2<-left_join(sumpart,tots, by = c("Sex"))
tot2$total<-tot2$totalsampled.y

tot2$Proportion<-tot2$totalsampled.x/tot2$total


datainf<-subset(tot2,OE.Y.N==1)
datainf$healthy=datainf$totalsampled.x
datainf$infected=datainf$totalsampled.y


OE.inf.all <- datainf %>% 
  group_by(Sex) %>% 
summarise(mean=mean(Proportion), n =sum(total))

OE.inf.all$se<-sqrt(OE.inf.all$mean*(1-OE.inf.all$mean)/OE.inf.all$n)
OE.inf.all

###########################################################################

# Adult sprore transfer data 

transfer<-read.csv("data/Adult_spore_transfer.csv")

table(transfer$Sex)
uniques<-unique(transfer$Animal.ID)

transfer$dates <- transfer$Released
betterDates <- as.Date(transfer$dates,
                       format = "%m/%d/%y")
transfer$months<-format(betterDates, format="%B")

hist(transfer$spores.acquired)


# Analysis of adult trasfer data
# Is spore transfer status influenced by month of collection and sex?
summary(model1.5<-glmer(got.spores ~  month + Sex + (1|Animal.ID), family = binomial,data = transfer, control = glmerControl(optimizer = "bobyqa"),na.action = na.omit))

# Is the number of spores transfered influenced by month of collection and sex?
summary(model1.5<-glmer.nb(spores.acquired ~  month + Sex + (1|Animal.ID), data = transfer, control = glmerControl(optimizer = "bobyqa"),na.action = na.omit))

 
# tbl = table(spores$got.spores, spores$Sex) 
# chisq.test(tbl)
# 
# tbl2 = table(spores$spores.acquired, spores$Sex) 
# chisq.test(tbl2)


###########################################################################

# Milkweed data 
Milkweed<-read.csv("data/Milkweed.csv")

# Analysis of milkweed data
summary(model1.6<-glmer(infections.status ~  Month + (1|Plant.ID), family = binomial,data = Milkweed, control = glmerControl(optimizer = "bobyqa"),na.action = na.omit))


