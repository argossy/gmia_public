# Author: Samuel Rosin (srosin@live.unc.edu)
# Date last modified: October 12, 2020

# if any of these packages are not installed, install them by running commands 
# install.packages("dplyr"), install.packages("plyr"), install.packages("here"), install.packages("ggplot")
library(dplyr)
library(plyr)
library(here)
library(ggplot2)

# 1. Load and clean data 
# the here package should "just work" when all files are in the same directory.
# if it doesn't, see https://github.com/jennybc/here_here for help. 
neo.filename <- here::here("neo_jan18.csv")
yr.filename <- here::here("yr_jan18.csv")

neo <- read.csv(neo.filename) %>%
  subset(Sub.ID!=2 & Sub.ID !=31) %>% 
  plyr::rename(c("SHANNON" = "Shannon.Entropy", "PD_WHOLE_TREE" = "Phylogenetic.Diversity",
                     "OBSERVED_SPECIES" = "Observed.Species", "CHAO" = "Chao1.Index"))

yr <- read.csv(yr.filename) %>%
  subset(Sub.ID!=2 & Sub.ID !=31) %>% 
  plyr::rename(c("SHANNON" = "Shannon.Entropy", "PD_WHOLE_TREE" = "Phylogenetic.Diversity",
                 "OBSERVED_SPECIES" = "Observed.Species", "CHAO" = "Chao1.Index"))

#grab relevant columns
neo.small <- neo[,1:7]
yr.small <- yr[,1:7]
neo.small$Time <- 1
yr.small$Time <- 2

#dataframe for spaghetti plots
spaghetti.dat <- rbind(neo.small,yr.small) %>% 
  group_by(Sub.ID) %>% 
  filter(n()>1)

#############
#
#  2. Make spaghetti plots
#

shannon.spaghetti <- ggplot(data=spaghetti.dat,
                            aes(x=Time,y=Shannon.Entropy,group=Sub.ID)) + 
  geom_point(shape=20, color = "grey53") + geom_line(size=.3, color = "grey53") + 
  ggtitle("Shannon Entropy Over Time") + 
  theme(plot.title = element_text(size=11,hjust=0.5)) + 
  scale_x_discrete(limits=factor(c(1,2)),
                   labels=c("Two-Week Neonate","One-Year Infant")) + 
  xlab("Visit Age") + 
  ylab("Shannon Entropy")
shannon.spaghetti

pd.spaghetti <- ggplot(data=spaghetti.dat,
                       aes(x=Time,y=Phylogenetic.Diversity,group=Sub.ID)) + 
  geom_point(shape=20, color = "grey53") + geom_line(size=.3, color = "grey53") + 
  ggtitle("Phylogenetic Diversity Over Time") + 
  theme(plot.title = element_text(size=11,hjust=0.5)) + 
  scale_x_discrete(limits=factor(c(1,2)),
                   labels=c("Two-Week Neonate","One-Year Infant")) + 
  xlab("Visit Age") + 
  ylab("Phylogenetic Diversity")
pd.spaghetti

obs.spaghetti <- ggplot(data=spaghetti.dat,
                        aes(x=Time,y=Observed.Species,group=Sub.ID)) + 
  geom_point(shape=20, color = "grey53") + geom_line(size=.3, color = "grey53") + 
  ggtitle("Observed Species Over Time") + 
  theme(plot.title = element_text(size=11,hjust=0.5)) + 
  scale_x_discrete(limits=factor(c(1,2)),
                   labels=c("Two-Week Neonate","One-Year Infant")) + 
  xlab("Visit Age") + 
  ylab("Observed Species")
obs.spaghetti

chao.spaghetti <- ggplot(data=spaghetti.dat,
                         aes(x=Time,y=Chao1.Index,group=Sub.ID)) + 
  geom_point(shape=20, color = "grey53") + geom_line(size=.3, color = "grey53") + 
  ggtitle("Chao1 Index Over Time") + 
  theme(plot.title = element_text(size=11,hjust=0.5)) + 
  scale_x_discrete(limits=factor(c(1,2)),
                   labels=c("Two-Week Neonate","One-Year Infant")) + 
  xlab("Visit Age") + 
  ylab("Chao1 Index")
chao.spaghetti
