
#Latent class analysis for the OIM window data.

#install.packages("tidyLPA")
#install.packages("ggplot2")

rm(list = ls())

library(tidyLPA)
library(ggplot2)


pfile <- 'X:/PhD/03-Original_OIM/02-analysis/2024-verification-dwgp_collation/normal/window_features_speech/summary_windows.csv'


datx <- read.csv(pfile)
datx< - data.frame(datx)

#estimate - c(1:3) is which variables are excluded. Exclude ref, emo, time.RT
c1 <- estimate_profiles(datx[,-c(1:4,9)],1)
c2 <- estimate_profiles(datx[,-c(1:4,9)],2)
c3 <- estimate_profiles(datx[,-c(1:4,9)],3)
c4 <- estimate_profiles(datx[,-c(1:4,9)],4)
c5 <- estimate_profiles(datx[,-c(1:4,9)],5)
c6 <- estimate_profiles(datx[,-c(1:4,9)],6)
c7 <- estimate_profiles(datx[,-c(1:4,9)],7)
 
# c8 <- estimate_profiles(datx[,-c(1:2)],8)


# Z score my vectors 

#fit
fits <- data.frame(rbind(t(data.frame(c1$model_1_class_1$fit)),
                         t(data.frame(c2$model_1_class_2$fit)),
                         t(data.frame(c3$model_1_class_3$fit)),
                         t(data.frame(c4$model_1_class_4$fit)),
                         t(data.frame(c5$model_1_class_5$fit)),
                         t(data.frame(c6$model_1_class_6$fit)),
                         t(data.frame(c7$model_1_class_7$fit))))
                         #t(data.frame(c8$model_1_class_8$fit))))

pal <- wes_palette("Darjeeling2", 5, type = "discrete")
br <- wes_palette("BottleRocket2", 5, type = "discrete")
cmPal <- c(br[3],br[2],br[1])
axT = 0.75
lT = 0.5

#plot fits.
p1 <- ggplot(data = fits) +
  geom_point(aes(x = Classes, y = AIC), color = cmPal[3], size = 2)+
  geom_point(aes(x = Classes, y = BIC), color = cmPal[2], size = 2) +
  geom_point(aes(x = Classes, y = CAIC), color = cmPal[1], size = 2) +
  geom_line(aes(x = Classes, y = AIC), color = cmPal[3], size = 1)+
  geom_line(aes(x = Classes, y = BIC), color = cmPal[2], size = 1) +
  geom_line(aes(x = Classes, y = CAIC), color = cmPal[1], size = 1) +
  theme_minimal() + xlab("K models of 1:7 classes") +
  ylab("AIC, BIC, and CIAC")

  
  
setwd('X:/PhD/03-Original_OIM/02-analysis/publication_figures_new/dwgp/LPA')

p1 <- p1+ theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                                axis.ticks.length = unit(.15,"cm"),
                                                text = element_text(size=10,family = "serif",colour = "black"),
                                                axis.ticks = element_line(size = axT,colour = 'black'),
                                                legend.position = "none") + 
  scale_linetype_discrete(guide = "none")+
  ggtitle(NULL)+
  scale_fill_jco() + 
  scale_color_jco()

ggsave("AIC_normal.svg",plot = p1,units = "cm", dpi = 500,width = 5,height = 5)




##class solution push through


fitdat <- data.frame(c5$model_1_class_5$dff)

fitdat$ID <- row.names(datx)
## Degreees of freedom and membership probs.
probs <-data.frame(c5$model_1_class_5$dff)

#  Plots
# # # # #

#Save and export for matlab.
datx$Class <- as.factor(fitdat$Class)
datx$CPROB1 <- as.factor(probs$CPROB1)
datx$CPROB2 <- as.factor(probs$CPROB2)
datx$CPROB3 <- as.factor(probs$CPROB3)
datx$CPROB4 <- as.factor(probs$CPROB4)
datx$CPROB5 <- as.factor(probs$CPROB5)
#datx$CPROB6 <- as.factor(probs$CPROB6)
#datx$CPROB7 <- as.factor(probs$CPROB7)





saveStr <- 'X:/PhD/03-Original_OIM/02-analysis/2024-verification-dwgp_collation/normal/window_features_speech/window_data_classes_5class.csv'

write.csv(datx,saveStr)


datX <- read.csv(saveStr)

#size by hz centre.
ggplot(data = datx) +
  geom_point(aes(x = windowSize, y = hzCentres, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = windowSize, y = hzCentres, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

##size by weight
ggplot(data = datx) +
  geom_point(aes(x = windowSize, y = windowWeights, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = windowSize, y = windowWeights, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)


#size by lengths
ggplot(data = datx) +
  geom_point(aes(x = windowSize, y = tLengths, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = windowSize, y = tLengths, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

#size by heightshz.
ggplot(data = datx) +
  geom_point(aes(x = windowSize, y = hzHeights, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = windowSize, y = hzHeights, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

##lengths first.

#length bt hzcentre
ggplot(data = datx) +
  geom_point(aes(x = tLengths, y = hzCentres, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = tLengths, y = hzCentres, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

#lengths by height
ggplot(data = datx) +
  geom_point(aes(x = tLengths, y = hzHeights, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = tLengths, y = hzHeights, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

#lengths by weights
ggplot(data = datx) +
  geom_point(aes(x = tLengths, y = windowWeights, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = tLengths, y = windowWeights, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

## Hz centre first

# centre by height
ggplot(data = datx) +
  geom_point(aes(x = hzHeights, y = hzCentres, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = hzHeights, y = hzCentres, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

# centre by weights
ggplot(data = datx) +
  geom_point(aes(x = windowWeights, y = hzCentres, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = windowWeights, y = hzCentres, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)

#weights by height
ggplot(data = datx) +
  geom_point(aes(x = windowWeights, y = hzHeights, color = emo), alpha = 0.4) +
  geom_smooth(aes(x = windowWeights, y = hzHeights, color = emo, group = emo), method = "lm") +
  facet_wrap(.~Class)





# Emotion IDx is on the rows, Class is on the columns

table(datx$emo, datx$Class)



#summary
c3$model_1_class_3$estimates









