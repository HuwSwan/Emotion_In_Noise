## For measuring how statistics of classes load across emotions.
# Final script in process. First run vin24_oin_LPA then paper_window_statistics_verification for both norm and scrambled.


rm(list = ls())

library(lme4)
library(afex)
library(nlme)
library(psych)
library(lsmeans)
library(ggplot2)
library(stringi)
library(stringr)
library(plyr)
library(ISLR)
library(multcomp)
library(sjPlot)
library(purrr)
#install.packages("ggeffects")
library(ggeffects)
library(wesanderson)
library(ggsci)
library(svglite)



normPath <- 'X:/PhD/03-Original_OIM/02-analysis/2024_dwgp_final_stage_R/lmeGlobalWeight/normal_lmeTable.csv'
scramPath <- 'X:/PhD/03-Original_OIM/02-analysis/2024_dwgp_final_stage_R/lmeGlobalWeight/scrambled_lmeTable.csv'

normDat  <- read.csv(normPath)
scramDat <- read.csv(scramPath)


## NORMAL

#normDat$hzCentres = as.factor(normDat$hzCentres)
normDat$part = as.factor(normDat$part) # Dummy variable for each individual 
normDat$emos = as.factor(normDat$emos)

##STANDARDISE VARIABLES
#normDat$GWeight <- scale(normDat$global_weights, center = TRUE, scale = TRUE)



#normDat$hzCentres = as.factor(normDat$hzCentres)
scramDat$part = as.factor(scramDat$part) # Dummy variable for each individual 
scramDat$emos = as.factor(scramDat$emos)

##STANDARDISE VARIABLES
#scramDat$GWeight <- scale(scramDat$global_weights, center = TRUE, scale = TRUE)



## New labels for graphing.
labels <- c("Neutral", "Angry", "Joyful")
normDat$Emotions <- factor(normDat$emos, levels = c(1, 2, 3), labels = labels)


labels <- c("Neutral", "Angry", "Joyful")
scramDat$Emotions <- factor(scramDat$emos, levels = c(1, 2, 3), labels = labels)


### THIS DOESNT INCLUDE CLASS SO IT'S GOING TO BE ONLY THE GLOBAL ANALYSIS
nrm <- rep(1, nrow(normDat))
scr <- rep(1, nrow(scramDat))+1

scramDat$Stimuli <- as.factor(scr)
normDat$Stimuli <- as.factor(nrm)

boundData <- rbind(normDat,scramDat)

boundData$GWeight <-boundData$global_weights

#Scale
boundData$GWeight <- scale(boundData$global_weights, center = TRUE, scale = TRUE)
boundData$snr <- scale(boundData$snr, center = TRUE, scale = TRUE)
boundData$snr_rt <- scale(boundData$snr_rt, center = TRUE, scale = TRUE)

labels <- c("Normal", "Scrambled")
boundData$Stimuli <- factor(boundData$Stimuli, levels = c(1, 2), labels = labels)

# palettes
pal <- wes_palette("Darjeeling2", 5, type = "discrete")
br <- wes_palette("BottleRocket2", 5, type = "discrete")
cmPal <- c(br[3],br[2],br[1])
axT = 0.75
lT = 0.5

setwd('X:/PhD/03-Original_OIM/02-analysis/publication_figures_new/COMBINED_GLM')

# -----------------------GLMS TO TELL STORY--------------------------------###

#Scaled because different experiments....

#Local weights
gweightModel <- lmer(snr ~   GWeight  * Stimuli * Emotions  + (1|part), data=boundData)

gweightONLYModel <- lmer(snr ~   GWeight  * Stimuli  + (1|part), data=boundData)

emoModel <- lmer(snr ~    Emotions * Stimuli + (1|part), data=boundData)


plot_model(gweightModel, type = "int")
plot_model(emoModel, type = "int")


performance::r2(emoModel)
performance::r2(gweightModel)
performance::r2(gweightONLYModel)


weightSummary <- summary(gweightModel)
emoSummary <- summary(emoModel)
onlySummary <- summary(gweightONLYModel)

modComp <- anova(gweightModel,emoModel,gweightONLYModel) #Model comparison. 

write.csv((modComp), "model_comparison.csv") # export

##Multicomp


emMComp <- emmeans(gweightModel, list(pairwise ~ Emotions), adjust = "tukey")


# Save models
coeffs <- coef(summary(gweightModel)) # get estimates, etc...
p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2 # add the much disputed p-values

coeffsp <- cbind(coeffs, "p value" = round(p,3)) # combine it into one object
write.csv(coeffsp, "full_LMER.csv") # export
# EMO ONLY
coeffs <- coef(summary(emoModel)) # get estimates, etc...
p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2 # add the much disputed p-values

coeffsp <- cbind(coeffs, "p value" = round(p,3)) # combine it into one object
write.csv(coeffsp, "emo_half_LMER.csv") # export
# Weight Only
coeffs <- coef(summary(gweightONLYModel)) # get estimates, etc...
p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2 # add the much disputed p-values

coeffsp <- cbind(coeffs, "p value" = round(p,3)) # combine it into one object
write.csv(coeffsp, "weight_half_LMER.csv") # export


## PLOTS
normPlots <- plot_model(
  gweightModel,
  type = "int",
  alpha = 0.4,
  dodge = 0.7,
  dot.size = 1
) %>% 
  map(function(plot) {
    # Apply color palette for Class distinction
    plot + scale_color_brewer(palette = "Set1")
    # Ensure uniform linetype, solid by default
  })

# All interaction
p1 <- normPlots[[4]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                               axis.ticks.length = unit(.15,"cm"),
                                               text = element_text(size=10,family = "serif",colour = "black"),
                                               axis.ticks = element_line(size = axT,colour = 'black'),
                                               legend.position = "none") + 
  scale_linetype_discrete(guide = "none")+
  xlab("STEP Distortion") + ylab("SNR")+
  ggtitle(NULL)+
  scale_fill_jco() + 
  scale_color_jco()
print(p1)

#Class
p2 <- normPlots[[1]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                               axis.ticks.length = unit(.15,"cm"),
                                               text = element_text(size=10,family = "serif",colour = "black"),
                                               axis.ticks = element_line(size = axT,colour = 'black'),
                                               legend.position = "none") + 
  scale_linetype_discrete(guide = "none")+
  xlab("STEP Distortion") + ylab("SNR")+
  ggtitle(NULL)+
  scale_fill_jco() + 
  scale_color_jco()

print(p2)

# Emotion Interaciton
p3 <- normPlots[[2]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                               axis.ticks.length = unit(.15,"cm"),
                                               text = element_text(size=10,family = "serif",colour = "black"),
                                               axis.ticks = element_line(size = axT,colour = 'black'),
                                               legend.position = "none") + 
  scale_linetype_discrete(guide = "none")+
  xlab("STEP Distortion") + ylab("SNR")+
  ggtitle(NULL)+
  scale_fill_manual(values = cmPal) + 
  scale_color_manual(values = cmPal)

print(p3)

## Plot Compare means.
p4 <- normPlots[[3]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                               axis.ticks.length = unit(.15,"cm"),
                                               text = element_text(size=10,family = "serif",colour = "black"),
                                               axis.ticks = element_line(size = axT,colour = 'black'),
                                               legend.position = "none") + 
  scale_linetype_discrete(guide = "none")+
  xlab("Stimuli") + ylab("SNR")+
  ggtitle(NULL)+ 
  scale_fill_manual(values = cmPal) + 
  scale_color_manual(values = cmPal)

print(p4)

## --------------------------------SAVE ALL PLOTS-------------------------------
ggsave("Full_LMER_interaction_all.svg",plot = p1,units = "cm", dpi = 500,width = 10,height = 5)
ggsave("Full_LMER_interaction_class.svg",plot = p2,units = "cm", dpi = 500,width = 5,height = 5)
ggsave("Full_LMER_interaction_emotion.svg",plot = p3,units = "cm", dpi = 500,width = 5,height = 5)
ggsave("Full_LMER_pairwise.svg",plot = p4,units = "cm", dpi = 500,width = 5,height = 5)

## PRINT LEGENDS
classLeg <- normPlots[[4]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                                     axis.ticks.length = unit(.15,"cm"),
                                                     text = element_text(size=10,family = "serif",colour = "black"),
                                                     axis.ticks = element_line(size = axT,colour = 'black'),
                                                     legend.position = "bottom") + 
  scale_linetype_discrete(guide = "none")+
  xlab("STEP Distortion") + ylab("SNR")+
  ggtitle(NULL)+
  scale_fill_jco() + 
  scale_color_jco()

EmoLeg <- normPlots[[3]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                                   axis.ticks.length = unit(.15,"cm"),
                                                   text = element_text(size=10,family = "serif",colour = "black"),
                                                   axis.ticks = element_line(size = axT,colour = 'black'),
                                                   legend.position = "bottom") + 
  scale_linetype_discrete(guide = "none")+
  xlab("STEP Distortion") + ylab("SNR")+
  ggtitle(NULL)+
  scale_fill_manual(values = cmPal) + 
  scale_color_manual(values = cmPal)

## SAVE LEGENDS
ggsave("Full_LMER_legend.svg",plot = classLeg,units = "cm", dpi = 500,width = 10,height = 7)
ggsave("Full_LMER_pair_legend.svg",plot = EmoLeg,units = "cm", dpi = 500,width = 10,height = 7)


## -----------------------Pairwise Comparisons---------------------------
EmoComp <- emmeans(gweightModel, list(pairwise ~ Emotions), adjust = "tukey")
StimuliComp <- emmeans(gweightModel, list(pairwise ~ Stimuli), adjust = "tukey")
AllComp <- emmeans(gweightModel, list(pairwise ~ Emotions * Stimuli), adjust = "tukey")


allCompSummary <- summary(AllComp)

# For comparisons within class for how they differ by emotion.

write.csv(allCompSummary$`emmeans of Emotions, Stimuli`, "Full_LMER_CombinedPairwise_emmeans.csv") # export
write.csv(allCompSummary$`pairwise differences of Emotions, Stimuli`, "Full_LMER_CombinedPairwise.csv") # export


write.csv(summary(StimuliComp$`emmeans of Stimuli`), "Full_LMER_StimuliPairwise_emmeans.csv") # export
write.csv(summary(StimuliComp$`pairwise differences of Stimuli`), "Full_LMER_StimuliPairwise.csv") # export

write.csv(summary(EmoComp$`emmeans of Emotions`), "Full_LMER_EmoPairwise_emmeans.csv") # export
write.csv(summary(EmoComp$`pairwise differences of Emotions`), "Full_LMER_EmoPairwise.csv") # export


######################## PLOT THE EMO MODEL TO DEMONSTRATE HOW MUCH THE EFFECT SHRINKS WHEN YOU REMOVE STEP DISTORTION



AllCompEMMOD <- emmeans(emoModel, list(pairwise ~ Emotions * Stimuli), adjust = "tukey")




write.csv(summary(AllCompEMMOD$`emmeans of Emotions, Stimuli`), "emonly_LMER_Pairwise_emmeans.csv") # export
write.csv(summary(AllCompEMMOD$`pairwise differences of Emotions, Stimuli`), "emonly_LMER_Pairwise.csv") # export


emoPlot <- list(plot_model(
  emoModel,
  type = "eff",
  terms = c("Stimuli","Emotions"),
  alpha = 0.4,
  dodge = 0.7,
  dot.size = 1
)) %>% 
  map(function(plot) {
    # Apply color palette for Class distinction
    plot + scale_color_brewer(palette = "Set1")
    # Ensure uniform linetype, solid by default
  })

# Emotion Interaciton
eP1 <- emoPlot[[1]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                               axis.ticks.length = unit(.15,"cm"),
                                               text = element_text(size=10,family = "serif",colour = "black"),
                                               axis.ticks = element_line(size = axT,colour = 'black'),
                                               legend.position = "none") + 
  scale_linetype_discrete(guide = "none")+
  xlab("STEP Distortion") + ylab("SNR")+
  ggtitle(NULL)+
  scale_fill_manual(values = cmPal) + 
  scale_color_manual(values = cmPal)


print(eP1)

## 
ggsave("emo_Model_Pairs.svg",plot = eP1,units = "cm", dpi = 500,width = 5,height = 5)


##################### Descriptives

nEmDesc <- aov_4(snr_rt ~ Emotions * (Emotions|part), data=normDat)



nEmo <- emmeans(nEmDesc, list(pairwise ~ Emotions), adjust = "tukey")

mean(normDat$RT_slope[normDat$emos==1])
mean(normDat$RT_slope[normDat$emos==1])


#################### DEMOSNTRATE INCREASED STEP PRESENCE BY CONDITON

 gweightModelSTEP <- lmer(GWeight ~ Stimuli * Emotions  + (1|part), data=boundData)
 
 summary(gweightModelSTEP)
 performance::r2((gweightModelSTEP))
 
 # Save model.
 coeffs <- coef(summary(gweightModel)) # get estimates, etc...
 p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2 # add the much disputed p-values
 
 coeffsp <- cbind(coeffs, "p value" = round(p,3)) # combine it into one object
 write.csv(coeffsp, "distortio_model.csv") # export
 
 
 distPlot <- list(plot_model(
   gweightModelSTEP,
   type = "int",
   terms = c("Stimuli","Emotions"),
   alpha = 0.4,
   dodge = 0.7,
   dot.size = 1
 )) %>% 
   map(function(plot) {
     # Apply color palette for Class distinction
     plot + scale_color_brewer(palette = "Set1")
     # Ensure uniform linetype, solid by default
   })
 
 # Emotion Interaciton
 dP1 <- distPlot[[1]] + theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = axT),
                                               axis.ticks.length = unit(.15,"cm"),
                                               text = element_text(size=10,family = "serif",colour = "black"),
                                               axis.ticks = element_line(size = axT,colour = 'black'),
                                               legend.position = "none") + 
   scale_linetype_discrete(guide = "none")+
   xlab("Stimuli") + ylab("STEP Distortion")+
   ggtitle(NULL)+
   scale_fill_manual(values = cmPal) + 
   scale_color_manual(values = cmPal)

 ggsave("distortion_Model_Pairs.svg",plot = dP1,units = "cm", dpi = 500,width = 5,height = 5)
 
 
 
 
 ## -----------------------Pairwise Comparisons---------------------------
 EmoDComp <- emmeans(gweightModelSTEP, list(pairwise ~ Emotions), adjust = "tukey")
 StimuliDComp <- emmeans(gweightModelSTEP, list(pairwise ~ Stimuli), adjust = "tukey")
 AllDComp <- emmeans(gweightModelSTEP, list(pairwise ~ Emotions * Stimuli), adjust = "tukey")
 
 
 # For comparisons within class for how they differ by emotion.
 #all
 write.csv(summary(AllDComp$`emmeans of Emotions, Stimuli`), "distortion_CombinedPairwise_emmeans.csv") # export
 write.csv(summary(AllDComp$`pairwise differences of Emotions, Stimuli`), "distortion_CombinedPairwise.csv") # export
 
 #stim
 write.csv(summary(StimuliDComp$`emmeans of Stimuli`), "distortion_StimuliPairwise_emmeans.csv") # export
 write.csv(summary(StimuliDComp$`pairwise differences of Stimuli`), "distortion_StimuliPairwise.csv") # export
 #emo
 write.csv(summary(EmoDComp$`emmeans of Emotions`), "distortion_EmoPairwise_emmeans.csv") # export
 write.csv(summary(EmoDComp$`pairwise differences of Emotions`), "distortion_EmoPairwise.csv") # export
 
 


