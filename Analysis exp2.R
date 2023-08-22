library(ggplot2)
library (lme4)
library(ptmixed)
library (MASS)
library(car)
library(fitdistrplus)
library (brms)
library(rstan)
library(agricolae)
library(magrittr)
library(multcomp)
library(emmeans)
library(repeated)
library(GLMMadaptive)
library(rcompanion)


setwd("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/1_Food calls Thyroptera/Analysis")


calls <- read.csv('Calls_produced.csv', header = TRUE, sep = ";", 
                  colClasses = c("Order"="character", "Bat_ID"="character", "Sex"="character", 
                                 "Age"="character", "Trial"="character", "Audio"="character", 
                                 "N_calls"="numeric", "Audio_length"="numeric", "Naïve"="character"))

data <- read.csv('E1_data.csv', header = TRUE, sep = ";", 
                 colClasses = c("id"="character", "trial"="character", "stimulus"="character", 
                                "occ"="numeric", "dur"="numeric", "first"="numeric", "first_graph"="character"))

data$stimulus <- factor(data$stimulus, levels=c('Silence', 'Pink noise', 'Chewing sounds', 'Food calls', "Chewing sounds + Food calls"))

data_2 <- read.csv('E2_data.csv', header = TRUE, sep = ";", 
                   colClasses = c("order"="numeric", "id"="character", "obs"="numeric", "trial"="character", "stimulus"="character", 
                                  "stim_order"="character", "pre_calls"="character", "pre_chew"="character", "pre_pn"="character",
                                  "playlist"="character", "dist_75_cm"="numeric", "dist_50_cm"="numeric"))

data_2$stimulus <- factor(data_2$stimulus, levels=c('Silence', 'Pink noise', 'Chewing sounds', 'Food calls', "Chewing sounds + Food calls"))


##__________________________________________________________________________________________________________##
#Calls produced

l <- ggplot(calls, aes(x=Sex, y=N_calls, fill=Naive)) + 
  geom_boxplot(alpha=0.6) +
  theme_classic()+
  labs(x="Sex", y = "Number of calls emitted")+ 
  facet_grid(vars(Trial))+
  guides(y.sec = guide_none("Trial")) +
  scale_fill_manual(values = c("#E89E16", "#59B29D"))
l

l + geom_point(colour="white", shape=21, size = 2, aes(x=Sex, y=N_calls, fill = factor(Naive)), position=position_jitterdodge()) 



ggsave("Figure 3_Number of calls emitted.tiff", 
       width = 5, height = 4, dpi = 600)



#Model to determine if sex, naive and trial affect call emitted

require(MuMIn)

options(na.action = "na.fail")


#Create global model with 4-way interactions
globalmodel_calls <- glmer(N_calls ~ (Sex + Trial + Naive)^3 + (1|Bat_ID), family = poisson, data=calls)


#Generate a table of models that include every possible combination of predictor variables
combinations_calls <- dredge(globalmodel_calls)


#Keep the best models based on delta values < 5
best_comb_calls <- subset(combinations_calls, delta < 5)

print(best_comb_calls)


#Estimate of the relative importance of each predictor variable 
imp <- sw(best_comb_calls)

imp


#Based on latter, the best model
best_calls <- glmer(N_calls ~ Naive + Sex + Trial + Naive*Trial + Sex*Trial + Naive*Sex + (1|Bat_ID), 
                    family = poisson, data=calls)



#Test effect of predictive variables
car::Anova(best_calls, type=3)


#Pre_chew non significant, drop from model for post-hoc tests
final_calls <- glmer(N_calls ~ Sex + Naive * Trial + (1|Bat_ID), 
                family = poisson, data=calls)



#Post-hoc tests

lssumm_calls <- summary(lsmeans(final_calls,
                           pairwise ~ Trial*Naive*Sex,
                           adjust="tukey"))

lssumm_calls 


lssumm_calls_1 <- summary(lsmeans(final_calls, pairwise ~ Naive*Sex, 
                            at = list(Trial = "1")))

lssumm_calls_2 <- summary(lsmeans(final_calls, pairwise ~ Naive*Sex, 
                                  at = list(Trial = "2")))

lssumm_calls_3 <- summary(lsmeans(final_calls, pairwise ~ Naive*Sex, 
                                  at = list(Trial = "3")))
                                          



##__________________________________________________________________________________________________________##

#Experiment 1

#Graph showing which stimulus prompted eating most often

a <- ggplot(data, aes(stimulus, fill = first_graph)) + 
  geom_bar(alpha=0.6, color = "black")+
  theme_classic()+
  labs(x="Stimulus", y = "Number of trials", fill = "Onset")+ 
  scale_fill_manual(values = c("#59B29D", "#E89E16"))


a +
  annotate("text", x=1, y=15.5, label= "B") + 
  annotate("text", x=2, y=15.5, label = "C")+
  annotate("text", x=3, y=15.5, label= "C") + 
  annotate("text", x=4, y=15.5, label = "AB")+
  annotate("text", x=5, y=15.5, label= "A")

ggsave("Figure 4_Number of trials in which a given stimulus prompted the onset of feeding.tiff", 
       width = 5, height = 3, dpi = 600)



#Model to test if stimulus affected when bats started eating
mod2 <- glm(first ~ stimulus, family = binomial, 
           data = data)

anova(mod2, test = "Chisq")

summary(mod2)


#To conduct a chi squared test we first converted data to % and then created data frame. Values 
#of 0 were increased to 1 (and 99) to avoid errors

perc = data.frame(stimulus = c("Calls", "Chewing", "Chewing + Calls", "Pink noise", "Silence"), 
                  yes = c(33, 1, 47, 1, 20), no = c(67, 99, 53, 99, 80))
perc.glm = glm(cbind(yes, no) ~ stimulus, family = binomial(), data = perc)
anova(perc.glm)
perc.lsm = lsmeans(perc.glm, "stimulus")
plot(perc.lsm)
pairs(perc.lsm)


##__________________________________________________________________________________________________________##

#Experiment 2, Model 1 to determine if stimulus affects time spent in the 75 cm area

require(MuMIn)

options(na.action = "na.fail")


#Create global model with 4-way interactions
globalmodel <- lmer(dist_75_cm ~ (stimulus + pre_calls + pre_chew + stim_order + trial)^4 + (1|id), data=data_2, REML=FALSE)


#Generate a table of models that include every possible combination of predictor variables
combinations <- dredge(globalmodel)


#Keep the best models based on delta values < 5
best_comb <- subset(combinations, delta < 5)

print(best_comb)


#Estimate of the relative importance of each predictor variable 
imp <- sw(best_comb)

imp


#Based on latter, the best model
best <- lmer(dist_75_cm ~ stimulus + pre_calls + stim_order + trial 
             + (1|id), data=data_2, REML=FALSE)


#Test effect of random variable
library(lmerTest)

ranova(best)


#Test effect of predictive variables
car::Anova(best, type=3)


#Pre_calls non significant, drop from model for post-hoc tests
best_1 <- lmer(dist_75_cm ~ stimulus + stim_order + trial + (1|id), data=data_2, REML=FALSE)
best_1_add <- lmer(dist_75_cm ~ stimulus + stim_order + trial + (1|id), data=data_2, REML=FALSE)
best_1_inter <- lmer(dist_75_cm ~ stimulus * stim_order * trial + (1|id), data=data_2, REML=FALSE)


#Post-hoc tests
leastsquare = lsmeans(best_1,
                      pairwise ~ stimulus*stim_order*trial,
                      adjust="tukey")

lssumm <- summary(lsmeans(best_1,
                          pairwise ~ stimulus*stim_order*trial,
                          adjust="tukey"))

#Create constrasts within trials and stim_order
ls_sub_1 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "1", stim_order = "1")))
ls_sub_2 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "1", stim_order = "2")))
ls_sub_3 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "1", stim_order = "3")))
ls_sub_4 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "1", stim_order = "4")))
ls_sub_5 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "1", stim_order = "5")))

ls_sub_6 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "2", stim_order = "1")))
ls_sub_7 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "2", stim_order = "2")))
ls_sub_8 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "2", stim_order = "3")))
ls_sub_9 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                            at = list(trial = "2", stim_order = "4")))
ls_sub_10 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                             at = list(trial = "2", stim_order = "5")))

ls_sub_11 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                             at = list(trial = "3", stim_order = "1")))
ls_sub_12 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                             at = list(trial = "3", stim_order = "2")))
ls_sub_13 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                             at = list(trial = "3", stim_order = "3")))
ls_sub_14 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                             at = list(trial = "3", stim_order = "4")))
ls_sub_15 <- summary(lsmeans(best_1, pairwise ~ stimulus, 
                             at = list(trial = "3", stim_order = "5")))


#Crear df de lista para contrasts
contrasts <- do.call(cbind.data.frame, lssumm$contrasts)

library(openxlsx)

# for writing a data.frame or list of data.frames to an xlsx file
write.xlsx(contrasts, 'contrasts.xlsx')



#Graph based on model results 

plot(leastsquare)

d <- summary(lsmeans(best_1_add, ~stimulus + stim_order + trial))
n <- summary(lsmeans(best_1_inter, ~stimulus + stim_order + trial))


#Graph for additive model
c <- ggplot(d, aes(x=stimulus, y=lsmean)) + 
  geom_line(aes(y = lsmean, group = stim_order)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, alpha = 0.5) +
  geom_point(aes(y = lsmean, fill = stim_order), size = 4, shape = 21, alpha = 0.5) +
  theme_classic()+
  labs(x="Stimulus", y = "Least squares mean", fill = "Stimulus order")+ 
  scale_fill_brewer(palette="Spectral") + facet_grid(vars(trial)) + ylim(-1, 5.5)

c 

#Graph for interaction model
c <- ggplot(n, aes(x=stimulus, y=lsmean)) + 
  geom_line(aes(y = lsmean, group = stim_order)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, alpha = 0.5) +
  geom_point(aes(y = lsmean, fill = stim_order), size = 4, shape = 21, alpha = 0.5) +
  theme_classic()+
  labs(x="Stimulus", y = "Least squares mean", fill = "Stimulus order")+ 
  scale_fill_brewer(palette="Spectral") + facet_grid(vars(trial)) + ylim(-1, 5.5)

c 



ggsave("Figure 9_model.tiff", width = 6, height = 6, dpi = 600)



#Graph based on real data

# Function to calculate the mean and the standard deviation for each group (in this case stim_order)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]]/sqrt(length(x)), na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#Summarize the data

df2 <- data_summary(data_2, varname="dist_75_cm", 
                    groupnames=c("stimulus", "stim_order", "trial"))

df2$stim_order=as.factor(df2$stim_order)

df2$trial=as.factor(df2$trial)

#Make graph
d <- ggplot(df2, aes(x=stimulus, y=dist_75_cm)) + 
  geom_line(aes(y = dist_75_cm, group = stim_order), alpha = 0.5) +
  geom_errorbar(aes(ymin=dist_75_cm-se, ymax=dist_75_cm+se), width = 0.2, alpha = 0.5) +
  geom_point(aes(y = dist_75_cm, fill = stim_order), size = 4, shape = 21, alpha = 0.5) +
  theme_classic()+
  labs(x="Stimulus", y = "Time spent in the 75 cm area (s)", fill = "Stimulus order")+ 
  scale_fill_brewer(palette="Spectral") + facet_grid(vars(trial))

d 

ggsave("Figure 9_real.tiff", width = 6, height = 6, dpi = 600)

##__________________________________________________________________________________________________________##

#Experiment 2, Model 2 to determine if stimulus affects number of times in the 50 cm area

require(MuMIn)

options(na.action = "na.fail")


#Create global model with 4-way interactions
globalmodel2 <- glmer(dist_50_cm ~ (stimulus + pre_calls + pre_chew + stim_order + trial)^4 + (1|id), family = poisson, data=data_2)


#Generate a table of models that include every possible combination of predictor variables
combinations2 <- dredge(globalmodel2)


#Keep the best models based on delta values < 5
best_comb2 <- subset(combinations2, delta < 5)

print(best_comb2)


#Estimate of the relative importance of each predictor variable 
imp <- sw(best_comb2)

imp


#Based on latter, the best model
best2 <- glmer(dist_50_cm ~ stimulus + trial + pre_calls + pre_calls * stimulus
               + pre_calls * trial + stim_order + stimulus * trial
               + pre_chew + (1|id), family = poisson, data=data_2)


#Test effect of random variable (does not work with glmer)
#library(lmerTest)

#ranova(best2)


#Test effect of predictive variables
car::Anova(best2, type=3)


#Pre_chew non significant, drop from model for post-hoc tests
best_2 <- glmer(dist_50_cm ~ stimulus + trial + pre_calls + pre_calls * stimulus
                + pre_calls * trial + stim_order + stimulus * trial
                + (1|id), family = poisson, data=data_2)


#Post-hoc tests

leastsquare2 = lsmeans(best_2,
                      pairwise ~ stimulus*trial*pre_calls+stim_order,
                      adjust="tukey")

lssumm2 <- summary(lsmeans(best_2,
                           pairwise ~ stimulus*trial+pre_calls+stim_order,
                           adjust="tukey"))


#Create constrasts within trials and stim_order
ls_sub_21 <- summary(emmeans(best_2, pairwise ~ "stimulus", "pre_calls",
                            at = list(trial = "1", stim_order = "1")))
ls_sub_22 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "1", stim_order = "2")))
ls_sub_23 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "1", stim_order = "3")))
ls_sub_24 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "1", stim_order = "4")))
ls_sub_25 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "1", stim_order = "5")))

ls_sub_26 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "2", stim_order = "1")))
ls_sub_27 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "2", stim_order = "2")))
ls_sub_28 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "2", stim_order = "3")))
ls_sub_29 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "2", stim_order = "4")))
ls_sub_30 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "2", stim_order = "5")))

ls_sub_31 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "3", stim_order = "1")))
ls_sub_32 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "3", stim_order = "2")))
ls_sub_33 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "3", stim_order = "3")))
ls_sub_34 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "3", stim_order = "4")))
ls_sub_35 <- summary(lsmeans(best_2, pairwise ~ "stimulus", "pre_calls", 
                            at = list(trial = "3", stim_order = "5")))




#Graph based on real data

# Function to calculate the mean and the standard deviation for each group (in this case stim_order)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]]/sqrt(length(x)), na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#Summarize the data

df3 <- data_summary(data_2, varname="dist_50_cm", 
                    groupnames=c("stimulus", "stim_order", "trial", "pre_calls"))


#Make graph
e <- ggplot(df3, aes(x=stimulus, y=dist_50_cm)) + 
  geom_line(aes(y = dist_50_cm, group = stim_order), alpha = 0.5) +
  geom_errorbar(aes(ymin=dist_50_cm-se, ymax=dist_50_cm+se), width = 0.2, alpha = 0.5) +
  geom_point(aes(y = dist_50_cm, fill = stim_order), size = 4, shape = 21, alpha = 0.5) +
  theme_classic()+
  labs(x="Stimulus", y = "Times bats visit the 50 cm area (s)", fill = "Stimulus order")+ 
  scale_fill_brewer(palette="Spectral") + facet_grid(vars(trial), vars(pre_calls))

e 

ggsave("Figure 10_real.tiff", width = 9, height = 6, dpi = 600)




##__________________________________________________________________________________________________________##

#OTHER GRAPHS (Supplementary)


#Graph showing how stimulus were distributed among trials

f <- ggplot(data_2, aes(stimulus, fill = stim_order)) + 
  geom_bar(alpha = 0.5, color = "black")+
  theme_classic()+
  labs(x="Stimulus", y = "Number of times a stimulus was used in a trial")+ 
  scale_fill_brewer(palette="Spectral")

f


#Graph showing how much time bats spent at 75 cm distance from speaker based on stimulus

g <- ggplot(data_2, aes(x=stimulus, y=dist_75_cm, fill=stimulus)) + 
  geom_boxplot(alpha = 0.5) +
  theme_classic()+
  theme(legend.position="none")+
  labs(x="Stimulus", y = "Time spent in the 75 cm area")+ 
  scale_fill_brewer(palette="Spectral")
g

g + geom_jitter(shape=16, position=position_jitter(0.2))




#Graph showing how much time bats spent at 75 cm distance from speaker based on stim_order

h <- ggplot(data_2, aes(x=stim_order, y=dist_75_cm, fill=stim_order)) + 
  geom_boxplot(alpha = 0.5) +
  theme_classic()+
  theme(legend.position="none")+
  labs(x="Stimulus order", y = "Time spent in the 75 cm area")+ 
  scale_fill_brewer(palette="Spectral")
h

h + geom_jitter(shape=16, position=position_jitter(0.2))




#Graph showing how much time bats spent at 75 cm distance from speaker based on stimulus and stim_order

i <- ggplot(data_2, aes(x=stimulus, y=dist_75_cm, fill=stim_order)) + 
  geom_boxplot(alpha = 0.5) +
  theme_classic()+
  labs(x="Stimulus", y = "Time spent in the 75 cm area")+ 
  scale_fill_brewer(palette="Spectral")
  
i  




#Graph showing number of times bats spent at 50 cm distance from speaker based on stimulus

j <- ggplot(data_2, aes(x=stimulus, y=dist_50_cm, fill=stimulus)) + 
  geom_boxplot(alpha = 0.5) +
  theme_classic()+
  theme(legend.position="none")+
  labs(x="Stimulus", y = "Times bat entered the 50 cm area")+ 
  scale_fill_brewer(palette="Spectral")
j

j + geom_jitter(shape=16, position=position_jitter(0.2))

ggsave("Figure 10.tiff", width = 6, height = 4, dpi = 600)


#Graph showing number of times bats spent at 50 cm distance from speaker based on stimulus order

k <- ggplot(data_2, aes(x=stim_order, y=dist_50_cm, fill=stim_order)) + 
  geom_boxplot(alpha = 0.5) +
  theme_classic()+
  theme(legend.position="none")+
  labs(x="Stimulus", y = "Times bat entered the 50 cm area")+ 
  scale_fill_brewer(palette="Spectral")
k

k + geom_jitter(shape=16, position=position_jitter(0.2))

ggsave("Figure 11.tiff", width = 6, height = 4, dpi = 600)






