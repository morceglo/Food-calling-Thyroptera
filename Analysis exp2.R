library(ggplot2)
library (lme4)
library(ptmixed)
library (MASS)
library(car)
library(fitdistrplus)
library(agricolae)
library(magrittr)
library(multcomp)
library(emmeans)
library(repeated)
library(GLMMadaptive)
library(rcompanion)

#Pacckages required to plot coefficients
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(jtools)
library(broom.mixed)
library(ggpubr)

setwd("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/1_Food calls Thyroptera/BEAS/Version 2/Analysis")

calls <- read.csv('Calls_produced.csv', header = TRUE, sep = ";", 
                  colClasses = c("Order"="character", "Bat_ID"="character", "Sex"="character", 
                                 "Age"="character", "Trial"="character", "Audio"="character", 
                                 "N_calls"="numeric", "Call"="numeric", "Audio_length"="numeric", "Naive"="character"))

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
#Calls produced using robust estimation of linear mixed-effects models
#From https://cran.r-project.org/web/packages/robustlmm/vignettes/rlmer.pdf


library(robustlmm)
library(modelsummary)

#Compare model that includes all predictors with a null model
models <- list(

"mixed_robust_model" <- rlmer(N_calls ~ Sex * Trial * Naive + (1|Bat_ID), data = calls),
"mixed_robust_model_null" <- rlmer(N_calls ~ (1|Bat_ID), data = calls)
)

modelsummary(models)

#Best model includes all factors
mixed_robust_model <- rlmer(N_calls ~ Sex * Trial * Naive + (1|Bat_ID), data = calls)

# Print summary of the model
summary(mixed_robust_model)

library(pandoc)

model1 <- modelsummary(mixed_robust_model, gof_omit = ".*",
             statistic = c("conf.int",
                           "s.e. = {std.error}", 
                           "t = {statistic}",
                           "p = {p.value}"), shape = term ~ statistic, output = "Table S2.xlsx")


#Generate contrasts for all predictors in the model (https://marginaleffects.com/vignettes/get_started.html) and save them in a table

library(marginaleffects)

comparisons(mixed_robust_model)

m <- avg_comparisons(mixed_robust_model, by = c("Sex", "Trial", "Naive"))

library("writexl")
write_xlsx(m, "Table 2.xlsx")



#Plot figure 3 that shows the predicted number of calls according to the combination of predictors

png("Figure 3_Predicted number of calls.tiff", res = 300, units = "cm", width = 12, height = 5)

d <- plot_predictions(
  mixed_robust_model, 
  by = c("Naive", "Sex", "Trial")) +
  ylab("Number of calls") + 
  coord_cartesian(ylim=c(-2.5,9)) +  
  set_theme(base = theme_sjplot(), 
            axis.title.size = 1.0,  
            axis.textsize.x = 0.8,  
            axis.textsize.y = 0.8,
            title.size = 0.8,
            title.align = "center",
            legend.size = 0.6,
            legend.title.size = 0.7,
            legend.title.face = "plain") +
  ggtitle("Trial number") +
  labs(x = "Naïve")

colors <- c("#59B29D", "#E89E16")
options(ggplot2.discrete.colour = colors)

dat_text <- data.frame(
  label = c("ns", "ns", "ns", "***", "ns", "ns"),
  Trial   = c(1, 2, 3),
  x     = c(1, 2, 1, 2, 1, 2),
  y     = c(8.75, 8.75, 8.75, 8.75, 8.75, 8.75)
)

d + geom_text(
  data    = dat_text,
  size = 2,
  mapping = aes(x = x, y = y, label = label)
)


dev.off()



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
mod2 <- glm(first ~ stimulus*trial, family = binomial, 
           data = data)

r <- anova(mod2, test = "Chisq")

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

#Create plots that show estimates of the model for time spent near the 75 cm region
plot_75 <- plot_model(
  best, 
  terms = c("stimulusChewing", "stimulusChewing + Calls", "stimulusPink noise", "stimulusSilence"),
  axis.labels = c("Silence", "Pink noise", "Chewing + Calls", "Chewing"),
  show.values = TRUE, value.offset = 0.3,
  title = "75 cm",
  colors = "bw",
  vline.color = "black", 
  width = 0.5, grid =TRUE, show.p = TRUE
)

plot_75 + set_theme(base = theme_sjplot(), 
              axis.title.size = 1.2,  
              axis.textsize.x = 1.0,  
              axis.textsize.y = 1.0,
              title.size = 1.2,
              title.align = "left")


#Supplementary graph (coefficients for other factors in the model)

supp_75cm <- plot_model(
  best, 
  rm.terms = c("stimulusChewing", "stimulusChewing + Calls", "stimulusPink noise", "stimulusSilence"),
  axis.labels = c("Trial 3", "Trial 2",  "Stimulus order 5", "Stimulus order 4", "Stimulus order 3", 
                  "Stimulus order 2", "Pre-calls (yes)"),
  show.values = TRUE, value.offset = 0.3,
  title = "75 cm",
  colors = "bw",
  vline.color = "black", 
  width = 0.5, grid =TRUE, show.p = TRUE
)

supp_75cm + set_theme(base = theme_sjplot(), 
              axis.title.size = 1.2,  
              axis.textsize.x = 1.0,  
              axis.textsize.y = 1.0,
              title.size = 1.2,
              title.align = "left") 



#Test effect of random variable
library(lmerTest)
ranova(best)


#Test effect of predictive variables
car::Anova(best, type=3)

lsm <- ls_means(best)
summary(best, ddf="lme4")


##__________________________________________________________________________________________________________##

#Experiment 2, Model 2 to determine if stimulus affects number of times in the 50 cm area

#Model that predicts time spent in the 50 cm area
best2 <- glmer(dist_50_cm ~ stimulus + trial + stim_order + pre_calls + pre_chew +
                 (1|id), family = poisson, data=data_2, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

best2_2 <- glmer(dist_50_cm ~ stimulus + trial + stim_order + pre_calls + 
                   (1|id), family = poisson, data=data_2, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

fullnull2 <- glmer(dist_50_cm ~ 1 + 
                 (1|id), family = poisson, data=data_2, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

null2 <- glmer(dist_50_cm ~ trial + stim_order + pre_calls + pre_chew +
                       (1|id), family = poisson, data=data_2, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))


anova (best2, fullnull2)
anova (best2, null2)
anova (best2_2, fullnull2)
anova (best2_2, null2)
anova (null2, fullnull2)
anova (fullnull2, null2)


summary(best2)


#Create plots that show estimates of the model for times visiting the 50 cm region
plot_50 <- plot_model(
  best2, 
  terms = c("stimulusChewing", "stimulusChewing + Calls", "stimulusPink noise", "stimulusSilence"),
  axis.labels = c("Silence", "Pink noise", "Chewing + Calls", "Chewing"),
  show.values = TRUE, value.offset = 0.3,
  title = "50 cm",
  colors = "bw",
  vline.color = "black", 
  width = 0.5, grid =TRUE, show.p = TRUE)

plot_50 + set_theme(base = theme_sjplot(), 
              axis.title.size = 1.2,  
              axis.textsize.x = 1.0,  
              axis.textsize.y = 1.0,
              title.size = 1.2,
              title.align = "left")

library(cowplot)
plot_grid(
  plot_grid(plot_75, nrow = 1) +
    theme(plot.background = element_rect(color = "black")),
  plot_grid(plot_50, nrow = 1) +
    theme(plot.background = element_rect(color = "black")), 
  nrow = 2)

ggsave("Figure 5_Effect of stimulus type on social recruitment.tiff", 
       width = 5, height = 5, dpi = 600)


#Supplementary graph (coefficients for other factors in the model)
supp_50cm <- plot_model(
  best2, 
  rm.terms = c("stimulusChewing", "stimulusChewing + Calls", "stimulusPink noise", "stimulusSilence"),
  axis.labels = c("Trial 3", "Trial 2",  "Stimulus order 5", "Stimulus order 4", "Stimulus order 3", 
                  "Stimulus order 2", "Pre-calls (yes)"),
  show.values = TRUE, value.offset = 0.3,
  title = "50 cm",
  colors = "bw",
  vline.color = "black", 
  width = 0.5, grid =TRUE, show.p = TRUE
)

supp_50cm + set_theme(base = theme_sjplot(), 
              axis.title.size = 1.2,  
              axis.textsize.x = 1.0,  
              axis.textsize.y = 1.0,
              title.size = 1.2,
              title.align = "left") 


plot_grid(
  plot_grid(supp_75cm, nrow = 1) +
    theme(plot.background = element_rect(color = "black")),
  plot_grid(supp_50cm, nrow = 1) +
    theme(plot.background = element_rect(color = "black")), 
  nrow = 2)

ggsave("Figure 2S_Effect of other variables on social recruitment.tiff", 
       width = 5, height = 7, dpi = 600)



#Test effect of predictive variables
car::Anova(best2, type=3)



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






