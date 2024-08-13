library(car)
library(lme4)
library(dplyr)
library(multcomp)
library(multcompView)
library(emmeans)
library(tidyverse)
library(nortest)
library(MuMIn)
library(ggplot2)
library(nlme)
library(DHARMa) ##for testing overdispersion
library(glmmTMB)
library(effects)
library(moments)
library(modelsummary)
library(piecewiseSEM)
library(performance)
library(sjPlot)



###########model to test thte effects of tmp, aridity and soil nutrients on percentage cover
############### 
cover_nutnet_mean  = read.csv("../data/nutnet.climate_cover/nutnet.cover_climate.csv")
names(cover_nutnet_mean)

#calculate mean of maximum.cover
#calculate_mean <- function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else first(x)
#cover_nutnet_mean <- cover_nutnet %>%
 # group_by(site_code, year, block, trt, plot) %>%
 # summarise_all(calculate_mean)
#nrow(cover_nutnet_mean) #27906
#head(cover_nutnet_mean)

#write.csv(cover_nutnet_mean, "cover_nutnet_mean.csv")


#create column for each trt
#cover_nutnet_mean 
cover_nutnet_mean$Ntrt = 0 ## creation of a column with Ntrt= 0 
cover_nutnet_mean$Ntrt[cover_nutnet_mean$trt == 'N' | cover_nutnet_mean$trt == 
                         'NP' | cover_nutnet_mean$trt == 'NK' | 
                         cover_nutnet_mean$trt == 'NPK' | cover_nutnet_mean$trt == 
                         'NPK+Fence'] = 1 ## set Ntrt=1 if any point receives N
### repeat the same for the other treatments 
cover_nutnet_mean$Ptrt = 0
cover_nutnet_mean$Ptrt[cover_nutnet_mean$trt == 'P' | cover_nutnet_mean$trt == 
                         'NP' | cover_nutnet_mean$trt == 'PK' | 
                         cover_nutnet_mean$trt == 'NPK' | cover_nutnet_mean$trt == 
                         'NPK+Fence'] = 1
cover_nutnet_mean$Ktrt = 0
cover_nutnet_mean$Ktrt[cover_nutnet_mean$trt == 'K' | cover_nutnet_mean$trt == 
                         'NK' | cover_nutnet_mean$trt == 'PK' |
                         cover_nutnet_mean$trt == 'NPK' | cover_nutnet_mean$trt == 
                         'NPK+Fence'] = 1


#converting trt factor
cover_nutnet_mean$Ntrt_fac <- as.factor(cover_nutnet_mean$Ntrt)
cover_nutnet_mean$Ptrt_fac <- as.factor(cover_nutnet_mean$Ptrt)
cover_nutnet_mean$Ktrt_fac <- as.factor(cover_nutnet_mean$Ktrt)
cover_nutnet_mean$block_fac <- as.factor(cover_nutnet_mean$block)
cover_nutnet_mean$trt_fac <- as.factor(cover_nutnet_mean$trt)

cover_nutnet_mean$plot <- as.factor(cover_nutnet_mean$plot)


#plotting aridity against max.cover
plot(cover_nutnet_mean$awi_pm_sr_yr, cover_nutnet_mean$c4_percent)
plot(cover_nutnet_mean$tmp, cover_nutnet_mean$c4_percent)
plot(cover_nutnet_mean$awi_pm_sr_yr, cover_nutnet_mean$tmp)
hist(cover_nutnet_mean$c4_percent)

#table(cover_nutnet_mean$ps_path) # C3- 14571, C4- 4063

#cover_subset = subset(cover_nutnet_mean, ps_path== "C3"| (ps_path== "C4"))
#nrow(cover_subset)

### add lifeform and include in the model??

str(cover_nutnet_mean)

#####################################
###considered zero infaltion####
#if variance is greater than the mean
####################################
summary(cover_nutnet_mean$c4_percent)
var(cover_nutnet_mean$c4_percent)

###converting percentage to porportion and saving as c4_proportion
#transform your vegetation cover data in a proportional cover data, 
#limiting the range of variation between 0 and 1. Then, you can model 
#your response variable using a beta distribution, which accounts for proportional data.
cover_nutnet_mean$c4_proportion <- (cover_nutnet_mean$c4_percent)/100



#####categorize aridity index by threshold
cover_nutnet_mean$aridity.cateogry <- with(cover_nutnet_mean,
                                           ifelse(awi_pm_sr_yr < 0.03, "HyperArid",
                                                  ifelse(awi_pm_sr_yr >= 0.03 & awi_pm_sr_yr < 0.20, "Arid",
                                                         ifelse(awi_pm_sr_yr >= 0.20 & awi_pm_sr_yr < 0.50, "SemiArid",
                                                                ifelse(awi_pm_sr_yr >= 0.50 & awi_pm_sr_yr < 0.65, "SubHumid",
                                                                       "Humid")))))
###combine arid and hyperarid
cover_nutnet_mean$aridity.cateogry <- with(cover_nutnet_mean,
                                                  ifelse(awi_pm_sr_yr <= 0.03 & awi_pm_sr_yr < 0.20, "Arid",
                                                         ifelse(awi_pm_sr_yr >= 0.20 & awi_pm_sr_yr < 0.50, "SemiArid",
                                                                ifelse(awi_pm_sr_yr >= 0.50 & awi_pm_sr_yr < 0.65, "SubHumid",
                                                                       "Humid"))))

cover_nutnet_mean$aridity.cateogry <- as.factor(cover_nutnet_mean$aridity.cateogry)
table(cover_nutnet_mean$aridity.cateogry)
cover_nutnet_mean[cover_nutnet_mean$aridity.cateogry == "SemiArid", "site_code"]


###############################################################################
#####c4.cover.model 
#######################################################################
#C4.percent.model_lmer <- lmer(sqrt(c4_percent) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp +
  #                               awi_pm_sr_yr  + (1|site_code) + (1|site_code:block_fac) + 
  #                               (1|site_code:block_fac:plot),data = cover_nutnet_mean)
                              

### using beta-distribution because the beta distribution inherently accommodates
#skewness and heteroskedasticity while also been a good distribution to use for proportion dataset
#With this function, the dependent variable varies between 0 and 1, 
#but no observation can equal exactly zero or exactly one. 
#The model assumes that the data follow a beta distribution.
#############################orderbeta(might not be the right family to use)
##orderbeta regression considers continuous distribution with lower and upper bounds
###Ordered beta regression from Kubinec (2022); fits continuous (e.g. proportion) data in the closed interval [0,1].

######trying out different glmmTMB
###ordbeta (the C4_proportion used here is with 0 and 1)
C4.percent.model_orderbeta <- glmmTMB(c4_proportion ~ tmp + aridity.cateogry + Ntrt_fac * Ptrt_fac * Ktrt_fac  +
                                       (1|site_code) +  
                                        (1|site_code:block_fac:plot),data = cover_nutnet_mean, 
                                      ziformula  =   ~ tmp + aridity.cateogry + Ntrt_fac * Ptrt_fac * Ktrt_fac  +
                                        (1|site_code) +
                                        (1|site_code:block_fac:plot),
                                      family = ordbeta(), na.action = na.omit)

###plotting model prediction against actual data
cover_nutnet_mean$predicted <- predict(C4.percent.model_orderbeta)
ggplot(cover_nutnet_mean, aes(c4_proportion, predicted)) +
  geom_point() + theme_minimal()

####model results
Anova(C4.percent.model_orderbeta)
summary(C4.percent.model_orderbeta)
fixef(C4.percent.model_orderbeta)



#plotting model for model coefficients
tab_model(C4.percent.model_orderbeta, show.ci = TRUE, show.se = TRUE, title = "model coefficients")
plot_model(C4.percent.model_orderbeta, grid = FALSE)
plot_model(C4.percent.model_orderbeta,  type = "pred", terms = "tmp[all]")


#######finding out position of the 1's in the dataset and changing them to 0.999
## to be able to use the proportion data in the model cos glmmTMB ....
cover_nutnet_mean$c4_proportion2 <- cover_nutnet_mean$c4_proportion
position <- which(cover_nutnet_mean$c4_proportion2==1) 
###insert changing from 1 to 0.999
cover_nutnet_mean$c4_proportion2[position] <- 0.999


#####beta_family
C4.percent.model_betafamily <-glmmTMB(c4_proportion2 ~ tmp + aridity.cateogry + Ntrt_fac * Ptrt_fac * Ktrt_fac  + 
                                        (1|site_code) + 
                                        (1|site_code:block_fac:plot),data = cover_nutnet_mean, 
                                      ziformula  =  ~ tmp + aridity.cateogry + Ntrt_fac * Ptrt_fac * Ktrt_fac  +
                                        (1|site_code) +
                                        (1|site_code:block_fac:plot),
                                      family = beta_family(),
                                      control = glmmTMBControl(optimizer = optim, profile = TRUE, 
                                                               optArgs = list(method = "BFGS")))



Anova(C4.percent.model_betafamily)
summary(C4.percent.model_betafamily)


cover_nutnet_mean$predicted <- predict(C4.percent.model_betafamily)
ggplot(cover_nutnet_mean, aes(c4_proportion, predicted)) +
  geom_point() + theme_minimal()
r2_nakagawa(C4.percent.model_betafamily)

plot_model(C4.percent.model_betafamily, grid = FALSE)
plot_model(C4.percent.model_betafamily,  type = "pred", terms = "tmp[all]", ci.lvl = 0.95)


##over-dispersion and zero inflation test
###################################
testDispersion(C4.percent.model_betafamily)
testZeroInflation(C4.percent.model_betafamily)

########model summary
###AIC - -3664.3, RMSE - 0.12
modelsummary(C4.percent.model_betafamily, statistic = "conf.int", conf_level = 0.95, metrics = "all")
plot_model(C4.percent.model_betafamily, transform = NULL, type = "est", grid = TRUE)
check_collinearity(C4.percent.model_betaordbeta2, component = "all")
tab_model(C4.percent.model_betafamily)

#####betaordebeta2( used for  ESA talk)
C4.percent.model_betaordbeta2 <- glmmTMB(c4_proportion2 ~ tmp + aridity.level + Ntrt_fac * Ptrt_fac * Ktrt_fac  + 
                                            (1|site_code) + 
                                            (1|site_code:block_fac:plot),data = cover_nutnet_mean, 
                                          ziformula  =  ~ tmp + aridity.cateogry + Ntrt_fac * Ptrt_fac * Ktrt_fac  +
                                         (1|site_code) + (1|site_code:block_fac) + 
                                           (1|site_code:block_fac:plot),
                                          family = beta_family())
                          
                                                                  


Anova(C4.percent.model_betaordbeta2)
summary(C4.percent.model_betaordbeta2)


cover_nutnet_mean$predicted <- predict(C4.percent.model_betaordbeta2)
ggplot(cover_nutnet_mean, aes(c4_proportion2, predicted)) +
  geom_point() + theme_minimal()

plot_model(C4.percent.model_betaordbeta2, grid = FALSE)
plot_model(C4.percent.model_betaordbeta2,  type = "pred", terms = "tmp[all]")


# Check model assumptions
#####################################
#r.squaredGLMM(C4.percent.model_betaordbeta2)
r2(C4.percent.model_betaordbeta2)
r2_nakagawa(C4.percent.model_betaordbeta2)

##over-dispersion and zero inflation test
###################################
testDispersion(C4.percent.model_betaordbeta2)
testZeroInflation(C4.percent.model_betaordbeta2)

########model summary
###AIC -3664.3, RMSE - 0.12
modelsummary(C4.percent.model_betaordbeta2, statistic = "conf.int", metrics = "all")
plot_model(C4.percent.model_betaordbeta2, transform = NULL, type = "est", grid = TRUE)
check_collinearity(C4.percent.model_betaordbeta2, component = "all")
tab_model(C4.percent.model_betaordbeta2)



C4.percent.model_betaordbeta2$sdr$pdHess


######################################post-hoc test
nb1 <- simulateResiduals(C4.percent.model_betaordbeta2)
plot(nb1)



N <- as.data.frame(emmeans(C4.percent.model_betaordbeta2, ~Ntrt_fac, type = 'response'))
emmeans(C4.percent.model_betaordbeta2, pairwise~Ntrt_fac, type = 'response')
cld(emmeans(C4.percent.model_orderbeta, pairwise~Ntrt_fac:Ktrt_fac))
cld(emmeans(C4.percent.model_betaordbeta2, pairwise~Ntrt_fac, type = "response"))
cld(emmeans(C4.percent.model_betaordbeta2, ~Ntrt_fac, type = 'response'))



cover_nutnet_mean$aridity <- as.numeric(as.factor(cover_nutnet_mean$aridity.cateogry))
aridity <- as.data.frame(emmeans(C4.percent.model_betaordbeta2, ~aridity.level, type = 'response'))


##################################################
######figs for ESA################################
figtheme <- theme_linedraw(base_size = 20) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(size = 1.5, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        #legend.title = element_text(face = "bold"),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

##PLOTS
cover_nutnet_mean$trt <- factor(cover_nutnet_mean$trt, levels = c("Control", "N", 
                                                                          "P", "K", "NP", "NK",
                                                                          "PK", "NPK"))

cover_nutnet_mean$aridity.level <- factor(cover_nutnet_mean$aridity.cateogry, 
                                             levels = c("SemiArid", "SubHumid", "Humid"))
                                                                                   
                                                                  
                                                               



####### all cover plot by trt and by growth form
cover_trt_fig <- ggplot(data = na.omit(cover_nutnet_mean), 
                      aes(trt, c4_proportion, fill = trt)) + 
  geom_violin() +
  scale_fill_manual(values = c('#2a9d8f', '#FF9200', '#abd17dff', "#d60404", "#2c46a5", "#d20e0f","#a63253", "purple")) +
  labs(x = "Treatment") +   
  figtheme +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(color = "black"))


###nitrogen availiability ##############
Ntrt.plot <- ggplot() +
  geom_violin(data = cover_nutnet_mean, aes(x = Ntrt_fac, y = c4_proportion2), 
              fill = "#FF9200", color = "black", alpha = 0.7) +
  geom_errorbar(data = N, 
                aes(x = Ntrt_fac, y = response, ymin=asymp.LCL, ymax = asymp.UCL),
                width = 0.06, linewidth = 1.3, color = "black", size = 0.5, alpha = 1) +
  geom_line(data = N, aes(Ntrt_fac, y = response, group = 1), 
            color = "black", linewidth = 1, alpha = 1.5, linetype = "solid") +
  geom_point(data = N, position = position_dodge(width= 0.75),
             aes(x = Ntrt_fac, y = response),
             color = "black", size = 3) + 
  labs(x = expression(bolditalic("Nitrogen")),
       y = expression(bolditalic("Relative C"["4"]*" cover")))  +
  figtheme +
  scale_x_discrete(labels = c("ambient", "addition")) +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))


Ntrt.plot


####################################
# Generate predictions
emm <- emmeans(C4.percent.model_betaordbeta2, specs = ~ tmp, 
               at = list(tmp = seq(3, 30, 5)), 
               type = "response")

# Convert to data frame
pred_df <- as.data.frame(emm)


temperature_values <- seq(10, 35, by = 5)

####################temperature regline plot
tmpregline <- ggplot() +
  geom_jitter(data = cover_nutnet_mean, aes(x = tmp, y = c4_proportion2), 
             alpha = 0.3, color = "black") +
  # Add predicted regression line
  geom_line(data = pred_df, aes(x = tmp, y = response), 
            color = "red", size = 1) +
  # Add confidence interval ribbon
  geom_ribbon(data = pred_df, aes(x = tmp, ymin = asymp.LCL, ymax = asymp.UCL), 
              alpha = 0.2, fill = "red") +
  labs(x = expression(bolditalic("Temperature")),
       y = expression(bolditalic("Relative C"["4"]*" cover"))) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 30, 35)) +
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))


tmpregline

#########################Aridity plot
aridity.plot <-  ggplot(data = na.omit(cover_nutnet_mean), 
                       aes(aridity.level, c4_proportion, fill = aridity.level)) + 
  geom_violin() +
  scale_fill_manual(values = c( '#FF9200',  "#d60404","#2c46a5")) +
  geom_errorbar(data = aridity , 
                aes(x = aridity.level, y = response, ymin=asymp.LCL, ymax = asymp.UCL),
                width = 0.06, linewidth = 1.3, color = "black", size = 0.5, alpha = 1) +
  geom_line(data = aridity, aes(aridity.level, y = response, group = 1), 
            color = "black", linewidth = 1, alpha = 1.5, linetype = "solid") +
  geom_point(data = aridity, position = position_dodge(width= 0.75),
             aes(x = aridity.level, y = response),
             color = "black", size = 3) + 
  labs(y = expression(bolditalic("Relative C"["4"]*" cover")))  +
  figtheme +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(color = "black"))

aridity.plot













###################phosphorous availiability 
Ptrt <- ggplot() +
  geom_violin(data = cover_nutnet_mean, aes(x = Ntrt_fac, y = c4_proportion), 
              fill = "#2c46a5", color = "black", alpha = 0.7) +
  geom_errorbar(data = N, 
                aes(x = Ntrt_fac, y = response, ymin=asymp.LCL, ymax = asymp.UCL),
                width = 0.06, linewidth = 1.3, color = "black", size = 0.5, alpha = 0.4) +
  geom_point(data = N, position = position_dodge(width= 0.75),
             aes(x = Ntrt_fac, y = response),
             color = "black", size = 3) +
  labs(x = expression(bolditalic("Nitrogen Availability")),
       y = expression(bolditalic("C4 cover proportion"))) +
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))
  
  
  
  


ggplot(N) +
  geom_violin(aes(x=Ntrt_fac, y=response), fill="blue", alpha=0.5, width = 0.2) +
  geom_errorbar( aes(x=Ntrt_fac, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1, colour="black", alpha=0.9, size=0.5) +
  labs(x = expression(bold("Phosphorous")),
       y = expression(bold(italic("max.cover"))),
       color = expression(bold("ps_path"))) +
  theme_minimal() +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(color = "black"))
  

ggplot(N, aes(x=aridity.cateogry, y=response)) +
  geom_point(position = position_dodge(width= 0.75), size=5) +
  geom_errorbar(position = position_dodge(width= 0.75),aes(ymin=asymp.LCL, 
                                                           ymax=asymp.UCL), width=.2, linewidth = 1.5) 






  




ggplot(nitrogen.phosphorous, aes(x=interaction(Ntrt_fac, Ktrt_fac), y=emmean, 
                                 color=ps_path)) +
  geom_point(position = position_dodge(width= 0.75), size=5) +
  geom_errorbar(position = position_dodge(width= 0.75),aes(ymin=asymp.LCL, 
                                ymax=asymp.UCL), width=.2, linewidth = 1.5) +
  scale_color_manual(values=c("C3"="blue", "C4"="red")) +
  labs(x = expression(bold("Nitrogen:potassium")),
       y = expression(bold(italic("max.cover"))),
       color = expression(bold("ps_path"))) +
  theme_minimal() +
  theme(text = element_text(size=12))


ggplot(temperature, aes(x=tmp, y=emmean, color=ps_path)) +
  geom_point(position = position_dodge(width= 0.75), size=5) +
  geom_errorbar(position = position_dodge(width= 0.75),aes(ymin=asymp.LCL,
                                ymax=asymp.UCL), width=.2, linewidth = 1.5) +
  scale_color_manual(values=c("C3"="blue", "C4"="red")) +
  labs(x = expression(bold("Temp")),
       y = expression(bold(italic("max.cover"))),
       color = expression(bold("ps_path"))) +
  theme_minimal() +
  theme(text = element_text(size=12))



###############################running SEM (figure this out after ESA)
SEModel.C4 <- psem(
  
  lme(tmp~awi_pm_sr_yr, random = ~1|site_code,
      data = cover_nutnet_mean, na.action = na.omit),
  
  glmmTMB(c4_proportion2 ~ tmp + awi_pm_sr_yr + Ntrt_fac * Ptrt_fac * Ktrt_fac +
            (1|site_code) + (1|site_code:block_fac) + 
            (1|site_code:block_fac:plot),data = cover_nutnet_mean, 
          ziformula  = ~1, family = beta_family(), na.action = na.omit)
  
)

summary(SEModel.C4)

plot(SEModel.C4)

summary(SEModel.C4, direction = c("awi_pm_sr_yr<-tmp"))

fisherC(SEModel.C4)







###############Code works, keep it for the publication subsititution result###
####Using Random forest to train and valid####################
# Install and load necessary packages
#install.packages(c("randomForestSRC", "caret"))
library(randomForestSRC) # For generating random forest model
library(caret)  #classification and regression training and to make prediction.

# Read in the vegetation dataset
# Replace 'your_file.csv' with your actual file name
data <- read.csv("your_file.csv")

# Check the structure of your data
str(cover_nutnet_mean)

# Prepare the data
# Ensure your response and predictor variables are in the correct format
# Ensure vegetation cover is between 0 and 100
cover_nutnet_mean$c4_percent <- pmin(pmax(cover_nutnet_mean$c4_percent, 0), 100)

#data$temperature <- as.numeric(data$temperature)
#data$aridity_index <- as.numeric(data$aridity_index)
#data$nutrient_availability <- as.factor(data$nutrient_availability)

# Convert random effect variables to factors if they're not already
cover_nutnet_mean$sitecode <- as.factor(cover_nutnet_mean$site_code)
#data$block <- as.factor(data$block)
#data$plot <- as.factor(data$plot)

# Split the data into training and testing sets
set.seed(123)  # for reproducibility
train_index <- createDataPartition(cover_nutnet_mean$c4_percent, p = 0.7, list = FALSE)
train_data <- cover_nutnet_mean[train_index, ]
test_data <- cover_nutnet_mean[-train_index, ]

# Train the Random Forest model with random effects
rf_model <- rfsrc(c4_percent ~ tmp + awi_pm_sr_yr + Ntrt_fac * Ptrt_fac * Ktrt_fac +
                    (1|sitecode) + (1|sitecode:block_fac) + 
                    (1|sitecode:block_fac:plot),
                  data = train_data,
                  ntree = 500,
                  importance = "permute")

# Print model summary
print(rf_model)

# Variable importance
print(vimp(rf_model, joint = TRUE))

# Plot variable importance
plot(vimp(rf_model))

# Make predictions on the test set
predictions <- predict(rf_model, newdata = test_data)$predicted

# Ensure predictions are between 0 and 100
predictions <- pmin(pmax(predictions, 0), 100)

# Calculate RMSE
rmse <- sqrt(mean((test_data$c4_percent - predictions)^2))
print(paste("RMSE:", rmse))

# Calculate R-squared
r_squared <- 1 - sum((test_data$c4_percent - predictions)^2) / 
  sum((test_data$c4_percent - mean(test_data$c4_percent))^2)
print(paste("R-squared:", r_squared))

# Calculate Mean Absolute Error (MAE)
mae <- mean(abs(test_data$c4_percent - predictions))
print(paste("MAE:", mae))

# Plot predicted vs observed
ggplot(data.frame(predicted = predictions, observed = test_data$c4_percent), 
       aes(x = observed, y = predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Observed Vegetation Cover (%)", y = "Predicted Vegetation Cover (%)") +
  ggtitle("Predicted vs Observed Vegetation Cover") +
  xlim(0, 100) + ylim(0, 100)

# Partial dependence plots
plot.variable(rf_model, xvar.names = c("temperature", "aridity_index", "nutrient_availability"))

# Calculate percent error
percent_error <- 100 * abs(test_data$c4_percent - predictions) / test_data$c4_percent
mean_percent_error <- mean(percent_error, na.rm = TRUE)
print(paste("Mean Percent Error:", mean_percent_error, "%"))

# Plot histogram of percent errors
ggplot(data.frame(percent_error = percent_error), aes(x = percent_error)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.7) +
  labs(x = "Percent Error", y = "Count") +
  ggtitle("Distribution of Percent Errors")



###training testing and validating
# Load necessary libraries
library(ranger)
library(caret)
library(ggplot2)
library(dplyr)
library(randomForest)

###Using proportion data


str(cover_nutnet_mean)



# Ensure vegetation cover is between 0 and 1 (proportion)
cover_nutnet_mean$c4_proportion <- pmin(pmax(cover_nutnet_mean$c4_proportion, 0), 1)


# Create interaction terms
data$site_block <- interaction(data$sitecode, data$block)
data$site_plot <- interaction(data$sitecode, data$plot)
data$block_plot <- interaction(data$block, data$plot)
data$site_block_plot <- interaction(data$sitecode, data$block, data$plot)

# Create interaction terms
cover_nutnet_mean$np_interaction <- interaction(cover_nutnet_mean$Ntrt_fac, cover_nutnet_mean$Ptrt_fac)
cover_nutnet_mean$nk_interaction <- interaction(cover_nutnet_mean$Ntrt_fac, cover_nutnet_mean$Ktrt_fac)
cover_nutnet_mean$pk_interaction <- interaction(cover_nutnet_mean$Ptrt_fac, cover_nutnet_mean$Ktrt_fac)
cover_nutnet_mean$npk_interaction <- interaction(cover_nutnet_mean$Ntrt_fac, cover_nutnet_mean$Ptrt_fac, cover_nutnet_mean$Ktrt_fac)
cover_nutnet_mean$plot_order <- as.ordered(cover_nutnet_mean$plot)


set.seed(123)
train_index <- createDataPartition(cover_nutnet_mean$c4_proportion, p = 0.6, list = FALSE)
temp_data <- cover_nutnet_mean[-train_index,]
valid_index <- createDataPartition(temp_data$c4_proportion, p = 0.5, list = FALSE)

train_data <- cover_nutnet_mean[train_index,]
valid_data <- temp_data[valid_index,]
test_data <- temp_data[-valid_index,]

# Define all variables to be used in the model
all_vars <- c("tmp", "awi_pm_sr_yr", "Ntrt_fac", "Ptrt_fac", "Ktrt_fac", 
              "np_interaction", "nk_interaction", "pk_interaction", "npk_interaction",
              "sitecode", "block_fac", "plot")

# Create the formula
formula <- as.formula(paste("c4_proportion ~", paste(all_vars, collapse = " + ")))

# Train the random forest model
rf_model <- ranger(
  formula = formula,
  data = train_data,
  num.trees = 500,
  importance = "permutation",
  keep.inbag = TRUE
)

# Print model summary
print(rf_model)

# Make predictions on the validation set
valid_predictions <- predict(rf_model, data = valid_data)$predictions

# Calculate RMSE for validation set
valid_rmse <- sqrt(mean((valid_data$c4_proportion - valid_predictions)^2))
print(paste("Validation RMSE:", valid_rmse))

# Calculate R-squared for validation set
valid_r_squared <- 1 - sum((valid_data$c4_proportion  - valid_predictions)^2) / 
  sum((valid_data$c4_proportion  - mean(valid_data$c4_proportion ))^2)
print(paste("Validation R-squared:", valid_r_squared))

# Make predictions on the test set
test_predictions <- predict(rf_model, data = test_data)$predictions

# Calculate RMSE for test set
test_rmse <- sqrt(mean((test_data$c4_proportion- test_predictions)^2))
print(paste("Test RMSE:", test_rmse))

# Calculate R-squared for test set
test_r_squared <- 1 - sum((test_data$c4_proportion - test_predictions)^2) / 
  sum((test_data$c4_proportion- mean(test_data$c4_proportion))^2)
print(paste("Test R-squared:", test_r_squared))

# Plot predicted vs observed for validation set
ggplot(data.frame(predicted = valid_predictions, observed = valid_data$c4_proportion), 
       aes(x = observed, y = predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Observed Vegetation Cover (proportion)", y = "Predicted Vegetation Cover (proportion)") +
  ggtitle("Predicted vs Observed Vegetation Cover (Validation Set)") +
  xlim(0, 1) + ylim(0, 1)

# Plot predicted vs observed for test set
ggplot(data.frame(predicted = test_predictions, observed = test_data$c4_proportion), 
       aes(x = observed, y = predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Observed Vegetation Cover (proportion)", y = "Predicted Vegetation Cover (proportion)") +
  ggtitle("Predicted vs Observed Vegetation Cover (Test Set)") +
  xlim(0, 1) + ylim(0, 1)

# Variable importance
var_importance <- importance(rf_model)
var_importance_df <- data.frame(
  variable = names(var_importance),
  importance = var_importance
)
var_importance_df <- var_importance_df[order(-var_importance_df$importance),]
print(var_importance_df)

# Plot variable importance
ggplot(var_importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Variables", y = "Importance") +
  theme_minimal() +
  ggtitle("Variable Importance")

# Create confusion matrix for validation set
bin_vegetation <- function(x) {
  cut(x, breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("Very Low", "Low", "Medium", "High"))
}

valid_actual_bins <- bin_vegetation(valid_data$c4_proportion)
valid_predicted_bins <- bin_vegetation(valid_predictions)

valid_conf_matrix <- confusionMatrix(factor(valid_predicted_bins), factor(valid_actual_bins))
print("Validation Set Confusion Matrix:")
print(valid_conf_matrix)

# Create confusion matrix for test set
test_actual_bins <- bin_vegetation(test_data$c4_proportion)
test_predicted_bins <- bin_vegetation(test_predictions)

test_conf_matrix <- confusionMatrix(factor(test_predicted_bins), factor(test_actual_bins))
print("Test Set Confusion Matrix:")
print(test_conf_matrix)



#Using randomForestSRC with proportion data

# Read in the vegetation dataset
# Replace 'your_file.csv' with your actual file name
data <- read.csv("your_file.csv")

# Print column names to see what variables are available
print("Available columns in the dataset:")
print(colnames(cover_nutnet_mean))

# Convert variables to appropriate types
#data$c4_percent <- as.numeric(data$c4_percent)
#data$tmp <- as.numeric(data$tmp)
#data$awi_pm_sr_yr <- as.numeric(data$awi_pm_sr_yr)
#data$Ntrt_fac <- as.factor(data$Ntrt_fac)
#data$Ptrt_fac <- as.factor(data$Ptrt_fac)
#data$Ktrt_fac <- as.factor(data$Ktrt_fac)
##data$sitecode <- as.factor(data$sitecode)
#data$block_fac <- as.factor(data$block_fac)
cover_nutnet_mean$plot_order <- as.factor(cover_nutnet_mean$plot)

# Create nested random effect variables
cover_nutnet_mean$site_block <- interaction(cover_nutnet_mean$sitecode, cover_nutnet_mean$block_fac)
cover_nutnet_mean$site_block_plot <- interaction(data$sitecode, cover_nutnet_mean$block_fac, cover_nutnet_mean$plot)

# Create nested random effect variables
cover_nutnet_mean[["site_block"]] <- interaction(cover_nutnet_mean[["sitecode"]], cover_nutnet_mean[["block_fac"]])
cover_nutnet_mean[["site_block_plot"]] <- interaction(cover_nutnet_mean[["sitecode"]], cover_nutnet_mean[["block_fac"]], cover_nutnet_mean[["plot"]])



# Split data into training (60%), validation (20%), and test (20%) sets
set.seed(123)
train_index <- createDataPartition(cover_nutnet_mean$c4_proportion, p = 0.6, list = FALSE)
temp_data <- cover_nutnet_mean[-train_index,]
valid_index <- createDataPartition(temp_data$c4_proportion, p = 0.5, list = FALSE)

train_data <- cover_nutnet_mean[train_index,]
valid_data <- temp_data[valid_index,]
test_data <- temp_data[-valid_index,]

# Create the formula
formula <- c4_proportion ~ tmp + awi_pm_sr_yr + Ntrt_fac * Ptrt_fac * Ktrt_fac +
  sitecode + site_block + site_block_plot

# Train the random forest model
rf_model <- rfsrc(formula,
                  data = train_data,
                  ntree = 500,
                  importance = "permute")

# Print model summary
print(rf_model)

# Make predictions on the validation set
valid_predictions <- predict(rf_model, newdata = valid_data)$predicted

# Calculate RMSE for validation set
valid_rmse <- sqrt(mean((valid_data$c4_percent - valid_predictions)^2))
print(paste("Validation RMSE:", valid_rmse))

# Calculate R-squared for validation set
valid_r_squared <- 1 - sum((valid_data$c4_percent - valid_predictions)^2) / 
  sum((valid_data$c4_percent - mean(valid_data$c4_percent))^2)
print(paste("Validation R-squared:", valid_r_squared))

# Make predictions on the test set
test_predictions <- predict(rf_model, newdata = test_data)$predicted

# Calculate RMSE for test set
test_rmse <- sqrt(mean((test_data$c4_percent - test_predictions)^2))
print(paste("Test RMSE:", test_rmse))

# Calculate R-squared for test set
test_r_squared <- 1 - sum((test_data$c4_percent - test_predictions)^2) / 
  sum((test_data$c4_percent - mean(test_data$c4_percent))^2)
print(paste("Test R-squared:", test_r_squared))

# Plot predicted vs observed for validation set
ggplot(data.frame(predicted = valid_predictions, observed = valid_data$c4_percent), 
       aes(x = observed, y = predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Observed C4 Percent", y = "Predicted C4 Percent") +
  ggtitle("Predicted vs Observed C4 Percent (Validation Set)") +
  xlim(0, 100) + ylim(0, 100)

# Plot predicted vs observed for test set
ggplot(data.frame(predicted = test_predictions, observed = test_data$c4_percent), 
       aes(x = observed, y = predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Observed C4 Percent", y = "Predicted C4 Percent") +
  ggtitle("Predicted vs Observed C4 Percent (Test Set)") +
  xlim(0, 100) + ylim(0, 100)

# Variable importance
var_importance <- rf_model$importance
var_importance_df <- data.frame(
  variable = names(var_importance),
  importance = var_importance
)
var_importance_df <- var_importance_df[order(-var_importance_df$importance),]
print(var_importance_df)

# Plot variable importance
ggplot(var_importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Variables", y = "Importance") +
  theme_minimal() +
  ggtitle("Variable Importance")

# Partial dependence plots
par(mfrow = c(2, 3))
plot.variable(rf_model, xvar.names = c("tmp", "awi_pm_sr_yr", "Ntrt_fac", "Ptrt_fac", "Ktrt_fac"), 
              partial = TRUE)
par(mfrow = c(1, 1))

# Interaction plots
plot.variable(rf_model, xvar.names = c("Ntrt_fac", "Ptrt_fac"), partial = TRUE)
plot.variable(rf_model, xvar.names = c("Ntrt_fac", "Ktrt_fac"), partial = TRUE)
plot.variable(rf_model, xvar.names = c("Ptrt_fac", "Ktrt_fac"), partial = TRUE)






#######using bayeasian brms####
library(brms)

model <- brm(
  bf(c4_proportion2 ~ tmp + awi_pm_sr_yr + Ntrt_fac * Ptrt_fac * Ktrt_fac + (1|site_code/block_fac/plot),
     phi ~ tmp + awi_pm_sr_yr + Ntrt_fac * Ptrt_fac * Ktrt_fac  + (1|site_code/block_fac/plot),
     zi ~ tmp + awi_pm_sr_yr + Ntrt_fac * Ptrt_fac * Ktrt_fac  + (1|site_code/block_fac/plot)),
  family = zero_inflated_beta(),
  data = cover_nutnet_mean
)




#####Remove from the final code###########
C3.C4.max_cover <- lmer(log(max_cover)~ Ntrt_fac * Ptrt_fac * Ktrt_fac + 
                          ps_path + tmp + awi_pm_sr_yr + ps_path*tmp + 
                          ps_path*awi_pm_sr_yr + Ntrt_fac*ps_path + 
                          Ptrt_fac*ps_path +  Ktrt_fac*ps_path +
                          (1|plot) + (1|Taxon:site_code) + 
                          (1|Taxon:site_code:block_fac), data = cover_subset)



plot(resid(C3.C4.max_cover) ~ fitted(C3.C4.max_cover))
qqnorm(residuals(C3.C4.max_cover))
qqline(residuals(C3.C4.max_cover))
densityPlot(residuals(C3.C4.max_cover))
#ks.test(residuals(C3.C4.max_cover))
outlierTest(C3.C4.max_cover)


###################model result#########################
summary(C3.C4.max_cover)
Anova(C3.C4.max_cover)
r.squaredGLMM(C3.C4.max_cover) ## 0.73

##########################################################################
# percentage increase of max cover in C3 compared to C4
emmeans(C3.C4.max_cover, ~ps_path)
temperature <- as.data.frame(cld(emmeans(C3.C4.max_cover, ~ps_path*tmp)))
emmeans(C3.C4.max_cover, pairwise~ps_path:awi_pm_sr_yr)
nitrogen <- as.data.frame(emmeans(C3.C4.max_cover, ~Ntrt_fac:ps_path))
nitrogen.phosphorous <- as.data.frame(emmeans(C3.C4.max_cover, ~Ntrt_fac:Ktrt_fac, 
                                              "ps_path"))

cld(emmeans(C3.C4.max_cover, ~Ntrt_fac:Ptrt_fac))


emmeans(C3.C4.max_cover, pairwise~Ntrt_fac:Ktrt_fac, "ps_path")


emmeans(C3.C4.max_cover, pairwise~ps_path, type = "response")
#AIC(C3.C4.max_cover) # 113.7041

###following Ben Bolker suggestion of having seperate different models
######zero-inflated data
C4.percent.model.zero.nonzero <- glmmTMB((c4_proportion==0) ~ tmp +
                                           aridity.cateogry + Ntrt_fac * Ptrt_fac * Ktrt_fac +
                                           (1|site_code) + (1|site_code:block_fac) + 
                                           (1|site_code:block_fac:plot),data = cover_nutnet_mean, 
                                         family = binomial(), na.action = na.omit)


Anova(C4.percent.model.zero.nonzero)
fixef(C4.percent.model.zero.nonzero)
summary(C4.percent.model.zero.nonzero)

C4.percent.model.nonzero.nonone <- glmmTMB(c4_proportion ~ tmp +
                                             aridity.cateogry + Ntrt_fac * Ptrt_fac * Ktrt_fac + 
                                             (1|site_code) + (1|site_code:block_fac) + 
                                             (1|site_code:block_fac:plot),
                                           data = subset(cover_nutnet_mean, 0<c4_proportion & c4_proportion<1),
                                           family = beta_family(), na.action = na.omit)


Anova(C4.percent.model.nonzero.nonone)
fixef(C4.percent.model.nonzero.nonone)
summary(C4.percent.model.nonzero.nonone)

C4.percent.model.nonone.one <- glmmTMB((c4_proportion==1) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp +
                                         aridity.cateogry + (1|site_code) + (1|site_code:block_fac) + 
                                         (1|site_code:block_fac:plot),
                                       data = cover_nutnet_mean,
                                       family = binomial(), na.action = na.omit)


Anova(C4.percent.model.nonone.one)
fixef(C4.percent.model.nonone.one)
summary(C4.percent.model.nonone.one)



