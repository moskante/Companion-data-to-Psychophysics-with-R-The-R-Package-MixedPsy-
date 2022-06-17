## ---- setup
# install.packages(c("MixedPsy", "tidyverse"))
library(MixedPsy)
library(lme4)
library(tidyverse)

theme_set(theme_classic() + theme(text=element_text(size=11), 
                                  legend.title = element_text(size = 11),
                                  legend.text=element_text(size=11),
                                  legend.key.size = unit(3,"line"),
                                  strip.text=element_text(size=9),
                                  panel.grid.minor.x = element_blank(),
                                  panel.grid.minor.y = element_blank(),
                                  strip.background = element_blank(),
                                  axis.title=element_text(size=11), 
                                  axis.text=element_text(size=9))
)


# Dataset ----
# ---- import
response <- read_csv("exampledata.csv")
# ---- preprocess
response.summary <- response %>%
  group_by(X) %>%
  summarise(Total = n(),
            Longer = sum(resp == 1),
            Shorter = sum(resp == 0))
# ---- 

# Results: Fitting GLMs with R ----
## The Psychometric Function ----
## ---- select
data.S1 <- filter(simul_data, Subject == "S1")

## ---- fitGLM
glm.S1 <- glm(formula = cbind(Longer, Total - Longer) ~ X, 
              family = binomial(link = "probit"), data = data.S1)

## ---- GLMcoeff
summary(glm.S1)$coefficients

## ---- PsychDelta
PsychDelta(glm.S1)

## ---- plot1
plotP1 <- PsychPlot(glm.S1, showData = TRUE, ps.lab = "glm S1") 
## ----

# Alternatively, estimate parameters without explicit gall to glm with PsychFunction()
data.S2 <- subset(simul_data, Subject == "S2")
psych.S2 <- PsychFunction(ps.formula = cbind(Longer, Total - Longer) ~ X, 
                ps.link = "probit", ps.data = data.S2)
psych.S2$estimate
# Add response of second participant in previous plot
plotP2 <- PsychPlot(psych.S2$model, addTo = plotP1, ps.lab = "S2")
## ----

## Second-level analysis ----
## ---- 2nd_level_fun
fit.glm <- function(data){
  mod <- glm(formula = cbind(Longer, Total - Longer) ~ X, 
             family = binomial(link = "probit"), data = data)
  psejnd <- PsychDelta(mod)
  return(list(pse = psejnd["pse","Estimate"], jnd = psejnd["jnd","Estimate"]))
  }
## ---- 2nd_level
data_bySub <- simul_data %>%
  group_by(Subject) %>%
  nest() 
data_bySub$psejnd <- map(data_bySub$data, fit.glm)
data_bySub <- unnest_wider(data_bySub, psejnd) %>%
  ungroup()
## ---- 

# This is a more compact form of the code above
psejnd_simul <- simul_data %>%
  group_by(Subject) %>%
  nest() %>%
  mutate(psejnd = map(data, fit.glm)) %>%
  unnest_wider(psejnd) %>%
  ungroup()
# For a dataset with a factorial predictor, the syntax is similar
estim_multi <- function(x){
  psejnd <- PsychFunction(ps.formula = cbind(faster, slower) ~ speed, ps.link = "probit", ps.data = x)$estimate
  return(list(pse = psejnd["pse","Estimate"], jnd = psejnd["jnd","Estimate"]))}
# Note the multiple grouping factor
table_exp <- vibro_exp3 %>%
  group_by(subject, vibration) %>%
  nest() %>%
  mutate(psejnd = map(data, estim_multi)) %>%
  unnest_wider(psejnd)
## ---- 

## ---- t-test
data_summary <- summarise(data_bySub, 
                           mean.pse = mean(pse),
                           sd.pse = sd(pse))
pse.t_test <- t.test(data_bySub$pse, mu = 80)
## ----


# Results: Fitting GLMMs with R ----
## Univariable GLMM ----
## ---- rnd.all 
glmm.uni <- glmer(formula = cbind(Longer, Total - Longer) ~ X + (1 + X| Subject),
                           family = binomial(link = "probit"), data = simul_data)

## ---- rnd.interc
uni.rndInterc <- glmer(formula = cbind(Longer, Total - Longer) ~ X + (1 | Subject),
                      family = binomial(link = "probit"), data = simul_data)

## ---- model.selection
anova(glmm.uni, uni.rndInterc, test = "chisq")
## ----

## Similarly, we could evaluate a model with only random slope
uni.rndSlope <- glmer(formula = cbind(Longer, Total - Longer) ~ X + (0 + X | Subject),
                     family = binomial(link = "probit"), data = simul_data)
anova(glmm.uni, uni.rndSlope, test = "chisq")
## ----

# Plot and estimates
## ---- xplode.uni
xplode.uni <- xplode(model = uni.rndInterc, name.cont = "X")

## ---- Mixplot.uni
MixPlot(xplode.uni) 

## ---- MixDelta.uni
MixDelta(xplode.uni)
## ----

start_time <- Sys.time()
## ---- boot.uni
BootEstim.uni <- pseMer(uni.rndInterc, B = 500)
## ----
end_time <- Sys.time()
time.boot1 <- end_time - start_time
## ---- boot.uni.summary
BootEstim.uni$summary
## ----

## Multivariable ----
## ---- glmm.vibro
glmm.multi <- glmer(cbind(faster, slower) ~ speed + vibration + speed:vibration + 
                      (1+speed | subject),
                    family = binomial(link = "probit"), 
                    data = vibro_exp3)

## ---- summary.vibro
summary(glmm.multi)$coefficients

## ---- plot.vibro
#create xplode object
xplode.multi <- xplode(glmm.multi, name.cont = "speed", name.factor = "vibration")
#plot
MixPlot(xplode.multi)
#delta method estimation
MixDelta(xplode.multi)
## ----
# bootstrap estimate of parameters
## ---- mod2fun
fun2mod <- function(mer.obj){
  #allocate space
  jndpse = vector(mode = "numeric", length = 6)
  #name parameters
  names(jndpse) = c("PSE 0", "PSE 32","JND 0", "JND 32", "diff PSE", "diff JND")
  #define statistics of interest
  jndpse[1] = -fixef(mer.obj)[1]/fixef(mer.obj)[2] #pse_0
  jndpse[2] = -(fixef(mer.obj)[1]+fixef(mer.obj)[3])/(fixef(mer.obj)[2]+ fixef(mer.obj)[4]) #pse_32
  jndpse[3] = qnorm(0.75)/fixef(mer.obj)[2] #jnd_0
  jndpse[4] = qnorm(0.75)/(fixef(mer.obj)[2]+ fixef(mer.obj)[4]) #jnd_32
  jndpse[5] = jndpse[2] - jndpse[1] # diff PSE
  jndpse[6] = jndpse[4] - jndpse[3] # diff JND
  # return
  return(jndpse)
}
## ----

start_time <- Sys.time()
## ---- boot.multi
BootEstim.multi <- pseMer(glmm.multi, B = 500, FUN = fun2mod)
## ----
end_time <- Sys.time()
time.boot2 <- end_time - start_time

## ---- boot.multi.summary
BootEstim.multi$summary
## ----

##### APPENDIX ##### 
#################### 

## ---- plotsGLM

ggplot(simul_data, aes(X, Longer/Total, color = Subject)) + geom_point() + 
  geom_smooth(method = "glm", method.args=list(family=binomial(link="probit")), se = FALSE)

## ---- plotsGLMM

formula.multi <- cbind(faster, slower) ~ speed * vibration  + (1 + speed| subject)
multi.mod <- glmer(formula = formula.multi, family = binomial(link = "probit"), 
                   data = vibro_exp3)

#create a new dataframe: names should match those of the originary dataframe
longData = expand.grid(speed = unique(c(pretty(vibro_exp3$speed, 1000), vibro_exp3$speed)),
                       subject = levels(vibro_exp3$subject), 
                       vibration = levels(vibro_exp3$vibration))

#add preditions
longData$predict = predict(multi.mod, type = "response", newdata = longData)

# alternatively, with modelr
library(modelr)
longData <- data_grid(vibro_exp3, speed = seq_range(speed, 200), .model = multi.mod) %>%
  add_predictions(multi.mod, var = "predict", type = "response") 

ggplot(longData, aes(x = speed, y = predict, color = subject)) +
  geom_line() +
  geom_point(data = vibro_exp3, aes(y = faster/(faster+slower)))+
  facet_wrap(~ vibration, ncol = 3) 


