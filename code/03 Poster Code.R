#' Code for generating analysis for poster project
#' Wesley Zuidema and Ellen Ahlness
#' MLE Poster Project

# Load library packages
library(ggplot2)
library(knitr)
library(MASS)
library(nlme)
library(boot)            # For cv.glm()
library(separationplot)  # For separation plot
library(pscl)            # Alternative PCP code
library(verification)    # For ROC area
library(tile)            # For some graphics; used by plot.binPredict()
library(RColorBrewer)    # For nice colors
source("binaryGOF.R")    # Percent correctly predicted and concordance indexes
source("binPredict.R")   # Code for making predicted vs actual plots
library(simcf)
library(xtable)
library(VGAM)
library(dplyr)
select <- dplyr::select


# Load data
load("data/wvs_data.rdata")
wvs_data <- WVS_Longitudinal_1981_2014_R_v2015_04_18
rm(WVS_Longitudinal_1981_2014_R_v2015_04_18)

# Rename variables we are using to make them easier to read
wvs_data$immigrationpref <- wvs_data$E143
wvs_data$age <- wvs_data$X003
wvs_data$incomelevel <- wvs_data$X047
wvs_data$female <- wvs_data$X001
wvs_data$maritalstatus <- wvs_data$X007
wvs_data$education <- wvs_data$X025R

wvs_data$religious <- wvs_data$F034
wvs_data$worldcitizen <- wvs_data$G019
wvs_data$nationalpride <- wvs_data$G006

wvs_data$infofriends <- wvs_data$E254

wvs_data$wave <- as.factor(wvs_data$S002)
wvs_data$country <- as.factor(wvs_data$S003)

# Have data for 20 - Andorra, 36 - Australia, 348 - Hungary, 360 - Indonesia, 
# 484 - Mexico, 840 - USA, and 858 - Uraguay
wvs_data$partypref <- wvs_data$E256
#wvs_filter <- filter(wvs_data, country == 840, partypref > 7)
#unique(wvs_filter$partypref)
#ggplot(wvs_filter) + geom_jitter(aes(x = country, y = partypref, col = country))

# Create binary variable to indicate support of far right party
wvs_filter <- filter(wvs_data, partypref > 7)
wvs_filter$farrightpref <- 0
wvs_filter <- wvs_filter %>%
  mutate(farrightpref = ifelse(partypref == 36003 | partypref == 36006 | partypref == 36008 | partypref == 36009 | partypref == 348001 | partypref == 348004 | partypref == 360007, 1, 0))
wvs_filter <- wvs_filter %>%
  filter(country != 20) %>%
  filter(country != 484) %>%
  filter(country != 840) %>%
  filter(country != 858)

# Recode data to remove non-responses

wvs_recode <- wvs_filter %>%
  filter(immigrationpref > 0) %>%
  filter(age > 0) %>%
  filter(incomelevel > 0) %>%
  filter(female > 0) %>%
  filter(maritalstatus > 0) %>%
  filter(education > 0) %>%
  filter(religious > 0) %>%
  filter(worldcitizen > 0) %>%
  filter(nationalpride > 0) %>%
  filter(infofriends > 0) %>%
  mutate(married = ifelse(maritalstatus == 1, 1, 0),
         domesticpartner = ifelse(maritalstatus == 2, 1, 0),
         divorced = ifelse(maritalstatus == 3, 1, 0),
         separated = ifelse(maritalstatus == 4, 1, 0),
         widowed = ifelse(maritalstatus == 5, 1, 0),
         nevermarried = ifelse(maritalstatus == 6, 1, 0),
         australia = ifelse(country == 36, 1, 0),
         hungary = ifelse(country == 348, 1, 0))
  

wvs_final <- wvs_recode %>%
  select(wave, country, farrightpref, immigrationpref, age, incomelevel, female, maritalstatus, education, religious, worldcitizen, nationalpride, infofriends, married, domesticpartner, divorced, separated, widowed, nevermarried, australia, hungary)


# Build model and estimate with least squares first for reference
# All models include country fixed effects, Indonesia is baseline
model_a <-  farrightpref ~ immigrationpref + australia + hungary #Binary Model
model_b <-  farrightpref ~ immigrationpref + religious + worldcitizen + nationalpride + australia + hungary#Preferences
model_c <- farrightpref ~ age + incomelevel + female + education + married + domesticpartner + divorced + separated + widowed + australia + hungary #Characteristics
model_d <- farrightpref ~ immigrationpref + age + incomelevel + female + education + religious + worldcitizen + nationalpride + married + domesticpartner + divorced + separated + widowed + australia + hungary #Full Model

model_state_specific <- model_d <- farrightpref ~ immigrationpref + age + incomelevel + female + education + religious + worldcitizen + nationalpride + married + domesticpartner + divorced + separated + widowed

model_optim <- farrightpref ~ immigrationpref + age + incomelevel + female + religious + worldcitizen + nationalpride + married + domesticpartner + divorced + separated + widowed + australia + hungary #Full Model

ls.result <- lm(model_optim, wvs_final)


# Binary logit model

# Define logit function
llk.logit <- function(param,y,xcovariates) {
  x <- as.matrix(xcovariates)
  y <- as.matrix(y)
  os <- rep(1,nrow(x))
  x <- cbind(os,x)
  b <- param[ 1 : ncol(x) ]
  xb <- x%*%b # systematic components for the pi
  
  -sum((-y * log(1 + exp(-xb))) - ((1-y) * log(1 + exp(xb))))
  
}

# Define the parameters we will use for optim
y <-  wvs_final$farrightpref
xcovariates <- cbind(wvs_final$immigrationpref, wvs_final$age, wvs_final$incomelevel, wvs_final$female, wvs_final$religious, wvs_final$worldcitizen, wvs_final$nationalpride, wvs_final$married, wvs_final$domesticpartner, wvs_final$divorced, wvs_final$separated, wvs_final$widowed, wvs_final$australia, wvs_final$hungary)
stval <- c(coef(ls.result))

# Run optim

# Use optim directly to get MLE
logit.farright <- optim(stval, llk.logit, method="BFGS", x=xcovariates, y=y, hessian=TRUE)
pe <- logit.farright$par                # point estimates
vc <- solve(logit.farright$hessian)     # var-cov matrix
se <- sqrt(diag(vc))                 # standard errors
ll <- -logit.farright$value             # likelihood at maximum

# Create table of results
logit.table <- rbind(pe, se)
logit.table <- t(logit.table)
colnames(logit.table) <- c("Estimate", "Standard Error")
rownames(logit.table) <- c("(Intercept)", "Immigration Preference", "Age", "Income Level", "Female", "Religious", "World Citizen", "National Pride", "Married", "Domestic Partner", "Divorced", "Separated", "Widowed", "Australia", "Hungary")

# GLM for comparing models
glm_a <- glm(model_a, family = binomial, data = wvs_final)
glm_b <- glm(model_b, family = binomial, data = wvs_final)
glm_c <- glm(model_c, family = binomial, data = wvs_final)
glm_d <- glm(model_d, family = binomial, data = wvs_final)

glm_indonesia <- glm(model_state_specific, family = binomial, data = filter(wvs_final, country == 360))
glm_hungary <- glm(model_state_specific, family = binomial, data = filter(wvs_final, country == 348))
glm_australia <- glm(model_state_specific, family = binomial, data = filter(wvs_final, country == 36))

summary(glm_d)

summary(glm_indonesia)
summary(glm_hungary)
summary(glm_australia)

col <- brewer.pal(5, "Set1")
blue <- col[2]
orange <- col[5]

## Goodness of Fit Tests
# Length of each model
k.m1 <- length(glm_a$coefficients)
k.m2 <- length(glm_b$coefficients)
k.m3 <- length(glm_c$coefficients)
k.m4 <- length(glm_d$coefficients)

# Percent Correctly Predicted
pcp.null <- pcp.glm(glm_a, wvs_final$farrightpref, type="null")
pcp.m1 <- pcp.glm(glm_a, wvs_final$farrightpref, type="model")
pcp.m2 <- pcp.glm(glm_d, wvs_final$farrightpref, type="model")
pcpi.m1 <- pcp.glm(glm_a, wvs_final$farrightpref, type="improve")
pcpi.m2 <- pcp.glm(glm_d, wvs_final$farrightpref, type="improve")

# Regression Table
glm_a_sum <- summary(glm_a)
glm_a_coef <- glm_a_sum$coefficients
table_a <- glm_a_coef[, 1:2]
glm_d_sum <- summary(glm_d)
glm_d_coef <- glm_d_sum$coefficients
table_d <- glm_d_coef[, 1:2]

       
## Simulations
# Model parameters
pe.glm <- glm_d$coefficients  # point estimates
vc.glm <- vcov(glm_d)         # var-cov matrix

# Simulate parameter distributions
sims <- 10000
simbetas <- mvrnorm(sims, pe.glm, vc.glm)

# Create example counterfactuals
xhyp <- cfMake(model_d, wvs_final, nscen=10)

xhyp <- cfName(xhyp,"Opposes Immigration", scen=1)
xhyp <- cfChange(xhyp, "immigrationpref", x=4, xpre=1, scen=1)

xhyp <- cfName(xhyp, "National Pride", scen=2)
xhyp <- cfChange(xhyp, "nationalpride", x=1, xpre=4, scen=2)

xhyp <- cfName(xhyp, "World Citizen", scen=3)
xhyp <- cfChange(xhyp, "worldcitizen", x=1, xpre=4, scen=3)

xhyp <- cfName(xhyp, "Religious", scen=4)
xhyp <- cfChange(xhyp, "religious", x=1, xpre=3, scen=4)

xhyp <- cfName(xhyp, "Age + 1sd = 56", scen=5)
xhyp <- cfChange(xhyp, "age",
                 x=mean(na.omit(wvs_final$age))+sd(na.omit(wvs_final$age)),
                 xpre=mean(na.omit(wvs_final$age)),
                 scen=5)

xhyp <- cfName(xhyp, "Age - 1sd = 24", scen=6)
xhyp <- cfChange(xhyp, "age",
                 x=mean(na.omit(wvs_final$age))-sd(na.omit(wvs_final$age)),
                 xpre=mean(na.omit(wvs_final$age)),
                 scen=6)

xhyp <- cfName(xhyp, "Male", scen=7)
xhyp <- cfChange(xhyp, "female", x=1, xpre=2, scen=7)

xhyp <- cfName(xhyp, "Income (1st to 3rd QT)", scen=8)
xhyp <- cfChange(xhyp, "incomelevel", x=7, xpre=3, scen=8)

xhyp <- cfName(xhyp, "Widowed", scen=9)
xhyp <- cfChange(xhyp, "widowed", x=1, xpre=0, scen=9)

xhyp <- cfName(xhyp, "Divorced", scen=10)
xhyp <- cfChange(xhyp, "female", x=1, xpre=0, scen=10)

# Simulate expected probabilities for all scenarios
logitsims <- logitsimev(xhyp, simbetas, ci=0.95)
logitfds <- logitsimfd(xhyp, simbetas, ci=0.95)

# Get 3 nice colors for traces
col <- brewer.pal(3,"Dark2")

# Plot predicted probabilities for all four categories, sorted by size
sorted <- order(logitsims$pe)
scenNames <- row.names(xhyp$x)

trace1 <- ropeladder(x = logitsims$pe[sorted],
                     lower = logitsims$lower[sorted],
                     upper = logitsims$upper[sorted],
                     labels = scenNames[sorted],
                     size=0.5,
                     lex=1.5,
                     lineend="square",
                     plot=1
)

tile(trace1,
     limits = c(0,0.5),
     gridlines = list(type="xt"),
     topaxis=list(add=TRUE, at=c(0,0.1,0.2,0.3,0.4,0.5)),
     plottitle=list(labels="Expected Values"),
     xaxistitle=list(labels="Probability of far right support"),
     width=list(spacer=3),
     height = list(plottitle=3,xaxistitle=3.5,topaxistitle=3.5),
     output = list(file="figures/expected_values")
)

# Plot First Difference
sorted <- order(logitfds$pe)
scenNames <- row.names(xhyp$x)
# Set up lineplot trace of relative risk
trace2 <- ropeladder(x = logitfds$pe[sorted],
                     lower = logitfds$lower[sorted],
                     upper = logitfds$upper[sorted],
                     labels = row.names(xhyp$x[sorted]),
                     size=0.5,
                     lex=1.5,
                     lineend="square",
                     baseline=list(at=0),
                     plot=1
)


# Plot traces using tile
tile(trace2,
     limits = c(-.40,0.25),
     gridlines = list(type="xt"),
     topaxis=list(add=TRUE, at=c(-.3,-.2,-.1,0,0.1,0.2)),
     plottitle=list(labels="First Differences"),
     xaxistitle=list(labels="Probability of far right support"),
     width=list(spacer=3),
     height = list(plottitle=3,xaxistitle=3.5,topaxistitle=3.5),
     output=list(file="figures/firstdifferences")
)

### Cross Validation
loocv <- function (obj) {
  data <- obj$data
  m <- dim(data)[1]
  form <- formula(obj)
  fam <- obj$family$family
  loo <- rep(NA, m)
  for (i in 1:m) {
    i.glm <- glm(form, data = data[-i, ], family = fam)
    loo[i] <- predict(i.glm, newdata = data[i,], family = fam, type = "response")
  }
  loo
}

# LOOCV for models 1 and 2
predCV_a <- loocv(glm_a)
predCV_b <- loocv(glm_b)
predCV_c <- loocv(glm_c)
predCV_d <- loocv(glm_d)

# Make cross-validated AVP and ROC plots; note use of newpred input in binPredict
binnedM1cv <- binPredict(glm_a, newpred=predCV_a, col=orange, label="M1: Binary Model", quantiles=TRUE, bins = 5)
binnedM4cv <- binPredict(glm_d, newpred=predCV_d, col=orange, label="M4: Full Model", quantiles=TRUE, bins = 5)

binnedM2cv <- binPredict(glm_b, newpred=predCV_b, col=blue, label="M2: Preferences Only", quantiles=TRUE, bins = 5)
binnedM3cv <- binPredict(glm_c, newpred=predCV_c, col=blue, label="M3: Characteristics Only", quantiles=TRUE, bins = 5)

plot(binnedM2cv, binnedM4cv, display=c("roc"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/ROC1 - Preferences vs FULL")

plot(binnedM2cv, binnedM4cv, display=c("avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/AVP1 - Preferences vs FULL")

plot(binnedM3cv, binnedM4cv, display=c("roc"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/ROC2 - Characateristics vs FULL")

plot(binnedM3cv, binnedM4cv, display=c("avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/AVP2 - Characateristics vs FULL")

plot(binnedM2cv, binnedM4cv, display=c("roc", "avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/Preferences vs FULL")
plot(binnedM3cv, binnedM4cv, display=c("roc", "avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/Characateristics vs FULL")


## Testing for heterogenous age effects with vgam
model_vgam <- farrightpref ~ immigrationpref + s(age) + incomelevel + female + education + religious + worldcitizen + nationalpride + married + domesticpartner + divorced + separated + widowed + australia + hungary #Full Model
fit <- vgam(model_vgam, binomialff, data = wvs_final)
summary(fit)
pfit <- predict(fit, type = "terms", raw = TRUE, se = TRUE)
head(pfit$fitted)
head(pfit$se.fit)
pfit$df
pfit$sigma
